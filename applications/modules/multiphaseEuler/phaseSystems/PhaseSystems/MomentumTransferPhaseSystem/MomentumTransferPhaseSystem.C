/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "MomentumTransferPhaseSystem.H"

#include "dragModel.H"
#include "virtualMassModel.H"
#include "liftModel.H"
#include "wallLubricationModel.H"
#include "turbulentDispersionModel.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"

#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcFlux.H"
#include "fvcSnGrad.H"
#include "fvcMeshPhi.H"
#include "fvcReconstruct.H"

#include "pimpleNoLoopControl.H"

// * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * //

template<class BasePhaseSystem>
void Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::addDmdtUfs
(
    const phaseSystem::dmdtfTable& dmdtfs,
    phaseSystem::momentumTransferTable& eqns
)
{
    forAllConstIter(phaseSystem::dmdtfTable, dmdtfs, dmdtfIter)
    {
        const phaseInterface interface(*this, dmdtfIter.key());

        const volScalarField& dmdtf = *dmdtfIter();
        const volScalarField dmdtf21(posPart(dmdtf));
        const volScalarField dmdtf12(negPart(dmdtf));

        phaseModel& phase1 = this->phases()[interface.phase1().name()];
        phaseModel& phase2 = this->phases()[interface.phase2().name()];

        if (!phase1.stationary())
        {
            *eqns[phase1.name()] +=
                dmdtf21*phase2.U() + fvm::Sp(dmdtf12, phase1.URef());
        }

        if (!phase2.stationary())
        {
            *eqns[phase2.name()] -=
                dmdtf12*phase1.U() + fvm::Sp(dmdtf21, phase2.URef());
        }
    }
}


template<class BasePhaseSystem>
void Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::addTmpField
(
    tmp<surfaceScalarField>& result,
    const tmp<surfaceScalarField>& field
) const
{
    if (result.valid())
    {
        result.ref() += field;
    }
    else
    {
        if (field.isTmp())
        {
            result = field;
        }
        else
        {
            result = field().clone();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::
MomentumTransferPhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh)
{
    this->generateInterfacialModels(dragModels_);
    this->generateInterfacialModels(virtualMassModels_);
    this->generateInterfacialModels(liftModels_);
    this->generateInterfacialModels(wallLubricationModels_);
    this->generateInterfacialModels(turbulentDispersionModels_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::
~MomentumTransferPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::momentumTransfer()
{
    // Create a momentum transfer matrix for each phase
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr
    (
        new phaseSystem::momentumTransferTable()
    );

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    forAll(this->movingPhases(), movingPhasei)
    {
        const phaseModel& phase = this->movingPhases()[movingPhasei];

        eqns.insert
        (
            phase.name(),
            new fvVectorMatrix(phase.U(), dimMass*dimVelocity/dimTime)
        );
    }

    // Initialise Kds table if required
    if (!Kds_.size())
    {
        forAllConstIter
        (
            dragModelTable,
            dragModels_,
            dragModelIter
        )
        {
            const phaseInterface& interface = dragModelIter()->interface();

            Kds_.insert
            (
                dragModelIter.key(),
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName("Kd", interface.name()),
                        this->mesh().time().name(),
                        this->mesh()
                    ),
                    this->mesh(),
                    dimensionedScalar(dragModel::dimK, 0)
                )
            );
        }
    }

    // Update the drag coefficients
    forAllConstIter
    (
        dragModelTable,
        dragModels_,
        dragModelIter
    )
    {
        *Kds_[dragModelIter.key()] = dragModelIter()->K();

        // Zero-gradient the drag coefficient to boundaries with fixed velocity
        const phaseInterface& interface = dragModelIter()->interface();
        volScalarField& K = *Kds_[dragModelIter.key()];

        forAll(K.boundaryField(), patchi)
        {
            if
            (
                (
                    !interface.phase1().stationary()
                 && interface.phase1().U()()
                   .boundaryField()[patchi].fixesValue()
                )
             && (
                    !interface.phase2().stationary()
                 && interface.phase2().U()()
                   .boundaryField()[patchi].fixesValue()
                )
            )
            {
                K.boundaryFieldRef()[patchi] =
                    K.boundaryField()[patchi].patchInternalField();
            }
        }
    }

    // Initialise Vms table if required
    if (!Vms_.size())
    {
        forAllConstIter
        (
            virtualMassModelTable,
            virtualMassModels_,
            virtualMassModelIter
        )
        {
            const phaseInterface& interface =
                virtualMassModelIter()->interface();

            Vms_.insert
            (
                interface,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName("Vm", interface.name()),
                        this->mesh().time().name(),
                        this->mesh()
                    ),
                    this->mesh(),
                    dimensionedScalar(virtualMassModel::dimK, 0)
                )
            );
        }
    }

    // Update the virtual mass coefficients
    forAllConstIter
    (
        virtualMassModelTable,
        virtualMassModels_,
        virtualMassModelIter
    )
    {
        *Vms_[virtualMassModelIter.key()] = virtualMassModelIter()->K();
    }

    // Add the virtual mass force
    forAllConstIter(VmTable, Vms_, VmIter)
    {
        const volScalarField& Vm(*VmIter());
        const phaseInterface interface(*this, VmIter.key());

        forAllConstIter(phaseInterface, interface, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            if (!phase.stationary())
            {
                fvVectorMatrix& eqn = *eqns[phase.name()];

                const volVectorField& U = eqn.psi();

                const surfaceScalarField& phi = phase.phiRef();
                const tmp<surfaceScalarField> taphi(fvc::absolute(phi, U));
                const surfaceScalarField& aphi(taphi());

                const volScalarField VmPhase
                (
                    (otherPhase/max(otherPhase, otherPhase.residualAlpha()))
                   *Vm
                );

                eqn -=
                    VmPhase
                   *(
                        fvm::ddt(U)
                      + fvm::div(aphi, U) - fvm::Sp(fvc::div(aphi), U)
                      - otherPhase.DUDt()
                    )
                  + this->MRF_.DDt(VmPhase, U - otherPhase.U());
            }
        }
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::momentumTransferf()
{
    // Create a momentum transfer matrix for each phase
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr
    (
        new phaseSystem::momentumTransferTable()
    );

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    forAll(this->movingPhases(), movingPhasei)
    {
        const phaseModel& phase = this->movingPhases()[movingPhasei];

        eqns.insert
        (
            phase.name(),
            new fvVectorMatrix(phase.U(), dimMass*dimVelocity/dimTime)
        );
    }

    // Initialise Kdfs table if required
    if (!Kdfs_.size())
    {
        forAllConstIter
        (
            dragModelTable,
            dragModels_,
            dragModelIter
        )
        {
            const phaseInterface& interface = dragModelIter()->interface();

            Kdfs_.insert
            (
                dragModelIter.key(),
                new surfaceScalarField
                (
                    IOobject
                    (
                        IOobject::groupName("Kdf", interface.name()),
                        this->mesh().time().name(),
                        this->mesh()
                    ),
                    this->mesh(),
                    dimensionedScalar(dragModel::dimK, 0)
                )
            );
        }
    }

    // Update the drag coefficients
    forAllConstIter
    (
        dragModelTable,
        dragModels_,
        dragModelIter
    )
    {
        *Kdfs_[dragModelIter.key()] = dragModelIter()->Kf();
    }

    // Initialise Vms table if required
    if (!Vms_.size())
    {
        forAllConstIter
        (
            virtualMassModelTable,
            virtualMassModels_,
            virtualMassModelIter
        )
        {
            const phaseInterface& interface =
                virtualMassModelIter()->interface();

            Vms_.insert
            (
                interface,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName("Vm", interface.name()),
                        this->mesh().time().name(),
                        this->mesh()
                    ),
                    this->mesh(),
                    dimensionedScalar(virtualMassModel::dimK, 0)
                )
            );
        }
    }

    // Update the virtual mass coefficients
    forAllConstIter
    (
        virtualMassModelTable,
        virtualMassModels_,
        virtualMassModelIter
    )
    {
        *Vms_[virtualMassModelIter.key()] = virtualMassModelIter()->K();
    }

    // Create U & grad(U) fields
    PtrList<fvVectorMatrix> UgradUs(this->phaseModels_.size());
    forAll(this->phaseModels_, phasei)
    {
        const phaseModel& phase = this->phaseModels_[phasei];

        if (!phase.stationary())
        {
            const volVectorField& U = phase.URef();
            const surfaceScalarField& phi = phase.phiRef();

            const tmp<surfaceScalarField> taphi(fvc::absolute(phi, U));
            const surfaceScalarField& aphi(taphi());

            UgradUs.set
            (
                phasei,
                new fvVectorMatrix
                (
                    fvm::div(aphi, U) - fvm::Sp(fvc::div(aphi), U)
                  + this->MRF().DDt(U)
                )
            );
        }
    }

    // Add the virtual mass force
    forAllConstIter(VmTable, Vms_, VmIter)
    {
        const volScalarField& Vm(*VmIter());
        const phaseInterface interface(*this, VmIter.key());

        forAllConstIter(phaseInterface, interface, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            if (!phase.stationary())
            {
                const volScalarField VmPhase
                (
                    (otherPhase/max(otherPhase, otherPhase.residualAlpha()))
                   *Vm
                );

                *eqns[phase.name()] -=
                    VmPhase
                   *(
                        UgradUs[phase.index()]
                      - (UgradUs[otherPhase.index()] & otherPhase.U())
                    );
            }
        }
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::surfaceScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::KdVmfs() const
{
    PtrList<surfaceScalarField> KdVmfs(this->phaseModels_.size());

    // Add the implicit part of the drag force
    forAllConstIter(KdfTable, Kdfs_, KdfIter)
    {
        const surfaceScalarField& Kf(*KdfIter());
        const phaseInterface interface(*this, KdfIter.key());

        forAllConstIter(phaseInterface, interface, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            const surfaceScalarField alphaf
            (
                fvc::interpolate(otherPhase)
            );

            addField
            (
                phase,
                "KdVmf",
                (alphaf/max(alphaf, otherPhase.residualAlpha()))
               *Kf,
                KdVmfs
            );
        }
    }

    // Add the implicit part of the virtual mass force
    forAllConstIter(VmTable, Vms_, VmIter)
    {
        const volScalarField& Vm(*VmIter());
        const phaseInterface interface(*this, VmIter.key());

        forAllConstIter(phaseInterface, interface, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            const volScalarField VmPhase
            (
                (otherPhase/max(otherPhase, otherPhase.residualAlpha()))
               *Vm
            );

            addField(phase, "KdVmf", byDt(fvc::interpolate(VmPhase)), KdVmfs);
        }
    }

    this->fillFields("KdVmf", dimDensity/dimTime, KdVmfs);

    return KdVmfs;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::surfaceScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::Fs() const
{
    PtrList<surfaceScalarField> Fs(this->phaseModels_.size());

    // Add the lift force
    forAllConstIter
    (
        liftModelTable,
        liftModels_,
        liftModelIter
    )
    {
        const phaseInterface& interface = liftModelIter()->interface();

        const volVectorField F(liftModelIter()->F());

        addField
        (
            interface.phase1(),
            "F",
            fvc::flux(F),
            Fs
        );
        addField
        (
            interface.phase2(),
            "F",
           -fvc::flux(F),
            Fs
        );
    }

    // Add the wall lubrication force
    forAllConstIter
    (
        wallLubricationModelTable,
        wallLubricationModels_,
        wallLubricationModelIter
    )
    {
        const phaseInterface& interface =
            wallLubricationModelIter()->interface();

        const volVectorField F(wallLubricationModelIter()->F());

        addField
        (
            interface.phase1(),
            "F",
            fvc::flux(F),
            Fs
        );
        addField
        (
            interface.phase2(),
            "F",
           -fvc::flux(F),
            Fs
        );
    }

    // Add the phase pressure
    forAll(this->phaseModels_, phasei)
    {
        const phaseModel& phase = this->phaseModels_[phasei];

        const tmp<volScalarField> pPrime(phase.pPrime());

        addField
        (
            phase,
            "F",
            fvc::interpolate(pPrime(), pPrime().name())
           *fvc::snGrad(phase)*this->mesh_.magSf(),
            Fs
        );
    }

    // Add the turbulent dispersion force
    forAllConstIter
    (
        turbulentDispersionModelTable,
        turbulentDispersionModels_,
        turbulentDispersionModelIter
    )
    {
        const phaseInterface& interface =
            turbulentDispersionModelIter()->interface();

        const surfaceScalarField Df
        (
            fvc::interpolate(turbulentDispersionModelIter()->D())
        );

        const volScalarField alpha12(interface.phase1() + interface.phase2());
        const surfaceScalarField snGradAlpha1By12
        (
            fvc::snGrad
            (
                interface.phase1()
               /max(alpha12, interface.phase1().residualAlpha())
            )*this->mesh_.magSf()
        );
        const surfaceScalarField snGradAlpha2By12
        (
            fvc::snGrad
            (
                interface.phase2()
               /max(alpha12, interface.phase2().residualAlpha())
            )*this->mesh_.magSf()
        );

        addField(interface.phase1(), "F", Df*snGradAlpha1By12, Fs);
        addField(interface.phase2(), "F", Df*snGradAlpha2By12, Fs);
    }

    if (this->fillFields_)
    {
        this->fillFields("F", dimArea*dimDensity*dimAcceleration, Fs);
    }

    return Fs;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::surfaceScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::Ffs() const
{
    PtrList<surfaceScalarField> Ffs(this->phaseModels_.size());

    // Add the explicit part of the virtual mass force
    forAllConstIter(VmTable, Vms_, VmIter)
    {
        const volScalarField& Vm(*VmIter());
        const phaseInterface interface(*this, VmIter.key());

        forAllConstIter(phaseInterface, interface, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            const volScalarField VmPhase
            (
                (otherPhase/max(otherPhase, otherPhase.residualAlpha()))
               *Vm
            );

            addField
            (
                phase,
                "Ff",
               -fvc::interpolate(VmPhase)
               *(
                   byDt
                   (
                       fvc::absolute
                       (
                           this->MRF().absolute(iter().phi()().oldTime()),
                           iter().U()
                       )
                   )
                 + otherPhase.DUDtf()
                ),
                Ffs
            );
        }
    }

    // Add the lift force
    forAllConstIter
    (
        liftModelTable,
        liftModels_,
        liftModelIter
    )
    {
        const phaseInterface& interface = liftModelIter()->interface();

        const surfaceScalarField Ff(liftModelIter()->Ff());

        addField
        (
            interface.phase1(),
            "Ff",
            Ff,
            Ffs
        );
        addField
        (
            interface.phase2(),
            "Ff",
           -Ff,
            Ffs
        );
    }

    // Add the wall lubrication force
    forAllConstIter
    (
        wallLubricationModelTable,
        wallLubricationModels_,
        wallLubricationModelIter
    )
    {
        const phaseInterface& interface =
            wallLubricationModelIter()->interface();

        const surfaceScalarField Ff(wallLubricationModelIter()->Ff());

        addField
        (
            interface.phase1(),
            "Ff",
            Ff,
            Ffs
        );
        addField
        (
            interface.phase2(),
            "Ff",
           -Ff,
            Ffs
        );
    }

    // Add the phase pressure
    forAll(this->phaseModels_, phasei)
    {
        const phaseModel& phase = this->phaseModels_[phasei];

        addField
        (
            phase,
            "Ff",
            fvc::interpolate(phase.pPrime())
           *fvc::snGrad(phase)*this->mesh_.magSf(),
            Ffs
        );
    }

    // Add the turbulent dispersion force
    forAllConstIter
    (
        turbulentDispersionModelTable,
        turbulentDispersionModels_,
        turbulentDispersionModelIter
    )
    {
        const phaseInterface& interface =
            turbulentDispersionModelIter()->interface();

        const surfaceScalarField Df
        (
            fvc::interpolate(turbulentDispersionModelIter()->D())
        );

        const volScalarField alpha12(interface.phase1() + interface.phase2());
        const surfaceScalarField snGradAlpha1By12
        (
            fvc::snGrad
            (
                interface.phase1()
               /max(alpha12, interface.phase1().residualAlpha())
            )*this->mesh_.magSf()
        );
        const surfaceScalarField snGradAlpha2By12
        (
            fvc::snGrad
            (
                interface.phase2()
               /max(alpha12, interface.phase2().residualAlpha())
            )*this->mesh_.magSf()
        );

        addField(interface.phase1(), "F", Df*snGradAlpha1By12, Ffs);
        addField(interface.phase2(), "F", Df*snGradAlpha2By12, Ffs);
    }

    if (this->fillFields_)
    {
        this->fillFields("Ff", dimArea*dimDensity*dimAcceleration, Ffs);
    }

    return Ffs;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::surfaceScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::KdPhis() const
{
    PtrList<surfaceScalarField> KdPhis(this->phaseModels_.size());

    // Add the explicit part of the drag force
    forAllConstIter(KdTable, Kds_, KdIter)
    {
        const volScalarField& K(*KdIter());
        const phaseInterface interface(*this, KdIter.key());

        forAllConstIter(phaseInterface, interface, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            addField
            (
                phase,
                "KdPhi",
                fvc::interpolate
                (
                  -(otherPhase/max(otherPhase, otherPhase.residualAlpha()))
                  *K
                )
               *fvc::absolute
                (
                    this->MRF().absolute(otherPhase.phi()),
                    otherPhase.U()
                ),
                KdPhis
            );
        }
    }

    if (this->fillFields_)
    {
        this->fillFields
        (
            "KdPhi",
            dimArea*dimDensity*dimAcceleration,
            KdPhis
        );
    }

    return KdPhis;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::surfaceScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::KdPhifs() const
{
    PtrList<surfaceScalarField> KdPhifs(this->phaseModels_.size());

    // Add the explicit part of the drag force
    forAllConstIter(KdfTable, Kdfs_, KdfIter)
    {
        const surfaceScalarField& Kf(*KdfIter());
        const phaseInterface interface(*this, KdfIter.key());

        forAllConstIter(phaseInterface, interface, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            const surfaceScalarField alphaf
            (
                fvc::interpolate(otherPhase)
            );

            addField
            (
                phase,
                "KdPhif",
               -(alphaf/max(alphaf, otherPhase.residualAlpha()))
               *Kf
               *fvc::absolute
                (
                    this->MRF().absolute(otherPhase.phi()),
                    otherPhase.U()
                ),
                KdPhifs
            );
        }
    }

    if (this->fillFields_)
    {
        this->fillFields
        (
            "KdPhif",
            dimArea*dimDensity*dimAcceleration,
            KdPhifs
        );
    }

    return KdPhifs;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::Kds() const
{
    PtrList<volScalarField> Kds(this->phaseModels_.size());

    forAllConstIter(KdTable, Kds_, KdIter)
    {
        const volScalarField& K(*KdIter());
        const phaseInterface interface(*this, KdIter.key());

        forAllConstIter(phaseInterface, interface, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            addField
            (
                phase,
                "Kd",
                (otherPhase/max(otherPhase, otherPhase.residualAlpha()))
               *K,
                Kds
            );
        }
    }

    if (this->fillFields_)
    {
        this->fillFields("Kd", dimDensity/dimTime, Kds);
    }

    return Kds;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::volVectorField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::KdUs() const
{
    PtrList<volVectorField> KdUs(this->phaseModels_.size());

    // Add the explicit part of the drag force
    forAllConstIter(KdTable, Kds_, KdIter)
    {
        const volScalarField& K(*KdIter());
        const phaseInterface interface(*this, KdIter.key());

        forAllConstIter(phaseInterface, interface, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            addField
            (
                phase,
                "KdU",
               -(otherPhase/max(otherPhase, otherPhase.residualAlpha()))
               *K*otherPhase.U(),
                KdUs
            );
        }
    }

    if (this->fillFields_)
    {
        this->fillFields("KdU", dimDensity*dimAcceleration, KdUs);
    }

    return KdUs;
}


template<class BasePhaseSystem>
bool Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::
implicitPhasePressure(const phaseModel& phase) const
{
    return
        this->mesh_.solution().solverDict(phase.volScalarField::name()).
        template lookupOrDefault<Switch>
        (
            "implicitPhasePressure",
            false
        );
}


template<class BasePhaseSystem>
bool Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::
implicitPhasePressure() const
{
    bool implicitPressure = false;

    forAll(this->phaseModels_, phasei)
    {
        const phaseModel& phase = this->phaseModels_[phasei];

        implicitPressure = implicitPressure || implicitPhasePressure(phase);
    }

    return implicitPressure;
}


template<class BasePhaseSystem>
Foam::tmp<Foam::surfaceScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::alphaDByAf
(
    const PtrList<volScalarField>& rAUs,
    const PtrList<surfaceScalarField>& rAUfs
) const
{
    tmp<surfaceScalarField> alphaDByAf;

    // Add the phase pressure
    forAll(this->phaseModels_, phasei)
    {
        const phaseModel& phase = this->phaseModels_[phasei];

        const tmp<volScalarField> pPrime(phase.pPrime());

        addTmpField
        (
            alphaDByAf,
            fvc::interpolate(max(phase, scalar(0)))
           *(
                rAUfs.size()

                // Face-momentum form
              ? rAUfs[phasei]*fvc::interpolate(pPrime())

                // Cell-momentum form
              : fvc::interpolate(rAUs[phasei]*pPrime(), pPrime().name())
            )
        );
    }

    // Add the turbulent dispersion
    forAllConstIter
    (
        turbulentDispersionModelTable,
        turbulentDispersionModels_,
        turbulentDispersionModelIter
    )
    {
        const phaseInterface& interface =
            turbulentDispersionModelIter()->interface();

        const surfaceScalarField alpha1f
        (
            fvc::interpolate(max(interface.phase1(), scalar(0)))
        );

        const surfaceScalarField alpha2f
        (
            fvc::interpolate(max(interface.phase2(), scalar(0)))
        );

        addTmpField
        (
            alphaDByAf,
            alpha1f*alpha2f
           /max(alpha1f + alpha2f, interface.phase1().residualAlpha())
           *(
                rAUfs.size()

                // Face-momentum form
              ? max
                (
                    rAUfs[interface.phase1().index()],
                    rAUfs[interface.phase2().index()]
                )
               *fvc::interpolate(turbulentDispersionModelIter()->D())

                // Cell-momentum form
              : fvc::interpolate
                (
                    max
                    (
                        rAUs[interface.phase1().index()],
                        rAUs[interface.phase2().index()]
                    )
                   *turbulentDispersionModelIter()->D()
                )
            )
        );
    }

    return alphaDByAf;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::surfaceScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::ddtCorrs() const
{
    PtrList<surfaceScalarField> ddtCorrs(this->phaseModels_.size());

    // Add correction
    forAll(this->movingPhases(), movingPhasei)
    {
        const phaseModel& phase = this->movingPhases()[movingPhasei];

        addField
        (
            phase,
            "ddtCorr",
            fvc::ddtCorr
            (
                phase,
                phase.rho(),
                phase.U()(),
                phase.phi()(),
                phase.Uf()
            ),
            ddtCorrs
        );
    }


    const pimpleNoLoopControl& pimple = this->pimple();
    const Switch VmDdtCorr
    (
        pimple.dict().lookupOrDefault<Switch>("VmDdtCorrection", false)
    );

    // Optionally ddd virtual mass correction
    if (VmDdtCorr)
    {
        PtrList<volScalarField> VmDdtCoeffs(this->phaseModels_.size());
        PtrList<surfaceScalarField> VmDdtCorrs(this->phaseModels_.size());

        forAll(this->movingPhases(), movingPhasei)
        {
            const phaseModel& phase = this->movingPhases()[movingPhasei];
            const label phasei = phase.index();

            VmDdtCoeffs.set
            (
                phasei,
                fvm::ddt(phase.U()())().A()
            );

            VmDdtCorrs.set
            (
                phasei,
                fvc::ddtCorr
                (
                    phase.U()(),
                    phase.phi()(),
                    phase.Uf()
                )
            );
        }

        forAllConstIter(VmTable, Vms_, VmIter)
        {
            const volScalarField& Vm(*VmIter());
            const phaseInterface interface(*this, VmIter.key());

            forAllConstIter(phaseInterface, interface, iter)
            {
                const phaseModel& phase = iter();
                const label phasei = phase.index();
                const phaseModel& otherPhase = iter.otherPhase();
                const label otherPhasei = otherPhase.index();

                const volScalarField VmPhase
                (
                    (otherPhase/max(otherPhase, otherPhase.residualAlpha()))
                   *Vm
                );

                addField
                (
                    phase,
                    "ddtCorr",
                    fvc::interpolate(VmPhase)
                   *(
                        VmDdtCorrs[phasei]
                      + (
                            fvc::interpolate(VmDdtCoeffs[otherPhasei])
                           *(
                               otherPhase.Uf().valid()
                             ? (this->mesh_.Sf() & otherPhase.Uf()())()
                             : otherPhase.phi()()
                            )
                          - fvc::flux(VmDdtCoeffs[otherPhasei]*otherPhase.U())
                        )
                      - VmDdtCorrs[otherPhase.index()]
                    ),
                    ddtCorrs
                );
            }
        }
    }

    return ddtCorrs;
}


template<class BasePhaseSystem>
void Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::dragCorrs
(
    PtrList<volVectorField>& dragCorrs,
    PtrList<surfaceScalarField>& dragCorrfs
) const
{
    const phaseSystem::phaseModelList& phases = this->phaseModels_;

    PtrList<volVectorField> Uphis(phases.size());

    forAll(phases, i)
    {
        if (!phases[i].stationary())
        {
            Uphis.set
            (
                i,
                fvc::reconstruct(phases[i].phi())
            );
        }
    }

    forAllConstIter(KdTable, Kds_, KdIter)
    {
        const volScalarField& K(*KdIter());
        const phaseInterface interface(*this, KdIter.key());

        const phaseModel& phase1 = interface.phase1();
        const phaseModel& phase2 = interface.phase2();

        const label phase1i = phase1.index();
        const label phase2i = phase2.index();

        if (!phase1.stationary())
        {
            const volScalarField K1
            (
                (phase2/max(phase2, phase2.residualAlpha()))*K
            );

            addField
            (
                phase1,
                "dragCorr",
                K1
               *(
                    phase2.stationary()
                 ? -Uphis[phase1i]
                 :  (Uphis[phase2i] - Uphis[phase1i])
                ),
                dragCorrs
            );

            addField
            (
                phase1,
                "dragCorrf",
                fvc::interpolate(K1)
               *(
                    phase2.stationary()
                 ? -phase1.phi()
                 :  (phase2.phi() - phase1.phi())
                ),
                dragCorrfs
            );
        }

        if (!phase2.stationary())
        {
            const volScalarField K2
            (
                (phase1/max(phase1, phase1.residualAlpha()))*K
            );

            addField
            (
                phase2,
                "dragCorr",
                K2
               *(
                    phase1.stationary()
                 ? -Uphis[phase2i]
                 : (Uphis[phase1i] - Uphis[phase2i])
                ),
                dragCorrs
            );

            addField
            (
                phase2,
                "dragCorrf",
                fvc::interpolate(K2)
               *(
                    phase1.stationary()
                 ? -phase2.phi()
                 :  (phase1.phi() - phase2.phi())
                ),
                dragCorrfs
            );
        }
    }
}


template<class BasePhaseSystem>
void Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::partialElimination
(
    const PtrList<volScalarField>& rAUs,
    const PtrList<volVectorField>& KdUs,
    const PtrList<surfaceScalarField>& alphafs,
    const PtrList<surfaceScalarField>& rAUfs,
    const PtrList<surfaceScalarField>& KdPhis
)
{
    Info<< "Inverting drag systems: ";

    phaseSystem::phaseModelList& phases = this->phaseModels_;

    // Calculate the mean velocity from the current velocity
    // of the moving phases
    volVectorField Um(this->movingPhases()[0]*this->movingPhases()[0].U());

    for
    (
        label movingPhasei=1;
        movingPhasei<this->movingPhases().size();
        movingPhasei++
    )
    {
        Um +=
            this->movingPhases()[movingPhasei]
           *this->movingPhases()[movingPhasei].U();
    }

    // Remove the drag contributions from the velocity and flux of the phases
    // in preparation for the partial elimination of these terms
    forAll(phases, i)
    {
        if (!phases[i].stationary())
        {
            phases[i].URef() += rAUs[i]*KdUs[i];
            phases[i].phiRef() += rAUfs[i]*KdPhis[i];
        }
    }

    {
        // Create drag coefficient matrices
        PtrList<PtrList<volScalarField>> KdByAs(phases.size());
        PtrList<PtrList<surfaceScalarField>> KdByAfs(phases.size());

        forAll(phases, i)
        {
            KdByAs.set
            (
                i,
                new PtrList<volScalarField>(phases.size())
            );

            KdByAfs.set
            (
                i,
                new PtrList<surfaceScalarField>(phases.size())
            );
        }

        forAllConstIter(KdTable, Kds_, KdIter)
        {
            const volScalarField& K(*KdIter());
            const phaseInterface interface(*this, KdIter.key());

            const label phase1i = interface.phase1().index();
            const label phase2i = interface.phase2().index();

            const volScalarField K1
            (
                interface.phase2()
               /max(interface.phase2(), interface.phase2().residualAlpha())
               *K
            );

            const volScalarField K2
            (
                interface.phase1()
               /max(interface.phase1(), interface.phase1().residualAlpha())
               *K
            );

            addField
            (
                interface.phase2(),
                "KdByA",
               -rAUs[phase1i]*K1,
                KdByAs[phase1i]
            );
            addField
            (
                interface.phase1(),
                "KdByA",
               -rAUs[phase2i]*K2,
                KdByAs[phase2i]
            );

            addField
            (
                interface.phase2(),
                "KdByAf",
               -rAUfs[phase1i]*fvc::interpolate(K1),
                KdByAfs[phase1i]
            );
            addField
            (
                interface.phase1(),
                "KdByAf",
               -rAUfs[phase2i]*fvc::interpolate(K2),
                KdByAfs[phase2i]
            );
        }

        forAll(phases, i)
        {
            this->fillFields("KdByAs", dimless, KdByAs[i]);
            this->fillFields("KdByAfs", dimless, KdByAfs[i]);

            KdByAs[i][i] = 1;
            KdByAfs[i][i] = 1;
        }

        // Decompose
        for (label i = 0; i < phases.size(); i++)
        {
            for (label j = i + 1; j < phases.size(); j++)
            {
                KdByAs[i][j] /= KdByAs[i][i];
                KdByAfs[i][j] /= KdByAfs[i][i];
                for (label k = i + 1; k < phases.size(); ++ k)
                {
                    KdByAs[j][k] -= KdByAs[j][i]*KdByAs[i][k];
                    KdByAfs[j][k] -= KdByAfs[j][i]*KdByAfs[i][k];
                }
            }
        }

        {
            volScalarField detKdByAs(KdByAs[0][0]);
            surfaceScalarField detPhiKdfs(KdByAfs[0][0]);

            for (label i = 1; i < phases.size(); i++)
            {
                detKdByAs *= KdByAs[i][i];
                detPhiKdfs *= KdByAfs[i][i];
            }

            Info<< "Min cell/face det = " << gMin(detKdByAs.primitiveField())
                << "/" << gMin(detPhiKdfs.primitiveField()) << endl;
        }

        // Solve for the velocities and fluxes
        for (label i = 1; i < phases.size(); i++)
        {
            if (!phases[i].stationary())
            {
                for (label j = 0; j < i; j ++)
                {
                    phases[i].URef() -= KdByAs[i][j]*phases[j].U();
                    phases[i].phiRef() -= KdByAfs[i][j]*phases[j].phi();
                }
            }
        }
        for (label i = phases.size() - 1; i >= 0; i--)
        {
            if (!phases[i].stationary())
            {
                for (label j = phases.size() - 1; j > i; j--)
                {
                    phases[i].URef() -= KdByAs[i][j]*phases[j].U();
                    phases[i].phiRef() -= KdByAfs[i][j]*phases[j].phi();
                }
                phases[i].URef() /= KdByAs[i][i];
                phases[i].phiRef() /= KdByAfs[i][i];
            }
        }
    }

    this->setMixtureU(Um);
    this->setMixturePhi(alphafs, this->phi());
}


template<class BasePhaseSystem>
void Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::partialEliminationf
(
    const PtrList<surfaceScalarField>& rAUfs,
    const PtrList<surfaceScalarField>& alphafs,
    const PtrList<surfaceScalarField>& KdPhifs
)
{
    Info<< "Inverting drag system: ";

    phaseSystem::phaseModelList& phases = this->phaseModels_;

    // Remove the drag contributions from the flux of the phases
    // in preparation for the partial elimination of these terms
    forAll(phases, i)
    {
        if (!phases[i].stationary())
        {
            phases[i].phiRef() += rAUfs[i]*KdPhifs[i];
        }
    }

    {
        // Create drag coefficient matrix
        PtrList<PtrList<surfaceScalarField>> phiKdfs(phases.size());

        forAll(phases, phasei)
        {
            phiKdfs.set
            (
                phasei,
                new PtrList<surfaceScalarField>(phases.size())
            );
        }

        forAllConstIter(KdfTable, Kdfs_, KdfIter)
        {
            const surfaceScalarField& Kf(*KdfIter());
            const phaseInterface interface(*this, KdfIter.key());

            const label phase1i = interface.phase1().index();
            const label phase2i = interface.phase2().index();

            const surfaceScalarField alpha1f
            (
                fvc::interpolate(interface.phase1())
            );

            const surfaceScalarField alpha2f
            (
                fvc::interpolate(interface.phase2())
            );

            addField
            (
                interface.phase2(),
                "phiKdf",
               -rAUfs[phase1i]
               *alpha2f/max(alpha2f, interface.phase2().residualAlpha())
               *Kf,
                phiKdfs[phase1i]
            );
            addField
            (
                interface.phase1(),
                "phiKdf",
               -rAUfs[phase2i]
               *alpha1f/max(alpha1f, interface.phase1().residualAlpha())
               *Kf,
                phiKdfs[phase2i]
            );
        }

        forAll(phases, phasei)
        {
            this->fillFields("phiKdf", dimless, phiKdfs[phasei]);

            phiKdfs[phasei][phasei] = 1;
        }

        // Decompose
        for (label i = 0; i < phases.size(); i++)
        {
            for (label j = i + 1; j < phases.size(); j++)
            {
                phiKdfs[i][j] /= phiKdfs[i][i];
                for (label k = i + 1; k < phases.size(); ++ k)
                {
                    phiKdfs[j][k] -= phiKdfs[j][i]*phiKdfs[i][k];
                }
            }
        }

        {
            surfaceScalarField detPhiKdfs(phiKdfs[0][0]);

            for (label i = 1; i < phases.size(); i++)
            {
                detPhiKdfs *= phiKdfs[i][i];
            }

            Info<< "Min face det = "
                << gMin(detPhiKdfs.primitiveField()) << endl;
        }

        // Solve for the fluxes
        for (label i = 1; i < phases.size(); i++)
        {
            if (!phases[i].stationary())
            {
                for (label j = 0; j < i; j ++)
                {
                    phases[i].phiRef() -= phiKdfs[i][j]*phases[j].phi();
                }
            }
        }
        for (label i = phases.size() - 1; i >= 0; i--)
        {
            if (!phases[i].stationary())
            {
                for (label j = phases.size() - 1; j > i; j--)
                {
                    phases[i].phiRef() -= phiKdfs[i][j]*phases[j].phi();
                }
                phases[i].phiRef() /= phiKdfs[i][i];
            }
        }
    }

    this->setMixturePhi(alphafs, this->phi());
}


template<class BasePhaseSystem>
bool Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::read()
{
    if (BasePhaseSystem::read())
    {
        bool readOK = true;

        // Read models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
