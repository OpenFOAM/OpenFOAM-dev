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

#include "scalarMatrices.H"

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
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::Vmfs() const
{
    PtrList<surfaceScalarField> Vmfs(this->phaseModels_.size());

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

            addField(phase, "Vmf", byDt(fvc::interpolate(VmPhase)), Vmfs);
        }
    }

    return Vmfs;
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

    return Ffs;
}


template<class BasePhaseSystem>
void Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::invADs
(
    List<UPtrList<scalarField>>& ADs
) const
{
    const label n = ADs.size();

    scalarSquareMatrix AD(n);
    scalarField source(n);
    labelList pivotIndices(n);

    forAll(ADs[0][0], ci)
    {
        for (label i=0; i<n; i++)
        {
            for (label j=0; j<n; j++)
            {
                AD(i, j) = ADs[i][j][ci];
            }
        }

        // Calculate the inverse of AD using LD decomposition
        // and back-substitution
        LUDecompose(AD, pivotIndices);

        for (label j=0; j<n; j++)
        {
            source = Zero;
            source[j] = 1;

            LUBacksubstitute(AD, pivotIndices, source);

            for (label i=0; i<n; i++)
            {
                ADs[i][j][ci] = source[i];
            }
        }
    }
}


template<class BasePhaseSystem>
template<template<class> class PatchField, class GeoMesh>
void Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::invADs
(
    PtrList<PtrList<GeometricField<scalar, PatchField, GeoMesh>>>& ADs
) const
{
    const label n = ADs.size();

    List<UPtrList<scalarField>> ADps(n);

    for (label i=0; i<n; i++)
    {
        ADps[i].setSize(n);

        for (label j=0; j<n; j++)
        {
            ADps[i].set(j, &ADs[i][j]);

            ADs[i][j].dimensions().reset(dimless/ADs[i][j].dimensions());
        }
    }

    invADs(ADps);

    forAll(ADs[0][0].boundaryField(), patchi)
    {
        for (label i=0; i<n; i++)
        {
            for (label j=0; j<n; j++)
            {
                ADps[i].set(j, &ADs[i][j].boundaryFieldRef()[patchi]);
            }
        }

        invADs(ADps);
    }
}


template<class BasePhaseSystem>
void Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::invADs
(
    const PtrList<volScalarField>& As,
    PtrList<PtrList<volScalarField>>& invADs,
    PtrList<PtrList<surfaceScalarField>>& invADfs
) const
{
    const label n = As.size();

    invADs.setSize(n);
    invADfs.setSize(n);

    forAll(invADs, i)
    {
        invADs.set(i, new PtrList<volScalarField>(n));
        invADs[i].set(i, As[i].clone());

        invADfs.set(i, new PtrList<surfaceScalarField>(n));
        invADfs[i].set(i, fvc::interpolate(As[i]));
    }

    labelList movingPhases(this->phases().size(), -1);

    forAll(this->movingPhases(), movingPhasei)
    {
        movingPhases[this->movingPhases()[movingPhasei].index()] = movingPhasei;
    }

    forAllConstIter(KdTable, Kds_, KdIter)
    {
        const volScalarField& K(*KdIter());
        const phaseInterface interface(*this, KdIter.key());

        forAllConstIter(phaseInterface, interface, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            const label i = movingPhases[phase.index()];

            if (i != -1)
            {
                const volScalarField Kij
                (
                    (otherPhase/max(otherPhase, otherPhase.residualAlpha()))*K
                );

                const surfaceScalarField Kijf(fvc::interpolate(Kij));

                invADs[i][i] += Kij;
                invADfs[i][i] += Kijf;

                const label j = movingPhases[otherPhase.index()];

                if (j != -1)
                {
                    invADs[i].set(j, -Kij);
                    invADfs[i].set(j, -Kijf);
                }
            }
        }
    }

    for (label i=0; i<n; i++)
    {
        for (label j=0; j<n; j++)
        {
            if (!invADs[i].set(j))
            {
                invADs[i].set
                (
                    j,
                    volScalarField::New
                    (
                        "0",
                        this->mesh(),
                        dimensionedScalar(As[0].dimensions(), 0)
                    )
                );

                invADfs[i].set
                (
                    j,
                    surfaceScalarField::New
                    (
                        "0",
                        this->mesh(),
                        dimensionedScalar(As[0].dimensions(), 0)
                    )
                );
            }
        }
    }

    MomentumTransferPhaseSystem<BasePhaseSystem>::invADs(invADs);
    MomentumTransferPhaseSystem<BasePhaseSystem>::invADs(invADfs);
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::PtrList<Foam::surfaceScalarField>>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::invADfs
(
    const PtrList<surfaceScalarField>& As
) const
{
    const label n = As.size();

    PtrList<PtrList<surfaceScalarField>> invADfs(n);

    forAll(invADfs, i)
    {
        invADfs.set(i, new PtrList<surfaceScalarField>(n));
        invADfs[i].set(i, As[i].clone());
    }

    labelList movingPhases(this->phases().size(), -1);

    forAll(this->movingPhases(), movingPhasei)
    {
        movingPhases[this->movingPhases()[movingPhasei].index()] = movingPhasei;
    }

    forAllConstIter(KdfTable, Kdfs_, KdfIter)
    {
        const surfaceScalarField& Kf(*KdfIter());
        const phaseInterface interface(*this, KdfIter.key());

        forAllConstIter(phaseInterface, interface, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            const label i = movingPhases[phase.index()];

            if (i != -1)
            {
                const surfaceScalarField alphaf
                (
                    fvc::interpolate(otherPhase)
                );

                const surfaceScalarField Kfij
                (
                    (alphaf/max(alphaf, otherPhase.residualAlpha()))*Kf
                );

                invADfs[i][i] += Kfij;

                const label j = movingPhases[otherPhase.index()];

                if (j != -1)
                {
                    invADfs[i].set(j, -Kfij);
                }
            }
        }
    }

    for (label i=0; i<n; i++)
    {
        for (label j=0; j<n; j++)
        {
            if (!invADfs[i].set(j))
            {
                invADfs[i].set
                (
                    j,
                    surfaceScalarField::New
                    (
                        "0",
                        this->mesh(),
                        dimensionedScalar(As[0].dimensions(), 0)
                    )
                );
            }
        }
    }

    invADs(invADfs);

    return invADfs;
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
    const PtrList<volScalarField>& rAs
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
           *fvc::interpolate(rAs[phasei]*pPrime(), pPrime().name())
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
           *fvc::interpolate
            (
                max
                (
                    rAs[interface.phase1().index()],
                    rAs[interface.phase2().index()]
                )
               *turbulentDispersionModelIter()->D()
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
    labelList movingPhases(this->phases().size(), -1);
    PtrList<volVectorField> Uphis(this->movingPhases().size());

    forAll(this->movingPhases(), movingPhasei)
    {
        movingPhases[this->movingPhases()[movingPhasei].index()] = movingPhasei;

        Uphis.set
        (
            movingPhasei,
            fvc::reconstruct(this->movingPhases()[movingPhasei].phi())
        );
    }

    forAllConstIter(KdTable, Kds_, KdIter)
    {
        const volScalarField& K(*KdIter());
        const phaseInterface interface(*this, KdIter.key());

        forAllConstIter(phaseInterface, interface, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            const label i = movingPhases[phase.index()];
            const label j = movingPhases[otherPhase.index()];

            if (i != -1)
            {
                const volScalarField K1
                (
                    (otherPhase/max(otherPhase, otherPhase.residualAlpha()))*K
                );

                addField
                (
                    i,
                    "dragCorr",
                    K1*(j == -1 ? -Uphis[i] : (Uphis[j] - Uphis[i])),
                    dragCorrs
                );

                addField
                (
                    i,
                    "dragCorrf",
                    fvc::interpolate(K1)
                   *(j == -1 ? -phase.phi() : (otherPhase.phi() - phase.phi())),
                    dragCorrfs
                );
            }
        }
    }
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
