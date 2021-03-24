/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
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

#include "BlendedInterfacialModel.H"
#include "dragModel.H"
#include "virtualMassModel.H"
#include "liftModel.H"
#include "wallLubricationModel.H"
#include "turbulentDispersionModel.H"

#include "HashPtrTable.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvcAverage.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcFlux.H"
#include "fvcSnGrad.H"
#include "fvcMeshPhi.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::Kd
(
    const phasePairKey& key
) const
{
    if (dragModels_.found(key))
    {
        return dragModels_[key]->K();
    }
    else
    {
        return volScalarField::New
        (
            dragModel::typeName + ":K",
            this->mesh_,
            dimensionedScalar(dragModel::dimK, 0)
        );
    }
}


template<class BasePhaseSystem>
Foam::tmp<Foam::surfaceScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::Kdf
(
    const phasePairKey& key
) const
{
    if (dragModels_.found(key))
    {
        return dragModels_[key]->Kf();
    }
    else
    {
        return surfaceScalarField::New
        (
            dragModel::typeName + ":K",
            this->mesh_,
            dimensionedScalar(dragModel::dimK, 0)
        );
    }
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::Vm
(
    const phasePairKey& key
) const
{
    if (virtualMassModels_.found(key))
    {
        return virtualMassModels_[key]->K();
    }
    else
    {
        return volScalarField::New
        (
            virtualMassModel::typeName + ":K",
            this->mesh_,
            dimensionedScalar(virtualMassModel::dimK, 0)
        );
    }
}


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
        const phasePairKey& key = dmdtfIter.key();
        const phasePair& pair(this->phasePairs_[key]);

        const volScalarField dmdtf(Pair<word>::compare(pair, key)**dmdtfIter());

        phaseModel& phase1 = this->phases()[pair.phase1().name()];
        phaseModel& phase2 = this->phases()[pair.phase2().name()];

        const volScalarField dmdtf21(posPart(dmdtf));
        const volScalarField dmdtf12(negPart(dmdtf));

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
    this->generatePairsAndSubModels
    (
        "drag",
        dragModels_
    );

    this->generatePairsAndSubModels
    (
        "virtualMass",
        virtualMassModels_
    );

    this->generatePairsAndSubModels
    (
        "lift",
        liftModels_
    );

    this->generatePairsAndSubModels
    (
        "wallLubrication",
        wallLubricationModels_
    );

    this->generatePairsAndSubModels
    (
        "turbulentDispersion",
        turbulentDispersionModels_
    );

    forAllConstIter
    (
        dragModelTable,
        dragModels_,
        dragModelIter
    )
    {
        const phasePair& pair(this->phasePairs_[dragModelIter.key()]);

        Kds_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("Kd", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh()
                ),
                this->mesh(),
                dimensionedScalar(dragModel::dimK, 0)
            )
        );

        Kdfs_.insert
        (
            pair,
            new surfaceScalarField
            (
                IOobject
                (
                    IOobject::groupName("Kdf", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh()
                ),
                this->mesh(),
                dimensionedScalar(dragModel::dimK, 0)
            )
        );
    }

    forAllConstIter
    (
        virtualMassModelTable,
        virtualMassModels_,
        virtualMassModelIter
    )
    {
        const phasePair& pair(this->phasePairs_[virtualMassModelIter.key()]);

        Vms_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("Vm", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh()
                ),
                this->mesh(),
                dimensionedScalar(virtualMassModel::dimK, 0)
            )
        );

        Vmfs_.insert
        (
            pair,
            new surfaceScalarField
            (
                IOobject
                (
                    IOobject::groupName("Vmf", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh()
                ),
                this->mesh(),
                dimensionedScalar(virtualMassModel::dimK, 0)
            )
        );
    }
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

    // Update the drag coefficients
    forAllConstIter
    (
        dragModelTable,
        dragModels_,
        dragModelIter
    )
    {
        *Kds_[dragModelIter.key()] = dragModelIter()->K();
        *Kdfs_[dragModelIter.key()] = dragModelIter()->Kf();
    }

    // Add the implicit part of the drag force
    forAllConstIter(KdTable, Kds_, KdIter)
    {
        const volScalarField& K(*KdIter());
        const phasePair& pair(this->phasePairs_[KdIter.key()]);

        forAllConstIter(phasePair, pair, iter)
        {
            if (!iter().stationary())
            {
                fvVectorMatrix& eqn = *eqns[iter().name()];

                eqn -= fvm::Sp(K, eqn.psi());
            }
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
        *Vmfs_[virtualMassModelIter.key()] = virtualMassModelIter()->Kf();
    }

    // Add the virtual mass force
    forAllConstIter(VmTable, Vms_, VmIter)
    {
        const volScalarField& Vm(*VmIter());
        const phasePair& pair(this->phasePairs_[VmIter.key()]);

        forAllConstIter(phasePair, pair, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            if (!phase.stationary())
            {
                fvVectorMatrix& eqn = *eqns[phase.name()];

                const volVectorField& U = eqn.psi();
                const surfaceScalarField& phi = phase.phi();
                const tmp<surfaceScalarField> taphi(fvc::absolute(phi, U));
                const surfaceScalarField& aphi(taphi());

                eqn -=
                    Vm
                   *(
                        fvm::ddt(U)
                      + fvm::div(aphi, U) - fvm::Sp(fvc::div(aphi), U)
                      - otherPhase.DUDt()
                    )
                  + this->MRF_.DDt(Vm, U - otherPhase.U());
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

    // Create U & grad(U) fields
    PtrList<fvVectorMatrix> UgradUs(this->phaseModels_.size());
    forAll(this->phaseModels_, phasei)
    {
        const phaseModel& phase = this->phaseModels_[phasei];

        if (!phase.stationary())
        {
            const volVectorField& U = phase.U();
            const surfaceScalarField& phi = phase.phi();
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
        const phasePair& pair(this->phasePairs_[VmIter.key()]);

        forAllConstIter(phasePair, pair, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            if (!phase.stationary())
            {
                *eqns[phase.name()] -=
                    Vm
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
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::AFfs() const
{
    PtrList<surfaceScalarField> AFfs(this->phaseModels_.size());

    // Add the implicit part of the drag force
    forAllConstIter(KdfTable, Kdfs_, KdfIter)
    {
        const surfaceScalarField& Kf(*KdfIter());
        const phasePair& pair(this->phasePairs_[KdfIter.key()]);

        forAllConstIter(phasePair, pair, iter)
        {
            addField(iter(), "AFf", Kf, AFfs);
        }
    }

    // Add the implicit part of the virtual mass force
    forAllConstIter(VmfTable, Vmfs_, VmfIter)
    {
        const surfaceScalarField& Vmf(*VmfIter());
        const phasePair& pair(this->phasePairs_[VmfIter.key()]);

        forAllConstIter(phasePair, pair, iter)
        {
            addField(iter(), "AFf", byDt(Vmf), AFfs);
        }
    }

    this->fillFields("AFf", dimDensity/dimTime, AFfs);

    return AFfs;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::surfaceScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::phiFs
(
    const PtrList<volScalarField>& rAUs
)
{
    PtrList<surfaceScalarField> phiFs(this->phaseModels_.size());

    // Add the lift force
    forAllConstIter
    (
        liftModelTable,
        liftModels_,
        liftModelIter
    )
    {
        const volVectorField F(liftModelIter()->F<vector>());
        const phasePair& pair(this->phasePairs_[liftModelIter.key()]);

        addField
        (
            pair.phase1(),
            "phiF",
            fvc::flux(rAUs[pair.phase1().index()]*F),
            phiFs
        );
        addField
        (
            pair.phase2(),
            "phiF",
           -fvc::flux(rAUs[pair.phase2().index()]*F),
            phiFs
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
        const volVectorField F(wallLubricationModelIter()->F<vector>());
        const phasePair&
            pair(this->phasePairs_[wallLubricationModelIter.key()]);

        addField
        (
            pair.phase1(),
            "phiF",
            fvc::flux(rAUs[pair.phase1().index()]*F),
            phiFs
        );
        addField
        (
            pair.phase2(),
            "phiF",
           -fvc::flux(rAUs[pair.phase2().index()]*F),
            phiFs
        );
    }

    // Add the phase pressure
    forAll(this->phaseModels_, phasei)
    {
        const phaseModel& phase = this->phaseModels_[phasei];

        const surfaceScalarField pPrimeByAf
        (
            fvc::interpolate(rAUs[phasei]*phase.pPrime())
        );

        const surfaceScalarField snGradAlpha1
        (
            fvc::snGrad(phase)*this->mesh_.magSf()
        );

        addField(phase, "phiF", pPrimeByAf*snGradAlpha1, phiFs);
    }

    // Add the turbulent dispersion force
    forAllConstIter
    (
        turbulentDispersionModelTable,
        turbulentDispersionModels_,
        turbulentDispersionModelIter
    )
    {
        const phasePair&
            pair(this->phasePairs_[turbulentDispersionModelIter.key()]);

        const volScalarField D(turbulentDispersionModelIter()->D());

        const surfaceScalarField DByA1f
        (
            fvc::interpolate(rAUs[pair.phase1().index()]*D)
        );
        const surfaceScalarField DByA2f
        (
            fvc::interpolate(rAUs[pair.phase2().index()]*D)
        );

        const volScalarField alpha12(pair.phase1() + pair.phase2());
        const surfaceScalarField snGradAlpha1By12
        (
            fvc::snGrad
            (
                pair.phase1()/max(alpha12, pair.phase1().residualAlpha())
            )*this->mesh_.magSf()
        );
        const surfaceScalarField snGradAlpha2By12
        (
            fvc::snGrad
            (
                pair.phase2()/max(alpha12, pair.phase2().residualAlpha())
            )*this->mesh_.magSf()
        );

        addField(pair.phase1(), "phiF", DByA1f*snGradAlpha1By12, phiFs);
        addField(pair.phase2(), "phiF", DByA2f*snGradAlpha2By12, phiFs);
    }

    if (this->fillFields_)
    {
        this->fillFields("phiF", dimForce/dimDensity/dimVelocity, phiFs);
    }

    return phiFs;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::surfaceScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::phiFfs
(
    const PtrList<surfaceScalarField>& rAUfs
)
{
    PtrList<surfaceScalarField> phiFfs(this->phaseModels_.size());

    // Add the explicit part of the virtual mass force
    forAllConstIter(VmfTable, Vmfs_, VmfIter)
    {
        const surfaceScalarField& Vmf(*VmfIter());
        const phasePair& pair(this->phasePairs_[VmfIter.key()]);

        forAllConstIter(phasePair, pair, iter)
        {
            addField
            (
                iter(),
                "phiFf",
               -rAUfs[iter().index()]*Vmf
               *(
                   byDt
                   (
                       fvc::absolute
                       (
                           this->MRF().absolute(iter().phi()().oldTime()),
                           iter().U()
                       )
                   )
                 + iter.otherPhase().DUDtf()
                ),
                phiFfs
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
        const surfaceScalarField Ff(liftModelIter()->Ff());
        const phasePair& pair(this->phasePairs_[liftModelIter.key()]);

        addField
        (
            pair.phase1(),
            "phiFf",
            rAUfs[pair.phase1().index()]*Ff,
            phiFfs
        );
        addField
        (
            pair.phase2(),
            "phiFf",
           -rAUfs[pair.phase2().index()]*Ff,
            phiFfs
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
        const surfaceScalarField Ff(wallLubricationModelIter()->Ff());
        const phasePair&
            pair(this->phasePairs_[wallLubricationModelIter.key()]);

        addField
        (
            pair.phase1(),
            "phiFf",
            rAUfs[pair.phase1().index()]*Ff,
            phiFfs
        );
        addField
        (
            pair.phase2(),
            "phiFf",
           -rAUfs[pair.phase2().index()]*Ff,
            phiFfs
        );
    }

    // Add the phase pressure
    forAll(this->phaseModels_, phasei)
    {
        const phaseModel& phase = this->phaseModels_[phasei];

        const surfaceScalarField pPrimeByAf
        (
            rAUfs[phasei]*fvc::interpolate(phase.pPrime())
        );

        const surfaceScalarField snGradAlpha1
        (
            fvc::snGrad(phase)*this->mesh_.magSf()
        );

        addField(phase, "phiFf", pPrimeByAf*snGradAlpha1, phiFfs);
    }

    // Add the turbulent dispersion force
    forAllConstIter
    (
        turbulentDispersionModelTable,
        turbulentDispersionModels_,
        turbulentDispersionModelIter
    )
    {
        const surfaceScalarField Ff(turbulentDispersionModelIter()->Ff());
        const phasePair&
            pair(this->phasePairs_[turbulentDispersionModelIter.key()]);

        addField
        (
            pair.phase1(),
            "phiFf",
            rAUfs[pair.phase1().index()]*Ff,
            phiFfs
        );
        addField
        (
            pair.phase2(),
            "phiFf",
           -rAUfs[pair.phase2().index()]*Ff,
            phiFfs
        );
    }

    if (this->fillFields_)
    {
        this->fillFields("phiFf", dimForce/dimDensity/dimVelocity, phiFfs);
    }

    return phiFfs;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::surfaceScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::phiKdPhis
(
    const PtrList<volScalarField>& rAUs
) const
{
    PtrList<surfaceScalarField> phiKdPhis(this->phaseModels_.size());

    // Add the explicit part of the drag force
    forAllConstIter(KdTable, Kds_, KdIter)
    {
        const volScalarField& K(*KdIter());
        const phasePair& pair(this->phasePairs_[KdIter.key()]);

        forAllConstIter(phasePair, pair, iter)
        {
            addField
            (
                iter(),
                "phiKdPhi",
               -fvc::interpolate(rAUs[iter().index()]*K)
               *fvc::absolute
                (
                    this->MRF().absolute(iter.otherPhase().phi()),
                    iter.otherPhase().U()
                ),
                phiKdPhis
            );
        }
    }

    if (this->fillFields_)
    {
        this->fillFields
        (
            "phiKdPhi",
            dimForce/dimDensity/dimVelocity,
            phiKdPhis
        );
    }

    return phiKdPhis;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::surfaceScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::phiKdPhifs
(
    const PtrList<surfaceScalarField>& rAUfs
) const
{
    PtrList<surfaceScalarField> phiKdPhifs(this->phaseModels_.size());

    // Add the explicit part of the drag force
    forAllConstIter(KdfTable, Kdfs_, KdfIter)
    {
        const surfaceScalarField& Kf(*KdfIter());
        const phasePair& pair(this->phasePairs_[KdfIter.key()]);

        forAllConstIter(phasePair, pair, iter)
        {
            addField
            (
                iter(),
                "phiKdPhif",
               -rAUfs[iter().index()]*Kf
               *fvc::absolute
                (
                    this->MRF().absolute(iter.otherPhase().phi()),
                    iter.otherPhase().U()
                ),
                phiKdPhifs
            );
        }
    }

    if (this->fillFields_)
    {
        this->fillFields
        (
            "phiKdPhif",
            dimForce/dimDensity/dimVelocity,
            phiKdPhifs
        );
    }

    return phiKdPhifs;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::volVectorField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::KdUByAs
(
    const PtrList<volScalarField>& rAUs
) const
{
    PtrList<volVectorField> KdUByAs(this->phaseModels_.size());

    // Add the explicit part of the drag force
    forAllConstIter(KdTable, Kds_, KdIter)
    {
        const volScalarField& K(*KdIter());
        const phasePair& pair(this->phasePairs_[KdIter.key()]);

        forAllConstIter(phasePair, pair, iter)
        {
            addField
            (
                iter(),
                "KdUByA",
               -rAUs[iter().index()]*K*iter.otherPhase().U(),
                KdUByAs
            );
        }
    }

    if (this->fillFields_)
    {
        this->fillFields("KdUByA", dimVelocity, KdUByAs);
    }

    return KdUByAs;
}


template<class BasePhaseSystem>
bool Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::
implicitPhasePressure(const phaseModel& phase) const
{
    return
        this->mesh_.solverDict(phase.volScalarField::name()).
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
Foam::PtrList<Foam::surfaceScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::DByAfs
(
    const PtrList<volScalarField>& rAUs,
    const PtrList<surfaceScalarField>& rAUfs
) const
{
    PtrList<surfaceScalarField> DByAfs(this->phaseModels_.size());

    if (rAUfs.size())
    {
        // Add the phase pressure
        forAll(this->phaseModels_, phasei)
        {
            const phaseModel& phase = this->phaseModels_[phasei];

            const surfaceScalarField pPrimeByAf
            (
                rAUfs[phasei]*fvc::interpolate(phase.pPrime())
            );

            addField(phase, "DByAf", pPrimeByAf, DByAfs);

            forAll(this->phaseModels_, phasej)
            {
                if (phasej != phasei)
                {
                    const phaseModel& phase2 = this->phaseModels_[phasej];

                    addField
                    (
                        phase2,
                        "DByAf",
                        fvc::interpolate
                        (
                            phase2
                           /max(1 - phase, phase2.residualAlpha())
                        )*pPrimeByAf,
                        DByAfs
                    );
                }
            }
        }

        // Add the turbulent dispersion
        forAllConstIter
        (
            turbulentDispersionModelTable,
            turbulentDispersionModels_,
            turbulentDispersionModelIter
        )
        {
            const phasePair&
                pair(this->phasePairs_[turbulentDispersionModelIter.key()]);

            const surfaceScalarField Df
            (
                fvc::interpolate(turbulentDispersionModelIter()->D())
            );

            const surfaceScalarField alpha12f
            (
                fvc::interpolate(pair.phase1() + pair.phase2())
            );

            addField
            (
                pair.phase1(),
                "DByAf",
                rAUfs[pair.phase1().index()]*Df
               /max(alpha12f, pair.phase1().residualAlpha()),
                DByAfs
            );
            addField
            (
                pair.phase2(),
                "DByAf",
                rAUfs[pair.phase2().index()]*Df
               /max(alpha12f, pair.phase2().residualAlpha()),
                DByAfs
            );
        }
    }
    else
    {
        // Add the phase pressure
        forAll(this->phaseModels_, phasei)
        {
            const phaseModel& phase = this->phaseModels_[phasei];

            const surfaceScalarField pPrimeByAf
            (
                fvc::interpolate(rAUs[phasei]*phase.pPrime())
            );

            addField(phase, "DByAf", pPrimeByAf, DByAfs);

            forAll(this->phaseModels_, phasej)
            {
                if (phasej != phasei)
                {
                    const phaseModel& phase2 = this->phaseModels_[phasej];

                    addField
                    (
                        phase2,
                        "DByAf",
                        fvc::interpolate
                        (
                            phase2
                           /max(1 - phase, phase2.residualAlpha())
                        )*pPrimeByAf,
                        DByAfs
                    );
                }
            }
        }

        // Add the turbulent dispersion
        forAllConstIter
        (
            turbulentDispersionModelTable,
            turbulentDispersionModels_,
            turbulentDispersionModelIter
        )
        {
            const phasePair&
                pair(this->phasePairs_[turbulentDispersionModelIter.key()]);

            const volScalarField D(turbulentDispersionModelIter()->D());

            const volScalarField alpha12(pair.phase1() + pair.phase2());

            addField
            (
                pair.phase1(),
                "DByAf",
                fvc::interpolate
                (
                    rAUs[pair.phase1().index()]*D
                   /max(alpha12, pair.phase1().residualAlpha())
                ),
                DByAfs
            );
            addField
            (
                pair.phase2(),
                "DByAf",
                fvc::interpolate
                (
                    rAUs[pair.phase2().index()]*D
                   /max(alpha12, pair.phase2().residualAlpha())
                ),
                DByAfs
            );
        }
    }

    return DByAfs;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::surfaceScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::ddtCorrByAs
(
    const PtrList<volScalarField>& rAUs,
    const bool includeVirtualMass
) const
{
    PtrList<surfaceScalarField> ddtCorrByAs(this->phaseModels_.size());

    // Construct phi differences
    PtrList<surfaceScalarField> phiCorrs(this->phaseModels_.size());

    forAll(this->movingPhases(), movingPhasei)
    {
        const phaseModel& phase = this->movingPhases()[movingPhasei];
        const label phasei = phase.index();

        phiCorrs.set
        (
            phasei,
            this->MRF().zeroFilter
            (
                (
                    phase.Uf().valid()
                  ? (this->mesh_.Sf() & phase.Uf()().oldTime())()
                  : phase.phi()().oldTime()
                )
              - fvc::flux(phase.U()().oldTime())
            )()
        );
    }

    // Add correction
    forAll(this->movingPhases(), movingPhasei)
    {
        const phaseModel& phase = this->movingPhases()[movingPhasei];
        const label phasei = phase.index();
        const volScalarField& alpha = phase;

        // Apply ddtPhiCorr filter in pure(ish) phases
        surfaceScalarField alphafBar
        (
            fvc::interpolate(fvc::average(fvc::interpolate(alpha)))
        );

        tmp<surfaceScalarField> phiCorrCoeff = pos0(alphafBar - 0.99);

        surfaceScalarField::Boundary& phiCorrCoeffBf =
            phiCorrCoeff.ref().boundaryFieldRef();

        forAll(this->mesh_.boundary(), patchi)
        {
            // Set ddtPhiCorr to 0 on non-coupled boundaries
            if
            (
                !this->mesh_.boundary()[patchi].coupled()
             || isA<cyclicAMIFvPatch>(this->mesh_.boundary()[patchi])
            )
            {
                phiCorrCoeffBf[patchi] = 0;
            }
        }

        addField
        (
            phase,
            "ddtCorrByA",
           -phiCorrCoeff*phiCorrs[phasei]*fvc::interpolate
            (
                byDt(alpha.oldTime()*phase.rho()().oldTime()*rAUs[phasei])
            ),
            ddtCorrByAs
        );
    }

    // Add virtual mass correction
    if (includeVirtualMass)
    {
        forAllConstIter(VmTable, Vms_, VmIter)
        {
            const volScalarField& Vm(*VmIter());
            const phasePair& pair(this->phasePairs_[VmIter.key()]);

            forAllConstIter(phasePair, pair, iter)
            {
                const phaseModel& phase = iter();
                const phaseModel& otherPhase = iter.otherPhase();

                addField
                (
                    iter(),
                    "ddtCorrByA",
                   -fvc::interpolate(Vm*byDt(rAUs[phase.index()]))
                   *(
                       phiCorrs[phase.index()]
                     + fvc::absolute
                       (
                           this->MRF().absolute(otherPhase.phi()),
                           otherPhase.U()
                       )
                     - fvc::flux(otherPhase.U())
                     - phiCorrs[otherPhase.index()]
                    ),
                    ddtCorrByAs
                );
            }
        }
    }

    return ddtCorrByAs;
}


template<class BasePhaseSystem>
void Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::partialElimination
(
    const PtrList<volScalarField>& rAUs,
    const PtrList<volVectorField>& KdUByAs,
    const PtrList<surfaceScalarField>& alphafs,
    const PtrList<surfaceScalarField>& phiKdPhis
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
            phases[i].URef() += KdUByAs[i];
            phases[i].phiRef() += phiKdPhis[i];
        }
    }

    {
        // Create drag coefficient matrices
        PtrList<PtrList<volScalarField>> KdByAs(phases.size());
        PtrList<PtrList<surfaceScalarField>> phiKds(phases.size());

        forAll(phases, i)
        {
            KdByAs.set
            (
                i,
                new PtrList<volScalarField>(phases.size())
            );

            phiKds.set
            (
                i,
                new PtrList<surfaceScalarField>(phases.size())
            );
        }

        forAllConstIter(KdTable, Kds_, KdIter)
        {
            const volScalarField& K(*KdIter());
            const phasePair& pair(this->phasePairs_[KdIter.key()]);

            const label phase1i = pair.phase1().index();
            const label phase2i = pair.phase2().index();

            addField
            (
                pair.phase2(),
                "KdByA",
                -rAUs[phase1i]*K,
                KdByAs[phase1i]
            );
            addField
            (
                pair.phase1(),
                "KdByA",
                -rAUs[phase2i]*K,
                KdByAs[phase2i]
            );

            addField
            (
                pair.phase2(),
                "phiKd",
                fvc::interpolate(KdByAs[phase1i][phase2i]),
                phiKds[phase1i]
            );
            addField
            (
                pair.phase1(),
                "phiKd",
                fvc::interpolate(KdByAs[phase2i][phase1i]),
                phiKds[phase2i]
            );
        }

        forAll(phases, i)
        {
            this->fillFields("KdByAs", dimless, KdByAs[i]);
            this->fillFields("phiKds", dimless, phiKds[i]);

            KdByAs[i][i] = 1;
            phiKds[i][i] = 1;
        }

        // Decompose
        for (label i = 0; i < phases.size(); i++)
        {
            for (label j = i + 1; j < phases.size(); j++)
            {
                KdByAs[i][j] /= KdByAs[i][i];
                phiKds[i][j] /= phiKds[i][i];
                for (label k = i + 1; k < phases.size(); ++ k)
                {
                    KdByAs[j][k] -= KdByAs[j][i]*KdByAs[i][k];
                    phiKds[j][k] -= phiKds[j][i]*phiKds[i][k];
                }
            }
        }

        {
            volScalarField detKdByAs(KdByAs[0][0]);
            surfaceScalarField detPhiKdfs(phiKds[0][0]);

            for (label i = 1; i < phases.size(); i++)
            {
                detKdByAs *= KdByAs[i][i];
                detPhiKdfs *= phiKds[i][i];
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
                    phases[i].phiRef() -= phiKds[i][j]*phases[j].phi();
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
                    phases[i].phiRef() -= phiKds[i][j]*phases[j].phi();
                }
                phases[i].URef() /= KdByAs[i][i];
                phases[i].phiRef() /= phiKds[i][i];
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
    const PtrList<surfaceScalarField>& phiKdPhifs
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
            phases[i].phiRef() += phiKdPhifs[i];
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
            const surfaceScalarField& K(*KdfIter());
            const phasePair& pair(this->phasePairs_[KdfIter.key()]);

            const label phase1i = pair.phase1().index();
            const label phase2i = pair.phase2().index();

            addField
            (
                pair.phase2(),
                "phiKdf",
                -rAUfs[phase1i]*K,
                phiKdfs[phase1i]
            );
            addField
            (
                pair.phase1(),
                "phiKdf",
                -rAUfs[phase2i]*K,
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
