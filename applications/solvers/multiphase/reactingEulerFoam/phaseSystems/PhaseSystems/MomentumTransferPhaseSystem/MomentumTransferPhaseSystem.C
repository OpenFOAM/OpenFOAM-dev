/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
#include "fvcDiv.H"
#include "fvMatrix.H"

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
                IOobject::groupName("Kd", pair.name()),
                dragModelIter()->K()
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
                IOobject::groupName("Vm", pair.name()),
                virtualMassModelIter()->K()
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
Foam::tmp<Foam::volScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::Kd
(
    const phasePairKey& key
) const
{
    return dragModels_[key]->K();
}


template<class BasePhaseSystem>
Foam::tmp<Foam::surfaceScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::Kdf
(
    const phasePairKey& key
) const
{
    return dragModels_[key]->Kf();
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
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    virtualMassModel::typeName + ":K",
                    this->mesh_.time().timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                this->mesh_,
                dimensionedScalar("zero", virtualMassModel::dimK, 0)
            )
        );
    }
}


template<class BasePhaseSystem>
Foam::tmp<Foam::surfaceScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::Vmf
(
    const phasePairKey& key
) const
{
    if (virtualMassModels_.found(key))
    {
        return virtualMassModels_[key]->Kf();
    }
    else
    {
        return tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                IOobject
                (
                    virtualMassModel::typeName + ":Kf",
                    this->mesh_.time().timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                this->mesh_,
                dimensionedScalar("zero", virtualMassModel::dimK, 0)
            )
        );
    }
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volVectorField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::F
(
    const phasePairKey& key
) const
{
    if (liftModels_.found(key) && wallLubricationModels_.found(key))
    {
        return
            liftModels_[key]->template F<vector>()
          + wallLubricationModels_[key]->template F<vector>();
    }
    else if (liftModels_.found(key))
    {
        return liftModels_[key]->template F<vector>();
    }
    else if (wallLubricationModels_.found(key))
    {
        return wallLubricationModels_[key]->template F<vector>();
    }
    else
    {
        return tmp<volVectorField>
        (
            new volVectorField
            (
                IOobject
                (
                    liftModel::typeName + ":F",
                    this->mesh_.time().timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                this->mesh_,
                dimensionedVector("zero", liftModel::dimF, vector::zero)
            )
        );
    }
}


template<class BasePhaseSystem>
Foam::tmp<Foam::surfaceScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::Ff
(
    const phasePairKey& key
) const
{
    if (liftModels_.found(key) && wallLubricationModels_.found(key))
    {
        return
            liftModels_[key]->Ff()
          + wallLubricationModels_[key]->Ff();
    }
    else if (liftModels_.found(key))
    {
        return liftModels_[key]->Ff();
    }
    else if (wallLubricationModels_.found(key))
    {
        return wallLubricationModels_[key]->Ff();
    }
    else
    {
        return tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                IOobject
                (
                    liftModel::typeName + ":Ff",
                    this->mesh_.time().timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                this->mesh_,
                dimensionedScalar("zero", liftModel::dimF*dimArea, 0)
            )
        );
    }
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::D
(
    const phasePairKey& key
) const
{
    if (turbulentDispersionModels_.found(key))
    {
        return turbulentDispersionModels_[key]->D();
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    turbulentDispersionModel::typeName + ":D",
                    this->mesh_.time().timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                this->mesh_,
                dimensionedScalar("zero", turbulentDispersionModel::dimD, 0)
            )
        );
    }
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::momentumTransfer() const
{
    // Create a momentum transfer matrix for each phase
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr
    (
        new phaseSystem::momentumTransferTable()
    );

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    forAllConstIter
    (
        phaseSystem::phaseModelTable,
        this->phaseModels_,
        phaseModelIter
    )
    {
        const phaseModel& phase(phaseModelIter());

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
    }

    // Add the implicit part of the drag force
    forAllConstIter
    (
        phaseSystem::KdTable,
        Kds_,
        KdIter
    )
    {
        const volScalarField& K(*KdIter());

        const phasePair& pair(this->phasePairs_[KdIter.key()]);

        const phaseModel* phase = &pair.phase1();
        const phaseModel* otherPhase = &pair.phase2();

        forAllConstIter(phasePair, pair, iter)
        {
            const volVectorField& U(phase->U());

            *eqns[phase->name()] -= fvm::Sp(K, U);

            Swap(phase, otherPhase);
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
    forAllConstIter
    (
        phaseSystem::VmTable,
        Vms_,
        VmIter
    )
    {
        const volScalarField& Vm(*VmIter());

        const phasePair& pair(this->phasePairs_[VmIter.key()]);

        const phaseModel* phase = &pair.phase1();
        const phaseModel* otherPhase = &pair.phase2();

        forAllConstIter(phasePair, pair, iter)
        {
            const volVectorField& U(phase->U());
            const surfaceScalarField& phi(phase->phi());

            *eqns[phase->name()] +=
              - Vm
               *(
                    fvm::ddt(U)
                  + fvm::div(phi, U)
                  - fvm::Sp(fvc::div(phi), U)
                  - otherPhase->DUDt()
                )
              - this->MRF_.DDt(Vm, U - otherPhase->U());

            Swap(phase, otherPhase);
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
        const volVectorField F(liftModelIter()->F<vector>());

        const phasePair& pair(this->phasePairs_[liftModelIter.key()]);

        *eqns[pair.phase1().name()] += F;
        *eqns[pair.phase2().name()] -= F;
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

        *eqns[pair.phase1().name()] += F;
        *eqns[pair.phase2().name()] -= F;
    }

    // Add the turbulent dispersion force
    forAllConstIter
    (
        turbulentDispersionModelTable,
        turbulentDispersionModels_,
        turbulentDispersionModelIter
    )
    {
        const volVectorField F(turbulentDispersionModelIter()->F<vector>());

        const phasePair&
            pair(this->phasePairs_[turbulentDispersionModelIter.key()]);

        *eqns[pair.phase1().name()] += F;
        *eqns[pair.phase2().name()] -= F;
    }

    return eqnsPtr;
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
