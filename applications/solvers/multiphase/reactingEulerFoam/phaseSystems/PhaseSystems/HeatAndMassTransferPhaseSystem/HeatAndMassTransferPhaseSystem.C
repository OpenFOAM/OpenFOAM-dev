/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "HeatAndMassTransferPhaseSystem.H"

#include "BlendedInterfacialModel.H"
#include "heatTransferModel.H"
#include "massTransferModel.H"

#include "HashPtrTable.H"

#include "fvcDiv.H"
#include "fvmSup.H"
#include "fvMatrix.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::
HeatAndMassTransferPhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh)
{
    this->generatePairsAndSubModels
    (
        "heatTransfer",
        heatTransferModels_
    );

    this->generatePairsAndSubModels
    (
        "massTransfer",
        massTransferModels_
    );

    forAllConstIter
    (
        phaseSystem::phasePairTable,
        this->phasePairs_,
        phasePairIter
    )
    {
        const phasePair& pair(phasePairIter());

        if (pair.ordered())
        {
            continue;
        }

        if
        (
            heatTransferModels_.found(pair)
         && heatTransferModels_[pair][pair.first()].valid()
         && heatTransferModels_[pair][pair.second()].valid()
        )
        {

            volScalarField H1(heatTransferModels_[pair][pair.first()]->K());
            volScalarField H2(heatTransferModels_[pair][pair.second()]->K());

            Tf_.insert
            (
                pair,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName("Tf", pair.name()),
                        this->mesh().time().timeName(),
                        this->mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    (
                        H1*pair.phase1().thermo().T()
                      + H2*pair.phase2().thermo().T()
                    )
                   /max
                    (
                        H1 + H2,
                        dimensionedScalar
                        (
                            "small",
                            heatTransferModel::dimK,
                            small
                        )
                    ),
                    zeroGradientFvPatchScalarField::typeName
                )
            );

            Tf_[pair]->correctBoundaryConditions();
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::
~HeatAndMassTransferPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
bool
Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::transfersMass() const
{
    return true;
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::dmdt
(
    const phasePairKey& key
) const
{
    return tmp<Foam::volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "dmdt",
                this->mesh_.time().timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar("zero", dimDensity/dimTime, 0)
        )
    );
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::dmdt
(
    const Foam::phaseModel& phase
) const
{
    tmp<volScalarField> tDmdt
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("dmdt", phase.name()),
                this->mesh_.time().timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar("zero", dimDensity/dimTime, 0)
        )
    );

    return tDmdt;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::momentumTransfer() const
{
    autoPtr<phaseSystem::momentumTransferTable>
        eqnsPtr(BasePhaseSystem::momentumTransfer());

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::momentumTransferf() const
{
    autoPtr<phaseSystem::momentumTransferTable>
        eqnsPtr(BasePhaseSystem::momentumTransferf());

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::heatTransfer() const
{
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr
    (
        new phaseSystem::heatTransferTable()
    );

    phaseSystem::heatTransferTable& eqns = eqnsPtr();

    forAll(this->phaseModels_, phasei)
    {
        const phaseModel& phase = this->phaseModels_[phasei];

        eqns.insert
        (
            phase.name(),
            new fvScalarMatrix(phase.thermo().he(), dimEnergy/dimTime)
        );
    }

    // Heat transfer with the interface
    forAllConstIter
    (
        heatTransferModelTable,
        heatTransferModels_,
        heatTransferModelIter
    )
    {
        if
        (
            heatTransferModels_.found(heatTransferModelIter.key())
        )
        {
            const phasePair& pair
            (
                this->phasePairs_[heatTransferModelIter.key()]
            );

            const volScalarField& Tf(*Tf_[pair]);

            const volScalarField K1
            (
                heatTransferModelIter()[pair.first()]->K()
            );
            const volScalarField K2
            (
                heatTransferModelIter()[pair.second()]->K()
            );
            const volScalarField KEff
            (
                K1*K2
               /max
                (
                    K1 + K2,
                    dimensionedScalar("small", heatTransferModel::dimK, small)
                )
            );

            const volScalarField* K = &K1;
            const volScalarField* otherK = &K2;

            forAllConstIter(phasePair, pair, iter)
            {
                const phaseModel& phase = iter();

                const volScalarField& he(phase.thermo().he());
                volScalarField Cpv(phase.thermo().Cpv());

                *eqns[phase.name()] +=
                    (*K)*(Tf - phase.thermo().T())
                  + KEff/Cpv*he - fvm::Sp(KEff/Cpv, he);

                Swap(K, otherK);
            }
        }
    }

    return eqnsPtr;
}

template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::massTransferTable>
Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::massTransfer() const
{
    // Create a mass transfer matrix for each species of each phase
    autoPtr<phaseSystem::massTransferTable> eqnsPtr
    (
        new phaseSystem::massTransferTable()
    );

    phaseSystem::massTransferTable& eqns = eqnsPtr();

    forAll(this->phaseModels_, phasei)
    {
        const phaseModel& phase = this->phaseModels_[phasei];

        const PtrList<volScalarField>& Yi = phase.Y();

        forAll(Yi, i)
        {
            eqns.insert
            (
                Yi[i].name(),
                new fvScalarMatrix(Yi[i], dimMass/dimTime)
            );
        }
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
void Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::correctThermo()
{

    phaseSystem::correctThermo();

    forAllConstIter
    (
        phaseSystem::phasePairTable,
        this->phasePairs_,
        phasePairIter
    )
    {
        if
        (
            this->heatTransferModels_.found(phasePairIter.key())
        )
        {
            const phasePair& pair(phasePairIter());

            if (pair.ordered())
            {
                continue;
            }

            const phaseModel& phase1 = pair.phase1();
            const phaseModel& phase2 = pair.phase2();

            const volScalarField& T1(phase1.thermo().T());
            const volScalarField& T2(phase2.thermo().T());

            volScalarField& Tf(*this->Tf_[pair]);

            volScalarField H1
            (
                this->heatTransferModels_[pair][pair.first()]->K()
            );

            volScalarField H2
            (
                this->heatTransferModels_[pair][pair.second()]->K()
            );

            // Limit the H[12] to avoid /0
            H1.max(small);
            H2.max(small);

            Tf = (H1*T1 + H2*T2)/(H1 + H2);

            Info<< "Tf." << pair.name()
                << ": min = " << min(Tf.primitiveField())
                << ", mean = " << average(Tf.primitiveField())
                << ", max = " << max(Tf.primitiveField())
                << endl;
        }
    }
}


template<class BasePhaseSystem>
bool Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::read()
{
    if (BasePhaseSystem::read())
    {
        bool readOK = true;

        // Models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
