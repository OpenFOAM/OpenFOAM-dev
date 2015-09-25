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

#include "HeatAndMassTransferPhaseSystem.H"

#include "BlendedInterfacialModel.H"
#include "heatTransferModel.H"
#include "massTransferModel.H"

#include "HashPtrTable.H"

#include "fvcDiv.H"
#include "fvmSup.H"
#include "fvMatrix.H"

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

        // Initialy assume no mass transfer

        dmdt_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("dmdt", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar("zero", dimDensity/dimTime, 0)
            )
        );

        dmdtExplicit_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("dmdtExplicit", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh()
                ),
                this->mesh(),
                dimensionedScalar("zero", dimDensity/dimTime, 0)
            )
        );

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
                    dimensionedScalar("small", heatTransferModel::dimK, SMALL)
                ),
                zeroGradientFvPatchScalarField::typeName
            )
        );
        Tf_[pair]->correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::
~HeatAndMassTransferPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
bool Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::transfersMass
(
    const phaseModel& phase
) const
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
    const scalar dmdtSign(Pair<word>::compare(dmdt_.find(key).key(), key));

    return dmdtSign**dmdt_[key];
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::dmdt
(
    const Foam::phaseModel& phase
) const
{
    tmp<volScalarField> tdmdt
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

        const phaseModel* phase1 = &pair.phase1();
        const phaseModel* phase2 = &pair.phase2();

        forAllConstIter(phasePair, pair, iter)
        {
            if (phase1 == &phase)
            {
                tdmdt() += this->dmdt(pair);
            }

            Swap(phase1, phase2);
        }
    }

    return tdmdt;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::momentumTransfer() const
{
    autoPtr<phaseSystem::momentumTransferTable>
        eqnsPtr(BasePhaseSystem::momentumTransfer());

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    // Source term due to mass trasfer
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

        const volVectorField& U1(pair.phase1().U());
        const volVectorField& U2(pair.phase2().U());

        const volScalarField dmdt(this->dmdt(pair));
        const volScalarField dmdt21(posPart(dmdt));
        const volScalarField dmdt12(negPart(dmdt));

        *eqns[pair.phase1().name()] += dmdt21*U2 - fvm::Sp(dmdt21, U1);
        *eqns[pair.phase2().name()] -= dmdt12*U1 - fvm::Sp(dmdt12, U2);
    }

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
        const phasePair& pair
        (
            this->phasePairs_[heatTransferModelIter.key()]
        );

        const phaseModel* phase = &pair.phase1();
        const phaseModel* otherPhase = &pair.phase2();

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
                dimensionedScalar("small", heatTransferModel::dimK, SMALL)
            )
        );

        const volScalarField* K = &K1;
        const volScalarField* otherK = &K2;

        forAllConstIter(phasePair, pair, iter)
        {
            const volScalarField& he(phase->thermo().he());
            volScalarField Cpv(phase->thermo().Cpv());

            *eqns[phase->name()] +=
                (*K)*(Tf - phase->thermo().T())
              + KEff/Cpv*he - fvm::Sp(KEff/Cpv, he);

            Swap(phase, otherPhase);
            Swap(K, otherK);
        }
    }

    // Source term due to mass transfer
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

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        const volScalarField& he1(phase1.thermo().he());
        const volScalarField& he2(phase2.thermo().he());

        const volScalarField& K1(phase1.K());
        const volScalarField& K2(phase2.K());

        const volScalarField dmdt(this->dmdt(pair));
        const volScalarField dmdt21(posPart(dmdt));
        const volScalarField dmdt12(negPart(dmdt));
        const volScalarField& Tf(*Tf_[pair]);

        *eqns[phase1.name()] +=
            dmdt21*(phase1.thermo().he(phase1.thermo().p(), Tf))
          - fvm::Sp(dmdt21, he1)
          + dmdt21*(K2 - K1);

        *eqns[phase2.name()] -=
            dmdt12*(phase2.thermo().he(phase2.thermo().p(), Tf))
          - fvm::Sp(dmdt12, he2)
          + dmdt12*(K1 - K2);
    }

    return eqnsPtr;
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
