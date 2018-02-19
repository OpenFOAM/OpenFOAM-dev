/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "PopulationBalancePhaseSystem.H"


// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class BasePhaseSystem>
void Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::
addMomentumTransfer(phaseSystem::momentumTransferTable& eqns) const
{
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

        const volVectorField& U1(pair.phase1().U());
        const volVectorField& U2(pair.phase2().U());

        const volScalarField dmdt(this->pDmdt(pair));
        const volScalarField dmdt21(posPart(dmdt));
        const volScalarField dmdt12(negPart(dmdt));

        *eqns[pair.phase1().name()] += dmdt21*U2 - fvm::Sp(dmdt21, U1);
        *eqns[pair.phase2().name()] -= dmdt12*U1 - fvm::Sp(dmdt12, U2);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::
PopulationBalancePhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh),

    populationBalances_
    (
        this->lookup("populationBalances"),
        diameterModels::populationBalanceModel::iNew(*this, pDmdt_)
    )
{
    forAll(populationBalances_, i)
    {
        const Foam::diameterModels::populationBalanceModel& popBal =
            populationBalances_[i];

        forAllConstIter(phaseSystem::phasePairTable, popBal.phasePairs(), iter)
        {
            const phasePairKey& key = iter.key();

            if (!this->phasePairs_.found(key))
            {
                this->phasePairs_.insert
                (
                    key,
                    autoPtr<phasePair>
                    (
                        new phasePair
                        (
                            this->phaseModels_[key.first()],
                            this->phaseModels_[key.second()]
                        )
                    )
                );
            }
        }
    }

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

        // Initially assume no mass transfer
        pDmdt_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("pDmdt", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar("zero", dimDensity/dimTime, 0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::
~PopulationBalancePhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::pDmdt
(
    const phasePairKey& key
) const
{
    if (pDmdt_.found(key))
    {
        const scalar pDmdtSign
        (
            Pair<word>::compare(pDmdt_.find(key).key(), key)
        );

        return pDmdtSign**pDmdt_[key];
    }

    tmp<volScalarField> tpDmdt
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

    return tpDmdt;
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::dmdt
(
    const phaseModel& phase
) const
{
    tmp<volScalarField> tDmdtFromKey
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("dmdtFromKey", phase.name()),
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

        if (pair.contains(phase))
        {
            tDmdtFromKey.ref() += this->dmdt
            (
                phasePairKey(phase.name(), pair.otherPhase(phase).name(), false)
            );
        }
    }

    return tDmdtFromKey;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::momentumTransfer() const
{
    autoPtr<phaseSystem::momentumTransferTable>
        eqnsPtr(BasePhaseSystem::momentumTransfer());

    addMomentumTransfer(eqnsPtr());

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::momentumTransferf() const
{
    autoPtr<phaseSystem::momentumTransferTable>
        eqnsPtr(BasePhaseSystem::momentumTransferf());

    addMomentumTransfer(eqnsPtr());

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::massTransferTable>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::heatTransfer() const
{
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr =
        BasePhaseSystem::heatTransfer();

    phaseSystem::heatTransferTable& eqns = eqnsPtr();

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

        const volScalarField& he1(pair.phase1().thermo().he());
        const volScalarField& he2(pair.phase2().thermo().he());

        const volScalarField dmdt(this->pDmdt(pair));
        const volScalarField dmdt21(posPart(dmdt));
        const volScalarField dmdt12(negPart(dmdt));

        *eqns[pair.phase1().name()] += dmdt21*he2 - fvm::Sp(dmdt21, he1);
        *eqns[pair.phase2().name()] -= dmdt12*he1 - fvm::Sp(dmdt12, he2);
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::massTransferTable>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::massTransfer() const
{
    autoPtr<phaseSystem::massTransferTable> eqnsPtr =
        BasePhaseSystem::massTransfer();

    phaseSystem::massTransferTable& eqns = eqnsPtr();

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

        const phaseModel& phase = pair.phase1();
        const phaseModel& otherPhase = pair.phase2();

        const volScalarField dmdt(this->pDmdt(pair));
        const volScalarField dmdt12(posPart(dmdt));
        const volScalarField dmdt21(negPart(dmdt));

        const PtrList<volScalarField>& Yi = phase.Y();

        forAll(Yi, i)
        {
            const word name
            (
                IOobject::groupName(Yi[i].member(), phase.name())
            );

            const word otherName
            (
                IOobject::groupName(Yi[i].member(), otherPhase.name())
            );

            *eqns[name] +=
                dmdt21*eqns[otherName]->psi()
              - fvm::Sp(dmdt21, eqns[name]->psi());

            *eqns[otherName] -=
                dmdt12*eqns[name]->psi()
              - fvm::Sp(dmdt12, eqns[otherName]->psi());
        }
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
bool Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::read()
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


template<class BasePhaseSystem>
void Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::solve()
{
    BasePhaseSystem::solve();

    forAll(populationBalances_, i)
    {
        populationBalances_[i].solve();
    }
}


// ************************************************************************* //
