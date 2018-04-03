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
Foam::tmp<Foam::volScalarField>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::pDmdt
(
    const phasePairKey& key
) const
{
    if (!pDmdt_.found(key))
    {
        return phaseSystem::dmdt(key);
    }

    const scalar pDmdtSign(Pair<word>::compare(pDmdt_.find(key).key(), key));

    return pDmdtSign**pDmdt_[key];
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
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::dmdt
(
    const phasePairKey& key
) const
{
    return BasePhaseSystem::dmdt(key) + this->pDmdt(key);
}


template<class BasePhaseSystem>
Foam::Xfer<Foam::PtrList<Foam::volScalarField>>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::dmdts() const
{
    PtrList<volScalarField> dmdts(BasePhaseSystem::dmdts());

    forAllConstIter(pDmdtTable, pDmdt_, pDmdtIter)
    {
        const phasePair& pair = this->phasePairs_[pDmdtIter.key()];
        const volScalarField& pDmdt = *pDmdtIter();

        this->addField(pair.phase1(), "dmdt", pDmdt, dmdts);
        this->addField(pair.phase2(), "dmdt", - pDmdt, dmdts);
    }

    return dmdts.xfer();
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
