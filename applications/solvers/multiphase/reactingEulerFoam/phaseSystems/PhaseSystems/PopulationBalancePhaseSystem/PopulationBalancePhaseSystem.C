/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2019 OpenFOAM Foundation
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
#include "rhoReactionThermo.H"


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

            pDmdt_.insert
            (
                key,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "pDmdt",
                            this->phasePairs_[key]->name()
                        ),
                        this->mesh().time().timeName(),
                        this->mesh(),
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    this->mesh(),
                    dimensionedScalar(dimDensity/dimTime, 0)
                )
            );
        }
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
Foam::PtrList<Foam::volScalarField>
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

    return dmdts;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::heatTransfer() const
{
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr =
        BasePhaseSystem::heatTransfer();

    phaseSystem::heatTransferTable& eqns = eqnsPtr();

    forAllConstIter
    (
        phaseSystem::phasePairTable,
        this->phasePairs_,
        phasePairIter
    )
    {
        const phasePair& pair(phasePairIter());

        if (pair.ordered() || !pDmdt_.found(pair))
        {
            continue;
        }

        forAll(populationBalances_, i)
        {
            const Foam::diameterModels::populationBalanceModel& popBal =
                populationBalances_[i];

            if (!popBal.phasePairs().found(pair))
            {
                continue;
            }

            // No correction needed for velocity group pairs
            if (popBal.isVelocityGroupPair(pair))
            {
                break;
            }

            const phaseModel& phase1 = pair.phase1();
            const phaseModel& phase2 = pair.phase2();

            const HashPtrTable<volScalarField>& sDmdt =
                popBal.speciesDmdt(pair);

            const rhoReactionThermo& thermo1 =
                refCast<const rhoReactionThermo>
                (
                    phase1.thermo()
                );

             const rhoReactionThermo& thermo2 =
                refCast<const rhoReactionThermo>
                (
                    phase2.thermo()
                );

            const volScalarField& p(thermo1.p());
            const volScalarField& T1(thermo1.T());
            const volScalarField& T2(thermo2.T());

            forAllConstIter
            (
                HashPtrTable<volScalarField>,
                sDmdt,
                sDmdtIter
            )
            {
                const label I1 =
                    thermo1.composition().species()[sDmdtIter.key()];

                const label I2 =
                    thermo2.composition().species()[sDmdtIter.key()];

                const volScalarField sDmdt12(negPart(**sDmdtIter));
                const volScalarField sDmdt21(posPart(**sDmdtIter));

                const volScalarField T
                (
                    neg0(**sDmdtIter)*T1 + pos(**sDmdtIter)*T2
                );

                const volScalarField L
                (
                    thermo1.composition().Ha(I1, p, T)
                  - thermo2.composition().Ha(I2, p, T)
                );

                *eqns[phase1.name()] -= sDmdt21*L;
                *eqns[phase2.name()] -= sDmdt12*L;
            }
        }
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

        if (pair.ordered() || !pDmdt_.found(pair))
        {
            continue;
        }

        const phaseModel& phase = pair.phase1();
        const phaseModel& otherPhase = pair.phase2();

        forAll(populationBalances_, i)
        {
            const Foam::diameterModels::populationBalanceModel& popBal =
                populationBalances_[i];

            if (popBal.phasePairs().found(pair))
            {
                // Both phases are velocity groups and belong to same population
                // balance -> transfer all species proportionally
                if (popBal.isVelocityGroupPair(pair))
                {
                    // Note that the phase YiEqn does not contain
                    // a continuity error term,
                    // so these additions represent the entire mass transfer

                    const volScalarField dmdt(this->pDmdt(pair));
                    const volScalarField dmdt12(negPart(dmdt));
                    const volScalarField dmdt21(posPart(dmdt));
                    const PtrList<volScalarField>& Yi = phase.Y();

                    forAll(Yi, i)
                    {
                        const word name
                        (
                            IOobject::groupName(Yi[i].member(), phase.name())
                        );

                        const word otherName
                        (
                            IOobject::groupName
                            (
                                Yi[i].member(),
                                otherPhase.name()
                            )
                        );

                        *eqns[name] +=
                            dmdt21*eqns[otherName]->psi()
                          + fvm::Sp(dmdt12, eqns[name]->psi());

                        *eqns[otherName] -=
                            dmdt12*eqns[name]->psi()
                          + fvm::Sp(dmdt21, eqns[otherName]->psi());
                    }
                }
                // The phases do not belong to the same population balance,
                // transfer each specie separately
                else
                {

                    const HashPtrTable<volScalarField>& sDmdt =
                        popBal.speciesDmdt(pair);

                    forAllConstIter
                    (
                        HashPtrTable<volScalarField>,
                        sDmdt,
                        sDmdtIter
                    )
                    {
                        const word name
                        (
                            IOobject::groupName(sDmdtIter.key(), phase.name())
                        );

                        const word otherName
                        (
                            IOobject::groupName
                            (
                                sDmdtIter.key(),
                                otherPhase.name()
                            )
                        );

                        const dimensionedScalar Yismall(dimless, rootVSmall);

                        const volScalarField& Y1 = eqns[name]->psi();
                        const volScalarField& Y2 = eqns[otherName]->psi();

                        const volScalarField sDmdt12(negPart(**sDmdtIter));
                        const volScalarField sDmdt21(posPart(**sDmdtIter));

                        *eqns[name] +=
                            sDmdt21
                          + fvm::Sp(sDmdt12/(Y1 + Yismall), Y1);

                        *eqns[otherName] -=
                            sDmdt12
                          + fvm::Sp(sDmdt21/(Y2 + Yismall), Y2);
                    }

                }
            }
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
