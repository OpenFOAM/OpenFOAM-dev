/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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

#include "PhaseTransferPhaseSystem.H"
#include "phaseTransferModel.H"
#include "fvmSup.H"

// * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * //

template<class BasePhaseSystem>
void Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::addDmdtY
(
    const phaseSystem::dmdtTable& dmdts,
    phaseSystem::massTransferTable& eqns
) const
{
    forAllConstIter(phaseSystem::dmdtTable, dmdts, dmdtIter)
    {
        const phasePairKey& key = dmdtIter.key();
        const phasePair& pair(this->phasePairs_[key]);

        const volScalarField dmdt(Pair<word>::compare(pair, key)**dmdtIter());
        const volScalarField dmdt12(negPart(dmdt));
        const volScalarField dmdt21(posPart(dmdt));

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        // Note that the phase YiEqn does not contain a continuity error term,
        // so the transfers below are complete.

        forAll(phase1.Y(), Yi1)
        {
            const volScalarField& Y1 = phase1.Y()[Yi1];
            const volScalarField& Y2 = phase2.Y(Y1.member());

            *eqns[Y1.name()] += dmdt21*Y2 + fvm::Sp(dmdt12, Y1);
            *eqns[Y2.name()] -= dmdt12*Y1 + fvm::Sp(dmdt21, Y2);
        }
    }
}


template<class BasePhaseSystem>
void Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::addDmidtY
(
    const phaseSystem::dmidtTable& dmidts,
    phaseSystem::massTransferTable& eqns
) const
{
    forAllConstIter(phaseSystem::dmidtTable, dmidts, dmidtIter)
    {
        const phasePairKey& key = dmidtIter.key();
        const phasePair& pair(this->phasePairs_[key]);

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        // Note that the phase YiEqn does not contain a continuity error term,
        // so the transfers below are complete.

        forAllConstIter(HashPtrTable<volScalarField>, *dmidtIter(), dmidtJter)
        {
            const word& member = dmidtJter.key();

            const volScalarField dmidt
            (
                Pair<word>::compare(pair, key)**dmidtJter()
            );

            if (!phase1.pure())
            {
                const volScalarField& Y1 = phase1.Y(member);
                *eqns[Y1.name()] += dmidt;
            }

            if (!phase2.pure())
            {
                const volScalarField& Y2 = phase2.Y(member);
                *eqns[Y2.name()] -= dmidt;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::PhaseTransferPhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh)
{
    this->generatePairsAndSubModels
    (
        "phaseTransfer",
        phaseTransferModels_,
        false
    );

    forAllConstIter
    (
        phaseTransferModelTable,
        phaseTransferModels_,
        phaseTransferModelIter
    )
    {
        this->rDmdt_.insert
        (
            phaseTransferModelIter.key(),
            phaseSystem::dmdt(phaseTransferModelIter.key()).ptr()
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::
~PhaseTransferPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::dmdt
(
    const phasePairKey& key
) const
{
    NotImplemented;

    return phaseSystem::dmdt(key);
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::dmdts() const
{
    PtrList<volScalarField> dmdts(BasePhaseSystem::dmdts());

    forAllConstIter(phaseSystem::dmdtTable, rDmdt_, rDmdtIter)
    {
        const phasePair& pair = this->phasePairs_[rDmdtIter.key()];
        const volScalarField& rDmdt = *rDmdtIter();

        this->addField(pair.phase1(), "dmdt", rDmdt, dmdts);
        this->addField(pair.phase2(), "dmdt", - rDmdt, dmdts);
    }

    return dmdts;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::momentumTransfer()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransfer();

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    this->addDmdtU(rDmdt_, eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::momentumTransferf()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransferf();

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    this->addDmdtU(rDmdt_, eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::heatTransfer() const
{
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr =
        BasePhaseSystem::heatTransfer();

    phaseSystem::heatTransferTable& eqns = eqnsPtr();

    this->addDmdtHe(rDmdt_, eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::massTransferTable>
Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::massTransfer() const
{
    autoPtr<phaseSystem::massTransferTable> eqnsPtr
    (
        new phaseSystem::massTransferTable()
    );

    phaseSystem::massTransferTable& eqns = eqnsPtr();

    // Create a mass transfer matrix for each species of each phase
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

    this->addDmdtY(rDmdt_, eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
void Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::correct()
{
    BasePhaseSystem::correct();

    forAllConstIter
    (
        phaseTransferModelTable,
        phaseTransferModels_,
        phaseTransferModelIter
    )
    {
        *rDmdt_[phaseTransferModelIter.key()] =
            dimensionedScalar(dimDensity/dimTime, 0);
    }

    forAllConstIter
    (
        phaseTransferModelTable,
        phaseTransferModels_,
        phaseTransferModelIter
    )
    {
        *rDmdt_[phaseTransferModelIter.key()] +=
            phaseTransferModelIter()->dmdt();
    }
}


template<class BasePhaseSystem>
bool Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::read()
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
