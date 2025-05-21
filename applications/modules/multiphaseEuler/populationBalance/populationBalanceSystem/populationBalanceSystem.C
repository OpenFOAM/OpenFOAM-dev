/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2025 OpenFOAM Foundation
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

#include "populationBalanceSystem.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(populationBalanceSystem, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::populationBalanceSystem::addDmdts
(
    const diameterModels::populationBalanceModel::dmdtfTable& dmdtfs,
    PtrList<volScalarField::Internal>& dmdts
) const
{
    forAll(populationBalances_, popBali)
    {
        forAllConstIter
        (
            diameterModels::populationBalanceModel::dmdtfTable,
            dmdtfs,
            dmdtfIter
        )
        {
            const phaseInterface interface(fluid_, dmdtfIter.key());

            addField(interface.phase1(), "dmdt", *dmdtfIter(), dmdts);
            addField(interface.phase2(), "dmdt", - *dmdtfIter(), dmdts);
        }
    }
}


void Foam::populationBalanceSystem::addDmdtUfs
(
    const diameterModels::populationBalanceModel::dmdtfTable& dmdtfs,
    HashPtrTable<fvVectorMatrix>& eqns
) const
{
    forAllConstIter
    (
        diameterModels::populationBalanceModel::dmdtfTable,
        dmdtfs,
        dmdtfIter
    )
    {
        const phaseInterface interface(fluid_, dmdtfIter.key());

        const volScalarField::Internal& dmdtf = *dmdtfIter();
        const volScalarField::Internal dmdtf21(posPart(dmdtf));
        const volScalarField::Internal dmdtf12(negPart(dmdtf));

        const phaseModel& phase1 = fluid_.phases()[interface.phase1().name()];
        const phaseModel& phase2 = fluid_.phases()[interface.phase2().name()];

        if (!phase1.stationary())
        {
            *eqns[phase1.name()] +=
                dmdtf21*phase2.U()()() + fvm::Sp(dmdtf12, phase1.URef());
        }

        if (!phase2.stationary())
        {
            *eqns[phase2.name()] -=
                dmdtf12*phase1.U()()() + fvm::Sp(dmdtf21, phase2.URef());
        }
    }
}


void Foam::populationBalanceSystem::addDmdtHefs
(
    const diameterModels::populationBalanceModel::dmdtfTable& dmdtfs,
    HashPtrTable<fvScalarMatrix>& eqns
) const
{
    // Loop the pairs
    forAllConstIter
    (
        diameterModels::populationBalanceModel::dmdtfTable,
        dmdtfs,
        dmdtfIter
    )
    {
        const phaseInterface interface(fluid_, dmdtfIter.key());

        const volScalarField::Internal& dmdtf = *dmdtfIter();
        const volScalarField::Internal dmdtf21(posPart(dmdtf));
        const volScalarField::Internal dmdtf12(negPart(dmdtf));

        const phaseModel& phase1 = interface.phase1();
        const phaseModel& phase2 = interface.phase2();
        const rhoFluidThermo& thermo1 = phase1.fluidThermo();
        const rhoFluidThermo& thermo2 = phase2.fluidThermo();
        const volScalarField& he1 = thermo1.he();
        const volScalarField& he2 = thermo2.he();
        const volScalarField::Internal hs1(thermo1.hs());
        const volScalarField::Internal hs2(thermo2.hs());
        const volScalarField::Internal K1(phase1.K());
        const volScalarField::Internal K2(phase2.K());

        // Transfer of sensible enthalpy within the phases
        *eqns[phase1.name()] +=
            dmdtf*hs1 + fvm::Sp(dmdtf12, he1) - dmdtf12*he1;
        *eqns[phase2.name()] -=
            dmdtf*hs2 + fvm::Sp(dmdtf21, he2) - dmdtf21*he2;

        // Transfer of sensible enthalpy between the phases
        *eqns[phase1.name()] += dmdtf21*(hs2 - hs1);
        *eqns[phase2.name()] -= dmdtf12*(hs1 - hs2);

        // Transfer of kinetic energy
        *eqns[phase1.name()] += dmdtf21*K2 + dmdtf12*K1;
        *eqns[phase2.name()] -= dmdtf12*K1 + dmdtf21*K2;
    }
}


void Foam::populationBalanceSystem::addDmdtYfs
(
    const diameterModels::populationBalanceModel::dmdtfTable& dmdtfs,
    HashPtrTable<fvScalarMatrix>& eqns
) const
{
    forAllConstIter
    (
        diameterModels::populationBalanceModel::dmdtfTable,
        dmdtfs,
        dmdtfIter
    )
    {
        const phaseInterface interface(fluid_, dmdtfIter.key());

        const volScalarField::Internal& dmdtf = *dmdtfIter();
        const volScalarField::Internal dmdtf12(negPart(dmdtf));
        const volScalarField::Internal dmdtf21(posPart(dmdtf));

        const phaseModel& phase1 = interface.phase1();
        const phaseModel& phase2 = interface.phase2();

        forAll(phase1.Y(), Yi1)
        {
            const volScalarField& Y1 = phase1.Y()[Yi1];
            const volScalarField& Y2 = phase2.Y(Y1.member());

            *eqns[Y1.name()] += dmdtf21*Y2 + fvm::Sp(dmdtf12, Y1);
            *eqns[Y2.name()] -= dmdtf12*Y1 + fvm::Sp(dmdtf21, Y2);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSystem::populationBalanceSystem
(
    const phaseSystem& fluid
)
:
    fluid_(fluid),
    populationBalances_()
{
    if (fluid_.found("populationBalances"))
    {
        PtrList<diameterModels::populationBalanceModel> popBals
        (
            fluid_.lookup("populationBalances"),
            diameterModels::populationBalanceModel::iNew(fluid_)
        );

        populationBalances_.transfer(popBals);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSystem::~populationBalanceSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::PtrList<Foam::volScalarField::Internal>
Foam::populationBalanceSystem::dmdts() const
{
    PtrList<volScalarField::Internal> dmdts(fluid_.phases().size());

    forAll(populationBalances_, popBali)
    {
        addDmdts(populationBalances_[popBali].dmdtfs(), dmdts);
        addDmdts(populationBalances_[popBali].expansionDmdtfs(), dmdts);
        addDmdts(populationBalances_[popBali].modelSourceDmdtfs(), dmdts);
    }

    return dmdts;
}


Foam::autoPtr<Foam::HashPtrTable<Foam::fvVectorMatrix>>
Foam::populationBalanceSystem::momentumTransfer()
{
    autoPtr<HashPtrTable<fvVectorMatrix>> eqnsPtr
    (
        new HashPtrTable<fvVectorMatrix>()
    );
    HashPtrTable<fvVectorMatrix>& eqns = eqnsPtr();

    forAll(fluid_.movingPhases(), movingPhasei)
    {
        const phaseModel& phase = fluid_.movingPhases()[movingPhasei];

        eqns.insert
        (
            phase.name(),
            new fvVectorMatrix(phase.U(), dimMass*dimVelocity/dimTime)
        );
    }

    forAll(populationBalances_, popBali)
    {
        addDmdtUfs(populationBalances_[popBali].dmdtfs(), eqns);
        addDmdtUfs(populationBalances_[popBali].expansionDmdtfs(), eqns);
        addDmdtUfs(populationBalances_[popBali].modelSourceDmdtfs(), eqns);
    }

    return eqnsPtr;
}


Foam::autoPtr<Foam::HashPtrTable<Foam::fvVectorMatrix>>
Foam::populationBalanceSystem::momentumTransferf()
{
    return momentumTransfer();
}


Foam::autoPtr<Foam::HashPtrTable<Foam::fvScalarMatrix>>
Foam::populationBalanceSystem::heatTransfer() const
{
    autoPtr<HashPtrTable<fvScalarMatrix>> eqnsPtr
    (
        new HashPtrTable<fvScalarMatrix>()
    );
    HashPtrTable<fvScalarMatrix>& eqns = eqnsPtr();

    forAll(fluid_.phases(), phasei)
    {
        const phaseModel& phase = fluid_.phases()[phasei];

        eqns.insert
        (
            phase.name(),
            new fvScalarMatrix(phase.thermo().he(), dimEnergy/dimTime)
        );
    }

    forAll(populationBalances_, popBali)
    {
        addDmdtHefs(populationBalances_[popBali].dmdtfs(), eqns);
        addDmdtHefs(populationBalances_[popBali].expansionDmdtfs(), eqns);
        addDmdtHefs(populationBalances_[popBali].modelSourceDmdtfs(), eqns);
    }

    return eqnsPtr;
}


Foam::autoPtr<Foam::HashPtrTable<Foam::fvScalarMatrix>>
Foam::populationBalanceSystem::specieTransfer() const
{
    autoPtr<HashPtrTable<fvScalarMatrix>> eqnsPtr
    (
        new HashPtrTable<fvScalarMatrix>()
    );
    HashPtrTable<fvScalarMatrix>& eqns = eqnsPtr();

    forAll(fluid_.multicomponentPhases(), multicomponentPhasei)
    {
        const phaseModel& phase =
            fluid_.multicomponentPhases()[multicomponentPhasei];

        const UPtrList<volScalarField>& Y = phase.Y();

        forAll(Y, i)
        {
            eqns.insert
            (
                Y[i].name(),
                new fvScalarMatrix(Y[i], dimMass/dimTime)
            );
        }
    }

    forAll(populationBalances_, popBali)
    {
        addDmdtYfs(populationBalances_[popBali].dmdtfs(), eqns);
        addDmdtYfs(populationBalances_[popBali].expansionDmdtfs(), eqns);
        addDmdtYfs(populationBalances_[popBali].modelSourceDmdtfs(), eqns);
    }

    return eqnsPtr;
}


void Foam::populationBalanceSystem::solve()
{
    forAll(populationBalances_, i)
    {
        populationBalances_[i].solve();
    }
}


void Foam::populationBalanceSystem::correct()
{
    forAll(populationBalances_, i)
    {
        populationBalances_[i].correct();
    }
}


// ************************************************************************* //
