/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2023 OpenFOAM Foundation
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
        diameterModels::populationBalanceModel::iNew(*this)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::
~PopulationBalancePhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::dmdtf
(
    const phaseInterfaceKey& key
) const
{
    tmp<volScalarField> tDmdtf = BasePhaseSystem::dmdtf(key);

    forAll(populationBalances_, popBali)
    {
        if (populationBalances_[popBali].dmdtfs().found(key))
        {
            tDmdtf.ref() += *populationBalances_[popBali].dmdtfs()[key];
        }
    }

    return tDmdtf;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::dmdts() const
{
    PtrList<volScalarField> dmdts(BasePhaseSystem::dmdts());

    forAll(populationBalances_, popBali)
    {
        forAllConstIter
        (
            phaseSystem::dmdtfTable,
            populationBalances_[popBali].dmdtfs(),
            dmdtfIter
        )
        {
            const phaseInterface interface(*this, dmdtfIter.key());

            addField(interface.phase1(), "dmdt", *dmdtfIter(), dmdts);
            addField(interface.phase2(), "dmdt", - *dmdtfIter(), dmdts);
        }
    }

    return dmdts;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::momentumTransfer()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransfer();

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    forAll(populationBalances_, popBali)
    {
        this->addDmdtUfs(populationBalances_[popBali].dmdtfs(), eqns);
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::momentumTransferf()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransferf();

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    forAll(populationBalances_, popBali)
    {
        this->addDmdtUfs(populationBalances_[popBali].dmdtfs(), eqns);
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::heatTransfer() const
{
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr =
        BasePhaseSystem::heatTransfer();

    phaseSystem::heatTransferTable& eqns = eqnsPtr();

    forAll(populationBalances_, popBali)
    {
        this->addDmdtHefs(populationBalances_[popBali].dmdtfs(), eqns);
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::specieTransferTable>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::specieTransfer() const
{
    autoPtr<phaseSystem::specieTransferTable> eqnsPtr =
        BasePhaseSystem::specieTransfer();

    phaseSystem::specieTransferTable& eqns = eqnsPtr();

    forAll(populationBalances_, popBali)
    {
        this->addDmdtYfs(populationBalances_[popBali].dmdtfs(), eqns);
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
void Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::solve
(
    const PtrList<volScalarField>& rAUs,
    const PtrList<surfaceScalarField>& rAUfs
)
{
    BasePhaseSystem::solve(rAUs, rAUfs);

    forAll(populationBalances_, i)
    {
        populationBalances_[i].solve();
    }
}


template<class BasePhaseSystem>
void Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::correct()
{
    BasePhaseSystem::correct();

    forAll(populationBalances_, i)
    {
        populationBalances_[i].correct();
    }
}


template<class BasePhaseSystem>
bool Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::read()
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
