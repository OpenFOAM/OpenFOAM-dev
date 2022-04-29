/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2022 OpenFOAM Foundation
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
        diameterModels::populationBalanceModel::iNew(*this, dmdtfs_)
    )
{
    forAll(populationBalances_, popBali)
    {
        const diameterModels::populationBalanceModel& popBal =
            populationBalances_[popBali];

        forAllConstIter
        (
            HashTable<const diameterModels::velocityGroup*>,
            popBal.velocityGroupPtrs(),
            iter1
        )
        {
            const diameterModels::velocityGroup& velGrp1 = *iter1();

            forAllConstIter
            (
                HashTable<const diameterModels::velocityGroup*>,
                popBal.velocityGroupPtrs(),
                iter2
            )
            {
                const diameterModels::velocityGroup& velGrp2 = *iter2();

                const phaseInterface interface
                (
                    velGrp1.phase(),
                    velGrp2.phase()
                );

                if
                (
                    &velGrp1 != &velGrp2
                    &&
                    !dmdtfs_.found(interface)
                )
                {
                    this->template validateMassTransfer
                        <diameterModels::populationBalanceModel>(interface);

                    dmdtfs_.insert
                    (
                        interface,
                        new volScalarField
                        (
                            IOobject
                            (
                                IOobject::groupName
                                (
                                    "populationBalance:dmdtf",
                                    interface.name()
                                ),
                                this->mesh().time().timeName(),
                                this->mesh()
                            ),
                            this->mesh(),
                            dimensionedScalar(dimDensity/dimTime, 0)
                        )
                    );
                }
            }
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
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::dmdtf
(
    const phaseInterfaceKey& key
) const
{
    tmp<volScalarField> tDmdtf = BasePhaseSystem::dmdtf(key);

    if (dmdtfs_.found(key))
    {
        tDmdtf.ref() += *dmdtfs_[key];
    }

    return tDmdtf;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::dmdts() const
{
    PtrList<volScalarField> dmdts(BasePhaseSystem::dmdts());

    forAllConstIter(phaseSystem::dmdtfTable, dmdtfs_, dmdtfIter)
    {
        const phaseInterface interface(*this, dmdtfIter.key());

        addField(interface.phase1(), "dmdt", *dmdtfIter(), dmdts);
        addField(interface.phase2(), "dmdt", - *dmdtfIter(), dmdts);
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

    this->addDmdtUfs(dmdtfs_, eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::momentumTransferf()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransferf();

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    this->addDmdtUfs(dmdtfs_, eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::heatTransfer() const
{
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr =
        BasePhaseSystem::heatTransfer();

    phaseSystem::heatTransferTable& eqns = eqnsPtr();

    this->addDmdtHefs(dmdtfs_, eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::specieTransferTable>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::specieTransfer() const
{
    autoPtr<phaseSystem::specieTransferTable> eqnsPtr =
        BasePhaseSystem::specieTransfer();

    phaseSystem::specieTransferTable& eqns = eqnsPtr();

    this->addDmdtYfs(dmdtfs_, eqns);

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
