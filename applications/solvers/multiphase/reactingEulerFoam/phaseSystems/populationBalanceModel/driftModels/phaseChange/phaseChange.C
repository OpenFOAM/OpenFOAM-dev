/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "phaseChange.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseSystem.H"
#include "phasePairKey.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace driftModels
{
    defineTypeNameAndDebug(phaseChange, 0);
    addToRunTimeSelectionTable(driftModel, phaseChange, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::driftModels::phaseChange::phaseChange
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    driftModel(popBal, dict),
    pairNames_(dict.lookup("pairNames")),
    iDmdt_
    (
        IOobject
        (
            "iDmdt",
            popBal.time().timeName(),
            popBal.mesh()
        ),
        popBal.mesh(),
        dimensionedScalar("Sui", dimDensity/dimTime, Zero)
    ),
    N_
    (
        IOobject
        (
            "N",
            popBal.mesh().time().timeName(),
            popBal.mesh()
        ),
        popBal.mesh(),
        dimensionedScalar("Sui", inv(dimVolume), Zero)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::diameterModels::driftModels::phaseChange::correct()
{
    iDmdt_ *= 0.0;

    forAll(pairNames_, i)
    {
        const word& pairName = pairNames_[i];

        iDmdt_ +=
            popBal_.mesh().lookupObject<volScalarField>
            (
                IOobject::groupName("iDmdt", pairName)
            );
    }

    N_ *= 0.0;

    forAll(popBal_.sizeGroups(), i)
    {
        const sizeGroup& fi = *popBal_.sizeGroups()[i];

        N_ += fi*max(fi.phase(), small)/fi.x();
    }
}


void Foam::diameterModels::driftModels::phaseChange::addToDriftRate
(
    volScalarField& driftRate,
    const label i
)
{
    const sizeGroup& fi = *popBal_.sizeGroups()[i];

    driftRate += iDmdt_/(N_*fi.phase().rho());
}


// ************************************************************************* //
