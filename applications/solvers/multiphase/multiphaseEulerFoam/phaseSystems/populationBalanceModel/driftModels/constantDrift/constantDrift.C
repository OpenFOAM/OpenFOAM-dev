/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
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

#include "constantDrift.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace driftModels
{
    defineTypeNameAndDebug(constantDrift, 0);
    addToRunTimeSelectionTable(driftModel, constantDrift, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::driftModels::constantDrift::constantDrift
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    driftModel(popBal, dict),
    N_
    (
        IOobject
        (
            "N",
            popBal.mesh().time().timeName(),
            popBal.mesh()
        ),
        popBal.mesh(),
        dimensionedScalar(inv(dimVolume), Zero)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::driftModels::constantDrift::correct()
{
    N_ = Zero;

    forAll(popBal_.sizeGroups(), i)
    {
        const sizeGroup& fi = popBal_.sizeGroups()[i];

        N_ += fi*fi.phase()/fi.x();
    }
}


void Foam::diameterModels::driftModels::constantDrift::addToDriftRate
(
    volScalarField& driftRate,
    const label i
)
{
    const sizeGroup& fi = popBal_.sizeGroups()[i];
    phaseModel& phase = const_cast<phaseModel&>(fi.phase());
    volScalarField& rho = phase.thermoRef().rho();

    driftRate += (popBal_.fluid().fvOptions()(phase, rho)&rho)/(N_*rho);
}


// ************************************************************************* //
