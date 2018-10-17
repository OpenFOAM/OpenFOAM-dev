/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

#include "constantNucleation.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace nucleationModels
{
    defineTypeNameAndDebug(constantNucleation, 0);
    addToRunTimeSelectionTable
    (
        nucleationModel,
        constantNucleation,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::nucleationModels::constantNucleation::
constantNucleation
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    nucleationModel(popBal, dict),
    d_("departureDiameter", dimLength, dict),
    velGroup_
    (
        refCast<const velocityGroup>
        (
            popBal.mesh().lookupObject<phaseModel>
            (
                IOobject::groupName
                (
                    "alpha",
                    dict.lookup("velocityGroup")
                )
            ).dPtr()()
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::nucleationModels::constantNucleation::correct()
{
    if
    (
        d_.value() < velGroup_.sizeGroups().first().d().value()
     || d_.value() > velGroup_.sizeGroups().last().d().value()
    )
    {
        WarningInFunction
            << "Departure diameter " << d_.value() << " m outside of range ["
            << velGroup_.sizeGroups().first().d().value() << ", "
            << velGroup_.sizeGroups().last().d().value() << "] m" << endl
            << "    The nucleation rate is set to zero." << endl
            << "    Adjust discretization over property space to suppress this"
            << " warning."
            << endl;
    }
}


void
Foam::diameterModels::nucleationModels::constantNucleation::addToNucleationRate
(
    volScalarField& nucleationRate,
    const label i
)
{
    const sizeGroup& fi = popBal_.sizeGroups()[i];
    phaseModel& phase = const_cast<phaseModel&>(fi.phase());
    volScalarField& rho = phase.thermoRef().rho();

    nucleationRate +=
        popBal_.gamma(i, velGroup_.formFactor()*pow3(d_))
       *(popBal_.fluid().fvOptions()(phase, rho)&rho)/rho/fi.x();
}


// ************************************************************************* //
