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

#include "LaakkonenAlopaeusAittamaa.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseCompressibleTurbulenceModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace breakupModels
{
    defineTypeNameAndDebug(LaakkonenAlopaeusAittamaa, 0);
    addToRunTimeSelectionTable
    (
        breakupModel,
        LaakkonenAlopaeusAittamaa,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::breakupModels::LaakkonenAlopaeusAittamaa::
LaakkonenAlopaeusAittamaa
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    breakupModel(popBal, dict),
    C1_
    (
        "C1",
        dimensionSet(0, -2.0/3.0, 0, 0, 0),
        dict.lookupOrDefault<scalar>("C1", 6.0)
    ),
    C2_("C2", dimless, dict.lookupOrDefault<scalar>("C2", 0.04)),
    C3_("C3", dimless, dict.lookupOrDefault<scalar>("C3", 0.01))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::diameterModels::breakupModels::LaakkonenAlopaeusAittamaa::setBreakupRate
(
    volScalarField& breakupRate,
    const label i
)
{
    const phaseModel& continuousPhase = popBal_.continuousPhase();
    const sizeGroup& fi = *popBal_.sizeGroups()[i];
    const volScalarField sigma
    (
        popBal_.fluid().sigma(phasePair(fi.phase(), continuousPhase))
    );

    breakupRate =
        C1_*cbrt(continuousTurbulence().epsilon())
       *erfc
        (
            sqrt
            (
                C2_*sigma
               /(
                    continuousPhase.rho()*pow(fi.d(), 5.0/3.0)
                   *pow(continuousTurbulence().epsilon(), 2.0/3.0)
                )
              + C3_*continuousPhase.mu()
               /(
                    sqrt(continuousPhase.rho()*fi.phase().rho())
                   *cbrt(continuousTurbulence().epsilon())*pow(fi.d(), 4.0/3.0)
                )
            )
        );
}


// ************************************************************************* //
