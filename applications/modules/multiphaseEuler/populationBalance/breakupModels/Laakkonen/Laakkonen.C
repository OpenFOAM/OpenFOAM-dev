/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2024 OpenFOAM Foundation
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

#include "Laakkonen.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace breakupModels
{
    defineTypeNameAndDebug(Laakkonen, 0);
    addToRunTimeSelectionTable
    (
        breakupModel,
        Laakkonen,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::breakupModels::Laakkonen::
Laakkonen
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    breakupModel(popBal, dict),
    C1_("C1", dimensionSet(0, -2.0/3.0, 0, 0, 0), dict, 2.25),
    C2_("C2", dimless, dict, 0.04),
    C3_("C3", dimless, dict, 0.01)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::diameterModels::breakupModels::Laakkonen::setBreakupRate
(
    volScalarField& breakupRate,
    const label i
)
{
    const phaseModel& continuousPhase = popBal_.continuousPhase();
    const sizeGroup& fi = popBal_.sizeGroups()[i];

    breakupRate =
        C1_*cbrt(popBal_.continuousTurbulence().epsilon())
       *erfc
        (
            sqrt
            (
                C2_*popBal_.sigmaWithContinuousPhase(fi.phase())
               /(
                    continuousPhase.rho()*pow(fi.dSph(), 5.0/3.0)
                   *pow(popBal_.continuousTurbulence().epsilon(), 2.0/3.0)
                )
              + C3_*continuousPhase.fluidThermo().mu()
               /(
                    sqrt(continuousPhase.rho()*fi.phase().rho())
                   *cbrt(popBal_.continuousTurbulence().epsilon())
                   *pow(fi.dSph(), 4.0/3.0)
                )
            )
        );
}


// ************************************************************************* //
