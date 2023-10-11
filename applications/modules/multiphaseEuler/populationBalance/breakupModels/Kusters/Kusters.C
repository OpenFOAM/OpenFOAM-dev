/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "Kusters.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace breakupModels
{
    defineTypeNameAndDebug(Kusters, 0);
    addToRunTimeSelectionTable
    (
        breakupModel,
        Kusters,
        dictionary
    );
}
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::breakupModels::Kusters::Kusters
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    breakupModel(popBal, dict),
    B_("B", dimensionSet(0, 3, -3, 0, 0), dict),
    dP_("dP", dimLength, dict),
    kc_(dimensionedScalar::lookupOrDefault("kc", dict, dimless, 1)),
    Df_("Df", dimless, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::breakupModels::Kusters::setBreakupRate
(
    volScalarField& breakupRate,
    const label i
)
{
    const sizeGroup& fi = popBal_.sizeGroups()[i];
    breakupRate =
        sqrt
        (
            4*popBal_.continuousTurbulence().epsilon()
           /(15*pi*popBal_.continuousPhase().fluidThermo().nu())
        )
       *exp
        (
          - B_/(dP_*0.5*pow(pow(fi.d()/dP_, Df_)/kc_, 1/Df_))
           /popBal_.continuousTurbulence().epsilon()
        );
}


// ************************************************************************* //
