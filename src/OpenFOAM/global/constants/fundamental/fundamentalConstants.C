/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Description
    Fundamental dimensioned constants

\*---------------------------------------------------------------------------*/

#include "fundamentalConstants.H"

#include "universalConstants.H"
#include "electromagneticConstants.H"
#include "atomicConstants.H"
#include "physicoChemicalConstants.H"

#include "dimensionedConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Universal constants

namespace Foam
{
namespace constant
{

defineDimensionedConstant
(
    universal::group,
    universal::c,
    constantuniversalc,
    "c"
);


defineDimensionedConstant
(
    universal::group,
    universal::G,
    constantuniversalG,
    "G"
);


defineDimensionedConstant
(
    universal::group,
    universal::h,
    constantuniversalh,
    "h"
);


// Electromagnetic

defineDimensionedConstant
(
    electromagnetic::group,
    electromagnetic::e,
    constantelectromagnetice,
    "e"
);


// Atomic

defineDimensionedConstant
(
    atomic::group,
    atomic::me,
    constantatomicme,
    "me"
);


defineDimensionedConstant
(
    atomic::group,
    atomic::mp,
    constantatomicmp,
    "mp"
);


// Physico-chemical

defineDimensionedConstant
(
    physicoChemical::group,
    physicoChemical::mu,
    constantphysicoChemicalmu,
    "mu"
);


// Note: cannot use dimless etc since not guaranteed to be constructed
defineDimensionedConstantWithDefault
(
    physicoChemical::group,
    physicoChemical::NA,
    Foam::dimensionedScalar
    (
        "NA",
        dimensionSet(0, 0, 0, 0, -1), //Foam::dimless/Foam::dimMoles,
        6.0221417930e+23
    ),
    constantphysicoChemicalNA,
    "NA"
);


defineDimensionedConstant
(
    physicoChemical::group,
    physicoChemical::k,
    constantphysicoChemicalk,
    "k"
);


// Standard

defineDimensionedConstant
(
    "standard",
    standard::Pstd,
    constantstandardPstd,
    "Pstd"
);


defineDimensionedConstant
(
    "standard",
    standard::Tstd,
    constantstandardTstd,
    "Tstd"
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace constant
} // End namespace Foam

// ************************************************************************* //
