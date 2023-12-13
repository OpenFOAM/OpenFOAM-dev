/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
#include "standardConstants.H"
#include "dimensionedConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Universal constants

namespace Foam
{
namespace constant
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Note: cannot use dimless etc. as they may not have been constructed yet


const Foam::dimensionedScalar universal::c
(
    dimensionedConstant(universal::group, "c", dimensionSet(0, 1, -1, 0, 0))
);

const Foam::dimensionedScalar universal::G
(
    dimensionedConstant(universal::group, "G", dimensionSet(-1, 3, -2, 0, 0))
);

const Foam::dimensionedScalar universal::h
(
    dimensionedConstant(universal::group, "h", dimensionSet(1, 2, -1, 0, 0))
);


// Electromagnetic

const Foam::dimensionedScalar electromagnetic::e
(
    dimensionedConstant
    (
        electromagnetic::group,
        "e",
        dimensionSet(0, 0, 1, 0, 0, 1, 0)
    )
);


// Atomic

const Foam::dimensionedScalar atomic::me
(
    dimensionedConstant(atomic::group, "me", dimensionSet(1, 0, 0, 0, 0))
);

const Foam::dimensionedScalar atomic::mp
(
    dimensionedConstant(atomic::group, "mp", dimensionSet(1, 0, 0, 0, 0))
);


// Physico-chemical

const Foam::dimensionedScalar physicoChemical::mu
(
    dimensionedConstant
    (
        physicoChemical::group,
        "mu",
        dimensionSet(1, 0, 0, 0, 0)
    )
);

const Foam::dimensionedScalar physicoChemical::NA
(
    dimensionedConstant
    (
        physicoChemical::group,
        "NA",
        dimensionSet(0, 0, 0, 0, -1)
    )
);

const Foam::dimensionedScalar physicoChemical::k
(
    dimensionedConstant
    (
        physicoChemical::group,
        "k",
        dimensionSet(1, 2, -2, -1, 0)
    )
);


// Standard

const Foam::dimensionedScalar standard::Pstd
(
    dimensionedConstant(standard::group, "Pstd", dimensionSet(1, -1, -2, 0, 0))
);

const Foam::dimensionedScalar standard::Tstd
(
    dimensionedConstant(standard::group, "Tstd", dimensionSet(0, 0, 0, 1, 0))
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace constant
} // End namespace Foam

// ************************************************************************* //
