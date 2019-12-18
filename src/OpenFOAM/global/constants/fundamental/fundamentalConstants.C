/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

defineDimensionedConstant(universal, c, dimensionSet(0, 1, -1, 0, 0));
defineDimensionedConstant(universal, G, dimensionSet(-1, 3, -2, 0, 0));
defineDimensionedConstant(universal, h, dimensionSet(1, 2, -1, 0, 0));

// Electromagnetic
defineDimensionedConstant
(
    electromagnetic,
    e,
    dimensionSet(0, 0, 1, 0, 0, 1, 0)
);

// Atomic
defineDimensionedConstant(atomic, me, dimensionSet(1, 0, 0, 0, 0));
defineDimensionedConstant(atomic, mp, dimensionSet(1, 0, 0, 0, 0));

// Physico-chemical
defineDimensionedConstant(physicoChemical, mu, dimensionSet(1, 0, 0, 0, 0));

// Note: cannot use dimless etc since not guaranteed to be constructed
defineDimensionedConstantWithDefault
(
    physicoChemical,
    NA,
    dimensionedScalar
    (
        dimensionSet(0, 0, 0, 0, -1), // dimless/dimMoles,
        6.0221417930e+23
    )
);

defineDimensionedConstant(physicoChemical, k, dimensionSet(1, 2, -2, -1, 0));

// Standard
defineDimensionedConstant(standard, Pstd, dimensionSet(1, -1, -2, 0, 0));
defineDimensionedConstant(standard, Tstd, dimensionSet(0, 0, 0, 1, 0));

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace constant
} // End namespace Foam

// ************************************************************************* //
