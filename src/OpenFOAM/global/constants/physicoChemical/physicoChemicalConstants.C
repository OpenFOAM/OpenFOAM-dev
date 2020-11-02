/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "mathematicalConstants.H"
#include "universalConstants.H"
#include "electromagneticConstants.H"
#include "physicoChemicalConstants.H"

#include "dimensionedConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constant
{

const char* const physicoChemical::group = "physicoChemical";

defineDimensionedConstantWithDefault
(
    physicoChemical,
    R,
    dimensionedScalar(physicoChemical::NA*physicoChemical::k)
);


// Note: the 1e3 converts from /mol to /kmol for consistency with the
// SI choice of kg rather than g for mass.
// This is not appropriate for USCS and will be changed to an entry in
// the DimensionedConstants dictionary in etc/controlDict
defineDimensionedConstantWithDefault
(
    physicoChemical,
    RR,
    dimensionedScalar(1000*physicoChemical::R)
);


defineDimensionedConstantWithDefault
(
    physicoChemical,
    F,
    dimensionedScalar(physicoChemical::NA*electromagnetic::e)
);


// Note: cannot use dimless etc. since not guaranteed to be constructed
defineDimensionedConstantWithDefault
(
    physicoChemical,
    sigma,
    dimensionedScalar
    (
        dimensionedScalar
        (
            "C",
            dimensionSet(0, 0, 0, 0, 0),
            sqr(mathematical::pi)/60.0
        )
       *pow4(physicoChemical::k)
       /(pow3(universal::hr)*sqr(universal::c))
    )
);


defineDimensionedConstantWithDefault
(
    physicoChemical,
    b,
    dimensionedScalar
    (
        (universal::h*universal::c/physicoChemical::k)
       /dimensionedScalar
        (
            "C",
            dimensionSet(0, 0, 0, 0, 0),
            4.965114231
        )
    )
);


defineDimensionedConstantWithDefault
(
    physicoChemical,
    c1,
    dimensionedScalar
    (
        dimensionedScalar
        (
            "C",
            dimensionSet(0, 0, 0, 0, 0),
            mathematical::twoPi
        )
       *universal::h*sqr(universal::c)
    )
);


defineDimensionedConstantWithDefault
(
    physicoChemical,
    c2,
    dimensionedScalar(universal::h*universal::c/physicoChemical::k)
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace constant
} // End namespace Foam

// ************************************************************************* //
