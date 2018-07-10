/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
#include "atomicConstants.H"

#include "dimensionedConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constant
{

const char* const atomic::group = "atomic";


// Note: cannot use dimless etc. since not guaranteed to be constructed
defineDimensionedConstantWithDefault
(
    atomic::group,
    atomic::alpha,
    dimensionedScalar
    (
        "alpha",
        sqr(electromagnetic::e)
       /(
            dimensionedScalar("C", dimensionSet(0, 0, 0, 0, 0), 2.0)
           *electromagnetic::epsilon0
           *universal::h
           *universal::c
        )
    ),
    constantatomicalpha,
    "alpha"
);


defineDimensionedConstantWithDefault
(
    atomic::group,
    atomic::Rinf,
    dimensionedScalar
    (
        "Rinf",
        sqr(atomic::alpha)
       *atomic::me
       *universal::c
       /(
            Foam::dimensionedScalar
            (
                "C",
                dimensionSet(0, 0, 0, 0, 0),
                2.0
            )
           *universal::h
        )
    ),
    constantatomicRinf,
    "Rinf"
);


defineDimensionedConstantWithDefault
(
    atomic::group,
    atomic::a0,
    dimensionedScalar
    (
        "a0",
        atomic::alpha
       /(
            Foam::dimensionedScalar
            (
                "C",
                dimensionSet(0, 0, 0, 0, 0),
                4.0*mathematical::pi
            )
           *atomic::Rinf
        )
    ),
    constantatomica0,
    "a0"
);


defineDimensionedConstantWithDefault
(
    atomic::group,
    atomic::re,
    dimensionedScalar
    (
        "re",
        Foam::sqr(electromagnetic::e)
       /(
            Foam::dimensionedScalar
            (
                "C",
                dimensionSet(0, 0, 0, 0, 0),
                4.0*mathematical::pi
            )
           *electromagnetic::epsilon0
           *atomic::me
           *Foam::sqr(universal::c)
        )
    ),
    constantatomicre,
    "re"
);


defineDimensionedConstantWithDefault
(
    atomic::group,
    atomic::Eh,
    dimensionedScalar
    (
        "Eh",
        Foam::dimensionedScalar("C", dimensionSet(0, 0, 0, 0, 0), 2.0)
       *atomic::Rinf*universal::h*universal::c
    ),
    constantatomicEh,
    "Eh"
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace constant
} // End namespace Foam

// ************************************************************************* //
