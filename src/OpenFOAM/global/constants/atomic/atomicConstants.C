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
    atomic,
    alpha,
    dimensionedScalar
    (
        sqr(electromagnetic::e)
       /(
            dimensionedScalar("C", dimensionSet(0, 0, 0, 0, 0), 2.0)
           *electromagnetic::epsilon0
           *universal::h
           *universal::c
        )
    )
);


defineDimensionedConstantWithDefault
(
    atomic,
    Rinf,
    dimensionedScalar
    (
        sqr(atomic::alpha)
       *atomic::me
       *universal::c
       /(
            dimensionedScalar
            (
                "C",
                dimensionSet(0, 0, 0, 0, 0),
                2.0
            )
           *universal::h
        )
    )
);


defineDimensionedConstantWithDefault
(
    atomic,
    a0,
    dimensionedScalar
    (
        atomic::alpha
       /(
            dimensionedScalar
            (
                "C",
                dimensionSet(0, 0, 0, 0, 0),
                4.0*mathematical::pi
            )
           *atomic::Rinf
        )
    )
);


defineDimensionedConstantWithDefault
(
    atomic,
    re,
    dimensionedScalar
    (
        sqr(electromagnetic::e)
       /(
            dimensionedScalar
            (
                "C",
                dimensionSet(0, 0, 0, 0, 0),
                4.0*mathematical::pi
            )
           *electromagnetic::epsilon0
           *atomic::me
           *sqr(universal::c)
        )
    )
);


defineDimensionedConstantWithDefault
(
    atomic,
    Eh,
    dimensionedScalar
    (
        dimensionedScalar("C", dimensionSet(0, 0, 0, 0, 0), 2.0)
       *atomic::Rinf*universal::h*universal::c
    )
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace constant
} // End namespace Foam

// ************************************************************************* //
