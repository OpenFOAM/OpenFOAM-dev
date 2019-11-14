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

const char* const electromagnetic::group = "electromagnetic";


defineDimensionedConstantWithDefault
(
    electromagnetic,
    mu0,
    dimensionedScalar
    (
        dimensionSet(1, 1, -2, 0, 0, -2, 0),
        4*mathematical::pi*1e-07
    )
);


// Note: cannot use dimless etc. since not guaranteed to be constructed
defineDimensionedConstantWithDefault
(
    electromagnetic,
    epsilon0,
    dimensionedScalar
    (
        dimensionedScalar
        (
            "C",
            dimensionSet(0, 0, 0, 0, 0),
            1
        )
       /(electromagnetic::mu0*sqr(universal::c))
    )
);


defineDimensionedConstantWithDefault
(
    electromagnetic,
    Z0,
    dimensionedScalar(electromagnetic::mu0*universal::c)
);


defineDimensionedConstantWithDefault
(
    electromagnetic,
    kappa,
    dimensionedScalar
    (
        dimensionedScalar
        (
            "C",
            dimensionSet(0, 0, 0, 0, 0),
            1/(4*mathematical::pi)
        )
       /electromagnetic::epsilon0
    )
);


defineDimensionedConstantWithDefault
(
    electromagnetic,
    G0,
    dimensionedScalar
    (
        dimensionedScalar("C", dimensionSet(0, 0, 0, 0, 0), 2)
       *sqr(electromagnetic::e)
       /universal::h
    )
);


defineDimensionedConstantWithDefault
(
    electromagnetic,
    KJ,
    dimensionedScalar
    (
        dimensionedScalar("C", dimensionSet(0, 0, 0, 0, 0), 2)
       *electromagnetic::e
       /universal::h
    )
);


defineDimensionedConstantWithDefault
(
    electromagnetic,
    phi0,
    dimensionedScalar
    (
        universal::h
       /(
            dimensionedScalar("C", dimensionSet(0, 0, 0, 0, 0), 2)
           *electromagnetic::e
        )
    )
);


defineDimensionedConstantWithDefault
(
    electromagnetic,
    RK,
    dimensionedScalar(universal::h/sqr(electromagnetic::e))
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace constant
} // End namespace Foam

// ************************************************************************* //
