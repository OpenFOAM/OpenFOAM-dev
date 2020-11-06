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
#include "atomicConstants.H"
#include "dimensionedConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constant
{

const char* const electromagnetic::group = "electromagnetic";


// Note: cannot use dimless etc. as they may not have been constructed yet

const Foam::dimensionedScalar electromagnetic::mu0
(
    dimensionedConstant
    (
        electromagnetic::group,
        "mu0",
        dimensionedScalar
        (
            dimensionSet(1, 1, -2, 0, 0, -2, 0),
            4*mathematical::pi*1e-07
        )
    )
);


const Foam::dimensionedScalar electromagnetic::epsilon0
(
    dimensionedConstant
    (
        electromagnetic::group,
        "epsilon0",
        1/(electromagnetic::mu0*sqr(universal::c))
    )
);


const Foam::dimensionedScalar electromagnetic::Z0
(
    dimensionedConstant
    (
        electromagnetic::group,
        "Z0",
        electromagnetic::mu0*universal::c
    )
);


const Foam::dimensionedScalar electromagnetic::kappa
(
    dimensionedConstant
    (
        electromagnetic::group,
        "kappa",
        (1/(4*mathematical::pi))/electromagnetic::epsilon0
    )
);


const Foam::dimensionedScalar electromagnetic::G0
(
    dimensionedConstant
    (
        electromagnetic::group,
        "G0",
        2*sqr(electromagnetic::e)/universal::h
    )
);


const Foam::dimensionedScalar electromagnetic::KJ
(
    dimensionedConstant
    (
        electromagnetic::group,
        "KJ",
        2*electromagnetic::e/universal::h
    )
);


const Foam::dimensionedScalar electromagnetic::phi0
(
    dimensionedConstant
    (
        electromagnetic::group,
        "phi0",
        universal::h/(2*electromagnetic::e)
    )
);


const Foam::dimensionedScalar electromagnetic::RK
(
    dimensionedConstant
    (
        electromagnetic::group,
        "RK",
        universal::h/sqr(electromagnetic::e)
    )
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace constant
} // End namespace Foam

// ************************************************************************* //
