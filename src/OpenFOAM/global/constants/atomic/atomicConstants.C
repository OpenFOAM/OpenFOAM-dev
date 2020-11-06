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

const char* const atomic::group = "atomic";


// Note: cannot use dimless etc. as they may not have been constructed yet

const Foam::dimensionedScalar atomic::alpha
(
    dimensionedConstant
    (
        atomic::group,
        "alpha",
        sqr(electromagnetic::e)
       /(2*electromagnetic::epsilon0 *universal::h *universal::c)
    )
);


const Foam::dimensionedScalar atomic::Rinf
(
    dimensionedConstant
    (
        atomic::group,
        "Rinf",
        sqr(atomic::alpha)*atomic::me*universal::c/(2*universal::h)
    )
);


const Foam::dimensionedScalar atomic::a0
(
    dimensionedConstant
    (
        atomic::group,
        "a0",
        atomic::alpha/(4*mathematical::pi*atomic::Rinf)
    )
);


const Foam::dimensionedScalar atomic::re
(
    dimensionedConstant
    (
        atomic::group,
        "re",
        sqr(electromagnetic::e)
       /(
            4*mathematical::pi
           *electromagnetic::epsilon0
           *atomic::me
           *sqr(universal::c)
        )
    )
);


const Foam::dimensionedScalar atomic::Eh
(
    dimensionedConstant
    (
        atomic::group,
        "Eh",
        2*atomic::Rinf*universal::h*universal::c
    )
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace constant
} // End namespace Foam

// ************************************************************************* //
