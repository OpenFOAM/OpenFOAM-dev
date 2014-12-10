/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Vector of floats.

\*---------------------------------------------------------------------------*/

#include "floatVector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const floatVector::typeName = "floatVector";

template<>
const char* floatVector::componentNames[] = {"x", "y", "z"};

template<>
const floatVector floatVector::zero(0, 0, 0);

template<>
const floatVector floatVector::one(1, 1, 1);

template<>
const floatVector floatVector::max
(
    floatScalarVGREAT,
    floatScalarVGREAT,
    floatScalarVGREAT
);

template<>
const floatVector floatVector::min
(
    -floatScalarVGREAT,
    -floatScalarVGREAT,
    -floatScalarVGREAT
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
