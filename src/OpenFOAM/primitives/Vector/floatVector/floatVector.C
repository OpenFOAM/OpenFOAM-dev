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
    Vector of floats.

\*---------------------------------------------------------------------------*/

#include "floatVector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::floatVector::vsType::typeName = "floatVector";

template<>
const char* const Foam::floatVector::vsType::componentNames[] =
{
    "x", "y", "z"
};

template<>
const Foam::floatVector Foam::floatVector::vsType::zero
(
    floatVector::uniform(0)
);

template<>
const Foam::floatVector Foam::floatVector::vsType::one
(
    floatVector::uniform(1)
);

template<>
const Foam::floatVector Foam::floatVector::vsType::max
(
    floatVector::uniform(floatScalarVGreat)
);

template<>
const Foam::floatVector Foam::floatVector::vsType::min
(
    floatVector::uniform(-floatScalarVGreat)
);

template<>
const Foam::floatVector Foam::floatVector::vsType::rootMax
(
    floatVector::uniform(floatScalarRootVGreat)
);

template<>
const Foam::floatVector Foam::floatVector::vsType::rootMin
(
    floatVector::uniform(-floatScalarRootVGreat)
);

template<>
const Foam::floatVector Foam::floatVector::vsType::nan
(
    floatVector::uniform(-floatScalarNaN)
);


// ************************************************************************* //
