/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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
    SpatialVector of scalars.

\*---------------------------------------------------------------------------*/

#include "spatialVector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::spatialVector::vsType::typeName = "spatialVector";

template<>
const char* const Foam::spatialVector::vsType::componentNames[] =
{
    "wx", "wy", "wz", "lx", "ly", "lz"
};

template<>
const Foam::spatialVector Foam::spatialVector::vsType::zero
(
    Foam::spatialVector::uniform(0)
);

template<>
const Foam::spatialVector Foam::spatialVector::vsType::one
(
    spatialVector::uniform(1)
);

template<>
const Foam::spatialVector Foam::spatialVector::vsType::max
(
    spatialVector::uniform(vGreat)
);

template<>
const Foam::spatialVector Foam::spatialVector::vsType::min
(
    spatialVector::uniform(-vGreat)
);

template<>
const Foam::spatialVector Foam::spatialVector::vsType::rootMax
(
    spatialVector::uniform(rootVGreat)
);

template<>
const Foam::spatialVector Foam::spatialVector::vsType::rootMin
(
    spatialVector::uniform(-rootVGreat)
);


// ************************************************************************* //
