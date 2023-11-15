/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2023 OpenFOAM Foundation
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

#include "floatTensor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::floatTensor::vsType::typeName = "floatTensor";

template<>
const char* const Foam::floatTensor::vsType::componentNames[] =
{
    "xx", "xy", "xz",
    "yx", "yy", "yz",
    "zx", "zy", "zz"
};

template<>
const Foam::floatTensor Foam::floatTensor::vsType::zero
(
    floatTensor::uniform(0)
);

template<>
const Foam::floatTensor Foam::floatTensor::vsType::one
(
    floatTensor::uniform(1)
);

template<>
const Foam::floatTensor Foam::floatTensor::vsType::max
(
    floatTensor::uniform(floatScalarVGreat)
);

template<>
const Foam::floatTensor Foam::floatTensor::vsType::min
(
    floatTensor::uniform(-floatScalarVGreat)
);

template<>
const Foam::floatTensor Foam::floatTensor::vsType::rootMax
(
    floatTensor::uniform(floatScalarRootVGreat)
);

template<>
const Foam::floatTensor Foam::floatTensor::vsType::rootMin
(
    floatTensor::uniform(-floatScalarRootVGreat)
);

template<>
const Foam::floatTensor Foam::floatTensor::vsType::nan
(
    floatTensor::uniform(-floatScalarNaN)
);

template<>
const Foam::floatTensor Foam::floatTensor::I
(
    1, 0, 0,
    0, 1, 0,
    0, 0, 1
);


// ************************************************************************* //
