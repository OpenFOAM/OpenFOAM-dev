/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "symmTensor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::symmTensor::vsType::typeName = "symmTensor";

template<>
const char* const Foam::symmTensor::vsType::componentNames[] =
{
    "xx", "xy", "xz",
          "yy", "yz",
                "zz"
};

template<>
const Foam::symmTensor Foam::symmTensor::vsType::vsType::zero
(
    symmTensor::uniform(0)
);

template<>
const Foam::symmTensor Foam::symmTensor::vsType::one
(
    symmTensor::uniform(1)
);

template<>
const Foam::symmTensor Foam::symmTensor::vsType::max
(
    symmTensor::uniform(VGREAT)
);

template<>
const Foam::symmTensor Foam::symmTensor::vsType::min
(
    symmTensor::uniform(-VGREAT)
);

template<>
const Foam::symmTensor Foam::symmTensor::vsType::rootMax
(
    symmTensor::uniform(ROOTVGREAT)
);

template<>
const Foam::symmTensor Foam::symmTensor::vsType::rootMin
(
    symmTensor::uniform(-ROOTVGREAT)
);

template<>
const Foam::symmTensor Foam::symmTensor::I
(
    1, 0, 0,
       1, 0,
          1
);


// ************************************************************************* //
