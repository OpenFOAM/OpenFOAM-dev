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

#include "labelSymmTensor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::labelSymmTensor::vsType::typeName = "labelSymmTensor";

template<>
const char* const Foam::labelSymmTensor::vsType::componentNames[] =
{
    "xx", "xy", "xz",
          "yy", "yz",
                "zz"
};

template<>
const Foam::labelSymmTensor Foam::labelSymmTensor::vsType::vsType::zero
(
    labelSymmTensor::uniform(0)
);

template<>
const Foam::labelSymmTensor Foam::labelSymmTensor::vsType::one
(
    labelSymmTensor::uniform(1)
);

template<>
const Foam::labelSymmTensor Foam::labelSymmTensor::vsType::max
(
    labelSymmTensor::uniform(labelMax)
);

template<>
const Foam::labelSymmTensor Foam::labelSymmTensor::vsType::min
(
    labelSymmTensor::uniform(-labelMax)
);

template<>
const Foam::labelSymmTensor Foam::labelSymmTensor::vsType::rootMax
(
    labelSymmTensor::uniform(sqrt(scalar(labelMax)))
);

template<>
const Foam::labelSymmTensor Foam::labelSymmTensor::vsType::rootMin
(
    labelSymmTensor::uniform(-sqrt(scalar(labelMax)))
);

template<>
const Foam::labelSymmTensor Foam::labelSymmTensor::I
(
    1, 0, 0,
       1, 0,
          1
);


// ************************************************************************* //
