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

#include "labelTensor.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::labelTensor::vsType::typeName = "labelTensor";

template<>
const char* const Foam::labelTensor::vsType::componentNames[] =
{
    "xx", "xy", "xz",
    "yx", "yy", "yz",
    "zx", "zy", "zz"
};

template<>
const Foam::labelTensor Foam::labelTensor::vsType::zero
(
    labelTensor::uniform(0)
);

template<>
const Foam::labelTensor Foam::labelTensor::vsType::one
(
    labelTensor::uniform(1)
);

template<>
const Foam::labelTensor Foam::labelTensor::vsType::max
(
    labelTensor::uniform(labelMax)
);

template<>
const Foam::labelTensor Foam::labelTensor::vsType::min
(
    labelTensor::uniform(-labelMax)
);

template<>
const Foam::labelTensor Foam::labelTensor::vsType::rootMax
(
    labelTensor::uniform(sqrt(scalar(labelMax)))
);

template<>
const Foam::labelTensor Foam::labelTensor::vsType::rootMin
(
    labelTensor::uniform(-sqrt(scalar(labelMax)))
);


// ************************************************************************* //
