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

\*---------------------------------------------------------------------------*/

#include "compactSpatialTensor.H"
#include "CompactSpatialTensorT.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::compactSpatialTensor::vsType::typeName =
    "compactSpatialTensor";

template<>
const char* const Foam::compactSpatialTensor::vsType::componentNames[] =
{
    "Exx",  "Exy",  "Exz",
    "Eyx",  "Eyy",  "Eyz",
    "Ezx",  "Ezy",  "Ezz",

    "Erxx", "Erxy", "Erxz",
    "Eryx", "Eryy", "Eryz",
    "Erzx", "Erzy", "Erzz",
};

template<>
const Foam::compactSpatialTensor Foam::compactSpatialTensor::vsType::zero
(
    Foam::compactSpatialTensor::uniform(0)
);

template<>
const Foam::compactSpatialTensor Foam::compactSpatialTensor::vsType::one
(
    compactSpatialTensor::uniform(1)
);

template<>
const Foam::compactSpatialTensor Foam::compactSpatialTensor::vsType::max
(
    compactSpatialTensor::uniform(vGreat)
);

template<>
const Foam::compactSpatialTensor Foam::compactSpatialTensor::vsType::min
(
    compactSpatialTensor::uniform(-vGreat)
);

template<>
const Foam::compactSpatialTensor Foam::compactSpatialTensor::vsType::rootMax
(
    compactSpatialTensor::uniform(rootVGreat)
);

template<>
const Foam::compactSpatialTensor Foam::compactSpatialTensor::vsType::rootMin
(
    compactSpatialTensor::uniform(-rootVGreat)
);


// ************************************************************************* //
