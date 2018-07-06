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

Description
    DiagTensor of scalars.

\*---------------------------------------------------------------------------*/

#include "diagTensor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::diagTensor::vsType::typeName = "diagTensor";

template<>
const char* const Foam::diagTensor::vsType::componentNames[] =
{
    "xx", "yy", "zz"
};

template<>
const Foam::diagTensor Foam::diagTensor::vsType::vsType::zero
(
    diagTensor::uniform(0)
);

template<>
const Foam::diagTensor Foam::diagTensor::vsType::one
(
    diagTensor::uniform(1)
);

template<>
const Foam::diagTensor Foam::diagTensor::vsType::max
(
    diagTensor::uniform(vGreat)
);

template<>
const Foam::diagTensor Foam::diagTensor::vsType::min
(
    diagTensor::uniform(-vGreat)
);

template<>
const Foam::diagTensor Foam::diagTensor::vsType::rootMax
(
    diagTensor::uniform(rootVGreat)
);

template<>
const Foam::diagTensor Foam::diagTensor::vsType::rootMin
(
    diagTensor::uniform(-rootVGreat)
);


// ************************************************************************* //
