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

\*---------------------------------------------------------------------------*/

#include "sphericalTensor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::sphericalTensor::vsType::typeName = "sphericalTensor";

template<>
const char* const Foam::sphericalTensor::vsType::componentNames[] = {"ii"};

template<>
const Foam::sphericalTensor Foam::sphericalTensor::vsType::zero
(
    sphericalTensor::uniform(0)
);

template<>
const Foam::sphericalTensor Foam::sphericalTensor::vsType::one
(
    sphericalTensor::uniform(1)
);

template<>
const Foam::sphericalTensor Foam::sphericalTensor::vsType::max
(
    sphericalTensor::uniform(vGreat)
);

template<>
const Foam::sphericalTensor Foam::sphericalTensor::vsType::min
(
    sphericalTensor::uniform(-vGreat)
);

template<>
const Foam::sphericalTensor Foam::sphericalTensor::vsType::rootMax
(
    sphericalTensor::uniform(rootVGreat)
);

template<>
const Foam::sphericalTensor Foam::sphericalTensor::vsType::rootMin
(
    sphericalTensor::uniform(-rootVGreat)
);

template<>
const Foam::sphericalTensor Foam::sphericalTensor::vsType::nan
(
    sphericalTensor::uniform(NaN)
);

template<>
const Foam::sphericalTensor Foam::sphericalTensor::I(1);

template<>
const Foam::sphericalTensor Foam::sphericalTensor::oneThirdI(1.0/3.0);

template<>
const Foam::sphericalTensor Foam::sphericalTensor::twoThirdsI(2.0/3.0);

// ************************************************************************* //
