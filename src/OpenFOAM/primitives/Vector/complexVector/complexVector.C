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
    Vector of complex numbers.

\*---------------------------------------------------------------------------*/

#include "complexVector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::complexVector::vsType::typeName = "complexVector";

template<>
const char* const Foam::complexVector::vsType::componentNames[] =
{
    "x", "y", "z"
};

template<>
const Foam::complexVector Foam::complexVector::vsType::zero
(
    complexVector::uniform(complex(0, 0))
);

template<>
const Foam::complexVector Foam::complexVector::vsType::one
(
    complexVector::uniform(complex(1, 1))
);

template<>
const Foam::complexVector Foam::complexVector::vsType::max
(
    complexVector::uniform(complex(vGreat, vGreat))
);

template<>
const Foam::complexVector Foam::complexVector::vsType::min
(
    complexVector::uniform(complex(-vGreat, -vGreat))
);

template<>
const Foam::complexVector Foam::complexVector::vsType::rootMax
(
    complexVector::uniform(complex(rootVGreat, rootVGreat))
);

template<>
const Foam::complexVector Foam::complexVector::vsType::rootMin
(
    complexVector::uniform(complex(-rootVGreat, -rootVGreat))
);

template<>
const Foam::complexVector Foam::complexVector::vsType::nan
(
    complexVector::uniform(complex(NaN, NaN))
);


// ************************************************************************* //
