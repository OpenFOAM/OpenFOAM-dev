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
    Vector2D of scalars.

\*---------------------------------------------------------------------------*/

#include "vector2D.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::vector2D::vsType::typeName = "vector2D";

template<>
const char* const Foam::vector2D::vsType::componentNames[] = {"x", "y"};

template<>
const Foam::vector2D Foam::vector2D::vsType::vsType::zero
(
    vector2D::uniform(0)
);

template<>
const Foam::vector2D Foam::vector2D::vsType::one
(
    vector2D::uniform(1)
);

template<>
const Foam::vector2D Foam::vector2D::vsType::max
(
    vector2D::uniform(vGreat)
);

template<>
const Foam::vector2D Foam::vector2D::vsType::min
(
    vector2D::uniform(-vGreat)
);

template<>
const Foam::vector2D Foam::vector2D::vsType::rootMax
(
    vector2D::uniform(rootVGreat)
);

template<>
const Foam::vector2D Foam::vector2D::vsType::rootMin
(
    vector2D::uniform(-rootVGreat)
);

template<>
const Foam::vector2D Foam::vector2D::vsType::nan
(
    vector2D::uniform(NaN)
);


// ************************************************************************* //
