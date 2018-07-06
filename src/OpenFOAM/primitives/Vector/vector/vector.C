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
    Vector of scalars.

\*---------------------------------------------------------------------------*/

#include "vector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::vector::vsType::typeName = "vector";

template<>
const char* const Foam::vector::vsType::componentNames[] = {"x", "y", "z"};

template<>
const Foam::vector Foam::vector::vsType::vsType::zero(vector::uniform(0));

template<>
const Foam::vector Foam::vector::vsType::one(vector::uniform(1));

template<>
const Foam::vector Foam::vector::vsType::max(vector::uniform(vGreat));

template<>
const Foam::vector Foam::vector::vsType::min(vector::uniform(-vGreat));

template<>
const Foam::vector Foam::vector::vsType::rootMax(vector::uniform(rootVGreat));

template<>
const Foam::vector Foam::vector::vsType::rootMin(vector::uniform(-rootVGreat));


// ************************************************************************* //
