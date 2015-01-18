/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "indexedVertexEnum.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char*
Foam::NamedEnum<Foam::indexedVertexEnum::vertexType, 15>::names[] =
{
    "Unassigned",
    "Internal",
    "InternalNearBoundary",
    "InternalSurface",
    "InternalSurfaceBaffle",
    "ExternalSurfaceBaffle",
    "InternalFeatureEdge",
    "InternalFeatureEdgeBaffle",
    "ExternalFeatureEdgeBaffle",
    "InternalFeaturePoint",
    "ExternalSurface",
    "ExternalFeatureEdge",
    "ExternalFeaturePoint",
    "Far",
    "Constrained"
};

const Foam::NamedEnum<Foam::indexedVertexEnum::vertexType, 15>
Foam::indexedVertexEnum::vertexTypeNames_;


template<>
const char*
Foam::NamedEnum<Foam::indexedVertexEnum::vertexMotion, 2>::names[] =
{
    "fixed",
    "movable"
};

const Foam::NamedEnum<Foam::indexedVertexEnum::vertexMotion, 2>
vertexMotionNames_;


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const Foam::indexedVertexEnum::vertexType& v
)
{
    os  << static_cast<int>(v);

    return os;
}

Foam::Istream& Foam::operator>>
(
    Istream& is,
    Foam::indexedVertexEnum::vertexType& v
)
{
    int type;
    is  >> type;

    v = static_cast<Foam::indexedVertexEnum::vertexType>(type);

    return is;
}

// ************************************************************************* //
