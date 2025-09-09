/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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

#include "generatedFaceZone.H"
#include "polyMesh.H"
#include "containsPoints.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::generatedFaceZone::generatedFaceZone
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::generatedFaceZone::~generatedFaceZone()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::generatedFaceZone::movePoints()
{
    faceZone_.movePoints();
}


void Foam::generatedFaceZone::topoChange(const polyTopoChangeMap& map)
{
    faceZone_.topoChange(map);
}


void Foam::generatedFaceZone::mapMesh(const polyMeshMap& map)
{
    faceZone_.mapMesh(map);
}


void Foam::generatedFaceZone::distribute(const polyDistributionMap& map)
{
    faceZone_.distribute(map);
}


bool Foam::generatedFaceZone::read(const dictionary& dict)
{
    faceZone_.read("faceZone", zoneTypes::face, mesh_, dict);

    return true;
}


// ************************************************************************* //
