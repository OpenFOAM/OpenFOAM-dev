/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "MeshedSurfaceIOAllocator.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MeshedSurfaceIOAllocator::MeshedSurfaceIOAllocator
(
    const IOobject& ioPoints,
    const IOobject& ioFaces,
    const IOobject& ioZones
)
:
    points_(ioPoints),
    faces_(ioFaces),
    zones_(ioZones)
{}


Foam::MeshedSurfaceIOAllocator::MeshedSurfaceIOAllocator
(
    const IOobject& ioPoints,
    const pointField& points,
    const IOobject& ioFaces,
    const faceList& faces,
    const IOobject& ioZones,
    const surfZoneList& zones
)
:
    points_(ioPoints, points),
    faces_(ioFaces, faces),
    zones_(ioZones, zones)
{}


Foam::MeshedSurfaceIOAllocator::MeshedSurfaceIOAllocator
(
    const IOobject& ioPoints,
    pointField&& points,
    const IOobject& ioFaces,
    faceList&& faces,
    const IOobject& ioZones,
    surfZoneList&& zones
)
:
    points_(ioPoints, move(points)),
    faces_(ioFaces, move(faces)),
    zones_(ioZones, move(zones))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::MeshedSurfaceIOAllocator::clear()
{
    points_.clear();
    faces_.clear();
    zones_.clear();
}


void Foam::MeshedSurfaceIOAllocator::resetFaces
(
    List<face>&& faces,
    surfZoneList&& zones
)
{
    if (notNull(faces))
    {
        faces_.transfer(faces);
    }

    if (notNull(zones))
    {
        zones_.transfer(zones);
    }
}


void Foam::MeshedSurfaceIOAllocator::reset
(
    pointField&& points,
    faceList&& faces,
    surfZoneList&& zones
)
{
    // Take over new primitive data.
    // Optimised to avoid overwriting data at all
    if (notNull(points))
    {
        points_.transfer(points);
    }

    resetFaces(move(faces), move(zones));
}


void Foam::MeshedSurfaceIOAllocator::reset
(
    List<point>&& points,
    faceList&& faces,
    surfZoneList&& zones
)
{
    // Take over new primitive data.
    // Optimised to avoid overwriting data at all
    if (notNull(points))
    {
        points_.transfer(points);
    }

    resetFaces(move(faces), move(zones));
}


// ************************************************************************* //
