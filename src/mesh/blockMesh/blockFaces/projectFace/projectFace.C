/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "projectFace.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blockFaces
{
    defineTypeNameAndDebug(projectFace, 0);
    addToRunTimeSelectionTable(blockFace, projectFace, Istream);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::searchableSurface& Foam::blockFaces::projectFace::lookupSurface
(
    const searchableSurfaces& geometry,
    Istream& is
) const
{
    word name(is);

    forAll(geometry, i)
    {
        if (geometry[i].name() == name)
        {
            return geometry[i];
        }
    }

    FatalIOErrorInFunction(is)
        << "Cannot find surface " << name << " in geometry"
        << exit(FatalIOError);

    return geometry[0];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockFaces::projectFace::projectFace
(
    const searchableSurfaces& geometry,
    Istream& is
)
:
    blockFace(is),
    surface_(lookupSurface(geometry, is))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blockFaces::projectFace::project(pointField& points) const
{
    List<pointIndexHit> hits;
    scalarField nearestDistSqr
    (
        points.size(),
        magSqr(points[0] - points[points.size()-1])
    );
    surface_.findNearest(points, nearestDistSqr, hits);

    forAll(hits, i)
    {
        if (hits[i].hit())
        {
            points[i] = hits[i].hitPoint();
        }
    }
}


// ************************************************************************* //
