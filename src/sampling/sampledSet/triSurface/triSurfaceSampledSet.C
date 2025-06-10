/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "triSurfaceSampledSet.H"
#include "meshSearch.H"
#include "polyMesh.H"
#include "searchableTriSurface.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace sampledSets
    {
        defineTypeNameAndDebug(triSurface, 0);

        addToRunTimeSelectionTable
        (
            sampledSet,
            triSurface,
            word
        );

        addBackwardCompatibleToRunTimeSelectionTable
        (
            sampledSet,
            triSurface,
            word,
            triSurface,
            "triSurface"
        );
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledSets::triSurface::calcSamples
(
    DynamicList<point>& samplingPositions,
    DynamicList<label>& samplingSegments,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces
) const
{
    forAll(points_, i)
    {
        const point& pt = points_[i];
        const label celli = searchEngine().findCell(pt);

        if (celli != -1)
        {
            samplingPositions.append(pt);
            samplingSegments.append(i);
            samplingCells.append(celli);
            samplingFaces.append(-1);
        }
    }
}


void Foam::sampledSets::triSurface::genSamples()
{
    DynamicList<point> samplingPositions;
    DynamicList<label> samplingSegments;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;

    calcSamples
    (
        samplingPositions,
        samplingSegments,
        samplingCells,
        samplingFaces
    );

    samplingPositions.shrink();
    samplingSegments.shrink();
    samplingCells.shrink();
    samplingFaces.shrink();

    setSamples
    (
        samplingPositions,
        samplingSegments,
        samplingCells,
        samplingFaces
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSets::triSurface::triSurface
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    surface_(dict.lookup("surface")),
    points_
    (
        mesh.time().foundObject<searchableSurfaces::triSurface>
        (
            surface_
        )
      ? mesh.time().lookupObject<searchableSurfaces::triSurface>
        (
            surface_
        ).points()
      : searchableSurfaces::triSurface
        (
            IOobject
            (
                surface_,
                mesh.time().constant(),
                searchableSurface::geometryDir(mesh.time()),
                mesh.time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).points()
    )
{
    genSamples();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSets::triSurface::~triSurface()
{}


// ************************************************************************* //
