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

#include "triSurfaceMeshSampledSet.H"
#include "meshSearch.H"
#include "DynamicList.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "triSurfaceMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSets
{
    defineTypeNameAndDebug(triSurfaceMeshSampledSet, 0);
    addToRunTimeSelectionTable(sampledSet, triSurfaceMeshSampledSet, word);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledSets::triSurfaceMeshSampledSet::calcSamples
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


void Foam::sampledSets::triSurfaceMeshSampledSet::genSamples()
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

Foam::sampledSets::triSurfaceMeshSampledSet::triSurfaceMeshSampledSet
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
        mesh.time().foundObject<triSurfaceMesh>(surface_)
      ? mesh.time().lookupObject<triSurfaceMesh>(surface_).points()
      : triSurfaceMesh
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

Foam::sampledSets::triSurfaceMeshSampledSet::~triSurfaceMeshSampledSet()
{}


// ************************************************************************* //
