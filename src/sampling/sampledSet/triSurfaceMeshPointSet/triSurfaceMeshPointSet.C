/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "triSurfaceMeshPointSet.H"
#include "meshSearch.H"
#include "DynamicList.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "triSurfaceMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(triSurfaceMeshPointSet, 0);
    addToRunTimeSelectionTable(sampledSet, triSurfaceMeshPointSet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::triSurfaceMeshPointSet::calcSamples
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
) const
{
    forAll(sampleCoords_, sampleI)
    {
        label cellI = searchEngine().findCell(sampleCoords_[sampleI]);

        if (cellI != -1)
        {
            samplingPts.append(sampleCoords_[sampleI]);
            samplingCells.append(cellI);
            samplingFaces.append(-1);
            samplingSegments.append(0);
            samplingCurveDist.append(1.0 * sampleI);
        }
    }
}


void Foam::triSurfaceMeshPointSet::genSamples()
{
    // Storage for sample points
    DynamicList<point> samplingPts;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;
    DynamicList<label> samplingSegments;
    DynamicList<scalar> samplingCurveDist;

    calcSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );

    samplingPts.shrink();
    samplingCells.shrink();
    samplingFaces.shrink();
    samplingSegments.shrink();
    samplingCurveDist.shrink();

    setSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurfaceMeshPointSet::triSurfaceMeshPointSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    surface_(dict.lookup("surface"))
{
    // Load surface.
    if (mesh.time().foundObject<triSurfaceMesh>(surface_))
    {
        // Note: should use localPoints() instead of points() but assume
        // trisurface is compact.
        sampleCoords_ = mesh.time().lookupObject<triSurfaceMesh>
        (
            surface_
        ).points();
    }
    else
    {
        sampleCoords_ = triSurfaceMesh
        (
            IOobject
            (
                surface_,
                mesh.time().constant(),     // instance
                "triSurface",               // local
                mesh.time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).points();
    }

    genSamples();

    if (debug)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::triSurfaceMeshPointSet::~triSurfaceMeshPointSet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::triSurfaceMeshPointSet::getRefPoint(const List<point>& pts)
 const
{
    if (pts.size())
    {
        // Use first samplePt as starting point
        return pts[0];
    }
    else
    {
        return vector::zero;
    }
}


// ************************************************************************* //
