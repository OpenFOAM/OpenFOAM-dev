/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2026 OpenFOAM Foundation
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

#include "layeredDisplacement_pointMeshMover.H"
#include "pointEdgeStructuredWalk.H"
#include "PointEdgeWave.H"
#include "pointConstraints.H"
#include "polyTopoChangeMap.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pointMeshMovers
{
    defineTypeNameAndDebug(layeredDisplacement, 0);

    addToRunTimeSelectionTable
    (
        pointMeshMover,
        layeredDisplacement,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pointMeshMovers::layeredDisplacement::walkLayers
(
    const polyPatch& startPatch,
    scalarField& distance,
    vectorField& displacement
) const
{
    // Get the start mesh points from the start patch
    const labelList& startMeshPoints = startPatch.meshPoints();

    // Set the displacement of the start points
    const vectorField startDisplacement
    (
        pointDisplacement_.boundaryField()[startPatch.index()]
       .patchInternalField()
    );

    // Set the wave info for the start patch
    List<pointEdgeStructuredWalk> startInfo(startMeshPoints.size());
    forAll(startMeshPoints, i)
    {
        startInfo[i] = pointEdgeStructuredWalk
        (
            points0()[startMeshPoints[i]],  // Location of start point
            points0()[startMeshPoints[i]],  // Previous location
            0,
            startDisplacement[i]            // Displacement of the start point
        );
    }

    // Initialise the wave info for all the points
    List<pointEdgeStructuredWalk> allPointInfo(poly().nPoints());
    forAll(allPointInfo, pointi)
    {
        allPointInfo[pointi] = pointEdgeStructuredWalk
        (
            points0()[pointi],  // Mesh point location
            vector::max,        // No valid previous location
            0,
            Zero                // Initial displacement = 0
        );
    }

    // Initialise the wave info for all edges
    List<pointEdgeStructuredWalk> allEdgeInfo(poly().nEdges());
    forAll(allEdgeInfo, edgei)
    {
        allEdgeInfo[edgei] = pointEdgeStructuredWalk
        (
            poly().edges()[edgei].centre(points0()), // Edge centre location
            vector::max,        // No valid previous location
            0,
            Zero                // Initial displacement = 0
        );
    }

    // Walk the distance and displacement from the startPatch
    PointEdgeWave<pointEdgeStructuredWalk> walk
    (
        poly(),
        startMeshPoints,
        startInfo,
        allPointInfo,
        allEdgeInfo,
        poly().globalData().nTotalPoints()  // Max number of iterations
    );

    // Extract distance and displacement from the wave info
    forAll(allPointInfo, pointi)
    {
        distance[pointi] = allPointInfo[pointi].dist();
        displacement[pointi] = allPointInfo[pointi].data();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointMeshMovers::layeredDisplacement::layeredDisplacement
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    displacement(mesh, dict, typeName),
    oppositePatchNames_(dict.lookup("oppositePatches")),
    oppositePatches_
    (
        mesh.boundary().findIndex(oppositePatchNames_.first()),
        mesh.boundary().findIndex(oppositePatchNames_.second())
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointMeshMovers::layeredDisplacement::~layeredDisplacement()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::pointMeshMovers::layeredDisplacement::newPoints()
{
    // The points have moved so before interpolation update the pointMeshMover
    movePoints(poly().points());

    // Update the displacement boundary conditions
    pointDisplacement_.boundaryFieldRef().updateCoeffs();

    // Walk the layers from patch0 to patch1
    const polyPatch& patch0 = poly().boundary()[oppositePatches_.first()];
    scalarField patchDist0(poly().nPoints());
    vectorField patchDisp0(pointDisplacement_);
    walkLayers(patch0, patchDist0, patchDisp0);

    // Walk the layers from patch1 to patch0
    const polyPatch& patch1 = poly().boundary()[oppositePatches_.second()];
    scalarField patchDist1(poly().nPoints());
    vectorField patchDisp1(pointDisplacement_);
    walkLayers(patch1, patchDist1, patchDisp1);

    // Calculate the interpolation factor from the distance to the
    // opposite patches
    const scalarField w(patchDist0/(patchDist0 + patchDist1 + small));

    // Linearly interpolate the displacements of the opposite patches
    pointDisplacement_.primitiveFieldRef() = (1 - w)*patchDisp0 + w*patchDisp1;

    // Constrain the pointDisplacement field
    pointConstraints::New(pointDisplacement_.mesh())
        .constrainDisplacement(pointDisplacement_, false);

    return points();
}


void Foam::pointMeshMovers::layeredDisplacement::mapMesh(const polyMeshMap& map)
{
    FatalErrorInFunction
        << "Mesh-to-mesh mapping in not implemented for displacement solvers"
        << nl
        << "    velocity based motion solvers are preferable for cases in which"
           " the mesh is reset periodically avoiding accumulation of error."
        << exit(FatalError);
}


// ************************************************************************* //
