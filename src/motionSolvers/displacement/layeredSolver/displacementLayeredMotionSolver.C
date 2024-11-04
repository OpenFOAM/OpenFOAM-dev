/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "displacementLayeredMotionSolver.H"
#include "pointEdgeStructuredWalk.H"
#include "PointEdgeWave.H"
#include "pointConstraints.H"
#include "polyTopoChangeMap.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(displacementLayeredMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        displacementLayeredMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::displacementLayeredMotionSolver::walkLayers
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
    List<pointEdgeStructuredWalk> allPointInfo(mesh().nPoints());
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
    List<pointEdgeStructuredWalk> allEdgeInfo(mesh().nEdges());
    forAll(allEdgeInfo, edgei)
    {
        allEdgeInfo[edgei] = pointEdgeStructuredWalk
        (
            mesh().edges()[edgei].centre(points0()), // Edge centre location
            vector::max,        // No valid previous location
            0,
            Zero                // Initial displacement = 0
        );
    }

    // Walk the distance and displacement from the startPatch
    PointEdgeWave<pointEdgeStructuredWalk> walk
    (
        mesh(),
        startMeshPoints,
        startInfo,
        allPointInfo,
        allEdgeInfo,
        mesh().globalData().nTotalPoints()  // Max number of iterations
    );

    // Extract distance and displacement from the wave info
    forAll(allPointInfo, pointi)
    {
        distance[pointi] = allPointInfo[pointi].dist();
        displacement[pointi] = allPointInfo[pointi].data();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementLayeredMotionSolver::displacementLayeredMotionSolver
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    displacementMotionSolver(name, mesh, dict, typeName),
    oppositePatchNames_(dict.lookup("oppositePatches")),
    oppositePatches_
    (
        mesh.boundaryMesh().findIndex(oppositePatchNames_.first()),
        mesh.boundaryMesh().findIndex(oppositePatchNames_.second())
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacementLayeredMotionSolver::~displacementLayeredMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::displacementLayeredMotionSolver::curPoints() const
{
    tmp<pointField> tcurPoints
    (
        points0() + pointDisplacement_.primitiveField()
    );

    return tcurPoints;
}


void Foam::displacementLayeredMotionSolver::solve()
{
    // The points have moved so before interpolation update the motionSolver
    movePoints(mesh().points());

    // Update the displacement boundary conditions
    pointDisplacement_.boundaryFieldRef().updateCoeffs();

    // Walk the layers from patch0 to patch1
    const polyPatch& patch0 = mesh().boundaryMesh()[oppositePatches_.first()];
    scalarField patchDist0(mesh().nPoints());
    vectorField patchDisp0(pointDisplacement_);
    walkLayers(patch0, patchDist0, patchDisp0);

    // Walk the layers from patch1 to patch0
    const polyPatch& patch1 = mesh().boundaryMesh()[oppositePatches_.second()];
    scalarField patchDist1(mesh().nPoints());
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
}


void Foam::displacementLayeredMotionSolver::topoChange
(
    const polyTopoChangeMap& map
)
{
    // Pending implementation of the inverse transformation of points0
    NotImplemented;

    // Alternatively using the current points instead of points0
    // would support topology change but the motion might be less robust
    // and potentially accumulate error
}


// ************************************************************************* //
