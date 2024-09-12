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

Application
    createEngineZones

Description
    Utility to create the cylinder head and piston bowl pointZones

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "multiValveEngine.H"
#include "pistonPointEdgeData.H"
#include "PointEdgeWave.H"

using namespace Foam;

void calcCylinderHeadPoints(const fvMeshMovers::multiValveEngine& mve)
{
    const fvMesh& mesh = mve.mesh();

    const pointField& points = mesh.points();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // Set the cylinderHead point flags
    boolList cylinderHeadPoints(points.size(), false);
    forAllConstIter(labelHashSet, mve.staticPatchSet, iter)
    {
        UIndirectList<bool>
        (
            cylinderHeadPoints,
            pbm[iter.key()].meshPoints()
        ) = true;
    }

    // Set the liner point flags
    boolList linerPoints(points.size(), false);
    forAllConstIter(labelHashSet, mve.linerPatchSet, iter)
    {
        UIndirectList<bool>
        (
            linerPoints,
            pbm[iter.key()].meshPoints()
        ) = true;
    }

    // Find the minimum axis_ coordinate of the cylinder head region
    // from the bottom of the intersections with the liner
    // for which both the cylinderHead and liner point flags are true
    // Assumes the piston moves in the positive axis direction

    scalar minZ = great;
    forAll(cylinderHeadPoints, pi)
    {
        if (cylinderHeadPoints[pi] && linerPoints[pi])
        {
            minZ = min(mve.piston.axis & points[pi], minZ);
        }
    }
    reduce(minZ, maxOp<scalar>());

    // Create the cylinderHead point set
    labelHashSet cylinderHeadPointSet;
    forAll(points, pi)
    {
        if ((mve.piston.axis & points[pi]) > minZ)
        {
            cylinderHeadPointSet.insert(pi);
        }
    }

    // Add the cylinderHead pointZone
    mesh.pointZones().append
    (
        new pointZone
        (
            mve.cylinderHeadName,
            cylinderHeadPointSet.toc(),
            mesh.pointZones()
        )
    );
}


void calcPistonBowlPoints(const fvMeshMovers::multiValveEngine& mve)
{
    const fvMesh& mesh = mve.mesh();
    const fvMeshMovers::multiValveEngine::pistonObject& piston = mve.piston;

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // Find the maximum axis coordinate of the piston patch-set
    // Assumes the piston moves in the positive axis direction
    scalar maxZ = -great;
    forAllConstIter(labelHashSet, piston.patchSet, iter)
    {
        const label patchi = iter.key();
        if (pbm[patchi].localPoints().size())
        {
            maxZ = max(maxZ, max(piston.axis & pbm[patchi].localPoints()));
        }
    }
    reduce(maxZ, maxOp<scalar>());


    // Search for points starting at the piston surface and stopping at maxZ
    // using a PointEdgeWave

    label nPistonPatchPoints = 0;
    forAllConstIter(labelHashSet, piston.patchSet, iter)
    {
        const label patchi = iter.key();
        nPistonPatchPoints += pbm[patchi].meshPoints().size();
    }


    const pointField& points = mesh.points();
    pistonPointEdgeData::trackingData td(points, piston.axis, maxZ);

    // Set initial changed points to all the patch points(if patch present)
    List<pistonPointEdgeData> pistonPatchPointData(nPistonPatchPoints);
    labelList pistonPatchPoints(nPistonPatchPoints);

    // Add the patch points to the pistonPatchPointData
    nPistonPatchPoints = 0;
    forAllConstIter(labelHashSet, piston.patchSet, iter)
    {
        const label patchi = iter.key();
        const labelList& mp = pbm[patchi].meshPoints();

        forAll(mp, ppi)
        {
            pistonPatchPoints[nPistonPatchPoints] = mp[ppi];
            pistonPatchPointData[nPistonPatchPoints] = pistonPointEdgeData
            (
                true
            );
            nPistonPatchPoints++;
        }
    }

    // Point data for wave
    List<pistonPointEdgeData> allPointData(mesh.nPoints());

    // Edge data for wave
    List<pistonPointEdgeData> allEdgeData(mesh.nEdges());

    PointEdgeWave
    <
        pistonPointEdgeData,
        pistonPointEdgeData::trackingData
    > patchCalc
    (
        mesh,
        pistonPatchPoints,
        pistonPatchPointData,

        allPointData,
        allEdgeData,
        mesh.globalData().nTotalPoints(), // max iterations
        td
    );

    // Create a labelHashSet of the point labels in the piston bowl
    labelHashSet pistonBowlPointSet;
    forAll(allPointData, pointi)
    {
        if (allPointData[pointi].inBowl())
        {
            pistonBowlPointSet.insert(pointi);
        }
    }

    // Convert the labelHashSet to a pointZone and add to the pointZones
    mesh.pointZones().append
    (
        new pointZone
        (
            piston.pistonBowlName,
            pistonBowlPointSet.toc(),
            mesh.pointZones()
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();

    #include "addMeshOption.H"
    #include "addRegionOption.H"

    argList::addBoolOption
    (
        "cylinderHead",
        "Create the cylinderHead pointZone"
    );

    argList::addBoolOption
    (
        "pistonBowl",
        "Create the piston pointZone"
    );

    #include "setRootCase.H"
    #include "createTimeNoFunctionObjects.H"

    const bool cylinderHeadZone = args.optionFound("cylinderHead");
    const bool pistonBowlZone = args.optionFound("pistonBowl");

    if (!cylinderHeadZone && !pistonBowlZone)
    {
        FatalErrorInFunction
            << "Neither cylinderHeadZone nor pistonBowl pointZones requested"
            << exit(FatalError);
    }

    #include "createSpecifiedMesh.H"

    const fvMeshMovers::multiValveEngine& mve
    (
        refCast<const fvMeshMovers::multiValveEngine>(mesh.mover())
    );

    if (cylinderHeadZone)
    {
        calcCylinderHeadPoints(mve);
    }

    if (pistonBowlZone)
    {
        calcPistonBowlPoints(mve);
    }

    // Take over refinement levels and write to new time directory.
    Info<< "Writing pointZones" << endl;
    mesh.pointZones().write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
