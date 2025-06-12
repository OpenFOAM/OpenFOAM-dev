/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "pistonBowlPoints.H"
#include "multiValveEngine.H"
#include "pistonPointEdgeData.H"
#include "PointEdgeWave.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(pistonBowlPoints, 0);
        addToRunTimeSelectionTable
        (
            zoneGenerator,
            pistonBowlPoints,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::pistonBowlPoints::pistonBowlPoints
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    zoneGenerator(name, mesh, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::pistonBowlPoints::~pistonBowlPoints()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::pistonBowlPoints::generate() const
{
    const fvMesh& mesh = refCast<const fvMesh>(mesh_);

    const fvMeshMovers::multiValveEngine& mve
    (
        refCast<const fvMeshMovers::multiValveEngine>(mesh.mover())
    );

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
    labelList pointIndices(mesh.nPoints());
    label zpi = 0;
    forAll(allPointData, pointi)
    {
        if (allPointData[pointi].inBowl())
        {
            pointIndices[zpi++] = pointi;
        }
    }
    pointIndices.setSize(zpi);

    return zoneSet
    (
        new pointZone
        (
            dict_.found("name")
              ? zoneName_
              : piston.pistonBowlName,
            pointIndices,
            mesh_.pointZones(),
            moveUpdate_,
            true
        )
    );
}


// ************************************************************************* //
