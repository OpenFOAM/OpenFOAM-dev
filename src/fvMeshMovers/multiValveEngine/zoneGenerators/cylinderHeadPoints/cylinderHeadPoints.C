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

#include "cylinderHeadPoints.H"
#include "multiValveEngine.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(cylinderHeadPoints, 0);
        addToRunTimeSelectionTable
        (
            zoneGenerator,
            cylinderHeadPoints,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::cylinderHeadPoints::cylinderHeadPoints
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    zoneGenerator(name, mesh, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::cylinderHeadPoints::~cylinderHeadPoints()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::cylinderHeadPoints::generate() const
{
    const fvMesh& mesh = refCast<const fvMesh>(mesh_);

    const fvMeshMovers::multiValveEngine& mve
    (
        refCast<const fvMeshMovers::multiValveEngine>(mesh.mover())
    );

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
    bool foundIntersection = false;
    forAll(cylinderHeadPoints, pi)
    {
        if (cylinderHeadPoints[pi] && linerPoints[pi])
        {
            minZ = min(mve.piston.axis & points[pi], minZ);
            foundIntersection = true;
        }
    }

    // If the intersection not found on this processor
    // set minZ to -great so that the reduction sets minZ
    if (!foundIntersection)
    {
        minZ = -great;
    }

    reduce(minZ, maxOp<scalar>());

    // Create the cylinderHead point set
    labelList pointIndices(points.size());
    label zpi = 0;
    forAll(points, pi)
    {
        if ((mve.piston.axis & points[pi]) > minZ)
        {
            pointIndices[zpi++] = pi;
        }
    }
    pointIndices.setSize(zpi);

    return zoneSet
    (
        new pointZone
        (
            dict_.found("name")
              ? zoneName_
              : mve.cylinderHeadName,
            pointIndices,
            mesh_.pointZones(),
            moveUpdate_,
            true
        )
    );
}


// ************************************************************************* //
