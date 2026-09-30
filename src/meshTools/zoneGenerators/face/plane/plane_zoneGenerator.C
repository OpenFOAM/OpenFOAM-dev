/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "plane_zoneGenerator.H"
#include "polyMesh.H"
#include "patchRegionSplit.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(plane, 0);
        addToRunTimeSelectionTable
        (
            zoneGenerator,
            plane,
            dictionary
        );
    }
}

const Foam::NamedEnum<Foam::zoneGenerators::plane::include, 2>
Foam::zoneGenerators::plane::includeNames
{
    "all",
    "closest"
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::plane::plane
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    zoneGenerator(name, mesh, dict),
    point_(dict.lookup<vector>("point", dimensions::length)),
    normal_(dict.lookup<vector>("normal", dimless)),
    include_(includeNames.lookupOrDefault("include", dict, include::all))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::plane::~plane()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::plane::generate() const
{
    // Mark all cells with centres above the plane
    boolList cellIsAbovePlane(mesh_.nCells());
    forAll(mesh_.cells(), celli)
    {
        cellIsAbovePlane[celli] =
            ((mesh_.cellCentres()[celli] - point_) & normal_) > 0;
    }

    // Mark all coupled neighbour cells with centres above the plane
    boolList bFaceNbrCellIsAbovePlane(mesh_.nFaces() - mesh_.nInternalFaces());
    {
        vectorField bFaceNbrCellCentres;
        syncTools::swapBoundaryCellPositions
        (
            mesh_,
            mesh_.cellCentres(),
            bFaceNbrCellCentres
        );
        forAll(bFaceNbrCellIsAbovePlane, bFacei)
        {
            bFaceNbrCellIsAbovePlane[bFacei] =
                ((bFaceNbrCellCentres[bFacei] - point_) & normal_) > 0;
        }
    }

    // Mark all faces that sit between cells above and below the plane
    boolList faceIsOnPlane(mesh_.nFaces(), false);
    forAll(mesh_.faceNeighbour(), facei)
    {
        faceIsOnPlane[facei] =
            cellIsAbovePlane[mesh_.faceOwner()[facei]]
         != cellIsAbovePlane[mesh_.faceNeighbour()[facei]];
    }
    forAll(mesh_.boundary(), patchi)
    {
        const polyPatch& patch = mesh_.boundary()[patchi];

        if (!patch.coupled()) continue;

        forAll(patch, patchFacei)
        {
            const label facei = patch.start() + patchFacei;
            faceIsOnPlane[facei] =
                cellIsAbovePlane[mesh_.faceOwner()[facei]]
             != bFaceNbrCellIsAbovePlane[facei - mesh_.nInternalFaces()];
        }
    }

    // Ensure consistency across couplings
    syncTools::syncFaceList(mesh_, faceIsOnPlane, orEqOp());

    // Convert marked faces to a list of indices
    labelList faceIndices(findIndices(faceIsOnPlane, true));

    // If constructing a single contiguous set, remove all faces except those
    // connected to the contiguous region closest to the specified point
    if (include_ == include::closest)
    {
        // Identify contiguous regions
        const patchRegionSplit prs(mesh_, faceIndices);

        // Find the smallest distance from each region to the point
        scalarField regionMinDistSqr(prs.nRegions(), vGreat);
        forAll(faceIndices, fi)
        {
            const label facei = faceIndices[fi];
            const label regioni = prs[fi];

            const scalar distSqr =
                magSqr
                (
                    mesh_.faces()[facei].nearestPoint
                    (
                        point_,
                        mesh_.points()
                    ).rawPoint()
                  - point_
                );

            regionMinDistSqr[regioni] = min(regionMinDistSqr[regioni], distSqr);
        }
        Pstream::listCombineGather(regionMinDistSqr, minEqOp());
        Pstream::listCombineScatter(regionMinDistSqr);

        // Choose the closest region
        const label selectedRegioni = findMin(regionMinDistSqr);

        // Remove faces in other regions from the list by shuffling up
        label fi0 = 0;
        forAll(faceIndices, fi)
        {
            faceIndices[fi0] = faceIndices[fi];
            if (prs[fi] == selectedRegioni) fi0 ++;
        }
        faceIndices.resize(fi0);
    }

    boolList flipMap(faceIndices.size());

    // Construct the flipMap
    forAll(faceIndices, fi)
    {
        flipMap[fi] = cellIsAbovePlane[mesh_.faceOwner()[faceIndices[fi]]];
    }

    return zoneSet
    (
        new faceZone
        (
            zoneName_,
            faceIndices,
            flipMap,
            mesh_.faceZones(),
            moveUpdate_,
            true
        )
    );
}


// ************************************************************************* //
