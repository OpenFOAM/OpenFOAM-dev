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

#include "planeZoneGenerator.H"
#include "polyMesh.H"
#include "uindirectPrimitivePatch.H"
#include "PatchTools.H"
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
    point_(dict.lookup<vector>("point", dimLength)),
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
    forAll(mesh_.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchi];

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
    syncTools::syncFaceList(mesh_, faceIsOnPlane, orEqOp<bool>());

    // Convert marked faces to a list of indices
    labelList faceIndices(findIndices(faceIsOnPlane, true));

    // If constructing a single contiguous set, remove all faces except those
    // connected to the contiguous region closest to the specified point
    if (include_ == include::closest)
    {
        // Step 1: Get locally contiguous regions for the new face set and the
        // total number of regions across all processors.
        labelList newSetFaceRegions(faceIndices.size(), -1);
        label nRegions = -1;
        {
            // Create a patch of the set faces
            const uindirectPrimitivePatch newSetPatch
            (
                UIndirectList<face>(mesh_.faces(), faceIndices),
                mesh_.points()
            );

            // Get the region ID-s and store the total number of regions on
            // each processor
            labelList procNRegions(Pstream::nProcs(), -1);
            procNRegions[Pstream::myProcNo()] =
                PatchTools::markZones
                (
                    newSetPatch,
                    boolList(newSetPatch.nEdges(), false),
                    newSetFaceRegions
                );
            Pstream::gatherList(procNRegions);
            Pstream::scatterList(procNRegions);

            // Cumulative sum the number of regions on each processor to get an
            // offset which makes the local region ID-s globally unique
            labelList procRegionOffset(Pstream::nProcs(), 0);
            for (label proci = 1; proci < Pstream::nProcs(); proci++)
            {
                procRegionOffset[proci] +=
                    procRegionOffset[proci - 1]
                  + procNRegions[proci - 1];
            }

            // Apply the offset
            forAll(faceIndices, fi)
            {
                newSetFaceRegions[fi] +=
                    procRegionOffset[Pstream::myProcNo()];
            }

            // Store the total number of regions across all processors
            nRegions = procRegionOffset.last() + procNRegions.last();
        }

        // Step 2: Create a region map which combines regions which are
        // connected across coupled interfaces
        labelList regionMap(identityMap(nRegions));
        {
            // Put region labels on connected boundary edges and synchronise to
            // create a list of all regions connected to a given edge
            labelListList meshEdgeRegions(mesh_.nEdges(), labelList());
            forAll(faceIndices, fi)
            {
                const label facei = faceIndices[fi];
                const label regioni = newSetFaceRegions[fi];
                forAll(mesh_.faceEdges()[facei], faceEdgei)
                {
                    const label edgei = mesh_.faceEdges()[facei][faceEdgei];
                    meshEdgeRegions[edgei] = labelList(1, regioni);
                }
            }
            syncTools::syncEdgeList
            (
                mesh_,
                meshEdgeRegions,
                globalMeshData::ListPlusEqOp<labelList>(),
                labelList()
            );

            // Combine edge regions to create a list of what regions a given
            // region is connected to
            List<labelHashSet> regionRegions(nRegions);
            forAll(faceIndices, fi)
            {
                const label facei = faceIndices[fi];
                const label regioni = newSetFaceRegions[fi];
                forAll(mesh_.faceEdges()[facei], faceEdgei)
                {
                    const label edgei = mesh_.faceEdges()[facei][faceEdgei];
                    forAll(meshEdgeRegions[edgei], edgeRegioni)
                    {
                        if (meshEdgeRegions[edgei][edgeRegioni] != regioni)
                        {
                            regionRegions[regioni].insert
                            (
                                meshEdgeRegions[edgei][edgeRegioni]
                            );
                        }
                    }
                }
            }
            Pstream::listCombineGather(regionRegions, plusEqOp<labelHashSet>());
            Pstream::listCombineScatter(regionRegions);

            // Collapse the region connections into a map between each region
            // and the lowest numbered region that it connects to
            forAll(regionRegions, regioni)
            {
                forAllConstIter(labelHashSet, regionRegions[regioni], iter)
                {
                    regionMap[iter.key()] =
                        min(regionMap[iter.key()], regionMap[regioni]);
                }
            }
        }

        // Step 3: Combine connected regions
        labelList regionNFaces;
        {
            // Remove duplicates from the region map
            label regioni0 = 0;
            forAll(regionMap, regioni)
            {
                if (regionMap[regioni] > regioni0)
                {
                    regioni0++;
                    regionMap[regioni] = regioni0;
                }
            }

            // Recompute the number of regions
            nRegions = regioni0 + 1;

            // Renumber the face region ID-s
            newSetFaceRegions =
                IndirectList<label>(regionMap, newSetFaceRegions);

            // Report the final number and size of the regions
            regionNFaces = labelList(nRegions, 0);
            forAll(faceIndices, fi)
            {
                regionNFaces[newSetFaceRegions[fi]] ++;
            }
            Pstream::listCombineGather(regionNFaces, plusEqOp<label>());
            Pstream::listCombineScatter(regionNFaces);
            Info<< "    Found " << nRegions << " contiguous regions with "
                << regionNFaces << " faces" << endl;
        }

        // Step 4: Choose the closest region to output
        label selectedRegioni = -1;
        {
            // Compute the region centres
            scalarField regionMagAreas(nRegions, 0);
            pointField regionCentres(nRegions, Zero);
            forAll(faceIndices, fi)
            {
                const label facei = faceIndices[fi];
                const label regioni = newSetFaceRegions[fi];

                const vector& a = mesh_.faceAreas()[facei];
                const point& c = mesh_.faceCentres()[facei];

                regionMagAreas[regioni] += mag(a);
                regionCentres[regioni] += mag(a)*c;
            }
            Pstream::listCombineGather(regionMagAreas, plusEqOp<scalar>());
            Pstream::listCombineGather(regionCentres, plusEqOp<point>());
            Pstream::listCombineScatter(regionMagAreas);
            Pstream::listCombineScatter(regionCentres);
            regionCentres /= regionMagAreas;

            // Find the region centroid closest to the reference point
            selectedRegioni = returnReduce
            (
                findMin(mag(regionCentres - point_)()),
                minOp<label>()
            );

            // Report the selection
            Info<< "    Selecting region " << selectedRegioni << " with "
                << regionNFaces[selectedRegioni]
                << " faces as the closest to point " << point_ << endl;
        }

        // Step 5: Remove any faces from the set list not in the selected region
        {
            // Remove faces from the list by shuffling up and resizing
            label fi0 = 0;
            forAll(faceIndices, fi)
            {
                faceIndices[fi0] = faceIndices[fi];

                if (newSetFaceRegions[fi] == selectedRegioni)
                {
                    fi0++;
                }
            }
            faceIndices.resize(fi0);
        }
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
