/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Description
    Post-processing mesh subset tool.  Given the original mesh and the
    list of selected cells, it creates the mesh consisting only of the
    desired cells, with the mapping list for points, faces, and cells.

\*---------------------------------------------------------------------------*/

#include "fvMeshSubset.H"
#include "boolList.H"
#include "Pstream.H"
#include "emptyPolyPatch.H"
#include "demandDrivenData.H"
#include "cyclicPolyPatch.H"
#include "removeCells.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::fvMeshSubset::checkCellSubset() const
{
    if (fvMeshSubsetPtr_.empty())
    {
        FatalErrorInFunction
            << "void setCellSubset(const labelHashSet& cellsToSubset)" << endl
            << "before attempting to access subset data"
            << abort(FatalError);

        return false;
    }
    else
    {
        return true;
    }
}


void Foam::fvMeshSubset::markPoints
(
    const labelList& curPoints,
    Map<label>& pointMap
)
{
    forAll(curPoints, pointi)
    {
        // Note: insert will only insert if not yet there.
        pointMap.insert(curPoints[pointi], 0);
    }
}


void Foam::fvMeshSubset::markPoints
(
    const labelList& curPoints,
    labelList& pointMap
)
{
    forAll(curPoints, pointi)
    {
        pointMap[curPoints[pointi]] = 0;
    }
}


void Foam::fvMeshSubset::doCoupledPatches
(
    const bool syncPar,
    Map<label>& facesToSubset,
    labelList& nCellsUsingFace
) const
{
    // Synchronize facesToSubset on both sides of coupled patches.
    // Marks faces that become 'uncoupled' with 3.

    const polyBoundaryMesh& oldPatches = baseMesh().boundaryMesh();

    label nUncoupled = 0;

    if (syncPar && Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send face usage across processor patches
        forAll(oldPatches, oldPatchi)
        {
            const polyPatch& pp = oldPatches[oldPatchi];

            if (isA<processorPolyPatch>(pp))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                UOPstream toNeighbour(procPatch.neighbProcNo(), pBufs);

                if (!facesToSubset.empty())
                {
                    DynamicList<label> patchFacesToSubset;
                    forAll(pp, i)
                    {
                        if
                        (
                            facesToSubset.found(pp.start()+i)
                         && facesToSubset[pp.start()+i] == 1
                        )
                        {
                            patchFacesToSubset.append(i);
                        }
                    }
                    toNeighbour << patchFacesToSubset;
                }
                else if (!nCellsUsingFace.empty())
                {
                    toNeighbour <<
                        SubList<label>(nCellsUsingFace, pp.size(), pp.start());
                }
                else
                {
                    toNeighbour << labelList();
                }
            }
        }

        pBufs.finishedSends();

        // Receive face usage count and check for faces that become uncoupled.
        forAll(oldPatches, oldPatchi)
        {
            const polyPatch& pp = oldPatches[oldPatchi];

            if (isA<processorPolyPatch>(pp))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                UIPstream fromNeighbour(procPatch.neighbProcNo(), pBufs);

                const labelList nbrList(fromNeighbour);

                // Combine with this side.

                if (!facesToSubset.empty())
                {
                    const labelHashSet nbrPatchFacesToSubset(nbrList);

                    forAll(pp, i)
                    {
                        if
                        (
                            facesToSubset.found(pp.start()+i)
                         && facesToSubset[pp.start()+i] == 1
                         && !nbrPatchFacesToSubset.found(i)
                        )
                        {
                            // Face's neighbour is no longer there. Mark face
                            // off as coupled
                            facesToSubset[pp.start()+i] = 3;
                            nUncoupled++;
                        }
                    }
                }
                else if (!nCellsUsingFace.empty())
                {
                    const labelList& nbrCellsUsingFace(nbrList);

                    // Combine with this side.

                    forAll(pp, i)
                    {
                        if
                        (
                            nCellsUsingFace[pp.start()+i] == 1
                         && nbrCellsUsingFace[i] == 0
                        )
                        {
                            // Face's neighbour is no longer there. Mark face
                            // off as coupled
                            nCellsUsingFace[pp.start()+i] = 3;
                            nUncoupled++;
                        }
                    }
                }
            }
        }
    }

    // Do same for cyclics.
    forAll(oldPatches, oldPatchi)
    {
        const polyPatch& pp = oldPatches[oldPatchi];

        if (isA<cyclicPolyPatch>(pp))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(pp);

            if (!facesToSubset.empty())
            {
                forAll(cycPatch, i)
                {
                    label thisFacei = cycPatch.start() + i;
                    label otherFacei = cycPatch.transformGlobalFace(thisFacei);

                    if
                    (
                        facesToSubset.found(thisFacei)
                     && facesToSubset[thisFacei] == 1
                     && !facesToSubset.found(otherFacei)
                    )
                    {
                        facesToSubset[thisFacei] = 3;
                        nUncoupled++;
                    }
                }
            }
            else if (!nCellsUsingFace.empty())
            {
                forAll(cycPatch, i)
                {
                    label thisFacei = cycPatch.start() + i;
                    label otherFacei = cycPatch.transformGlobalFace(thisFacei);

                    if
                    (
                        nCellsUsingFace[thisFacei] == 1
                     && nCellsUsingFace[otherFacei] == 0
                    )
                    {
                        nCellsUsingFace[thisFacei] = 3;
                        nUncoupled++;
                    }
                }
            }
        }
    }

    if (syncPar)
    {
        reduce(nUncoupled, sumOp<label>());
    }

    if (nUncoupled > 0)
    {
        Info<< "Uncoupled " << nUncoupled << " faces on coupled patches. "
            << "(processorPolyPatch, cyclicPolyPatch)" << endl;
    }
}


labelList Foam::fvMeshSubset::subset
(
    const label nElems,
    const labelList& selectedElements,
    const labelList& subsetMap
)
{
    // Mark selected elements.
    boolList selected(nElems, false);
    forAll(selectedElements, i)
    {
        selected[selectedElements[i]] = true;
    }

    // Count subset of selected elements
    label n = 0;
    forAll(subsetMap, i)
    {
        if (selected[subsetMap[i]])
        {
            n++;
        }
    }

    // Collect selected elements
    labelList subsettedElements(n);
    n = 0;

    forAll(subsetMap, i)
    {
        if (selected[subsetMap[i]])
        {
            subsettedElements[n++] = i;
        }
    }

    return subsettedElements;
}


void Foam::fvMeshSubset::subsetZones()
{
    // Keep all zones, even if zero size.

    const pointZoneMesh& pointZones = baseMesh().pointZones();

    // PointZones
    List<pointZone*> pZonePtrs(pointZones.size());

    forAll(pointZones, i)
    {
        const pointZone& pz = pointZones[i];

        pZonePtrs[i] = new pointZone
        (
            pz.name(),
            subset(baseMesh().nPoints(), pz, pointMap()),
            i,
            fvMeshSubsetPtr_().pointZones()
        );
    }


    // FaceZones

    const faceZoneMesh& faceZones = baseMesh().faceZones();


    // Do we need to remove zones where the side we're interested in
    // no longer exists? Guess not.
    List<faceZone*> fZonePtrs(faceZones.size());

    forAll(faceZones, i)
    {
        const faceZone& fz = faceZones[i];

        // Expand faceZone to full mesh
        // +1 : part of faceZone, flipped
        // -1 :    ,,           , unflipped
        //  0 : not part of faceZone
        labelList zone(baseMesh().nFaces(), 0);
        forAll(fz, j)
        {
            if (fz.flipMap()[j])
            {
                zone[fz[j]] = 1;
            }
            else
            {
                zone[fz[j]] = -1;
            }
        }

        // Select faces
        label nSub = 0;
        forAll(faceMap(), j)
        {
            if (zone[faceMap()[j]] != 0)
            {
                nSub++;
            }
        }
        labelList subAddressing(nSub);
        boolList subFlipStatus(nSub);
        nSub = 0;
        forAll(faceMap(), subFacei)
        {
            label meshFacei = faceMap()[subFacei];
            if (zone[meshFacei] != 0)
            {
                subAddressing[nSub] = subFacei;
                label subOwner = subMesh().faceOwner()[subFacei];
                label baseOwner = baseMesh().faceOwner()[meshFacei];
                // If subowner is the same cell as the base keep the flip status
                bool sameOwner = (cellMap()[subOwner] == baseOwner);
                bool flip = (zone[meshFacei] == 1);
                subFlipStatus[nSub] = (sameOwner == flip);

                nSub++;
            }
        }

        fZonePtrs[i] = new faceZone
        (
            fz.name(),
            subAddressing,
            subFlipStatus,
            i,
            fvMeshSubsetPtr_().faceZones()
        );
    }


    const cellZoneMesh& cellZones = baseMesh().cellZones();

    List<cellZone*> cZonePtrs(cellZones.size());

    forAll(cellZones, i)
    {
        const cellZone& cz = cellZones[i];

        cZonePtrs[i] = new cellZone
        (
            cz.name(),
            subset(baseMesh().nCells(), cz, cellMap()),
            i,
            fvMeshSubsetPtr_().cellZones()
        );
    }


    // Add the zones
    fvMeshSubsetPtr_().addZones(pZonePtrs, fZonePtrs, cZonePtrs);
}


Foam::labelList Foam::fvMeshSubset::getCellsToRemove
(
    const labelList& region,
    const label currentRegion
) const
{
    // Count
    label nKeep = 0;
    forAll(region, cellI)
    {
        if (region[cellI] == currentRegion)
        {
            nKeep++;
        }
    }

    // Collect cells to remove
    label nRemove = baseMesh().nCells() - nKeep;
    labelList cellsToRemove(nRemove);

    nRemove = 0;
    forAll(region, cellI)
    {
        if (region[cellI] != currentRegion)
        {
            cellsToRemove[nRemove++] = cellI;
        }
    }

    return cellsToRemove;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshSubset::fvMeshSubset(const fvMesh& baseMesh)
:
    baseMesh_(baseMesh),
    fvMeshSubsetPtr_(nullptr),
    pointMap_(0),
    faceMap_(0),
    cellMap_(0),
    patchMap_(0),
    faceFlipMapPtr_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshSubset::setCellSubset
(
    const labelHashSet& globalCellMap,
    const label patchID,
    const bool syncPar
)
{
    // Initial check on patches before doing anything time consuming.
    const polyBoundaryMesh& oldPatches = baseMesh().boundaryMesh();
    const cellList& oldCells = baseMesh().cells();
    const faceList& oldFaces = baseMesh().faces();
    const pointField& oldPoints = baseMesh().points();
    const labelList& oldOwner = baseMesh().faceOwner();
    const labelList& oldNeighbour = baseMesh().faceNeighbour();

    label wantedPatchID = patchID;

    if (wantedPatchID == -1)
    {
        // No explicit patch specified. Put in oldInternalFaces patch.
        // Check if patch with this name already exists.
        wantedPatchID = oldPatches.findPatchID("oldInternalFaces");
    }
    else if (wantedPatchID < 0 || wantedPatchID >= oldPatches.size())
    {
        FatalErrorInFunction
            << "Non-existing patch index " << wantedPatchID << endl
            << "Should be between 0 and " << oldPatches.size()-1
            << abort(FatalError);
    }


    // Clear demand driven data
    faceFlipMapPtr_.clear();


    cellMap_ = globalCellMap.toc();

    // Sort the cell map in the ascending order
    sort(cellMap_);

    // Approximate sizing parameters for face and point lists
    const label avgNFacesPerCell = 6;
    const label avgNPointsPerFace = 4;


    label nCellsInSet = cellMap_.size();

    // Mark all used faces

    Map<label> facesToSubset(avgNFacesPerCell*nCellsInSet);

    forAll(cellMap_, celli)
    {
        // Mark all faces from the cell
        const labelList& curFaces = oldCells[cellMap_[celli]];

        forAll(curFaces, facei)
        {
            if (!facesToSubset.found(curFaces[facei]))
            {
                facesToSubset.insert(curFaces[facei], 1);
            }
            else
            {
                facesToSubset[curFaces[facei]]++;
            }
        }
    }

    // Handle coupled faces. Modifies patch faces to be uncoupled to 3.
    labelList empty;
    doCoupledPatches(syncPar, facesToSubset, empty);

    // Mark all used points and make a global-to-local face map
    Map<label> globalFaceMap(facesToSubset.size());

    // Make a global-to-local point map
    Map<label> globalPointMap(avgNPointsPerFace*facesToSubset.size());

    // This is done in two goes, so that the boundary faces are last
    // in the list.  Because of this, I need to create the face map
    // along the way rather than just grab the table of contents.
    labelList facesToc = facesToSubset.toc();
    sort(facesToc);
    faceMap_.setSize(facesToc.size());

    // 1. Get all faces that will be internal to the submesh.
    forAll(facesToc, facei)
    {
        if (facesToSubset[facesToc[facei]] == 2)
        {
            // Mark face and increment number of points in set
            faceMap_[globalFaceMap.size()] = facesToc[facei];
            globalFaceMap.insert(facesToc[facei], globalFaceMap.size());

            // Mark all points from the face
            markPoints(oldFaces[facesToc[facei]], globalPointMap);
        }
    }

    // These are all the internal faces in the mesh.
    label nInternalFaces = globalFaceMap.size();


    // Where to insert old internal faces.
    label oldPatchStart = labelMax;
    if (wantedPatchID != -1)
    {
        oldPatchStart = oldPatches[wantedPatchID].start();
    }


    label facei = 0;

    // 2. Boundary faces up to where we want to insert old internal faces
    for (; facei< facesToc.size(); facei++)
    {
        if (facesToc[facei] >= oldPatchStart)
        {
            break;
        }
        if
        (
            !baseMesh().isInternalFace(facesToc[facei])
         && facesToSubset[facesToc[facei]] == 1
        )
        {
            // Mark face and increment number of points in set
            faceMap_[globalFaceMap.size()] = facesToc[facei];
            globalFaceMap.insert(facesToc[facei], globalFaceMap.size());

            // Mark all points from the face
            markPoints(oldFaces[facesToc[facei]], globalPointMap);
        }
    }

    // 3. old internal faces and uncoupled faces
    forAll(facesToc, intFacei)
    {
        if
        (
            (
                baseMesh().isInternalFace(facesToc[intFacei])
             && facesToSubset[facesToc[intFacei]] == 1
            )
         || (
                !baseMesh().isInternalFace(facesToc[intFacei])
             && facesToSubset[facesToc[intFacei]] == 3
            )
        )
        {
            // Mark face and increment number of points in set
            faceMap_[globalFaceMap.size()] = facesToc[intFacei];
            globalFaceMap.insert(facesToc[intFacei], globalFaceMap.size());

            // Mark all points from the face
            markPoints(oldFaces[facesToc[intFacei]], globalPointMap);
        }
    }

    // 4. Remaining boundary faces
    for (; facei< facesToc.size(); facei++)
    {
        if
        (
            !baseMesh().isInternalFace(facesToc[facei])
         && facesToSubset[facesToc[facei]] == 1
        )
        {
            // Mark face and increment number of points in set
            faceMap_[globalFaceMap.size()] = facesToc[facei];
            globalFaceMap.insert(facesToc[facei], globalFaceMap.size());

            // Mark all points from the face
            markPoints(oldFaces[facesToc[facei]], globalPointMap);
        }
    }



    // Grab the points map
    pointMap_ = globalPointMap.toc();
    sort(pointMap_);

    forAll(pointMap_, pointi)
    {
        globalPointMap[pointMap_[pointi]] = pointi;
    }

    // Pout<< "Number of cells in new mesh: " << nCellsInSet << endl;
    // Pout<< "Number of faces in new mesh: " << globalFaceMap.size() << endl;
    // Pout<< "Number of points in new mesh: " << globalPointMap.size() << endl;

    // Make a new mesh
    pointField newPoints(globalPointMap.size());

    label nNewPoints = 0;

    forAll(pointMap_, pointi)
    {
        newPoints[nNewPoints] = oldPoints[pointMap_[pointi]];
        nNewPoints++;
    }

    faceList newFaces(globalFaceMap.size());

    label nNewFaces = 0;

    // Make internal faces
    for (label facei = 0; facei < nInternalFaces; facei++)
    {
        const face& oldF = oldFaces[faceMap_[facei]];

        face newF(oldF.size());

        forAll(newF, i)
        {
            newF[i] = globalPointMap[oldF[i]];
        }

        newFaces[nNewFaces] = newF;
        nNewFaces++;
    }

    // Make boundary faces

    label nbSize = oldPatches.size();
    label oldInternalPatchID  = -1;

    if (wantedPatchID == -1)
    {
        // Create 'oldInternalFaces' patch at the end
        // and put all exposed internal faces in there.
        oldInternalPatchID = nbSize;
        nbSize++;

    }
    else
    {
        oldInternalPatchID = wantedPatchID;
    }


    // Grad size and start of each patch on the fly.  Because of the
    // structure of the underlying mesh, the patches will appear in the
    // ascending order
    labelList boundaryPatchSizes(nbSize, 0);

    // Assign boundary faces. Visited in order of faceMap_.
    for (label facei = nInternalFaces; facei < faceMap_.size(); facei++)
    {
        label oldFacei = faceMap_[facei];

        face oldF = oldFaces[oldFacei];

        // Turn the faces as necessary to point outwards
        if (baseMesh().isInternalFace(oldFacei))
        {
            // Internal face. Possibly turned the wrong way round
            if
            (
                !globalCellMap.found(oldOwner[oldFacei])
             && globalCellMap.found(oldNeighbour[oldFacei])
            )
            {
                oldF = oldFaces[oldFacei].reverseFace();
            }

            // Update count for patch
            boundaryPatchSizes[oldInternalPatchID]++;
        }
        else if (facesToSubset[oldFacei] == 3)
        {
            // Uncoupled face. Increment the old patch.
            boundaryPatchSizes[oldInternalPatchID]++;
        }
        else
        {
            // Boundary face. Increment the appropriate patch
            label patchOfFace = oldPatches.whichPatch(oldFacei);

            // Update count for patch
            boundaryPatchSizes[patchOfFace]++;
        }

        face newF(oldF.size());

        forAll(newF, i)
        {
            newF[i] = globalPointMap[oldF[i]];
        }

        newFaces[nNewFaces] = newF;
        nNewFaces++;
    }



    // Create cells
    cellList newCells(nCellsInSet);

    label nNewCells = 0;

    forAll(cellMap_, celli)
    {
        const labelList& oldC = oldCells[cellMap_[celli]];

        labelList newC(oldC.size());

        forAll(newC, i)
        {
            newC[i] = globalFaceMap[oldC[i]];
        }

        newCells[nNewCells] = cell(newC);
        nNewCells++;
    }


    // Delete any old mesh
    fvMeshSubsetPtr_.clear();
    // Make a new mesh
    fvMeshSubsetPtr_.reset
    (
        new fvMesh
        (
            IOobject
            (
                baseMesh().name() + "SubSet",
                baseMesh().time().timeName(),
                baseMesh().time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            xferMove(newPoints),
            xferMove(newFaces),
            xferMove(newCells)
        )
    );


    // Add old patches
    List<polyPatch*> newBoundary(nbSize);
    patchMap_.setSize(nbSize);
    label nNewPatches = 0;
    label patchStart = nInternalFaces;


    forAll(oldPatches, patchi)
    {
        if (boundaryPatchSizes[patchi] > 0)
        {
            // Patch still exists. Add it
            newBoundary[nNewPatches] = oldPatches[patchi].clone
            (
                fvMeshSubsetPtr_().boundaryMesh(),
                nNewPatches,
                boundaryPatchSizes[patchi],
                patchStart
            ).ptr();

            patchStart += boundaryPatchSizes[patchi];
            patchMap_[nNewPatches] = patchi;
            nNewPatches++;
        }
    }

    if (wantedPatchID == -1)
    {
        // Newly created patch so is at end. Check if any faces in it.
        if (boundaryPatchSizes[oldInternalPatchID] > 0)
        {
            newBoundary[nNewPatches] = new emptyPolyPatch
            (
                "oldInternalFaces",
                boundaryPatchSizes[oldInternalPatchID],
                patchStart,
                nNewPatches,
                fvMeshSubsetPtr_().boundaryMesh(),
                emptyPolyPatch::typeName
            );

            // The index for the first patch is -1 as it originates from
            // the internal faces
            patchMap_[nNewPatches] = -1;
            nNewPatches++;
        }
    }

    // Reset the patch lists
    newBoundary.setSize(nNewPatches);
    patchMap_.setSize(nNewPatches);

    // Add the fvPatches
    fvMeshSubsetPtr_().addFvPatches(newBoundary);

    // Subset and add any zones
    subsetZones();
}


void Foam::fvMeshSubset::setLargeCellSubset
(
    const labelList& region,
    const label currentRegion,
    const label patchID,
    const bool syncPar
)
{
    const cellList& oldCells = baseMesh().cells();
    const faceList& oldFaces = baseMesh().faces();
    const pointField& oldPoints = baseMesh().points();
    const labelList& oldOwner = baseMesh().faceOwner();
    const labelList& oldNeighbour = baseMesh().faceNeighbour();
    const polyBoundaryMesh& oldPatches = baseMesh().boundaryMesh();
    const label oldNInternalFaces = baseMesh().nInternalFaces();

    // Initial checks

    if (region.size() != oldCells.size())
    {
        FatalErrorInFunction
            << "Size of region " << region.size()
            << " is not equal to number of cells in mesh " << oldCells.size()
            << abort(FatalError);
    }


    label wantedPatchID = patchID;

    if (wantedPatchID == -1)
    {
        // No explicit patch specified. Put in oldInternalFaces patch.
        // Check if patch with this name already exists.
        wantedPatchID = oldPatches.findPatchID("oldInternalFaces");
    }
    else if (wantedPatchID < 0 || wantedPatchID >= oldPatches.size())
    {
        FatalErrorInFunction
            << "Non-existing patch index " << wantedPatchID << endl
            << "Should be between 0 and " << oldPatches.size()-1
            << abort(FatalError);
    }

    // Clear demand driven data
    faceFlipMapPtr_.clear();

    // Get the cells for the current region.
    cellMap_.setSize(oldCells.size());
    label nCellsInSet = 0;

    forAll(region, oldCelli)
    {
        if (region[oldCelli] == currentRegion)
        {
            cellMap_[nCellsInSet++] = oldCelli;
        }
    }
    cellMap_.setSize(nCellsInSet);


    // Mark all used faces. Count number of cells using them
    // 0: face not used anymore
    // 1: face used by one cell, face becomes/stays boundary face
    // 2: face still used and remains internal face
    // 3: face coupled and used by one cell only (so should become normal,
    //    non-coupled patch face)
    //
    // Note that this is not really necessary - but means we can size things
    // correctly. Also makes handling coupled faces much easier.

    labelList nCellsUsingFace(oldFaces.size(), 0);

    label nFacesInSet = 0;
    forAll(oldFaces, oldFacei)
    {
        bool faceUsed = false;

        if (region[oldOwner[oldFacei]] == currentRegion)
        {
            nCellsUsingFace[oldFacei]++;
            faceUsed = true;
        }

        if
        (
            baseMesh().isInternalFace(oldFacei)
         && (region[oldNeighbour[oldFacei]] == currentRegion)
        )
        {
            nCellsUsingFace[oldFacei]++;
            faceUsed = true;
        }

        if (faceUsed)
        {
            nFacesInSet++;
        }
    }
    faceMap_.setSize(nFacesInSet);

    // Handle coupled faces. Modifies patch faces to be uncoupled to 3.
    Map<label> empty;
    doCoupledPatches(syncPar, empty, nCellsUsingFace);


    // See which patch to use for exposed internal faces.
    label oldInternalPatchID = 0;

    // Insert faces before which patch
    label nextPatchID = oldPatches.size();

    // old to new patches
    labelList globalPatchMap(oldPatches.size());

    // New patch size
    label nbSize = oldPatches.size();

    if (wantedPatchID == -1)
    {
        // Create 'oldInternalFaces' patch at the end (or before
        // processorPatches)
        // and put all exposed internal faces in there.

        forAll(oldPatches, patchi)
        {
            if (isA<processorPolyPatch>(oldPatches[patchi]))
            {
                nextPatchID = patchi;
                break;
            }
            oldInternalPatchID++;
        }

        nbSize++;

        // adapt old to new patches for inserted patch
        for (label oldPatchi = 0; oldPatchi < nextPatchID; oldPatchi++)
        {
            globalPatchMap[oldPatchi] = oldPatchi;
        }
        for
        (
            label oldPatchi = nextPatchID;
            oldPatchi < oldPatches.size();
            oldPatchi++
        )
        {
            globalPatchMap[oldPatchi] = oldPatchi+1;
        }
    }
    else
    {
        oldInternalPatchID = wantedPatchID;
        nextPatchID = wantedPatchID+1;

        // old to new patches
        globalPatchMap = identity(oldPatches.size());
    }

    labelList boundaryPatchSizes(nbSize, 0);


    // Make a global-to-local point map
    labelList globalPointMap(oldPoints.size(), -1);

    labelList globalFaceMap(oldFaces.size(), -1);
    label facei = 0;

    // 1. Pick up all preserved internal faces.
    for (label oldFacei = 0; oldFacei < oldNInternalFaces; oldFacei++)
    {
        if (nCellsUsingFace[oldFacei] == 2)
        {
            globalFaceMap[oldFacei] = facei;
            faceMap_[facei++] = oldFacei;

            // Mark all points from the face
            markPoints(oldFaces[oldFacei], globalPointMap);
        }
    }

    // These are all the internal faces in the mesh.
    label nInternalFaces = facei;

    // 2. Boundary faces up to where we want to insert old internal faces
    for
    (
        label oldPatchi = 0;
        oldPatchi < oldPatches.size()
     && oldPatchi < nextPatchID;
        oldPatchi++
    )
    {
        const polyPatch& oldPatch = oldPatches[oldPatchi];

        label oldFacei = oldPatch.start();

        forAll(oldPatch, i)
        {
            if (nCellsUsingFace[oldFacei] == 1)
            {
                // Boundary face is kept.

                // Mark face and increment number of points in set
                globalFaceMap[oldFacei] = facei;
                faceMap_[facei++] = oldFacei;

                // Mark all points from the face
                markPoints(oldFaces[oldFacei], globalPointMap);

                // Increment number of patch faces
                boundaryPatchSizes[globalPatchMap[oldPatchi]]++;
            }
            oldFacei++;
        }
    }

    // 3a. old internal faces that have become exposed.
    for (label oldFacei = 0; oldFacei < oldNInternalFaces; oldFacei++)
    {
        if (nCellsUsingFace[oldFacei] == 1)
        {
            globalFaceMap[oldFacei] = facei;
            faceMap_[facei++] = oldFacei;

            // Mark all points from the face
            markPoints(oldFaces[oldFacei], globalPointMap);

            // Increment number of patch faces
            boundaryPatchSizes[oldInternalPatchID]++;
        }
    }

    // 3b. coupled patch faces that have become uncoupled.
    for
    (
        label oldFacei = oldNInternalFaces;
        oldFacei < oldFaces.size();
        oldFacei++
    )
    {
        if (nCellsUsingFace[oldFacei] == 3)
        {
            globalFaceMap[oldFacei] = facei;
            faceMap_[facei++] = oldFacei;

            // Mark all points from the face
            markPoints(oldFaces[oldFacei], globalPointMap);

            // Increment number of patch faces
            boundaryPatchSizes[oldInternalPatchID]++;
        }
    }

    // 4. Remaining boundary faces
    for
    (
        label oldPatchi = nextPatchID;
        oldPatchi < oldPatches.size();
        oldPatchi++
    )
    {
        const polyPatch& oldPatch = oldPatches[oldPatchi];

        label oldFacei = oldPatch.start();

        forAll(oldPatch, i)
        {
            if (nCellsUsingFace[oldFacei] == 1)
            {
                // Boundary face is kept.

                // Mark face and increment number of points in set
                globalFaceMap[oldFacei] = facei;
                faceMap_[facei++] = oldFacei;

                // Mark all points from the face
                markPoints(oldFaces[oldFacei], globalPointMap);

                // Increment number of patch faces
                boundaryPatchSizes[globalPatchMap[oldPatchi]]++;
            }
            oldFacei++;
        }
    }

    if (facei != nFacesInSet)
    {
        FatalErrorInFunction
            << "Problem" << abort(FatalError);
    }


    // Grab the points map
    label nPointsInSet = 0;

    forAll(globalPointMap, pointi)
    {
        if (globalPointMap[pointi] != -1)
        {
            nPointsInSet++;
        }
    }
    pointMap_.setSize(nPointsInSet);

    nPointsInSet = 0;

    forAll(globalPointMap, pointi)
    {
        if (globalPointMap[pointi] != -1)
        {
            pointMap_[nPointsInSet] = pointi;
            globalPointMap[pointi] = nPointsInSet;
            nPointsInSet++;
        }
    }

    // Pout<< "Number of cells in new mesh : " << cellMap_.size() << endl;
    // Pout<< "Number of faces in new mesh : " << faceMap_.size() << endl;
    // Pout<< "Number of points in new mesh: " << pointMap_.size() << endl;

    // Make a new mesh
    pointField newPoints(pointMap_.size());

    label nNewPoints = 0;

    forAll(pointMap_, pointi)
    {
        newPoints[nNewPoints] = oldPoints[pointMap_[pointi]];
        nNewPoints++;
    }

    faceList newFaces(faceMap_.size());

    label nNewFaces = 0;

    // Make internal faces
    for (label facei = 0; facei < nInternalFaces; facei++)
    {
        const face& oldF = oldFaces[faceMap_[facei]];

        face newF(oldF.size());

        forAll(newF, i)
        {
            newF[i] = globalPointMap[oldF[i]];
        }

        newFaces[nNewFaces] = newF;
        nNewFaces++;
    }


    // Make boundary faces. (different from internal since might need to be
    // flipped)
    for (label facei = nInternalFaces; facei < faceMap_.size(); facei++)
    {
        label oldFacei = faceMap_[facei];

        face oldF = oldFaces[oldFacei];

        // Turn the faces as necessary to point outwards
        if (baseMesh().isInternalFace(oldFacei))
        {
            // Was internal face. Possibly turned the wrong way round
            if
            (
                region[oldOwner[oldFacei]] != currentRegion
             && region[oldNeighbour[oldFacei]] == currentRegion
            )
            {
                oldF = oldFaces[oldFacei].reverseFace();
            }
        }

        // Relabel vertices of the (possibly turned) face.
        face newF(oldF.size());

        forAll(newF, i)
        {
            newF[i] = globalPointMap[oldF[i]];
        }

        newFaces[nNewFaces] = newF;
        nNewFaces++;
    }



    // Create cells
    cellList newCells(nCellsInSet);

    label nNewCells = 0;

    forAll(cellMap_, celli)
    {
        const labelList& oldC = oldCells[cellMap_[celli]];

        labelList newC(oldC.size());

        forAll(newC, i)
        {
            newC[i] = globalFaceMap[oldC[i]];
        }

        newCells[nNewCells] = cell(newC);
        nNewCells++;
    }


    // Delete any old one
    fvMeshSubsetPtr_.clear();

    // Make a new mesh
    // Note that mesh gets registered with same name as original mesh. This is
    // not proper but cannot be avoided since otherwise surfaceInterpolation
    // cannot find its fvSchemes (it will try to read e.g.
    // system/region0SubSet/fvSchemes)
    // Make a new mesh
    fvMeshSubsetPtr_.reset
    (
        new fvMesh
        (
            IOobject
            (
                baseMesh().name(),
                baseMesh().time().timeName(),
                baseMesh().time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            xferMove(newPoints),
            xferMove(newFaces),
            xferMove(newCells),
            syncPar           // parallel synchronisation
        )
    );

    // Add old patches
    List<polyPatch*> newBoundary(nbSize);
    patchMap_.setSize(nbSize);
    label nNewPatches = 0;
    label patchStart = nInternalFaces;


    // For parallel: only remove patch if none of the processors has it.
    // This only gets done for patches before the one being inserted
    // (so patches < nextPatchID)

    // Get sum of patch sizes. Zero if patch can be deleted.
    labelList globalPatchSizes(boundaryPatchSizes);
    globalPatchSizes.setSize(nextPatchID);

    if (syncPar && Pstream::parRun())
    {
        // Get patch names (up to nextPatchID)
        List<wordList> patchNames(Pstream::nProcs());
        patchNames[Pstream::myProcNo()] = oldPatches.names();
        patchNames[Pstream::myProcNo()].setSize(nextPatchID);
        Pstream::gatherList(patchNames);
        Pstream::scatterList(patchNames);

        // Get patch sizes (up to nextPatchID).
        // Note that up to nextPatchID the globalPatchMap is an identity so
        // no need to index through that.
        Pstream::listCombineGather(globalPatchSizes, plusEqOp<label>());
        Pstream::listCombineScatter(globalPatchSizes);

        // Now all processors have all the patchnames.
        // Decide: if all processors have the same patch names and size is zero
        // everywhere remove the patch.
        bool samePatches = true;

        for (label proci = 1; proci < patchNames.size(); proci++)
        {
            if (patchNames[proci] != patchNames[0])
            {
                samePatches = false;
                break;
            }
        }

        if (!samePatches)
        {
            // Patchnames not sync on all processors so disable removal of
            // zero sized patches.
            globalPatchSizes = labelMax;
        }
    }


    // Old patches

    for
    (
        label oldPatchi = 0;
        oldPatchi < oldPatches.size()
     && oldPatchi < nextPatchID;
        oldPatchi++
    )
    {
        label newSize = boundaryPatchSizes[globalPatchMap[oldPatchi]];

        // Clone (even if 0 size)
        newBoundary[nNewPatches] = oldPatches[oldPatchi].clone
        (
            fvMeshSubsetPtr_().boundaryMesh(),
            nNewPatches,
            newSize,
            patchStart
        ).ptr();

        patchStart += newSize;
        patchMap_[nNewPatches] = oldPatchi;    // compact patchMap
        nNewPatches++;
    }

    // Inserted patch

    if (wantedPatchID == -1)
    {
        label oldInternalSize = boundaryPatchSizes[oldInternalPatchID];

        if (syncPar)
        {
            reduce(oldInternalSize, sumOp<label>());
        }

        // Newly created patch so is at end. Check if any faces in it.
        if (oldInternalSize > 0)
        {
            newBoundary[nNewPatches] = new emptyPolyPatch
            (
                "oldInternalFaces",
                boundaryPatchSizes[oldInternalPatchID],
                patchStart,
                nNewPatches,
                fvMeshSubsetPtr_().boundaryMesh(),
                emptyPolyPatch::typeName
            );

            // Pout<< "    oldInternalFaces : "
            //    << boundaryPatchSizes[oldInternalPatchID] << endl;

            // The index for the first patch is -1 as it originates from
            // the internal faces
            patchStart += boundaryPatchSizes[oldInternalPatchID];
            patchMap_[nNewPatches] = -1;
            nNewPatches++;
        }
    }

    // Old patches

    for
    (
        label oldPatchi = nextPatchID;
        oldPatchi < oldPatches.size();
        oldPatchi++
    )
    {
        label newSize = boundaryPatchSizes[globalPatchMap[oldPatchi]];

        // Patch still exists. Add it
        newBoundary[nNewPatches] = oldPatches[oldPatchi].clone
        (
            fvMeshSubsetPtr_().boundaryMesh(),
            nNewPatches,
            newSize,
            patchStart
        ).ptr();

        // Pout<< "    " << oldPatches[oldPatchi].name() << " : "
        //    << newSize << endl;

        patchStart += newSize;
        patchMap_[nNewPatches] = oldPatchi;    // compact patchMap
        nNewPatches++;
    }


    // Reset the patch lists
    newBoundary.setSize(nNewPatches);
    patchMap_.setSize(nNewPatches);


    // Add the fvPatches
    fvMeshSubsetPtr_().addFvPatches(newBoundary, syncPar);

    // Subset and add any zones
    subsetZones();
}


void Foam::fvMeshSubset::setLargeCellSubset
(
    const labelHashSet& globalCellMap,
    const label patchID,
    const bool syncPar
)
{
    labelList region(baseMesh().nCells(), 0);

    forAllConstIter(labelHashSet, globalCellMap, iter)
    {
        region[iter.key()] = 1;
    }
    setLargeCellSubset(region, 1, patchID, syncPar);
}


Foam::labelList Foam::fvMeshSubset::getExposedFaces
(
    const labelList& region,
    const label currentRegion,
    const bool syncCouples
) const
{
    // Collect cells to remove
    labelList cellsToRemove(getCellsToRemove(region, currentRegion));

    return removeCells(baseMesh(), syncCouples).getExposedFaces(cellsToRemove);
}


void Foam::fvMeshSubset::setLargeCellSubset
(
    const labelList& region,
    const label currentRegion,
    const labelList& exposedFaces,
    const labelList& patchIDs,
    const bool syncCouples
)
{
    // Collect cells to remove
    labelList cellsToRemove(getCellsToRemove(region, currentRegion));

    // Mesh changing engine.
    polyTopoChange meshMod(baseMesh());

    removeCells cellRemover(baseMesh(), syncCouples);

    cellRemover.setRefinement
    (
        cellsToRemove,
        exposedFaces,
        patchIDs,
        meshMod
    );

    // Create mesh, return map from old to new mesh.
    autoPtr<mapPolyMesh> map = meshMod.makeMesh
    (
        fvMeshSubsetPtr_,
        IOobject
        (
            baseMesh().name(),
            baseMesh().time().timeName(),
            baseMesh().time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        baseMesh(),
        syncCouples
    );

    pointMap_ = map().pointMap();
    faceMap_ = map().faceMap();
    cellMap_ = map().cellMap();
    patchMap_ = identity(baseMesh().boundaryMesh().size());
}


bool Foam::fvMeshSubset::hasSubMesh() const
{
    return fvMeshSubsetPtr_.valid();
}


const fvMesh& Foam::fvMeshSubset::subMesh() const
{
    checkCellSubset();

    return fvMeshSubsetPtr_();
}


fvMesh& Foam::fvMeshSubset::subMesh()
{
    checkCellSubset();

    return fvMeshSubsetPtr_();
}


const labelList& Foam::fvMeshSubset::pointMap() const
{
    checkCellSubset();

    return pointMap_;
}


const labelList& Foam::fvMeshSubset::faceMap() const
{
    checkCellSubset();

    return faceMap_;
}


const labelList& Foam::fvMeshSubset::faceFlipMap() const
{
    if (!faceFlipMapPtr_.valid())
    {
        const labelList& subToBaseFace = faceMap();
        const labelList& subToBaseCell = cellMap();

        faceFlipMapPtr_.reset(new labelList(subToBaseFace.size()));
        labelList& faceFlipMap = faceFlipMapPtr_();

        // Only exposed internal faces might be flipped (since we don't do
        // any cell renumbering, just compacting)
        label subInt = subMesh().nInternalFaces();
        const labelList& subOwn = subMesh().faceOwner();
        const labelList& own = baseMesh_.faceOwner();

        for (label subFaceI = 0; subFaceI < subInt; subFaceI++)
        {
            faceFlipMap[subFaceI] = subToBaseFace[subFaceI]+1;
        }
        for (label subFaceI = subInt; subFaceI < subOwn.size(); subFaceI++)
        {
            label faceI = subToBaseFace[subFaceI];
            if (subToBaseCell[subOwn[subFaceI]] == own[faceI])
            {
                faceFlipMap[subFaceI] = faceI+1;
            }
            else
            {
                faceFlipMap[subFaceI] = -faceI-1;
            }
        }
    }

    return faceFlipMapPtr_();
}


const labelList& Foam::fvMeshSubset::cellMap() const
{
    checkCellSubset();

    return cellMap_;
}


const labelList& Foam::fvMeshSubset::patchMap() const
{
    checkCellSubset();

    return patchMap_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
