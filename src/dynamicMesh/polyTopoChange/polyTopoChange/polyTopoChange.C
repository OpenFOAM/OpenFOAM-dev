/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "polyTopoChange.H"
#include "SortableList.H"
#include "polyMesh.H"
#include "polyAddPoint.H"
#include "polyModifyPoint.H"
#include "polyRemovePoint.H"
#include "polyAddFace.H"
#include "polyModifyFace.H"
#include "polyRemoveFace.H"
#include "polyAddCell.H"
#include "polyModifyCell.H"
#include "polyRemoveCell.H"
#include "objectMap.H"
#include "processorPolyPatch.H"
#include "fvMesh.H"
#include "CompactListList.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyTopoChange, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Renumber with special handling for merged items (marked with <-1)
void Foam::polyTopoChange::renumberReverseMap
(
    const labelList& map,
    DynamicList<label>& elems
)
{
    forAll(elems, elemI)
    {
        label val = elems[elemI];

        if (val >= 0)
        {
            elems[elemI] = map[val];
        }
        else if (val < -1)
        {
            label mergedVal = -val-2;
            elems[elemI] = -map[mergedVal]-2;
        }
    }
}


void Foam::polyTopoChange::renumber
(
    const labelList& map,
    labelHashSet& elems
)
{
    labelHashSet newElems(elems.size());

    forAllConstIter(labelHashSet, elems, iter)
    {
        label newElem = map[iter.key()];

        if (newElem >= 0)
        {
            newElems.insert(newElem);
        }
    }

    elems.transfer(newElems);
}


// Renumber and remove -1 elements.
void Foam::polyTopoChange::renumberCompact
(
    const labelList& map,
    labelList& elems
)
{
    label newElemI = 0;

    forAll(elems, elemI)
    {
        label newVal = map[elems[elemI]];

        if (newVal != -1)
        {
            elems[newElemI++] = newVal;
        }
    }
    elems.setSize(newElemI);
}


void Foam::polyTopoChange::countMap
(
    const labelList& map,
    const labelList& reverseMap,
    label& nAdd,
    label& nInflate,
    label& nMerge,
    label& nRemove
)
{
    nAdd = 0;
    nInflate = 0;
    nMerge = 0;
    nRemove = 0;

    forAll(map, newCellI)
    {
        label oldCellI = map[newCellI];

        if (oldCellI >= 0)
        {
            if (reverseMap[oldCellI] == newCellI)
            {
                // unchanged
            }
            else
            {
                // Added (from another cell v.s. inflated from face/point)
                nAdd++;
            }
        }
        else if (oldCellI == -1)
        {
            // Created from nothing
            nInflate++;
        }
        else
        {
            FatalErrorIn("countMap") << "old:" << oldCellI
                << " new:" << newCellI << abort(FatalError);
        }
    }

    forAll(reverseMap, oldCellI)
    {
        label newCellI = reverseMap[oldCellI];

        if (newCellI >= 0)
        {
            // unchanged
        }
        else if (newCellI == -1)
        {
            // removed
            nRemove++;
        }
        else
        {
            // merged into -newCellI-2
            nMerge++;
        }
    }
}


Foam::labelHashSet Foam::polyTopoChange::getSetIndices
(
    const PackedBoolList& lst
)
{
    labelHashSet values(lst.count());
    forAll(lst, i)
    {
        if (lst[i])
        {
            values.insert(i);
        }
    }
    return values;
}


void Foam::polyTopoChange::writeMeshStats(const polyMesh& mesh, Ostream& os)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    labelList patchSizes(patches.size());
    labelList patchStarts(patches.size());
    forAll(patches, patchI)
    {
        patchSizes[patchI] = patches[patchI].size();
        patchStarts[patchI] = patches[patchI].start();
    }

    os  << "    Points      : " << mesh.nPoints() << nl
        << "    Faces       : " << mesh.nFaces() << nl
        << "    Cells       : " << mesh.nCells() << nl
        << "    PatchSizes  : " << patchSizes << nl
        << "    PatchStarts : " << patchStarts << nl
        << endl;
}


void Foam::polyTopoChange::getMergeSets
(
    const labelList& reverseCellMap,
    const labelList& cellMap,
    List<objectMap>& cellsFromCells
)
{
    // Per new cell the number of old cells that have been merged into it
    labelList nMerged(cellMap.size(), 1);

    forAll(reverseCellMap, oldCellI)
    {
        label newCellI = reverseCellMap[oldCellI];

        if (newCellI < -1)
        {
            label mergeCellI = -newCellI-2;

            nMerged[mergeCellI]++;
        }
    }

    // From merged cell to set index
    labelList cellToMergeSet(cellMap.size(), -1);

    label nSets = 0;

    forAll(nMerged, cellI)
    {
        if (nMerged[cellI] > 1)
        {
            cellToMergeSet[cellI] = nSets++;
        }
    }

    // Collect cell labels.
    // Each objectMap will have
    // - index : new mesh cell label
    // - masterObjects : list of old cells that have been merged. Element 0
    //                   will be the original destination cell label.

    cellsFromCells.setSize(nSets);

    forAll(reverseCellMap, oldCellI)
    {
        label newCellI = reverseCellMap[oldCellI];

        if (newCellI < -1)
        {
            label mergeCellI = -newCellI-2;

            // oldCellI was merged into mergeCellI

            label setI = cellToMergeSet[mergeCellI];

            objectMap& mergeSet = cellsFromCells[setI];

            if (mergeSet.masterObjects().empty())
            {
                // First occurrence of master cell mergeCellI

                mergeSet.index() = mergeCellI;
                mergeSet.masterObjects().setSize(nMerged[mergeCellI]);

                // old master label
                mergeSet.masterObjects()[0] = cellMap[mergeCellI];

                // old slave label
                mergeSet.masterObjects()[1] = oldCellI;

                nMerged[mergeCellI] = 2;
            }
            else
            {
                mergeSet.masterObjects()[nMerged[mergeCellI]++] = oldCellI;
            }
        }
    }
}


bool Foam::polyTopoChange::hasValidPoints(const face& f) const
{
    forAll(f, fp)
    {
        if (f[fp] < 0 || f[fp] >= points_.size())
        {
            return false;
        }
    }
    return true;
}


Foam::pointField Foam::polyTopoChange::facePoints(const face& f) const
{
    pointField points(f.size());
    forAll(f, fp)
    {
        if (f[fp] < 0 && f[fp] >= points_.size())
        {
            FatalErrorIn("polyTopoChange::facePoints(const face&) const")
                << "Problem." << abort(FatalError);
        }
        points[fp] = points_[f[fp]];
    }
    return points;
}


void Foam::polyTopoChange::checkFace
(
    const face& f,
    const label faceI,
    const label own,
    const label nei,
    const label patchI,
    const label zoneI
) const
{
    if (nei == -1)
    {
        if (own == -1 && zoneI != -1)
        {
            // retired face
        }
        else if (patchI == -1 || patchI >= nPatches_)
        {
            FatalErrorIn
            (
                "polyTopoChange::checkFace(const face&, const label"
                ", const label, const label, const label)"
            )   << "Face has no neighbour (so external) but does not have"
                << " a valid patch" << nl
                << "f:" << f
                << " faceI(-1 if added face):" << faceI
                << " own:" << own << " nei:" << nei
                << " patchI:" << patchI << nl;
            if (hasValidPoints(f))
            {
                FatalError
                        << "points (removed points marked with "
                        << vector::max << ") " << facePoints(f);
            }
            FatalError << abort(FatalError);
        }
    }
    else
    {
        if (patchI != -1)
        {
            FatalErrorIn
            (
                "polyTopoChange::checkFace(const face&, const label"
                ", const label, const label, const label)"
            )   << "Cannot both have valid patchI and neighbour" << nl
                << "f:" << f
                << " faceI(-1 if added face):" << faceI
                << " own:" << own << " nei:" << nei
                << " patchI:" << patchI << nl;
            if (hasValidPoints(f))
            {
                FatalError
                        << "points (removed points marked with "
                        << vector::max << ") : " << facePoints(f);
            }
            FatalError << abort(FatalError);
        }

        if (nei <= own)
        {
            FatalErrorIn
            (
                "polyTopoChange::checkFace(const face&, const label"
                ", const label, const label, const label)"
            )   << "Owner cell label should be less than neighbour cell label"
                << nl
                << "f:" << f
                << " faceI(-1 if added face):" << faceI
                << " own:" << own << " nei:" << nei
                << " patchI:" << patchI << nl;
            if (hasValidPoints(f))
            {
                FatalError
                        << "points (removed points marked with "
                        << vector::max << ") : " << facePoints(f);
            }
            FatalError << abort(FatalError);
        }
    }

    if (f.size() < 3 || findIndex(f, -1) != -1)
    {
        FatalErrorIn
        (
            "polyTopoChange::checkFace(const face&, const label"
            ", const label, const label, const label)"
        )   << "Illegal vertices in face"
            << nl
            << "f:" << f
            << " faceI(-1 if added face):" << faceI
            << " own:" << own << " nei:" << nei
            << " patchI:" << patchI << nl;
            if (hasValidPoints(f))
            {
                FatalError
                        << "points (removed points marked with "
                        << vector::max << ") : " << facePoints(f);
            }
            FatalError << abort(FatalError);
    }
    if (faceI >= 0 && faceI < faces_.size() && faceRemoved(faceI))
    {
        FatalErrorIn
        (
            "polyTopoChange::checkFace(const face&, const label"
            ", const label, const label, const label)"
        )   << "Face already marked for removal"
            << nl
            << "f:" << f
            << " faceI(-1 if added face):" << faceI
            << " own:" << own << " nei:" << nei
            << " patchI:" << patchI << nl;
            if (hasValidPoints(f))
            {
                FatalError
                        << "points (removed points marked with "
                        << vector::max << ") : " << facePoints(f);
            }
            FatalError << abort(FatalError);
    }
    forAll(f, fp)
    {
        if (f[fp] < points_.size() && pointRemoved(f[fp]))
        {
            FatalErrorIn
            (
                "polyTopoChange::checkFace(const face&, const label"
                ", const label, const label, const label)"
            )   << "Face uses removed vertices"
                << nl
                << "f:" << f
                << " faceI(-1 if added face):" << faceI
                << " own:" << own << " nei:" << nei
                << " patchI:" << patchI << nl;
            if (hasValidPoints(f))
            {
                FatalError
                        << "points (removed points marked with "
                        << vector::max << ") : " << facePoints(f);
            }
            FatalError << abort(FatalError);
        }
    }
}


void Foam::polyTopoChange::makeCells
(
    const label nActiveFaces,
    labelList& cellFaces,
    labelList& cellFaceOffsets
) const
{
    cellFaces.setSize(2*nActiveFaces);
    cellFaceOffsets.setSize(cellMap_.size() + 1);

    // Faces per cell
    labelList nNbrs(cellMap_.size(), 0);

    // 1. Count faces per cell

    for (label faceI = 0; faceI < nActiveFaces; faceI++)
    {
        if (faceOwner_[faceI] < 0)
        {
            FatalErrorIn
            (
                "polyTopoChange::makeCells\n"
                "(\n"
                "    const label,\n"
                "    labelList&,\n"
                "    labelList&\n"
                ") const\n"
            )   << "Face " << faceI << " is active but its owner has"
                << " been deleted. This is usually due to deleting cells"
                << " without modifying exposed faces to be boundary faces."
                << exit(FatalError);
        }
        nNbrs[faceOwner_[faceI]]++;
    }
    for (label faceI = 0; faceI < nActiveFaces; faceI++)
    {
        if (faceNeighbour_[faceI] >= 0)
        {
            nNbrs[faceNeighbour_[faceI]]++;
        }
    }

    // 2. Calculate offsets

    cellFaceOffsets[0] = 0;
    forAll(nNbrs, cellI)
    {
        cellFaceOffsets[cellI+1] = cellFaceOffsets[cellI] + nNbrs[cellI];
    }

    // 3. Fill faces per cell

    // reset the whole list to use as counter
    nNbrs = 0;

    for (label faceI = 0; faceI < nActiveFaces; faceI++)
    {
        label cellI = faceOwner_[faceI];

        cellFaces[cellFaceOffsets[cellI] + nNbrs[cellI]++] = faceI;
    }

    for (label faceI = 0; faceI < nActiveFaces; faceI++)
    {
        label cellI = faceNeighbour_[faceI];

        if (cellI >= 0)
        {
            cellFaces[cellFaceOffsets[cellI] + nNbrs[cellI]++] = faceI;
        }
    }

    // Last offset points to beyond end of cellFaces.
    cellFaces.setSize(cellFaceOffsets[cellMap_.size()]);
}


// Create cell-cell addressing. Called after compaction (but before ordering)
// of faces
void Foam::polyTopoChange::makeCellCells
(
    const label nActiveFaces,
    CompactListList<label>& cellCells
) const
{
    // Neighbours per cell
    labelList nNbrs(cellMap_.size(), 0);

    // 1. Count neighbours (through internal faces) per cell

    for (label faceI = 0; faceI < nActiveFaces; faceI++)
    {
        if (faceNeighbour_[faceI] >= 0)
        {
            nNbrs[faceOwner_[faceI]]++;
            nNbrs[faceNeighbour_[faceI]]++;
        }
    }

    // 2. Construct csr
    cellCells.setSize(nNbrs);


    // 3. Fill faces per cell

    // reset the whole list to use as counter
    nNbrs = 0;

    for (label faceI = 0; faceI < nActiveFaces; faceI++)
    {
        label nei = faceNeighbour_[faceI];

        if (nei >= 0)
        {
            label own = faceOwner_[faceI];
            cellCells.m()[cellCells.index(own, nNbrs[own]++)] = nei;
            cellCells.m()[cellCells.index(nei, nNbrs[nei]++)] = own;
        }
    }
}


// Cell ordering (based on bandCompression).
// Handles removed cells. Returns number of remaining cells.
Foam::label Foam::polyTopoChange::getCellOrder
(
    const CompactListList<label>& cellCellAddressing,
    labelList& oldToNew
) const
{
    labelList newOrder(cellCellAddressing.size());

    // Fifo buffer for string of cells
    SLList<label> nextCell;

    // Whether cell has been done already
    PackedBoolList visited(cellCellAddressing.size());

    label cellInOrder = 0;


    // Work arrays. Kept outside of loop to minimise allocations.
    // - neighbour cells
    DynamicList<label> nbrs;
    // - corresponding weights
    DynamicList<label> weights;

    // - ordering
    labelList order;


    while (true)
    {
        // For a disconnected region find the lowest connected cell.

        label currentCell = -1;
        label minWeight = labelMax;

        forAll(visited, cellI)
        {
            // find the lowest connected cell that has not been visited yet
            if (!cellRemoved(cellI) && !visited[cellI])
            {
                if (cellCellAddressing[cellI].size() < minWeight)
                {
                    minWeight = cellCellAddressing[cellI].size();
                    currentCell = cellI;
                }
            }
        }


        if (currentCell == -1)
        {
            break;
        }


        // Starting from currentCell walk breadth-first


        // use this cell as a start
        nextCell.append(currentCell);

        // loop through the nextCell list. Add the first cell into the
        // cell order if it has not already been visited and ask for its
        // neighbours. If the neighbour in question has not been visited,
        // add it to the end of the nextCell list

        while (nextCell.size())
        {
            currentCell = nextCell.removeHead();

            if (!visited[currentCell])
            {
                visited[currentCell] = 1;

                // add into cellOrder
                newOrder[cellInOrder] = currentCell;
                cellInOrder++;

                // find if the neighbours have been visited
                const labelUList neighbours = cellCellAddressing[currentCell];

                // Add in increasing order of connectivity

                // 1. Count neighbours of unvisited neighbours
                nbrs.clear();
                weights.clear();

                forAll(neighbours, nI)
                {
                    label nbr = neighbours[nI];
                    if (!cellRemoved(nbr) && !visited[nbr])
                    {
                        // not visited, add to the list
                        nbrs.append(nbr);
                        weights.append(cellCellAddressing[nbr].size());
                    }
                }
                // 2. Sort
                sortedOrder(weights, order);
                // 3. Add in sorted order
                forAll(order, i)
                {
                    nextCell.append(nbrs[i]);
                }
            }
        }
    }

    // Now we have new-to-old in newOrder.
    newOrder.setSize(cellInOrder);

    // Invert to get old-to-new. Make sure removed (i.e. unmapped) cells are -1.
    oldToNew = invert(cellCellAddressing.size(), newOrder);

    return cellInOrder;
}


// Determine order for faces:
// - upper-triangular order for internal faces
// - external faces after internal faces and in patch order.
void Foam::polyTopoChange::getFaceOrder
(
    const label nActiveFaces,
    const labelList& cellFaces,
    const labelList& cellFaceOffsets,

    labelList& oldToNew,
    labelList& patchSizes,
    labelList& patchStarts
) const
{
    oldToNew.setSize(faceOwner_.size());
    oldToNew = -1;

    // First unassigned face
    label newFaceI = 0;

    labelList nbr;
    labelList order;

    forAll(cellMap_, cellI)
    {
        label startOfCell = cellFaceOffsets[cellI];
        label nFaces = cellFaceOffsets[cellI+1] - startOfCell;

        // Neighbouring cells
        //SortableList<label> nbr(nFaces);
        nbr.setSize(nFaces);

        for (label i = 0; i < nFaces; i++)
        {
            label faceI = cellFaces[startOfCell + i];

            label nbrCellI = faceNeighbour_[faceI];

            if (faceI >= nActiveFaces)
            {
                // Retired face.
                nbr[i] = -1;
            }
            else if (nbrCellI != -1)
            {
                // Internal face. Get cell on other side.
                if (nbrCellI == cellI)
                {
                    nbrCellI = faceOwner_[faceI];
                }

                if (cellI < nbrCellI)
                {
                    // CellI is master
                    nbr[i] = nbrCellI;
                }
                else
                {
                    // nbrCell is master. Let it handle this face.
                    nbr[i] = -1;
                }
            }
            else
            {
                // External face. Do later.
                nbr[i] = -1;
            }
        }

        //nbr.sort();
        order.setSize(nFaces);
        sortedOrder(nbr, order);

        //forAll(nbr, i)
        //{
        //    if (nbr[i] != -1)
        //    {
        //        oldToNew[cellFaces[startOfCell + nbr.indices()[i]]] =
        //            newFaceI++;
        //    }
        //}
        forAll(order, i)
        {
            label index = order[i];
            if (nbr[index] != -1)
            {
                oldToNew[cellFaces[startOfCell + index]] = newFaceI++;
            }
        }
    }


    // Pick up all patch faces in patch face order.
    patchStarts.setSize(nPatches_);
    patchStarts = 0;
    patchSizes.setSize(nPatches_);
    patchSizes = 0;

    if (nPatches_ > 0)
    {
        patchStarts[0] = newFaceI;

        for (label faceI = 0; faceI < nActiveFaces; faceI++)
        {
            if (region_[faceI] >= 0)
            {
                patchSizes[region_[faceI]]++;
            }
        }

        label faceI = patchStarts[0];

        forAll(patchStarts, patchI)
        {
            patchStarts[patchI] = faceI;
            faceI += patchSizes[patchI];
        }
    }

    //if (debug)
    //{
    //    Pout<< "patchSizes:" << patchSizes << nl
    //        << "patchStarts:" << patchStarts << endl;
    //}

    labelList workPatchStarts(patchStarts);

    for (label faceI = 0; faceI < nActiveFaces; faceI++)
    {
        if (region_[faceI] >= 0)
        {
            oldToNew[faceI] = workPatchStarts[region_[faceI]]++;
        }
    }

    // Retired faces.
    for (label faceI = nActiveFaces; faceI < oldToNew.size(); faceI++)
    {
        oldToNew[faceI] = faceI;
    }

    // Check done all faces.
    forAll(oldToNew, faceI)
    {
        if (oldToNew[faceI] == -1)
        {
            FatalErrorIn
            (
                "polyTopoChange::getFaceOrder"
                "(const label, const labelList&, const labelList&)"
                " const"
            )   << "Did not determine new position"
                << " for face " << faceI
                << " owner " << faceOwner_[faceI]
                << " neighbour " << faceNeighbour_[faceI]
                << " region " << region_[faceI] << endl
                << "This is usually caused by not specifying a patch for"
                << " a boundary face." << nl
                << "Switch on the polyTopoChange::debug flag to catch"
                << " this error earlier." << nl;
            if (hasValidPoints(faces_[faceI]))
            {
                FatalError
                        << "points (removed points marked with "
                        << vector::max << ") " << facePoints(faces_[faceI]);
            }
            FatalError << abort(FatalError);
        }
    }
}


// Reorder and compact faces according to map.
void Foam::polyTopoChange::reorderCompactFaces
(
    const label newSize,
    const labelList& oldToNew
)
{
    reorder(oldToNew, faces_);
    faces_.setCapacity(newSize);

    reorder(oldToNew, region_);
    region_.setCapacity(newSize);

    reorder(oldToNew, faceOwner_);
    faceOwner_.setCapacity(newSize);

    reorder(oldToNew, faceNeighbour_);
    faceNeighbour_.setCapacity(newSize);

    // Update faceMaps.
    reorder(oldToNew, faceMap_);
    faceMap_.setCapacity(newSize);

    renumberReverseMap(oldToNew, reverseFaceMap_);

    renumberKey(oldToNew, faceFromPoint_);
    renumberKey(oldToNew, faceFromEdge_);
    inplaceReorder(oldToNew, flipFaceFlux_);
    flipFaceFlux_.setCapacity(newSize);
    renumberKey(oldToNew, faceZone_);
    inplaceReorder(oldToNew, faceZoneFlip_);
    faceZoneFlip_.setCapacity(newSize);
}


// Compact all and orders points and faces:
// - points into internal followed by external points
// - internalfaces upper-triangular
// - externalfaces after internal ones.
void Foam::polyTopoChange::compact
(
    const bool orderCells,
    const bool orderPoints,
    label& nInternalPoints,
    labelList& patchSizes,
    labelList& patchStarts
)
{
    points_.shrink();
    pointMap_.shrink();
    reversePointMap_.shrink();

    faces_.shrink();
    region_.shrink();
    faceOwner_.shrink();
    faceNeighbour_.shrink();
    faceMap_.shrink();
    reverseFaceMap_.shrink();

    cellMap_.shrink();
    reverseCellMap_.shrink();
    cellZone_.shrink();


    // Compact points
    label nActivePoints = 0;
    {
        labelList localPointMap(points_.size(), -1);
        label newPointI = 0;

        if (!orderPoints)
        {
            nInternalPoints = -1;

            forAll(points_, pointI)
            {
                if (!pointRemoved(pointI) && !retiredPoints_.found(pointI))
                {
                    localPointMap[pointI] = newPointI++;
                }
            }
            nActivePoints = newPointI;
        }
        else
        {
            forAll(points_, pointI)
            {
                if (!pointRemoved(pointI) && !retiredPoints_.found(pointI))
                {
                    nActivePoints++;
                }
            }

            // Mark boundary points
            forAll(faceOwner_, faceI)
            {
                if
                (
                   !faceRemoved(faceI)
                 && faceOwner_[faceI] >= 0
                 && faceNeighbour_[faceI] < 0
                )
                {
                    // Valid boundary face
                    const face& f = faces_[faceI];

                    forAll(f, fp)
                    {
                        label pointI = f[fp];

                        if (localPointMap[pointI] == -1)
                        {
                            if
                            (
                                pointRemoved(pointI)
                             || retiredPoints_.found(pointI)
                            )
                            {
                                FatalErrorIn("polyTopoChange::compact(..)")
                                    << "Removed or retired point " << pointI
                                    << " in face " << f
                                    << " at position " << faceI << endl
                                    << "Probably face has not been adapted for"
                                    << " removed points." << abort(FatalError);
                            }
                            localPointMap[pointI] = newPointI++;
                        }
                    }
                }
            }

            label nBoundaryPoints = newPointI;
            nInternalPoints = nActivePoints - nBoundaryPoints;

            // Move the boundary addressing up
            forAll(localPointMap, pointI)
            {
                if (localPointMap[pointI] != -1)
                {
                    localPointMap[pointI] += nInternalPoints;
                }
            }

            newPointI = 0;

            // Mark internal points
            forAll(faceOwner_, faceI)
            {
                if
                (
                   !faceRemoved(faceI)
                 && faceOwner_[faceI] >= 0
                 && faceNeighbour_[faceI] >= 0
                )
                {
                    // Valid internal face
                    const face& f = faces_[faceI];

                    forAll(f, fp)
                    {
                        label pointI = f[fp];

                        if (localPointMap[pointI] == -1)
                        {
                            if
                            (
                                pointRemoved(pointI)
                             || retiredPoints_.found(pointI)
                            )
                            {
                                FatalErrorIn("polyTopoChange::compact(..)")
                                    << "Removed or retired point " << pointI
                                    << " in face " << f
                                    << " at position " << faceI << endl
                                    << "Probably face has not been adapted for"
                                    << " removed points." << abort(FatalError);
                            }
                            localPointMap[pointI] = newPointI++;
                        }
                    }
                }
            }

            if (newPointI != nInternalPoints)
            {
                FatalErrorIn("polyTopoChange::compact(..)")
                    << "Problem." << abort(FatalError);
            }
            newPointI = nActivePoints;
        }

        forAllConstIter(labelHashSet, retiredPoints_, iter)
        {
            localPointMap[iter.key()] = newPointI++;
        }


        if (debug)
        {
            Pout<< "Points : active:" << nActivePoints
                << "  removed:" << points_.size()-newPointI << endl;
        }

        reorder(localPointMap, points_);
        points_.setCapacity(newPointI);

        // Update pointMaps
        reorder(localPointMap, pointMap_);
        pointMap_.setCapacity(newPointI);
        renumberReverseMap(localPointMap, reversePointMap_);

        renumberKey(localPointMap, pointZone_);
        renumber(localPointMap, retiredPoints_);

        // Use map to relabel face vertices
        forAll(faces_, faceI)
        {
            face& f = faces_[faceI];

            //labelList oldF(f);
            renumberCompact(localPointMap, f);

            if (!faceRemoved(faceI) && f.size() < 3)
            {
                FatalErrorIn("polyTopoChange::compact(..)")
                    << "Created illegal face " << f
                    //<< " from face " << oldF
                    << " at position:" << faceI
                    << " when filtering removed points"
                    << abort(FatalError);
            }
        }
    }


    // Compact faces.
    {
        labelList localFaceMap(faces_.size(), -1);
        label newFaceI = 0;

        forAll(faces_, faceI)
        {
            if (!faceRemoved(faceI) && faceOwner_[faceI] >= 0)
            {
                localFaceMap[faceI] = newFaceI++;
            }
        }
        nActiveFaces_ = newFaceI;

        forAll(faces_, faceI)
        {
            if (!faceRemoved(faceI) && faceOwner_[faceI] < 0)
            {
                // Retired face
                localFaceMap[faceI] = newFaceI++;
            }
        }

        if (debug)
        {
            Pout<< "Faces : active:" << nActiveFaces_
                << "  removed:" << faces_.size()-newFaceI << endl;
        }

        // Reorder faces.
        reorderCompactFaces(newFaceI, localFaceMap);
    }

    // Compact cells.
    {
        labelList localCellMap;
        label newCellI;

        if (orderCells)
        {
            // Construct cellCell addressing
            CompactListList<label> cellCells;
            makeCellCells(nActiveFaces_, cellCells);

            // Cell ordering (based on bandCompression). Handles removed cells.
            newCellI = getCellOrder(cellCells, localCellMap);
        }
        else
        {
            // Compact out removed cells
            localCellMap.setSize(cellMap_.size());
            localCellMap = -1;

            newCellI = 0;
            forAll(cellMap_, cellI)
            {
                if (!cellRemoved(cellI))
                {
                    localCellMap[cellI] = newCellI++;
                }
            }
        }

        if (debug)
        {
            Pout<< "Cells : active:" << newCellI
                << "  removed:" << cellMap_.size()-newCellI << endl;
        }

        // Renumber -if cells reordered or -if cells removed
        if (orderCells || (newCellI != cellMap_.size()))
        {
            reorder(localCellMap, cellMap_);
            cellMap_.setCapacity(newCellI);
            renumberReverseMap(localCellMap, reverseCellMap_);

            reorder(localCellMap, cellZone_);
            cellZone_.setCapacity(newCellI);

            renumberKey(localCellMap, cellFromPoint_);
            renumberKey(localCellMap, cellFromEdge_);
            renumberKey(localCellMap, cellFromFace_);

            // Renumber owner/neighbour. Take into account if neighbour suddenly
            // gets lower cell than owner.
            forAll(faceOwner_, faceI)
            {
                label own = faceOwner_[faceI];
                label nei = faceNeighbour_[faceI];

                if (own >= 0)
                {
                    // Update owner
                    faceOwner_[faceI] = localCellMap[own];

                    if (nei >= 0)
                    {
                        // Update neighbour.
                        faceNeighbour_[faceI] = localCellMap[nei];

                        // Check if face needs reversing.
                        if
                        (
                            faceNeighbour_[faceI] >= 0
                         && faceNeighbour_[faceI] < faceOwner_[faceI]
                        )
                        {
                            faces_[faceI].flip();
                            Swap(faceOwner_[faceI], faceNeighbour_[faceI]);
                            flipFaceFlux_[faceI] =
                            (
                                flipFaceFlux_[faceI]
                              ? 0
                              : 1
                            );
                            faceZoneFlip_[faceI] =
                            (
                                faceZoneFlip_[faceI]
                              ? 0
                              : 1
                            );
                        }
                    }
                }
                else if (nei >= 0)
                {
                    // Update neighbour.
                    faceNeighbour_[faceI] = localCellMap[nei];
                }
            }
        }
    }

    // Reorder faces into upper-triangular and patch ordering
    {
        // Create cells (packed storage)
        labelList cellFaces;
        labelList cellFaceOffsets;
        makeCells(nActiveFaces_, cellFaces, cellFaceOffsets);

        // Do upper triangular order and patch sorting
        labelList localFaceMap;
        getFaceOrder
        (
            nActiveFaces_,
            cellFaces,
            cellFaceOffsets,

            localFaceMap,
            patchSizes,
            patchStarts
        );

        // Reorder faces.
        reorderCompactFaces(localFaceMap.size(), localFaceMap);
    }
}


// Find faces to interpolate to create value for new face. Only used if
// face was inflated from edge or point. Internal faces should only be
// created from internal faces, external faces only from external faces
// (and ideally the same patch)
// Is bit problematic if there are no faces to select, i.e. in polyDualMesh
// an internal face can be created from a boundary edge with no internal
// faces connected to it.
Foam::labelList Foam::polyTopoChange::selectFaces
(
    const primitiveMesh& mesh,
    const labelList& faceLabels,
    const bool internalFacesOnly
)
{
    label nFaces = 0;

    forAll(faceLabels, i)
    {
        label faceI = faceLabels[i];

        if (internalFacesOnly == mesh.isInternalFace(faceI))
        {
            nFaces++;
        }
    }

    labelList collectedFaces;

    if (nFaces == 0)
    {
        // Did not find any faces of the correct type so just use any old
        // face.
        collectedFaces = faceLabels;
    }
    else
    {
        collectedFaces.setSize(nFaces);

        nFaces = 0;

        forAll(faceLabels, i)
        {
            label faceI = faceLabels[i];

            if (internalFacesOnly == mesh.isInternalFace(faceI))
            {
                collectedFaces[nFaces++] = faceI;
            }
        }
    }

    return collectedFaces;
}


// Calculate pointMap per patch (so from patch point label to old patch point
// label)
void Foam::polyTopoChange::calcPatchPointMap
(
    const List<Map<label> >& oldPatchMeshPointMaps,
    const polyBoundaryMesh& boundary,
    labelListList& patchPointMap
) const
{
    patchPointMap.setSize(boundary.size());

    forAll(boundary, patchI)
    {
        const labelList& meshPoints = boundary[patchI].meshPoints();

        const Map<label>& oldMeshPointMap = oldPatchMeshPointMaps[patchI];

        labelList& curPatchPointRnb = patchPointMap[patchI];

        curPatchPointRnb.setSize(meshPoints.size());

        forAll(meshPoints, i)
        {
            if (meshPoints[i] < pointMap_.size())
            {
                // Check if old point was part of same patch
                Map<label>::const_iterator ozmpmIter = oldMeshPointMap.find
                (
                    pointMap_[meshPoints[i]]
                );

                if (ozmpmIter != oldMeshPointMap.end())
                {
                    curPatchPointRnb[i] = ozmpmIter();
                }
                else
                {
                    curPatchPointRnb[i] = -1;
                }
            }
            else
            {
                curPatchPointRnb[i] = -1;
            }
        }
    }
}


void Foam::polyTopoChange::calcFaceInflationMaps
(
    const polyMesh& mesh,
    List<objectMap>& facesFromPoints,
    List<objectMap>& facesFromEdges,
    List<objectMap>& facesFromFaces
) const
{
    // Faces inflated from points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    facesFromPoints.setSize(faceFromPoint_.size());

    if (faceFromPoint_.size())
    {
        label nFacesFromPoints = 0;

        // Collect all still existing faces connected to this point.
        forAllConstIter(Map<label>, faceFromPoint_, iter)
        {
            label newFaceI = iter.key();

            if (region_[newFaceI] == -1)
            {
                // Get internal faces using point on old mesh
                facesFromPoints[nFacesFromPoints++] = objectMap
                (
                    newFaceI,
                    selectFaces
                    (
                        mesh,
                        mesh.pointFaces()[iter()],
                        true
                    )
                );
            }
            else
            {
                // Get patch faces using point on old mesh
                facesFromPoints[nFacesFromPoints++] = objectMap
                (
                    newFaceI,
                    selectFaces
                    (
                        mesh,
                        mesh.pointFaces()[iter()],
                        false
                    )
                );
            }
        }
    }


    // Faces inflated from edges
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    facesFromEdges.setSize(faceFromEdge_.size());

    if (faceFromEdge_.size())
    {
        label nFacesFromEdges = 0;

        // Collect all still existing faces connected to this edge.
        forAllConstIter(Map<label>, faceFromEdge_, iter)
        {
            label newFaceI = iter.key();

            if (region_[newFaceI] == -1)
            {
                // Get internal faces using edge on old mesh
                facesFromEdges[nFacesFromEdges++] = objectMap
                (
                    newFaceI,
                    selectFaces
                    (
                        mesh,
                        mesh.edgeFaces(iter()),
                        true
                    )
                );
            }
            else
            {
                // Get patch faces using edge on old mesh
                facesFromEdges[nFacesFromEdges++] = objectMap
                (
                    newFaceI,
                    selectFaces
                    (
                        mesh,
                        mesh.edgeFaces(iter()),
                        false
                    )
                );
            }
        }
    }


    // Faces from face merging
    // ~~~~~~~~~~~~~~~~~~~~~~~

    getMergeSets
    (
        reverseFaceMap_,
        faceMap_,
        facesFromFaces
    );
}


void Foam::polyTopoChange::calcCellInflationMaps
(
    const polyMesh& mesh,
    List<objectMap>& cellsFromPoints,
    List<objectMap>& cellsFromEdges,
    List<objectMap>& cellsFromFaces,
    List<objectMap>& cellsFromCells
) const
{
    cellsFromPoints.setSize(cellFromPoint_.size());

    if (cellFromPoint_.size())
    {
        label nCellsFromPoints = 0;

        // Collect all still existing faces connected to this point.
        forAllConstIter(Map<label>, cellFromPoint_, iter)
        {
            cellsFromPoints[nCellsFromPoints++] = objectMap
            (
                iter.key(),
                mesh.pointCells()[iter()]
            );
        }
    }


    cellsFromEdges.setSize(cellFromEdge_.size());

    if (cellFromEdge_.size())
    {
        label nCellsFromEdges = 0;

        // Collect all still existing faces connected to this point.
        forAllConstIter(Map<label>, cellFromEdge_, iter)
        {
            cellsFromEdges[nCellsFromEdges++] = objectMap
            (
                iter.key(),
                mesh.edgeCells()[iter()]
            );
        }
    }


    cellsFromFaces.setSize(cellFromFace_.size());

    if (cellFromFace_.size())
    {
        label nCellsFromFaces = 0;

        labelList twoCells(2);

        // Collect all still existing faces connected to this point.
        forAllConstIter(Map<label>, cellFromFace_, iter)
        {
            label oldFaceI = iter();

            if (mesh.isInternalFace(oldFaceI))
            {
                twoCells[0] = mesh.faceOwner()[oldFaceI];
                twoCells[1] = mesh.faceNeighbour()[oldFaceI];
                cellsFromFaces[nCellsFromFaces++] = objectMap
                (
                    iter.key(),
                    twoCells
                );
            }
            else
            {
                cellsFromFaces[nCellsFromFaces++] = objectMap
                (
                    iter.key(),
                    labelList(1, mesh.faceOwner()[oldFaceI])
                );
            }
        }
    }


    // Cells from cell merging
    // ~~~~~~~~~~~~~~~~~~~~~~~

    getMergeSets
    (
        reverseCellMap_,
        cellMap_,
        cellsFromCells
    );
}


void Foam::polyTopoChange::resetZones
(
    const polyMesh& mesh,
    polyMesh& newMesh,
    labelListList& pointZoneMap,
    labelListList& faceZoneFaceMap,
    labelListList& cellZoneMap
) const
{
    // pointZones
    // ~~~~~~~~~~

    pointZoneMap.setSize(mesh.pointZones().size());
    {
        const pointZoneMesh& pointZones = mesh.pointZones();

        // Count points per zone

        labelList nPoints(pointZones.size(), 0);

        forAllConstIter(Map<label>, pointZone_, iter)
        {
            label zoneI = iter();

            if (zoneI < 0 || zoneI >= pointZones.size())
            {
                FatalErrorIn
                (
                    "resetZones(const polyMesh&, polyMesh&, labelListList&"
                    "labelListList&, labelListList&)"
                )   << "Illegal zoneID " << zoneI << " for point "
                    << iter.key() << " coord " << mesh.points()[iter.key()]
                    << abort(FatalError);
            }
            nPoints[zoneI]++;
        }

        // Distribute points per zone

        labelListList addressing(pointZones.size());
        forAll(addressing, zoneI)
        {
            addressing[zoneI].setSize(nPoints[zoneI]);
        }
        nPoints = 0;

        forAllConstIter(Map<label>, pointZone_, iter)
        {
            label zoneI = iter();

            addressing[zoneI][nPoints[zoneI]++] = iter.key();
        }
        // Sort the addressing
        forAll(addressing, zoneI)
        {
            stableSort(addressing[zoneI]);
        }

        // So now we both have old zones and the new addressing.
        // Invert the addressing to get pointZoneMap.
        forAll(addressing, zoneI)
        {
            const pointZone& oldZone = pointZones[zoneI];
            const labelList& newZoneAddr = addressing[zoneI];

            labelList& curPzRnb = pointZoneMap[zoneI];
            curPzRnb.setSize(newZoneAddr.size());

            forAll(newZoneAddr, i)
            {
                if (newZoneAddr[i] < pointMap_.size())
                {
                    curPzRnb[i] = oldZone.whichPoint(pointMap_[newZoneAddr[i]]);
                }
                else
                {
                    curPzRnb[i] = -1;
                }
            }
        }

        // Reset the addresing on the zone
        newMesh.pointZones().clearAddressing();
        forAll(newMesh.pointZones(), zoneI)
        {
            if (debug)
            {
                Pout<< "pointZone:" << zoneI
                    << "  name:" << newMesh.pointZones()[zoneI].name()
                    << "  size:" << addressing[zoneI].size()
                    << endl;
            }

            newMesh.pointZones()[zoneI] = addressing[zoneI];
        }
    }


    // faceZones
    // ~~~~~~~~~

    faceZoneFaceMap.setSize(mesh.faceZones().size());
    {
        const faceZoneMesh& faceZones = mesh.faceZones();

        labelList nFaces(faceZones.size(), 0);

        forAllConstIter(Map<label>, faceZone_, iter)
        {
            label zoneI = iter();

            if (zoneI < 0 || zoneI >= faceZones.size())
            {
                FatalErrorIn
                (
                    "resetZones(const polyMesh&, polyMesh&, labelListList&"
                    "labelListList&, labelListList&)"
                )   << "Illegal zoneID " << zoneI << " for face "
                    << iter.key()
                    << abort(FatalError);
            }
            nFaces[zoneI]++;
        }

        labelListList addressing(faceZones.size());
        boolListList flipMode(faceZones.size());

        forAll(addressing, zoneI)
        {
            addressing[zoneI].setSize(nFaces[zoneI]);
            flipMode[zoneI].setSize(nFaces[zoneI]);
        }
        nFaces = 0;

        forAllConstIter(Map<label>, faceZone_, iter)
        {
            label zoneI = iter();
            label faceI = iter.key();

            label index = nFaces[zoneI]++;

            addressing[zoneI][index] = faceI;
            flipMode[zoneI][index] = faceZoneFlip_[faceI];
        }
        // Sort the addressing
        forAll(addressing, zoneI)
        {
            labelList newToOld;
            sortedOrder(addressing[zoneI], newToOld);
            {
                labelList newAddressing(addressing[zoneI].size());
                forAll(newAddressing, i)
                {
                    newAddressing[i] = addressing[zoneI][newToOld[i]];
                }
                addressing[zoneI].transfer(newAddressing);
            }
            {
                boolList newFlipMode(flipMode[zoneI].size());
                forAll(newFlipMode, i)
                {
                    newFlipMode[i] = flipMode[zoneI][newToOld[i]];
                }
                flipMode[zoneI].transfer(newFlipMode);
            }
        }

        // So now we both have old zones and the new addressing.
        // Invert the addressing to get faceZoneFaceMap.
        forAll(addressing, zoneI)
        {
            const faceZone& oldZone = faceZones[zoneI];
            const labelList& newZoneAddr = addressing[zoneI];

            labelList& curFzFaceRnb = faceZoneFaceMap[zoneI];

            curFzFaceRnb.setSize(newZoneAddr.size());

            forAll(newZoneAddr, i)
            {
                if (newZoneAddr[i] < faceMap_.size())
                {
                    curFzFaceRnb[i] =
                        oldZone.whichFace(faceMap_[newZoneAddr[i]]);
                }
                else
                {
                    curFzFaceRnb[i] = -1;
                }
            }
        }


        // Reset the addresing on the zone
        newMesh.faceZones().clearAddressing();
        forAll(newMesh.faceZones(), zoneI)
        {
            if (debug)
            {
                Pout<< "faceZone:" << zoneI
                    << "  name:" << newMesh.faceZones()[zoneI].name()
                    << "  size:" << addressing[zoneI].size()
                    << endl;
            }

            newMesh.faceZones()[zoneI].resetAddressing
            (
                addressing[zoneI],
                flipMode[zoneI]
            );
        }
    }


    // cellZones
    // ~~~~~~~~~

    cellZoneMap.setSize(mesh.cellZones().size());
    {
        const cellZoneMesh& cellZones = mesh.cellZones();

        labelList nCells(cellZones.size(), 0);

        forAll(cellZone_, cellI)
        {
            label zoneI = cellZone_[cellI];

            if (zoneI >= cellZones.size())
            {
                FatalErrorIn
                (
                    "resetZones(const polyMesh&, polyMesh&, labelListList&"
                    "labelListList&, labelListList&)"
                )   << "Illegal zoneID " << zoneI << " for cell "
                    << cellI << abort(FatalError);
            }

            if (zoneI >= 0)
            {
                nCells[zoneI]++;
            }
        }

        labelListList addressing(cellZones.size());
        forAll(addressing, zoneI)
        {
            addressing[zoneI].setSize(nCells[zoneI]);
        }
        nCells = 0;

        forAll(cellZone_, cellI)
        {
            label zoneI = cellZone_[cellI];

            if (zoneI >= 0)
            {
                addressing[zoneI][nCells[zoneI]++] = cellI;
            }
        }
        // Sort the addressing
        forAll(addressing, zoneI)
        {
            stableSort(addressing[zoneI]);
        }

        // So now we both have old zones and the new addressing.
        // Invert the addressing to get cellZoneMap.
        forAll(addressing, zoneI)
        {
            const cellZone& oldZone = cellZones[zoneI];
            const labelList& newZoneAddr = addressing[zoneI];

            labelList& curCellRnb = cellZoneMap[zoneI];

            curCellRnb.setSize(newZoneAddr.size());

            forAll(newZoneAddr, i)
            {
                if (newZoneAddr[i] < cellMap_.size())
                {
                    curCellRnb[i] =
                        oldZone.whichCell(cellMap_[newZoneAddr[i]]);
                }
                else
                {
                    curCellRnb[i] = -1;
                }
            }
        }

        // Reset the addresing on the zone
        newMesh.cellZones().clearAddressing();
        forAll(newMesh.cellZones(), zoneI)
        {
            if (debug)
            {
                Pout<< "cellZone:" << zoneI
                    << "  name:" << newMesh.cellZones()[zoneI].name()
                    << "  size:" << addressing[zoneI].size()
                    << endl;
            }

            newMesh.cellZones()[zoneI] = addressing[zoneI];
        }
    }
}


void Foam::polyTopoChange::calcFaceZonePointMap
(
    const polyMesh& mesh,
    const List<Map<label> >& oldFaceZoneMeshPointMaps,
    labelListList& faceZonePointMap
) const
{
    const faceZoneMesh& faceZones = mesh.faceZones();

    faceZonePointMap.setSize(faceZones.size());

    forAll(faceZones, zoneI)
    {
        const faceZone& newZone = faceZones[zoneI];

        const labelList& newZoneMeshPoints = newZone().meshPoints();

        const Map<label>& oldZoneMeshPointMap = oldFaceZoneMeshPointMaps[zoneI];

        labelList& curFzPointRnb = faceZonePointMap[zoneI];

        curFzPointRnb.setSize(newZoneMeshPoints.size());

        forAll(newZoneMeshPoints, pointI)
        {
            if (newZoneMeshPoints[pointI] < pointMap_.size())
            {
                Map<label>::const_iterator ozmpmIter =
                    oldZoneMeshPointMap.find
                    (
                        pointMap_[newZoneMeshPoints[pointI]]
                    );

                if (ozmpmIter != oldZoneMeshPointMap.end())
                {
                    curFzPointRnb[pointI] = ozmpmIter();
                }
                else
                {
                    curFzPointRnb[pointI] = -1;
                }
            }
            else
            {
                curFzPointRnb[pointI] = -1;
            }
        }
    }
}


void Foam::polyTopoChange::reorderCoupledFaces
(
    const bool syncParallel,
    const polyBoundaryMesh& boundary,
    const labelList& patchStarts,
    const labelList& patchSizes,
    const pointField& points
)
{
    // Mapping for faces (old to new). Extends over all mesh faces for
    // convenience (could be just the external faces)
    labelList oldToNew(identity(faces_.size()));

    // Rotation on new faces.
    labelList rotation(faces_.size(), 0);

    PstreamBuffers pBufs(Pstream::nonBlocking);

    // Send ordering
    forAll(boundary, patchI)
    {
        if (syncParallel || !isA<processorPolyPatch>(boundary[patchI]))
        {
            boundary[patchI].initOrder
            (
                pBufs,
                primitivePatch
                (
                    SubList<face>
                    (
                        faces_,
                        patchSizes[patchI],
                        patchStarts[patchI]
                    ),
                    points
                )
            );
        }
    }

    if (syncParallel)
    {
        pBufs.finishedSends();
    }

    // Receive and calculate ordering

    bool anyChanged = false;

    forAll(boundary, patchI)
    {
        if (syncParallel || !isA<processorPolyPatch>(boundary[patchI]))
        {
            labelList patchFaceMap(patchSizes[patchI], -1);
            labelList patchFaceRotation(patchSizes[patchI], 0);

            bool changed = boundary[patchI].order
            (
                pBufs,
                primitivePatch
                (
                    SubList<face>
                    (
                        faces_,
                        patchSizes[patchI],
                        patchStarts[patchI]
                    ),
                    points
                ),
                patchFaceMap,
                patchFaceRotation
            );

            if (changed)
            {
                // Merge patch face reordering into mesh face reordering table
                label start = patchStarts[patchI];

                forAll(patchFaceMap, patchFaceI)
                {
                    oldToNew[patchFaceI + start] =
                        start + patchFaceMap[patchFaceI];
                }

                forAll(patchFaceRotation, patchFaceI)
                {
                    rotation[patchFaceI + start] =
                        patchFaceRotation[patchFaceI];
                }

                anyChanged = true;
            }
        }
    }

    if (syncParallel)
    {
        reduce(anyChanged, orOp<bool>());
    }

    if (anyChanged)
    {
        // Reorder faces according to oldToNew.
        reorderCompactFaces(oldToNew.size(), oldToNew);

        // Rotate faces (rotation is already in new face indices).
        forAll(rotation, faceI)
        {
            if (rotation[faceI] != 0)
            {
                inplaceRotateList<List, label>(faces_[faceI], rotation[faceI]);
            }
        }
    }
}


void Foam::polyTopoChange::compactAndReorder
(
    const polyMesh& mesh,
    const bool syncParallel,
    const bool orderCells,
    const bool orderPoints,

    label& nInternalPoints,
    pointField& newPoints,
    labelList& patchSizes,
    labelList& patchStarts,
    List<objectMap>& pointsFromPoints,
    List<objectMap>& facesFromPoints,
    List<objectMap>& facesFromEdges,
    List<objectMap>& facesFromFaces,
    List<objectMap>& cellsFromPoints,
    List<objectMap>& cellsFromEdges,
    List<objectMap>& cellsFromFaces,
    List<objectMap>& cellsFromCells,
    List<Map<label> >& oldPatchMeshPointMaps,
    labelList& oldPatchNMeshPoints,
    labelList& oldPatchStarts,
    List<Map<label> >& oldFaceZoneMeshPointMaps
)
{
    if (mesh.boundaryMesh().size() != nPatches_)
    {
        FatalErrorIn("polyTopoChange::compactAndReorder(..)")
            << "polyTopoChange was constructed with a mesh with "
            << nPatches_ << " patches." << endl
            << "The mesh now provided has a different number of patches "
            << mesh.boundaryMesh().size()
            << " which is illegal" << endl
            << abort(FatalError);
    }

    // Remove any holes from points/faces/cells and sort faces.
    // Sets nActiveFaces_.
    compact(orderCells, orderPoints, nInternalPoints, patchSizes, patchStarts);

    // Transfer points to pointField. points_ are now cleared!
    // Only done since e.g. reorderCoupledFaces requires pointField.
    newPoints.transfer(points_);

    // Reorder any coupled faces
    reorderCoupledFaces
    (
        syncParallel,
        mesh.boundaryMesh(),
        patchStarts,
        patchSizes,
        newPoints
    );


    // Calculate inflation/merging maps
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // These are for the new face(/point/cell) the old faces whose value
    // needs to be
    // averaged/summed to get the new value. These old faces come either from
    // merged old faces (face remove into other face),
    // the old edgeFaces (inflate from edge) or the old pointFaces (inflate
    // from point). As an additional complexity will use only internal faces
    // to create new value for internal face and vice versa only patch
    // faces to to create patch face value.

    // For point only point merging
    getMergeSets
    (
        reversePointMap_,
        pointMap_,
        pointsFromPoints
    );

    calcFaceInflationMaps
    (
        mesh,
        facesFromPoints,
        facesFromEdges,
        facesFromFaces
    );

    calcCellInflationMaps
    (
        mesh,
        cellsFromPoints,
        cellsFromEdges,
        cellsFromFaces,
        cellsFromCells
    );

    // Clear inflation info
    {
        faceFromPoint_.clearStorage();
        faceFromEdge_.clearStorage();

        cellFromPoint_.clearStorage();
        cellFromEdge_.clearStorage();
        cellFromFace_.clearStorage();
    }


    const polyBoundaryMesh& boundary = mesh.boundaryMesh();

    // Grab patch mesh point maps
    oldPatchMeshPointMaps.setSize(boundary.size());
    oldPatchNMeshPoints.setSize(boundary.size());
    oldPatchStarts.setSize(boundary.size());

    forAll(boundary, patchI)
    {
        // Copy old face zone mesh point maps
        oldPatchMeshPointMaps[patchI] = boundary[patchI].meshPointMap();
        oldPatchNMeshPoints[patchI] = boundary[patchI].meshPoints().size();
        oldPatchStarts[patchI] = boundary[patchI].start();
    }

    // Grab old face zone mesh point maps.
    // These need to be saved before resetting the mesh and are used
    // later on to calculate the faceZone pointMaps.
    oldFaceZoneMeshPointMaps.setSize(mesh.faceZones().size());

    forAll(mesh.faceZones(), zoneI)
    {
        const faceZone& oldZone = mesh.faceZones()[zoneI];

        oldFaceZoneMeshPointMaps[zoneI] = oldZone().meshPointMap();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::polyTopoChange::polyTopoChange(const label nPatches, const bool strict)
:
    strict_(strict),
    nPatches_(nPatches),
    points_(0),
    pointMap_(0),
    reversePointMap_(0),
    pointZone_(0),
    retiredPoints_(0),
    faces_(0),
    region_(0),
    faceOwner_(0),
    faceNeighbour_(0),
    faceMap_(0),
    reverseFaceMap_(0),
    faceFromPoint_(0),
    faceFromEdge_(0),
    flipFaceFlux_(0),
    faceZone_(0),
    faceZoneFlip_(0),
    nActiveFaces_(0),
    cellMap_(0),
    reverseCellMap_(0),
    cellFromPoint_(0),
    cellFromEdge_(0),
    cellFromFace_(0),
    cellZone_(0)
{}


// Construct from components
Foam::polyTopoChange::polyTopoChange
(
    const polyMesh& mesh,
    const bool strict
)
:
    strict_(strict),
    nPatches_(0),
    points_(0),
    pointMap_(0),
    reversePointMap_(0),
    pointZone_(0),
    retiredPoints_(0),
    faces_(0),
    region_(0),
    faceOwner_(0),
    faceNeighbour_(0),
    faceMap_(0),
    reverseFaceMap_(0),
    faceFromPoint_(0),
    faceFromEdge_(0),
    flipFaceFlux_(0),
    faceZone_(0),
    faceZoneFlip_(0),
    nActiveFaces_(0),
    cellMap_(0),
    reverseCellMap_(0),
    cellFromPoint_(0),
    cellFromEdge_(0),
    cellFromFace_(0),
    cellZone_(0)
{
    addMesh
    (
        mesh,
        identity(mesh.boundaryMesh().size()),
        identity(mesh.pointZones().size()),
        identity(mesh.faceZones().size()),
        identity(mesh.cellZones().size())
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyTopoChange::clear()
{
    points_.clearStorage();
    pointMap_.clearStorage();
    reversePointMap_.clearStorage();
    pointZone_.clearStorage();
    retiredPoints_.clearStorage();

    faces_.clearStorage();
    region_.clearStorage();
    faceOwner_.clearStorage();
    faceNeighbour_.clearStorage();
    faceMap_.clearStorage();
    reverseFaceMap_.clearStorage();
    faceFromPoint_.clearStorage();
    faceFromEdge_.clearStorage();
    flipFaceFlux_.clearStorage();
    faceZone_.clearStorage();
    faceZoneFlip_.clearStorage();
    nActiveFaces_ = 0;

    cellMap_.clearStorage();
    reverseCellMap_.clearStorage();
    cellZone_.clearStorage();
    cellFromPoint_.clearStorage();
    cellFromEdge_.clearStorage();
    cellFromFace_.clearStorage();
}


void Foam::polyTopoChange::addMesh
(
    const polyMesh& mesh,
    const labelList& patchMap,
    const labelList& pointZoneMap,
    const labelList& faceZoneMap,
    const labelList& cellZoneMap
)
{
    label maxRegion = nPatches_ - 1;
    forAll(patchMap, i)
    {
        maxRegion = max(maxRegion, patchMap[i]);
    }
    nPatches_ = maxRegion + 1;


    // Add points
    {
        const pointField& points = mesh.points();
        const pointZoneMesh& pointZones = mesh.pointZones();

        // Extend
        points_.setCapacity(points_.size() + points.size());
        pointMap_.setCapacity(pointMap_.size() + points.size());
        reversePointMap_.setCapacity(reversePointMap_.size() + points.size());
        pointZone_.resize(pointZone_.size() + points.size()/100);

        // Precalc offset zones
        labelList newZoneID(points.size(), -1);

        forAll(pointZones, zoneI)
        {
            const labelList& pointLabels = pointZones[zoneI];

            forAll(pointLabels, j)
            {
                newZoneID[pointLabels[j]] = pointZoneMap[zoneI];
            }
        }

        // Add points in mesh order
        for (label pointI = 0; pointI < mesh.nPoints(); pointI++)
        {
            addPoint
            (
                points[pointI],
                pointI,
                newZoneID[pointI],
                true
            );
        }
    }

    // Add cells
    {
        const cellZoneMesh& cellZones = mesh.cellZones();

        // Resize

        // Note: polyMesh does not allow retired cells anymore. So allCells
        // always equals nCells
        label nAllCells = mesh.nCells();

        cellMap_.setCapacity(cellMap_.size() + nAllCells);
        reverseCellMap_.setCapacity(reverseCellMap_.size() + nAllCells);
        cellFromPoint_.resize(cellFromPoint_.size() + nAllCells/100);
        cellFromEdge_.resize(cellFromEdge_.size() + nAllCells/100);
        cellFromFace_.resize(cellFromFace_.size() + nAllCells/100);
        cellZone_.setCapacity(cellZone_.size() + nAllCells);


        // Precalc offset zones
        labelList newZoneID(nAllCells, -1);

        forAll(cellZones, zoneI)
        {
            const labelList& cellLabels = cellZones[zoneI];

            forAll(cellLabels, j)
            {
                label cellI = cellLabels[j];

                if (newZoneID[cellI] != -1)
                {
                    WarningIn
                    (
                        "polyTopoChange::addMesh"
                        "(const polyMesh&, const labelList&,"
                        "const labelList&, const labelList&,"
                        "const labelList&)"
                    )   << "Cell:" << cellI
                        << " centre:" << mesh.cellCentres()[cellI]
                        << " is in two zones:"
                        << cellZones[newZoneID[cellI]].name()
                        << " and " << cellZones[zoneI].name() << endl
                        << "    This is not supported."
                        << " Continuing with first zone only." << endl;
                }
                else
                {
                    newZoneID[cellI] = cellZoneMap[zoneI];
                }
            }
        }

        // Add cells in mesh order
        for (label cellI = 0; cellI < nAllCells; cellI++)
        {
            // Add cell from cell
            addCell(-1, -1, -1, cellI, newZoneID[cellI]);
        }
    }

    // Add faces
    {
        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        const faceList& faces = mesh.faces();
        const labelList& faceOwner = mesh.faceOwner();
        const labelList& faceNeighbour = mesh.faceNeighbour();
        const faceZoneMesh& faceZones = mesh.faceZones();

        // Resize
        label nAllFaces = mesh.faces().size();

        faces_.setCapacity(faces_.size() + nAllFaces);
        region_.setCapacity(region_.size() + nAllFaces);
        faceOwner_.setCapacity(faceOwner_.size() + nAllFaces);
        faceNeighbour_.setCapacity(faceNeighbour_.size() + nAllFaces);
        faceMap_.setCapacity(faceMap_.size() + nAllFaces);
        reverseFaceMap_.setCapacity(reverseFaceMap_.size() + nAllFaces);
        faceFromPoint_.resize(faceFromPoint_.size() + nAllFaces/100);
        faceFromEdge_.resize(faceFromEdge_.size() + nAllFaces/100);
        flipFaceFlux_.setCapacity(faces_.size() + nAllFaces);
        faceZone_.resize(faceZone_.size() + nAllFaces/100);
        faceZoneFlip_.setCapacity(faces_.size() + nAllFaces);


        // Precalc offset zones
        labelList newZoneID(nAllFaces, -1);
        boolList zoneFlip(nAllFaces, false);

        forAll(faceZones, zoneI)
        {
            const labelList& faceLabels = faceZones[zoneI];
            const boolList& flipMap = faceZones[zoneI].flipMap();

            forAll(faceLabels, j)
            {
                newZoneID[faceLabels[j]] = faceZoneMap[zoneI];
                zoneFlip[faceLabels[j]] = flipMap[j];
            }
        }

        // Add faces in mesh order

        // 1. Internal faces
        for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
        {
            addFace
            (
                faces[faceI],
                faceOwner[faceI],
                faceNeighbour[faceI],
                -1,                         // masterPointID
                -1,                         // masterEdgeID
                faceI,                      // masterFaceID
                false,                      // flipFaceFlux
                -1,                         // patchID
                newZoneID[faceI],           // zoneID
                zoneFlip[faceI]             // zoneFlip
            );
        }

        // 2. Patch faces
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (pp.start() != faces_.size())
            {
                FatalErrorIn
                (
                    "polyTopoChange::polyTopoChange"
                    "(const polyMesh& mesh, const bool strict)"
                )   << "Problem : "
                    << "Patch " << pp.name() << " starts at " << pp.start()
                    << endl
                    << "Current face counter at " << faces_.size() << endl
                    << "Are patches in incremental order?"
                    << abort(FatalError);
            }
            forAll(pp, patchFaceI)
            {
                label faceI = pp.start() + patchFaceI;

                addFace
                (
                    faces[faceI],
                    faceOwner[faceI],
                    -1,                         // neighbour
                    -1,                         // masterPointID
                    -1,                         // masterEdgeID
                    faceI,                      // masterFaceID
                    false,                      // flipFaceFlux
                    patchMap[patchI],           // patchID
                    newZoneID[faceI],           // zoneID
                    zoneFlip[faceI]             // zoneFlip
                );
            }
        }
    }
}


void Foam::polyTopoChange::setCapacity
(
    const label nPoints,
    const label nFaces,
    const label nCells
)
{
    points_.setCapacity(nPoints);
    pointMap_.setCapacity(nPoints);
    reversePointMap_.setCapacity(nPoints);
    pointZone_.resize(pointZone_.size() + nPoints/100);

    faces_.setCapacity(nFaces);
    region_.setCapacity(nFaces);
    faceOwner_.setCapacity(nFaces);
    faceNeighbour_.setCapacity(nFaces);
    faceMap_.setCapacity(nFaces);
    reverseFaceMap_.setCapacity(nFaces);
    faceFromPoint_.resize(faceFromPoint_.size() + nFaces/100);
    faceFromEdge_.resize(faceFromEdge_.size() + nFaces/100);
    flipFaceFlux_.setCapacity(nFaces);
    faceZone_.resize(faceZone_.size() + nFaces/100);
    faceZoneFlip_.setCapacity(nFaces);

    cellMap_.setCapacity(nCells);
    reverseCellMap_.setCapacity(nCells);
    cellFromPoint_.resize(cellFromPoint_.size() + nCells/100);
    cellFromEdge_.resize(cellFromEdge_.size() + nCells/100);
    cellFromFace_.resize(cellFromFace_.size() + nCells/100);
    cellZone_.setCapacity(nCells);
}


Foam::label Foam::polyTopoChange::setAction(const topoAction& action)
{
    if (isType<polyAddPoint>(action))
    {
        const polyAddPoint& pap = refCast<const polyAddPoint>(action);

        return addPoint
        (
            pap.newPoint(),
            pap.masterPointID(),
            pap.zoneID(),
            pap.inCell()
        );
    }
    else if (isType<polyModifyPoint>(action))
    {
        const polyModifyPoint& pmp = refCast<const polyModifyPoint>(action);

        modifyPoint
        (
            pmp.pointID(),
            pmp.newPoint(),
            pmp.zoneID(),
            pmp.inCell()
        );

        return -1;
    }
    else if (isType<polyRemovePoint>(action))
    {
        const polyRemovePoint& prp = refCast<const polyRemovePoint>(action);

        removePoint(prp.pointID(), prp.mergePointID());

        return -1;
    }
    else if (isType<polyAddFace>(action))
    {
        const polyAddFace& paf = refCast<const polyAddFace>(action);

        return addFace
        (
            paf.newFace(),
            paf.owner(),
            paf.neighbour(),
            paf.masterPointID(),
            paf.masterEdgeID(),
            paf.masterFaceID(),
            paf.flipFaceFlux(),
            paf.patchID(),
            paf.zoneID(),
            paf.zoneFlip()
        );
    }
    else if (isType<polyModifyFace>(action))
    {
        const polyModifyFace& pmf = refCast<const polyModifyFace>(action);

        modifyFace
        (
            pmf.newFace(),
            pmf.faceID(),
            pmf.owner(),
            pmf.neighbour(),
            pmf.flipFaceFlux(),
            pmf.patchID(),
            pmf.zoneID(),
            pmf.zoneFlip()
        );

        return -1;
    }
    else if (isType<polyRemoveFace>(action))
    {
        const polyRemoveFace& prf = refCast<const polyRemoveFace>(action);

        removeFace(prf.faceID(), prf.mergeFaceID());

        return -1;
    }
    else if (isType<polyAddCell>(action))
    {
        const polyAddCell& pac = refCast<const polyAddCell>(action);

        return addCell
        (
            pac.masterPointID(),
            pac.masterEdgeID(),
            pac.masterFaceID(),
            pac.masterCellID(),
            pac.zoneID()
        );
    }
    else if (isType<polyModifyCell>(action))
    {
        const polyModifyCell& pmc = refCast<const polyModifyCell>(action);

        if (pmc.removeFromZone())
        {
            modifyCell(pmc.cellID(), -1);
        }
        else
        {
            modifyCell(pmc.cellID(), pmc.zoneID());
        }

        return -1;
    }
    else if (isType<polyRemoveCell>(action))
    {
        const polyRemoveCell& prc = refCast<const polyRemoveCell>(action);

        removeCell(prc.cellID(), prc.mergeCellID());

        return -1;
    }
    else
    {
        FatalErrorIn
        (
            "label polyTopoChange::setAction(const topoAction& action)"
        )   << "Unknown type of topoChange: " << action.type()
            << abort(FatalError);

        // Dummy return to keep compiler happy
        return -1;
    }
}


Foam::label Foam::polyTopoChange::addPoint
(
    const point& pt,
    const label masterPointID,
    const label zoneID,
    const bool inCell
)
{
    label pointI = points_.size();

    points_.append(pt);
    pointMap_.append(masterPointID);
    reversePointMap_.append(pointI);

    if (zoneID >= 0)
    {
        pointZone_.insert(pointI, zoneID);
    }

    if (!inCell)
    {
        retiredPoints_.insert(pointI);
    }

    return pointI;
}


void Foam::polyTopoChange::modifyPoint
(
    const label pointI,
    const point& pt,
    const label newZoneID,
    const bool inCell
)
{
    if (pointI < 0 || pointI >= points_.size())
    {
        FatalErrorIn
        (
            "polyTopoChange::modifyPoint(const label, const point&)"
        )   << "illegal point label " << pointI << endl
            << "Valid point labels are 0 .. " << points_.size()-1
            << abort(FatalError);
    }
    if (pointRemoved(pointI) || pointMap_[pointI] == -1)
    {
        FatalErrorIn
        (
            "polyTopoChange::modifyPoint(const label, const point&)"
        )   << "point " << pointI << " already marked for removal"
            << abort(FatalError);
    }
    points_[pointI] = pt;

    Map<label>::iterator pointFnd = pointZone_.find(pointI);

    if (pointFnd != pointZone_.end())
    {
        if (newZoneID >= 0)
        {
            pointFnd() = newZoneID;
        }
        else
        {
            pointZone_.erase(pointFnd);
        }
    }
    else if (newZoneID >= 0)
    {
        pointZone_.insert(pointI, newZoneID);
    }

    if (inCell)
    {
        retiredPoints_.erase(pointI);
    }
    else
    {
        retiredPoints_.insert(pointI);
    }
}


void Foam::polyTopoChange::movePoints(const pointField& newPoints)
{
    if (newPoints.size() != points_.size())
    {
        FatalErrorIn("polyTopoChange::movePoints(const pointField&)")
            << "illegal pointField size." << endl
            << "Size:" << newPoints.size() << endl
            << "Points in mesh:" << points_.size()
            << abort(FatalError);
    }

    forAll(points_, pointI)
    {
        points_[pointI] = newPoints[pointI];
    }
}


void Foam::polyTopoChange::removePoint
(
    const label pointI,
    const label mergePointI
)
{
    if (pointI < 0 || pointI >= points_.size())
    {
        FatalErrorIn("polyTopoChange::removePoint(const label, const label)")
            << "illegal point label " << pointI << endl
            << "Valid point labels are 0 .. " << points_.size()-1
            << abort(FatalError);
    }

    if
    (
        strict_
     && (pointRemoved(pointI) || pointMap_[pointI] == -1)
    )
    {
        FatalErrorIn("polyTopoChange::removePoint(const label, const label)")
            << "point " << pointI << " already marked for removal" << nl
            << "Point:" << points_[pointI] << " pointMap:" << pointMap_[pointI]
            << abort(FatalError);
    }

    if (pointI == mergePointI)
    {
        FatalErrorIn("polyTopoChange::removePoint(const label, const label)")
            << "Cannot remove/merge point " << pointI << " onto itself."
            << abort(FatalError);
    }

    points_[pointI] = point::max;
    pointMap_[pointI] = -1;
    if (mergePointI >= 0)
    {
        reversePointMap_[pointI] = -mergePointI-2;
    }
    else
    {
        reversePointMap_[pointI] = -1;
    }
    pointZone_.erase(pointI);
    retiredPoints_.erase(pointI);
}


Foam::label Foam::polyTopoChange::addFace
(
    const face& f,
    const label own,
    const label nei,
    const label masterPointID,
    const label masterEdgeID,
    const label masterFaceID,
    const bool flipFaceFlux,
    const label patchID,
    const label zoneID,
    const bool zoneFlip
)
{
    // Check validity
    if (debug)
    {
        checkFace(f, -1, own, nei, patchID, zoneID);
    }

    label faceI = faces_.size();

    faces_.append(f);
    region_.append(patchID);
    faceOwner_.append(own);
    faceNeighbour_.append(nei);

    if (masterPointID >= 0)
    {
        faceMap_.append(-1);
        faceFromPoint_.insert(faceI, masterPointID);
    }
    else if (masterEdgeID >= 0)
    {
        faceMap_.append(-1);
        faceFromEdge_.insert(faceI, masterEdgeID);
    }
    else if (masterFaceID >= 0)
    {
        faceMap_.append(masterFaceID);
    }
    else
    {
        // Allow inflate-from-nothing?
        //FatalErrorIn("polyTopoChange::addFace")
        //    << "Need to specify a master point, edge or face"
        //    << "face:" << f << " own:" << own << " nei:" << nei
        //    << abort(FatalError);
        faceMap_.append(-1);
    }
    reverseFaceMap_.append(faceI);

    flipFaceFlux_[faceI] = (flipFaceFlux ? 1 : 0);

    if (zoneID >= 0)
    {
        faceZone_.insert(faceI, zoneID);
    }
    faceZoneFlip_[faceI] = (zoneFlip ? 1 : 0);

    return faceI;
}


void Foam::polyTopoChange::modifyFace
(
    const face& f,
    const label faceI,
    const label own,
    const label nei,
    const bool flipFaceFlux,
    const label patchID,
    const label zoneID,
    const bool zoneFlip
)
{
    // Check validity
    if (debug)
    {
        checkFace(f, faceI, own, nei, patchID, zoneID);
    }

    faces_[faceI] = f;
    faceOwner_[faceI] = own;
    faceNeighbour_[faceI] = nei;
    region_[faceI] = patchID;

    flipFaceFlux_[faceI] = (flipFaceFlux ? 1 : 0);

    Map<label>::iterator faceFnd = faceZone_.find(faceI);

    if (faceFnd != faceZone_.end())
    {
        if (zoneID >= 0)
        {
            faceFnd() = zoneID;
        }
        else
        {
            faceZone_.erase(faceFnd);
        }
    }
    else if (zoneID >= 0)
    {
        faceZone_.insert(faceI, zoneID);
    }
    faceZoneFlip_[faceI] = (zoneFlip ? 1 : 0);
}


void Foam::polyTopoChange::removeFace(const label faceI, const label mergeFaceI)
{
    if (faceI < 0 || faceI >= faces_.size())
    {
        FatalErrorIn("polyTopoChange::removeFace(const label, const label)")
            << "illegal face label " << faceI << endl
            << "Valid face labels are 0 .. " << faces_.size()-1
            << abort(FatalError);
    }

    if
    (
        strict_
     && (faceRemoved(faceI) || faceMap_[faceI] == -1)
    )
    {
        FatalErrorIn("polyTopoChange::removeFace(const label, const label)")
            << "face " << faceI
            << " already marked for removal"
            << abort(FatalError);
    }

    faces_[faceI].setSize(0);
    region_[faceI] = -1;
    faceOwner_[faceI] = -1;
    faceNeighbour_[faceI] = -1;
    faceMap_[faceI] = -1;
    if (mergeFaceI >= 0)
    {
        reverseFaceMap_[faceI] = -mergeFaceI-2;
    }
    else
    {
        reverseFaceMap_[faceI] = -1;
    }
    faceFromEdge_.erase(faceI);
    faceFromPoint_.erase(faceI);
    flipFaceFlux_[faceI] = 0;
    faceZone_.erase(faceI);
    faceZoneFlip_[faceI] = 0;
}


Foam::label Foam::polyTopoChange::addCell
(
    const label masterPointID,
    const label masterEdgeID,
    const label masterFaceID,
    const label masterCellID,
    const label zoneID
)
{
    label cellI = cellMap_.size();

    if (masterPointID >= 0)
    {
        cellMap_.append(-1);
        cellFromPoint_.insert(cellI, masterPointID);
    }
    else if (masterEdgeID >= 0)
    {
        cellMap_.append(-1);
        cellFromEdge_.insert(cellI, masterEdgeID);
    }
    else if (masterFaceID >= 0)
    {
        cellMap_.append(-1);
        cellFromFace_.insert(cellI, masterFaceID);
    }
    else
    {
        cellMap_.append(masterCellID);
    }
    reverseCellMap_.append(cellI);
    cellZone_.append(zoneID);

    return cellI;
}


void Foam::polyTopoChange::modifyCell
(
    const label cellI,
    const label zoneID
)
{
    cellZone_[cellI] = zoneID;
}


void Foam::polyTopoChange::removeCell(const label cellI, const label mergeCellI)
{
    if (cellI < 0 || cellI >= cellMap_.size())
    {
        FatalErrorIn("polyTopoChange::removeCell(const label, const label)")
            << "illegal cell label " << cellI << endl
            << "Valid cell labels are 0 .. " << cellMap_.size()-1
            << abort(FatalError);
    }

    if (strict_ && cellMap_[cellI] == -2)
    {
        FatalErrorIn("polyTopoChange::removeCell(const label, const label)")
            << "cell " << cellI
            << " already marked for removal"
            << abort(FatalError);
    }

    cellMap_[cellI] = -2;
    if (mergeCellI >= 0)
    {
        reverseCellMap_[cellI] = -mergeCellI-2;
    }
    else
    {
        reverseCellMap_[cellI] = -1;
    }
    cellFromPoint_.erase(cellI);
    cellFromEdge_.erase(cellI);
    cellFromFace_.erase(cellI);
    cellZone_[cellI] = -1;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::polyTopoChange::changeMesh
(
    polyMesh& mesh,
    const bool inflate,
    const bool syncParallel,
    const bool orderCells,
    const bool orderPoints
)
{
    if (debug)
    {
        Pout<< "polyTopoChange::changeMesh"
            << "(polyMesh&, const bool, const bool, const bool, const bool)"
            << endl;
    }

    if (debug)
    {
        Pout<< "Old mesh:" << nl;
        writeMeshStats(mesh, Pout);
    }

    // new mesh points
    pointField newPoints;
    // number of internal points
    label nInternalPoints;
    // patch slicing
    labelList patchSizes;
    labelList patchStarts;
    // inflate maps
    List<objectMap> pointsFromPoints;
    List<objectMap> facesFromPoints;
    List<objectMap> facesFromEdges;
    List<objectMap> facesFromFaces;
    List<objectMap> cellsFromPoints;
    List<objectMap> cellsFromEdges;
    List<objectMap> cellsFromFaces;
    List<objectMap> cellsFromCells;
    // old mesh info
    List<Map<label> > oldPatchMeshPointMaps;
    labelList oldPatchNMeshPoints;
    labelList oldPatchStarts;
    List<Map<label> > oldFaceZoneMeshPointMaps;

    // Compact, reorder patch faces and calculate mesh/patch maps.
    compactAndReorder
    (
        mesh,
        syncParallel,
        orderCells,
        orderPoints,

        nInternalPoints,
        newPoints,
        patchSizes,
        patchStarts,
        pointsFromPoints,
        facesFromPoints,
        facesFromEdges,
        facesFromFaces,
        cellsFromPoints,
        cellsFromEdges,
        cellsFromFaces,
        cellsFromCells,
        oldPatchMeshPointMaps,
        oldPatchNMeshPoints,
        oldPatchStarts,
        oldFaceZoneMeshPointMaps
    );

    const label nOldPoints(mesh.nPoints());
    const label nOldFaces(mesh.nFaces());
    const label nOldCells(mesh.nCells());
    autoPtr<scalarField> oldCellVolumes(new scalarField(mesh.cellVolumes()));


    // Change the mesh
    // ~~~~~~~~~~~~~~~
    // This will invalidate any addressing so better make sure you have
    // all the information you need!!!

    if (inflate)
    {
        // Keep (renumbered) mesh points, store new points in map for inflation
        // (appended points (i.e. from nowhere) get value zero)
        pointField renumberedMeshPoints(newPoints.size());

        forAll(pointMap_, newPointI)
        {
            label oldPointI = pointMap_[newPointI];

            if (oldPointI >= 0)
            {
                renumberedMeshPoints[newPointI] = mesh.points()[oldPointI];
            }
            else
            {
                renumberedMeshPoints[newPointI] = vector::zero;
            }
        }

        mesh.resetPrimitives
        (
            xferMove(renumberedMeshPoints),
            faces_.xfer(),
            faceOwner_.xfer(),
            faceNeighbour_.xfer(),
            patchSizes,
            patchStarts,
            syncParallel
        );

        mesh.topoChanging(true);
        // Note: could already set moving flag as well
        //       mesh.moving(true);
    }
    else
    {
        // Set new points.
        mesh.resetPrimitives
        (
            xferMove(newPoints),
            faces_.xfer(),
            faceOwner_.xfer(),
            faceNeighbour_.xfer(),
            patchSizes,
            patchStarts,
            syncParallel
        );
        mesh.topoChanging(true);
    }

    // Clear out primitives
    {
        retiredPoints_.clearStorage();
        region_.clearStorage();
    }


    if (debug)
    {
        // Some stats on changes
        label nAdd, nInflate, nMerge, nRemove;
        countMap(pointMap_, reversePointMap_, nAdd, nInflate, nMerge, nRemove);
        Pout<< "Points:"
            << "  added(from point):" << nAdd
            << "  added(from nothing):" << nInflate
            << "  merged(into other point):" << nMerge
            << "  removed:" << nRemove
            << nl;

        countMap(faceMap_, reverseFaceMap_, nAdd, nInflate, nMerge, nRemove);
        Pout<< "Faces:"
            << "  added(from face):" << nAdd
            << "  added(inflated):" << nInflate
            << "  merged(into other face):" << nMerge
            << "  removed:" << nRemove
            << nl;

        countMap(cellMap_, reverseCellMap_, nAdd, nInflate, nMerge, nRemove);
        Pout<< "Cells:"
            << "  added(from cell):" << nAdd
            << "  added(inflated):" << nInflate
            << "  merged(into other cell):" << nMerge
            << "  removed:" << nRemove
            << nl
            << endl;
    }

    if (debug)
    {
        Pout<< "New mesh:" << nl;
        writeMeshStats(mesh, Pout);
    }


    // Zones
    // ~~~~~

    // Inverse of point/face/cell zone addressing.
    // For every preserved point/face/cells in zone give the old position.
    // For added points, the index is set to -1
    labelListList pointZoneMap(mesh.pointZones().size());
    labelListList faceZoneFaceMap(mesh.faceZones().size());
    labelListList cellZoneMap(mesh.cellZones().size());

    resetZones(mesh, mesh, pointZoneMap, faceZoneFaceMap, cellZoneMap);

    // Clear zone info
    {
        pointZone_.clearStorage();
        faceZone_.clearStorage();
        faceZoneFlip_.clearStorage();
        cellZone_.clearStorage();
    }


    // Patch point renumbering
    // For every preserved point on a patch give the old position.
    // For added points, the index is set to -1
    labelListList patchPointMap(mesh.boundaryMesh().size());
    calcPatchPointMap
    (
        oldPatchMeshPointMaps,
        mesh.boundaryMesh(),
        patchPointMap
    );

    // Create the face zone mesh point renumbering
    labelListList faceZonePointMap(mesh.faceZones().size());
    calcFaceZonePointMap(mesh, oldFaceZoneMeshPointMaps, faceZonePointMap);

    labelHashSet flipFaceFluxSet(getSetIndices(flipFaceFlux_));

    return autoPtr<mapPolyMesh>
    (
        new mapPolyMesh
        (
            mesh,
            nOldPoints,
            nOldFaces,
            nOldCells,

            pointMap_,
            pointsFromPoints,

            faceMap_,
            facesFromPoints,
            facesFromEdges,
            facesFromFaces,

            cellMap_,
            cellsFromPoints,
            cellsFromEdges,
            cellsFromFaces,
            cellsFromCells,

            reversePointMap_,
            reverseFaceMap_,
            reverseCellMap_,

            flipFaceFluxSet,

            patchPointMap,

            pointZoneMap,

            faceZonePointMap,
            faceZoneFaceMap,
            cellZoneMap,

            newPoints,          // if empty signals no inflation.
            oldPatchStarts,
            oldPatchNMeshPoints,

            oldCellVolumes,

            true                // steal storage.
        )
    );

    // At this point all member DynamicList (pointMap_, cellMap_ etc.) will
    // be invalid.
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::polyTopoChange::makeMesh
(
    autoPtr<fvMesh>& newMeshPtr,
    const IOobject& io,
    const polyMesh& mesh,
    const bool syncParallel,
    const bool orderCells,
    const bool orderPoints
)
{
    if (debug)
    {
        Pout<< "polyTopoChange::changeMesh"
            << "(autoPtr<fvMesh>&, const IOobject&, const fvMesh&"
            << ", const bool, const bool, const bool)"
            << endl;
    }

    if (debug)
    {
        Pout<< "Old mesh:" << nl;
        writeMeshStats(mesh, Pout);
    }

    // new mesh points
    pointField newPoints;
    // number of internal points
    label nInternalPoints;
    // patch slicing
    labelList patchSizes;
    labelList patchStarts;
    // inflate maps
    List<objectMap> pointsFromPoints;
    List<objectMap> facesFromPoints;
    List<objectMap> facesFromEdges;
    List<objectMap> facesFromFaces;
    List<objectMap> cellsFromPoints;
    List<objectMap> cellsFromEdges;
    List<objectMap> cellsFromFaces;
    List<objectMap> cellsFromCells;

    // old mesh info
    List<Map<label> > oldPatchMeshPointMaps;
    labelList oldPatchNMeshPoints;
    labelList oldPatchStarts;
    List<Map<label> > oldFaceZoneMeshPointMaps;

    // Compact, reorder patch faces and calculate mesh/patch maps.
    compactAndReorder
    (
        mesh,
        syncParallel,
        orderCells,
        orderPoints,

        nInternalPoints,
        newPoints,
        patchSizes,
        patchStarts,
        pointsFromPoints,
        facesFromPoints,
        facesFromEdges,
        facesFromFaces,
        cellsFromPoints,
        cellsFromEdges,
        cellsFromFaces,
        cellsFromCells,
        oldPatchMeshPointMaps,
        oldPatchNMeshPoints,
        oldPatchStarts,
        oldFaceZoneMeshPointMaps
    );

    const label nOldPoints(mesh.nPoints());
    const label nOldFaces(mesh.nFaces());
    const label nOldCells(mesh.nCells());
    autoPtr<scalarField> oldCellVolumes(new scalarField(mesh.cellVolumes()));


    // Create the mesh
    // ~~~~~~~~~~~~~~~

    IOobject noReadIO(io);
    noReadIO.readOpt() = IOobject::NO_READ;
    newMeshPtr.reset
    (
        new fvMesh
        (
            noReadIO,
            xferMove(newPoints),
            faces_.xfer(),
            faceOwner_.xfer(),
            faceNeighbour_.xfer()
        )
    );
    fvMesh& newMesh = newMeshPtr();

    // Clear out primitives
    {
        retiredPoints_.clearStorage();
        region_.clearStorage();
    }


    if (debug)
    {
        // Some stats on changes
        label nAdd, nInflate, nMerge, nRemove;
        countMap(pointMap_, reversePointMap_, nAdd, nInflate, nMerge, nRemove);
        Pout<< "Points:"
            << "  added(from point):" << nAdd
            << "  added(from nothing):" << nInflate
            << "  merged(into other point):" << nMerge
            << "  removed:" << nRemove
            << nl;

        countMap(faceMap_, reverseFaceMap_, nAdd, nInflate, nMerge, nRemove);
        Pout<< "Faces:"
            << "  added(from face):" << nAdd
            << "  added(inflated):" << nInflate
            << "  merged(into other face):" << nMerge
            << "  removed:" << nRemove
            << nl;

        countMap(cellMap_, reverseCellMap_, nAdd, nInflate, nMerge, nRemove);
        Pout<< "Cells:"
            << "  added(from cell):" << nAdd
            << "  added(inflated):" << nInflate
            << "  merged(into other cell):" << nMerge
            << "  removed:" << nRemove
            << nl
            << endl;
    }


    {
        const polyBoundaryMesh& oldPatches = mesh.boundaryMesh();

        List<polyPatch*> newBoundary(oldPatches.size());

        forAll(oldPatches, patchI)
        {
            newBoundary[patchI] = oldPatches[patchI].clone
            (
                newMesh.boundaryMesh(),
                patchI,
                patchSizes[patchI],
                patchStarts[patchI]
            ).ptr();
        }
        newMesh.addFvPatches(newBoundary);
    }


    // Zones
    // ~~~~~

    // Start off from empty zones.
    const pointZoneMesh& oldPointZones = mesh.pointZones();
    List<pointZone*> pZonePtrs(oldPointZones.size());
    {
        forAll(oldPointZones, i)
        {
            pZonePtrs[i] = new pointZone
            (
                oldPointZones[i].name(),
                labelList(0),
                i,
                newMesh.pointZones()
            );
        }
    }

    const faceZoneMesh& oldFaceZones = mesh.faceZones();
    List<faceZone*> fZonePtrs(oldFaceZones.size());
    {
        forAll(oldFaceZones, i)
        {
            fZonePtrs[i] = new faceZone
            (
                oldFaceZones[i].name(),
                labelList(0),
                boolList(0),
                i,
                newMesh.faceZones()
            );
        }
    }

    const cellZoneMesh& oldCellZones = mesh.cellZones();
    List<cellZone*> cZonePtrs(oldCellZones.size());
    {
        forAll(oldCellZones, i)
        {
            cZonePtrs[i] = new cellZone
            (
                oldCellZones[i].name(),
                labelList(0),
                i,
                newMesh.cellZones()
            );
        }
    }

    newMesh.addZones(pZonePtrs, fZonePtrs, cZonePtrs);

    // Inverse of point/face/cell zone addressing.
    // For every preserved point/face/cells in zone give the old position.
    // For added points, the index is set to -1
    labelListList pointZoneMap(mesh.pointZones().size());
    labelListList faceZoneFaceMap(mesh.faceZones().size());
    labelListList cellZoneMap(mesh.cellZones().size());

    resetZones(mesh, newMesh, pointZoneMap, faceZoneFaceMap, cellZoneMap);

    // Clear zone info
    {
        pointZone_.clearStorage();
        faceZone_.clearStorage();
        faceZoneFlip_.clearStorage();
        cellZone_.clearStorage();
    }

    // Patch point renumbering
    // For every preserved point on a patch give the old position.
    // For added points, the index is set to -1
    labelListList patchPointMap(newMesh.boundaryMesh().size());
    calcPatchPointMap
    (
        oldPatchMeshPointMaps,
        newMesh.boundaryMesh(),
        patchPointMap
    );

    // Create the face zone mesh point renumbering
    labelListList faceZonePointMap(newMesh.faceZones().size());
    calcFaceZonePointMap(newMesh, oldFaceZoneMeshPointMaps, faceZonePointMap);

    if (debug)
    {
        Pout<< "New mesh:" << nl;
        writeMeshStats(mesh, Pout);
    }

    labelHashSet flipFaceFluxSet(getSetIndices(flipFaceFlux_));

    return autoPtr<mapPolyMesh>
    (
        new mapPolyMesh
        (
            newMesh,
            nOldPoints,
            nOldFaces,
            nOldCells,

            pointMap_,
            pointsFromPoints,

            faceMap_,
            facesFromPoints,
            facesFromEdges,
            facesFromFaces,

            cellMap_,
            cellsFromPoints,
            cellsFromEdges,
            cellsFromFaces,
            cellsFromCells,

            reversePointMap_,
            reverseFaceMap_,
            reverseCellMap_,

            flipFaceFluxSet,

            patchPointMap,

            pointZoneMap,

            faceZonePointMap,
            faceZoneFaceMap,
            cellZoneMap,

            newPoints,          // if empty signals no inflation.
            oldPatchStarts,
            oldPatchNMeshPoints,
            oldCellVolumes,
            true                // steal storage.
        )
    );

    // At this point all member DynamicList (pointMap_, cellMap_ etc.) will
    // be invalid.
}


// ************************************************************************* //
