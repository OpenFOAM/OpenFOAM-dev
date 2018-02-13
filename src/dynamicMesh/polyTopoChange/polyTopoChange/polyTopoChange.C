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

    forAll(map, newCelli)
    {
        label oldCelli = map[newCelli];

        if (oldCelli >= 0)
        {
            if (reverseMap[oldCelli] == newCelli)
            {
                // unchanged
            }
            else
            {
                // Added (from another cell v.s. inflated from face/point)
                nAdd++;
            }
        }
        else if (oldCelli == -1)
        {
            // Created from nothing
            nInflate++;
        }
        else
        {
            FatalErrorInFunction
                << " new:" << newCelli << abort(FatalError);
        }
    }

    forAll(reverseMap, oldCelli)
    {
        label newCelli = reverseMap[oldCelli];

        if (newCelli >= 0)
        {
            // unchanged
        }
        else if (newCelli == -1)
        {
            // removed
            nRemove++;
        }
        else
        {
            // merged into -newCelli-2
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
    forAll(patches, patchi)
    {
        patchSizes[patchi] = patches[patchi].size();
        patchStarts[patchi] = patches[patchi].start();
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

    forAll(reverseCellMap, oldCelli)
    {
        label newCelli = reverseCellMap[oldCelli];

        if (newCelli < -1)
        {
            label mergeCelli = -newCelli-2;

            nMerged[mergeCelli]++;
        }
    }

    // From merged cell to set index
    labelList cellToMergeSet(cellMap.size(), -1);

    label nSets = 0;

    forAll(nMerged, celli)
    {
        if (nMerged[celli] > 1)
        {
            cellToMergeSet[celli] = nSets++;
        }
    }

    // Collect cell labels.
    // Each objectMap will have
    // - index : new mesh cell label
    // - masterObjects : list of old cells that have been merged. Element 0
    //                   will be the original destination cell label.

    cellsFromCells.setSize(nSets);

    forAll(reverseCellMap, oldCelli)
    {
        label newCelli = reverseCellMap[oldCelli];

        if (newCelli < -1)
        {
            label mergeCelli = -newCelli-2;

            // oldCelli was merged into mergeCelli

            label setI = cellToMergeSet[mergeCelli];

            objectMap& mergeSet = cellsFromCells[setI];

            if (mergeSet.masterObjects().empty())
            {
                // First occurrence of master cell mergeCelli

                mergeSet.index() = mergeCelli;
                mergeSet.masterObjects().setSize(nMerged[mergeCelli]);

                // old master label
                mergeSet.masterObjects()[0] = cellMap[mergeCelli];

                // old slave label
                mergeSet.masterObjects()[1] = oldCelli;

                nMerged[mergeCelli] = 2;
            }
            else
            {
                mergeSet.masterObjects()[nMerged[mergeCelli]++] = oldCelli;
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
            FatalErrorInFunction
                << "Problem." << abort(FatalError);
        }
        points[fp] = points_[f[fp]];
    }
    return points;
}


void Foam::polyTopoChange::checkFace
(
    const face& f,
    const label facei,
    const label own,
    const label nei,
    const label patchi,
    const label zoneI
) const
{
    if (nei == -1)
    {
        if (own == -1 && zoneI != -1)
        {
            // retired face
        }
        else if (patchi == -1 || patchi >= nPatches_)
        {
            FatalErrorInFunction
                << "Face has no neighbour (so external) but does not have"
                << " a valid patch" << nl
                << "f:" << f
                << " facei(-1 if added face):" << facei
                << " own:" << own << " nei:" << nei
                << " patchi:" << patchi << nl;
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
        if (patchi != -1)
        {
            FatalErrorInFunction
                << "Cannot both have valid patchi and neighbour" << nl
                << "f:" << f
                << " facei(-1 if added face):" << facei
                << " own:" << own << " nei:" << nei
                << " patchi:" << patchi << nl;
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
            FatalErrorInFunction
                << "Owner cell label should be less than neighbour cell label"
                << nl
                << "f:" << f
                << " facei(-1 if added face):" << facei
                << " own:" << own << " nei:" << nei
                << " patchi:" << patchi << nl;
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
        FatalErrorInFunction
            << "Illegal vertices in face"
            << nl
            << "f:" << f
            << " facei(-1 if added face):" << facei
            << " own:" << own << " nei:" << nei
            << " patchi:" << patchi << nl;
            if (hasValidPoints(f))
            {
                FatalError
                        << "points (removed points marked with "
                        << vector::max << ") : " << facePoints(f);
            }
            FatalError << abort(FatalError);
    }
    if (facei >= 0 && facei < faces_.size() && faceRemoved(facei))
    {
        FatalErrorInFunction
            << "Face already marked for removal"
            << nl
            << "f:" << f
            << " facei(-1 if added face):" << facei
            << " own:" << own << " nei:" << nei
            << " patchi:" << patchi << nl;
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
            FatalErrorInFunction
                << "Face uses removed vertices"
                << nl
                << "f:" << f
                << " facei(-1 if added face):" << facei
                << " own:" << own << " nei:" << nei
                << " patchi:" << patchi << nl;
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

    for (label facei = 0; facei < nActiveFaces; facei++)
    {
        if (faceOwner_[facei] < 0)
        {
            FatalErrorInFunction
                << "Face " << facei << " is active but its owner has"
                << " been deleted. This is usually due to deleting cells"
                << " without modifying exposed faces to be boundary faces."
                << exit(FatalError);
        }
        nNbrs[faceOwner_[facei]]++;
    }
    for (label facei = 0; facei < nActiveFaces; facei++)
    {
        if (faceNeighbour_[facei] >= 0)
        {
            nNbrs[faceNeighbour_[facei]]++;
        }
    }

    // 2. Calculate offsets

    cellFaceOffsets[0] = 0;
    forAll(nNbrs, celli)
    {
        cellFaceOffsets[celli+1] = cellFaceOffsets[celli] + nNbrs[celli];
    }

    // 3. Fill faces per cell

    // reset the whole list to use as counter
    nNbrs = 0;

    for (label facei = 0; facei < nActiveFaces; facei++)
    {
        label celli = faceOwner_[facei];

        cellFaces[cellFaceOffsets[celli] + nNbrs[celli]++] = facei;
    }

    for (label facei = 0; facei < nActiveFaces; facei++)
    {
        label celli = faceNeighbour_[facei];

        if (celli >= 0)
        {
            cellFaces[cellFaceOffsets[celli] + nNbrs[celli]++] = facei;
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

    for (label facei = 0; facei < nActiveFaces; facei++)
    {
        if (faceNeighbour_[facei] >= 0)
        {
            nNbrs[faceOwner_[facei]]++;
            nNbrs[faceNeighbour_[facei]]++;
        }
    }

    // 2. Construct csr
    cellCells.setSize(nNbrs);


    // 3. Fill faces per cell

    // reset the whole list to use as counter
    nNbrs = 0;

    for (label facei = 0; facei < nActiveFaces; facei++)
    {
        label nei = faceNeighbour_[facei];

        if (nei >= 0)
        {
            label own = faceOwner_[facei];
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

        forAll(visited, celli)
        {
            // find the lowest connected cell that has not been visited yet
            if (!cellRemoved(celli) && !visited[celli])
            {
                if (cellCellAddressing[celli].size() < minWeight)
                {
                    minWeight = cellCellAddressing[celli].size();
                    currentCell = celli;
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
    label newFacei = 0;

    labelList nbr;
    labelList order;

    forAll(cellMap_, celli)
    {
        label startOfCell = cellFaceOffsets[celli];
        label nFaces = cellFaceOffsets[celli+1] - startOfCell;

        // Neighbouring cells
        //SortableList<label> nbr(nFaces);
        nbr.setSize(nFaces);

        for (label i = 0; i < nFaces; i++)
        {
            label facei = cellFaces[startOfCell + i];

            label nbrCelli = faceNeighbour_[facei];

            if (facei >= nActiveFaces)
            {
                // Retired face.
                nbr[i] = -1;
            }
            else if (nbrCelli != -1)
            {
                // Internal face. Get cell on other side.
                if (nbrCelli == celli)
                {
                    nbrCelli = faceOwner_[facei];
                }

                if (celli < nbrCelli)
                {
                    // Celli is master
                    nbr[i] = nbrCelli;
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
        //            newFacei++;
        //    }
        //}
        forAll(order, i)
        {
            label index = order[i];
            if (nbr[index] != -1)
            {
                oldToNew[cellFaces[startOfCell + index]] = newFacei++;
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
        patchStarts[0] = newFacei;

        for (label facei = 0; facei < nActiveFaces; facei++)
        {
            if (region_[facei] >= 0)
            {
                patchSizes[region_[facei]]++;
            }
        }

        label facei = patchStarts[0];

        forAll(patchStarts, patchi)
        {
            patchStarts[patchi] = facei;
            facei += patchSizes[patchi];
        }
    }

    //if (debug)
    //{
    //    Pout<< "patchSizes:" << patchSizes << nl
    //        << "patchStarts:" << patchStarts << endl;
    //}

    labelList workPatchStarts(patchStarts);

    for (label facei = 0; facei < nActiveFaces; facei++)
    {
        if (region_[facei] >= 0)
        {
            oldToNew[facei] = workPatchStarts[region_[facei]]++;
        }
    }

    // Retired faces.
    for (label facei = nActiveFaces; facei < oldToNew.size(); facei++)
    {
        oldToNew[facei] = facei;
    }

    // Check done all faces.
    forAll(oldToNew, facei)
    {
        if (oldToNew[facei] == -1)
        {
            FatalErrorInFunction
                << "Did not determine new position"
                << " for face " << facei
                << " owner " << faceOwner_[facei]
                << " neighbour " << faceNeighbour_[facei]
                << " region " << region_[facei] << endl
                << "This is usually caused by not specifying a patch for"
                << " a boundary face." << nl
                << "Switch on the polyTopoChange::debug flag to catch"
                << " this error earlier." << nl;
            if (hasValidPoints(faces_[facei]))
            {
                FatalError
                        << "points (removed points marked with "
                        << vector::max << ") " << facePoints(faces_[facei]);
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
        label newPointi = 0;

        if (!orderPoints)
        {
            nInternalPoints = -1;

            forAll(points_, pointi)
            {
                if (!pointRemoved(pointi) && !retiredPoints_.found(pointi))
                {
                    localPointMap[pointi] = newPointi++;
                }
            }
            nActivePoints = newPointi;
        }
        else
        {
            forAll(points_, pointi)
            {
                if (!pointRemoved(pointi) && !retiredPoints_.found(pointi))
                {
                    nActivePoints++;
                }
            }

            // Mark boundary points
            forAll(faceOwner_, facei)
            {
                if
                (
                   !faceRemoved(facei)
                 && faceOwner_[facei] >= 0
                 && faceNeighbour_[facei] < 0
                )
                {
                    // Valid boundary face
                    const face& f = faces_[facei];

                    forAll(f, fp)
                    {
                        label pointi = f[fp];

                        if (localPointMap[pointi] == -1)
                        {
                            if
                            (
                                pointRemoved(pointi)
                             || retiredPoints_.found(pointi)
                            )
                            {
                                FatalErrorInFunction
                                    << "Removed or retired point " << pointi
                                    << " in face " << f
                                    << " at position " << facei << endl
                                    << "Probably face has not been adapted for"
                                    << " removed points." << abort(FatalError);
                            }
                            localPointMap[pointi] = newPointi++;
                        }
                    }
                }
            }

            label nBoundaryPoints = newPointi;
            nInternalPoints = nActivePoints - nBoundaryPoints;

            // Move the boundary addressing up
            forAll(localPointMap, pointi)
            {
                if (localPointMap[pointi] != -1)
                {
                    localPointMap[pointi] += nInternalPoints;
                }
            }

            newPointi = 0;

            // Mark internal points
            forAll(faceOwner_, facei)
            {
                if
                (
                   !faceRemoved(facei)
                 && faceOwner_[facei] >= 0
                 && faceNeighbour_[facei] >= 0
                )
                {
                    // Valid internal face
                    const face& f = faces_[facei];

                    forAll(f, fp)
                    {
                        label pointi = f[fp];

                        if (localPointMap[pointi] == -1)
                        {
                            if
                            (
                                pointRemoved(pointi)
                             || retiredPoints_.found(pointi)
                            )
                            {
                                FatalErrorInFunction
                                    << "Removed or retired point " << pointi
                                    << " in face " << f
                                    << " at position " << facei << endl
                                    << "Probably face has not been adapted for"
                                    << " removed points." << abort(FatalError);
                            }
                            localPointMap[pointi] = newPointi++;
                        }
                    }
                }
            }

            if (newPointi != nInternalPoints)
            {
                FatalErrorInFunction
                    << "Problem." << abort(FatalError);
            }
            newPointi = nActivePoints;
        }

        forAllConstIter(labelHashSet, retiredPoints_, iter)
        {
            localPointMap[iter.key()] = newPointi++;
        }


        if (debug)
        {
            Pout<< "Points : active:" << nActivePoints
                << "  removed:" << points_.size()-newPointi << endl;
        }

        reorder(localPointMap, points_);
        points_.setCapacity(newPointi);

        // Update pointMaps
        reorder(localPointMap, pointMap_);
        pointMap_.setCapacity(newPointi);
        renumberReverseMap(localPointMap, reversePointMap_);

        renumberKey(localPointMap, pointZone_);
        renumber(localPointMap, retiredPoints_);

        // Use map to relabel face vertices
        forAll(faces_, facei)
        {
            face& f = faces_[facei];

            //labelList oldF(f);
            renumberCompact(localPointMap, f);

            if (!faceRemoved(facei) && f.size() < 3)
            {
                FatalErrorInFunction
                    << "Created illegal face " << f
                    //<< " from face " << oldF
                    << " at position:" << facei
                    << " when filtering removed points"
                    << abort(FatalError);
            }
        }
    }


    // Compact faces.
    {
        labelList localFaceMap(faces_.size(), -1);
        label newFacei = 0;

        forAll(faces_, facei)
        {
            if (!faceRemoved(facei) && faceOwner_[facei] >= 0)
            {
                localFaceMap[facei] = newFacei++;
            }
        }
        nActiveFaces_ = newFacei;

        forAll(faces_, facei)
        {
            if (!faceRemoved(facei) && faceOwner_[facei] < 0)
            {
                // Retired face
                localFaceMap[facei] = newFacei++;
            }
        }

        if (debug)
        {
            Pout<< "Faces : active:" << nActiveFaces_
                << "  removed:" << faces_.size()-newFacei << endl;
        }

        // Reorder faces.
        reorderCompactFaces(newFacei, localFaceMap);
    }

    // Compact cells.
    {
        labelList localCellMap;
        label newCelli;

        if (orderCells)
        {
            // Construct cellCell addressing
            CompactListList<label> cellCells;
            makeCellCells(nActiveFaces_, cellCells);

            // Cell ordering (based on bandCompression). Handles removed cells.
            newCelli = getCellOrder(cellCells, localCellMap);
        }
        else
        {
            // Compact out removed cells
            localCellMap.setSize(cellMap_.size());
            localCellMap = -1;

            newCelli = 0;
            forAll(cellMap_, celli)
            {
                if (!cellRemoved(celli))
                {
                    localCellMap[celli] = newCelli++;
                }
            }
        }

        if (debug)
        {
            Pout<< "Cells : active:" << newCelli
                << "  removed:" << cellMap_.size()-newCelli << endl;
        }

        // Renumber -if cells reordered or -if cells removed
        if (orderCells || (newCelli != cellMap_.size()))
        {
            reorder(localCellMap, cellMap_);
            cellMap_.setCapacity(newCelli);
            renumberReverseMap(localCellMap, reverseCellMap_);

            reorder(localCellMap, cellZone_);
            cellZone_.setCapacity(newCelli);

            renumberKey(localCellMap, cellFromPoint_);
            renumberKey(localCellMap, cellFromEdge_);
            renumberKey(localCellMap, cellFromFace_);

            // Renumber owner/neighbour. Take into account if neighbour suddenly
            // gets lower cell than owner.
            forAll(faceOwner_, facei)
            {
                label own = faceOwner_[facei];
                label nei = faceNeighbour_[facei];

                if (own >= 0)
                {
                    // Update owner
                    faceOwner_[facei] = localCellMap[own];

                    if (nei >= 0)
                    {
                        // Update neighbour.
                        faceNeighbour_[facei] = localCellMap[nei];

                        // Check if face needs reversing.
                        if
                        (
                            faceNeighbour_[facei] >= 0
                         && faceNeighbour_[facei] < faceOwner_[facei]
                        )
                        {
                            faces_[facei].flip();
                            Swap(faceOwner_[facei], faceNeighbour_[facei]);
                            flipFaceFlux_[facei] =
                            (
                                flipFaceFlux_[facei]
                              ? 0
                              : 1
                            );
                            faceZoneFlip_[facei] =
                            (
                                faceZoneFlip_[facei]
                              ? 0
                              : 1
                            );
                        }
                    }
                }
                else if (nei >= 0)
                {
                    // Update neighbour.
                    faceNeighbour_[facei] = localCellMap[nei];
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
        label facei = faceLabels[i];

        if (internalFacesOnly == mesh.isInternalFace(facei))
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
            label facei = faceLabels[i];

            if (internalFacesOnly == mesh.isInternalFace(facei))
            {
                collectedFaces[nFaces++] = facei;
            }
        }
    }

    return collectedFaces;
}


// Calculate pointMap per patch (so from patch point label to old patch point
// label)
void Foam::polyTopoChange::calcPatchPointMap
(
    const List<Map<label>>& oldPatchMeshPointMaps,
    const polyBoundaryMesh& boundary,
    labelListList& patchPointMap
) const
{
    patchPointMap.setSize(boundary.size());

    forAll(boundary, patchi)
    {
        const labelList& meshPoints = boundary[patchi].meshPoints();

        const Map<label>& oldMeshPointMap = oldPatchMeshPointMaps[patchi];

        labelList& curPatchPointRnb = patchPointMap[patchi];

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
            label newFacei = iter.key();

            if (region_[newFacei] == -1)
            {
                // Get internal faces using point on old mesh
                facesFromPoints[nFacesFromPoints++] = objectMap
                (
                    newFacei,
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
                    newFacei,
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
            label newFacei = iter.key();

            if (region_[newFacei] == -1)
            {
                // Get internal faces using edge on old mesh
                facesFromEdges[nFacesFromEdges++] = objectMap
                (
                    newFacei,
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
                    newFacei,
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
            label oldFacei = iter();

            if (mesh.isInternalFace(oldFacei))
            {
                twoCells[0] = mesh.faceOwner()[oldFacei];
                twoCells[1] = mesh.faceNeighbour()[oldFacei];
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
                    labelList(1, mesh.faceOwner()[oldFacei])
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
                FatalErrorInFunction
                    << "Illegal zoneID " << zoneI << " for point "
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

        // Reset the addressing on the zone
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
                FatalErrorInFunction
                    << "Illegal zoneID " << zoneI << " for face "
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
            label facei = iter.key();

            label index = nFaces[zoneI]++;

            addressing[zoneI][index] = facei;
            flipMode[zoneI][index] = faceZoneFlip_[facei];
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


        // Reset the addressing on the zone
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

        forAll(cellZone_, celli)
        {
            label zoneI = cellZone_[celli];

            if (zoneI >= cellZones.size())
            {
                FatalErrorInFunction
                    << "Illegal zoneID " << zoneI << " for cell "
                    << celli << abort(FatalError);
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

        forAll(cellZone_, celli)
        {
            label zoneI = cellZone_[celli];

            if (zoneI >= 0)
            {
                addressing[zoneI][nCells[zoneI]++] = celli;
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

        // Reset the addressing on the zone
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
    const List<Map<label>>& oldFaceZoneMeshPointMaps,
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

        forAll(newZoneMeshPoints, pointi)
        {
            if (newZoneMeshPoints[pointi] < pointMap_.size())
            {
                Map<label>::const_iterator ozmpmIter =
                    oldZoneMeshPointMap.find
                    (
                        pointMap_[newZoneMeshPoints[pointi]]
                    );

                if (ozmpmIter != oldZoneMeshPointMap.end())
                {
                    curFzPointRnb[pointi] = ozmpmIter();
                }
                else
                {
                    curFzPointRnb[pointi] = -1;
                }
            }
            else
            {
                curFzPointRnb[pointi] = -1;
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

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // Send ordering
    forAll(boundary, patchi)
    {
        if (syncParallel || !isA<processorPolyPatch>(boundary[patchi]))
        {
            boundary[patchi].initOrder
            (
                pBufs,
                primitivePatch
                (
                    SubList<face>
                    (
                        faces_,
                        patchSizes[patchi],
                        patchStarts[patchi]
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

    forAll(boundary, patchi)
    {
        if (syncParallel || !isA<processorPolyPatch>(boundary[patchi]))
        {
            labelList patchFaceMap(patchSizes[patchi], -1);
            labelList patchFaceRotation(patchSizes[patchi], 0);

            bool changed = boundary[patchi].order
            (
                pBufs,
                primitivePatch
                (
                    SubList<face>
                    (
                        faces_,
                        patchSizes[patchi],
                        patchStarts[patchi]
                    ),
                    points
                ),
                patchFaceMap,
                patchFaceRotation
            );

            if (changed)
            {
                // Merge patch face reordering into mesh face reordering table
                label start = patchStarts[patchi];

                forAll(patchFaceMap, patchFacei)
                {
                    oldToNew[patchFacei + start] =
                        start + patchFaceMap[patchFacei];
                }

                forAll(patchFaceRotation, patchFacei)
                {
                    rotation[patchFacei + start] =
                        patchFaceRotation[patchFacei];
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
        forAll(rotation, facei)
        {
            if (rotation[facei] != 0)
            {
                inplaceRotateList<List, label>(faces_[facei], rotation[facei]);
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
    List<Map<label>>& oldPatchMeshPointMaps,
    labelList& oldPatchNMeshPoints,
    labelList& oldPatchStarts,
    List<Map<label>>& oldFaceZoneMeshPointMaps
)
{
    if (mesh.boundaryMesh().size() != nPatches_)
    {
        FatalErrorInFunction
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

    forAll(boundary, patchi)
    {
        // Copy old face zone mesh point maps
        oldPatchMeshPointMaps[patchi] = boundary[patchi].meshPointMap();
        oldPatchNMeshPoints[patchi] = boundary[patchi].meshPoints().size();
        oldPatchStarts[patchi] = boundary[patchi].start();
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
        for (label pointi = 0; pointi < mesh.nPoints(); pointi++)
        {
            addPoint
            (
                points[pointi],
                pointi,
                newZoneID[pointi],
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
                label celli = cellLabels[j];

                if (newZoneID[celli] != -1)
                {
                    WarningInFunction
                        << "Cell:" << celli
                        << " centre:" << mesh.cellCentres()[celli]
                        << " is in two zones:"
                        << cellZones[newZoneID[celli]].name()
                        << " and " << cellZones[zoneI].name() << endl
                        << "    This is not supported."
                        << " Continuing with first zone only." << endl;
                }
                else
                {
                    newZoneID[celli] = cellZoneMap[zoneI];
                }
            }
        }

        // Add cells in mesh order
        for (label celli = 0; celli < nAllCells; celli++)
        {
            // Add cell from cell
            addCell(-1, -1, -1, celli, newZoneID[celli]);
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
        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            addFace
            (
                faces[facei],
                faceOwner[facei],
                faceNeighbour[facei],
                -1,                         // masterPointID
                -1,                         // masterEdgeID
                facei,                      // masterFaceID
                false,                      // flipFaceFlux
                -1,                         // patchID
                newZoneID[facei],           // zoneID
                zoneFlip[facei]             // zoneFlip
            );
        }

        // 2. Patch faces
        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];

            if (pp.start() != faces_.size())
            {
                FatalErrorInFunction
                    << "Problem : "
                    << "Patch " << pp.name() << " starts at " << pp.start()
                    << endl
                    << "Current face counter at " << faces_.size() << endl
                    << "Are patches in incremental order?"
                    << abort(FatalError);
            }
            forAll(pp, patchFacei)
            {
                label facei = pp.start() + patchFacei;

                addFace
                (
                    faces[facei],
                    faceOwner[facei],
                    -1,                         // neighbour
                    -1,                         // masterPointID
                    -1,                         // masterEdgeID
                    facei,                      // masterFaceID
                    false,                      // flipFaceFlux
                    patchMap[patchi],           // patchID
                    newZoneID[facei],           // zoneID
                    zoneFlip[facei]             // zoneFlip
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
        FatalErrorInFunction
            << "Unknown type of topoChange: " << action.type()
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
    label pointi = points_.size();

    points_.append(pt);
    pointMap_.append(masterPointID);
    reversePointMap_.append(pointi);

    if (zoneID >= 0)
    {
        pointZone_.insert(pointi, zoneID);
    }

    if (!inCell)
    {
        retiredPoints_.insert(pointi);
    }

    return pointi;
}


void Foam::polyTopoChange::modifyPoint
(
    const label pointi,
    const point& pt,
    const label newZoneID,
    const bool inCell
)
{
    if (pointi < 0 || pointi >= points_.size())
    {
        FatalErrorInFunction
            << "illegal point label " << pointi << endl
            << "Valid point labels are 0 .. " << points_.size()-1
            << abort(FatalError);
    }
    if (pointRemoved(pointi) || pointMap_[pointi] == -1)
    {
        FatalErrorInFunction
            << "point " << pointi << " already marked for removal"
            << abort(FatalError);
    }
    points_[pointi] = pt;

    Map<label>::iterator pointFnd = pointZone_.find(pointi);

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
        pointZone_.insert(pointi, newZoneID);
    }

    if (inCell)
    {
        retiredPoints_.erase(pointi);
    }
    else
    {
        retiredPoints_.insert(pointi);
    }
}


void Foam::polyTopoChange::movePoints(const pointField& newPoints)
{
    if (newPoints.size() != points_.size())
    {
        FatalErrorInFunction
            << "illegal pointField size." << endl
            << "Size:" << newPoints.size() << endl
            << "Points in mesh:" << points_.size()
            << abort(FatalError);
    }

    forAll(points_, pointi)
    {
        points_[pointi] = newPoints[pointi];
    }
}


void Foam::polyTopoChange::removePoint
(
    const label pointi,
    const label mergePointi
)
{
    if (pointi < 0 || pointi >= points_.size())
    {
        FatalErrorInFunction
            << "illegal point label " << pointi << endl
            << "Valid point labels are 0 .. " << points_.size()-1
            << abort(FatalError);
    }

    if
    (
        strict_
     && (pointRemoved(pointi) || pointMap_[pointi] == -1)
    )
    {
        FatalErrorInFunction
            << "point " << pointi << " already marked for removal" << nl
            << "Point:" << points_[pointi] << " pointMap:" << pointMap_[pointi]
            << abort(FatalError);
    }

    if (pointi == mergePointi)
    {
        FatalErrorInFunction
            << "Cannot remove/merge point " << pointi << " onto itself."
            << abort(FatalError);
    }

    points_[pointi] = point::max;
    pointMap_[pointi] = -1;
    if (mergePointi >= 0)
    {
        reversePointMap_[pointi] = -mergePointi-2;
    }
    else
    {
        reversePointMap_[pointi] = -1;
    }
    pointZone_.erase(pointi);
    retiredPoints_.erase(pointi);
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

    label facei = faces_.size();

    faces_.append(f);
    region_.append(patchID);
    faceOwner_.append(own);
    faceNeighbour_.append(nei);

    if (masterPointID >= 0)
    {
        faceMap_.append(-1);
        faceFromPoint_.insert(facei, masterPointID);
    }
    else if (masterEdgeID >= 0)
    {
        faceMap_.append(-1);
        faceFromEdge_.insert(facei, masterEdgeID);
    }
    else if (masterFaceID >= 0)
    {
        faceMap_.append(masterFaceID);
    }
    else
    {
        // Allow inflate-from-nothing?
        //FatalErrorInFunction
        //    << "Need to specify a master point, edge or face"
        //    << "face:" << f << " own:" << own << " nei:" << nei
        //    << abort(FatalError);
        faceMap_.append(-1);
    }
    reverseFaceMap_.append(facei);

    flipFaceFlux_[facei] = (flipFaceFlux ? 1 : 0);

    if (zoneID >= 0)
    {
        faceZone_.insert(facei, zoneID);
    }
    faceZoneFlip_[facei] = (zoneFlip ? 1 : 0);

    return facei;
}


void Foam::polyTopoChange::modifyFace
(
    const face& f,
    const label facei,
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
        checkFace(f, facei, own, nei, patchID, zoneID);
    }

    faces_[facei] = f;
    faceOwner_[facei] = own;
    faceNeighbour_[facei] = nei;
    region_[facei] = patchID;

    flipFaceFlux_[facei] = (flipFaceFlux ? 1 : 0);

    Map<label>::iterator faceFnd = faceZone_.find(facei);

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
        faceZone_.insert(facei, zoneID);
    }
    faceZoneFlip_[facei] = (zoneFlip ? 1 : 0);
}


void Foam::polyTopoChange::removeFace(const label facei, const label mergeFacei)
{
    if (facei < 0 || facei >= faces_.size())
    {
        FatalErrorInFunction
            << "illegal face label " << facei << endl
            << "Valid face labels are 0 .. " << faces_.size()-1
            << abort(FatalError);
    }

    if
    (
        strict_
     && (faceRemoved(facei) || faceMap_[facei] == -1)
    )
    {
        FatalErrorInFunction
            << "face " << facei
            << " already marked for removal"
            << abort(FatalError);
    }

    faces_[facei].setSize(0);
    region_[facei] = -1;
    faceOwner_[facei] = -1;
    faceNeighbour_[facei] = -1;
    faceMap_[facei] = -1;
    if (mergeFacei >= 0)
    {
        reverseFaceMap_[facei] = -mergeFacei-2;
    }
    else
    {
        reverseFaceMap_[facei] = -1;
    }
    faceFromEdge_.erase(facei);
    faceFromPoint_.erase(facei);
    flipFaceFlux_[facei] = 0;
    faceZone_.erase(facei);
    faceZoneFlip_[facei] = 0;
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
    label celli = cellMap_.size();

    if (masterPointID >= 0)
    {
        cellMap_.append(-1);
        cellFromPoint_.insert(celli, masterPointID);
    }
    else if (masterEdgeID >= 0)
    {
        cellMap_.append(-1);
        cellFromEdge_.insert(celli, masterEdgeID);
    }
    else if (masterFaceID >= 0)
    {
        cellMap_.append(-1);
        cellFromFace_.insert(celli, masterFaceID);
    }
    else
    {
        cellMap_.append(masterCellID);
    }
    reverseCellMap_.append(celli);
    cellZone_.append(zoneID);

    return celli;
}


void Foam::polyTopoChange::modifyCell
(
    const label celli,
    const label zoneID
)
{
    cellZone_[celli] = zoneID;
}


void Foam::polyTopoChange::removeCell(const label celli, const label mergeCelli)
{
    if (celli < 0 || celli >= cellMap_.size())
    {
        FatalErrorInFunction
            << "illegal cell label " << celli << endl
            << "Valid cell labels are 0 .. " << cellMap_.size()-1
            << abort(FatalError);
    }

    if (strict_ && cellMap_[celli] == -2)
    {
        FatalErrorInFunction
            << "cell " << celli
            << " already marked for removal"
            << abort(FatalError);
    }

    cellMap_[celli] = -2;
    if (mergeCelli >= 0)
    {
        reverseCellMap_[celli] = -mergeCelli-2;
    }
    else
    {
        reverseCellMap_[celli] = -1;
    }
    cellFromPoint_.erase(celli);
    cellFromEdge_.erase(celli);
    cellFromFace_.erase(celli);
    cellZone_[celli] = -1;
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
    List<Map<label>> oldPatchMeshPointMaps;
    labelList oldPatchNMeshPoints;
    labelList oldPatchStarts;
    List<Map<label>> oldFaceZoneMeshPointMaps;

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

        forAll(pointMap_, newPointi)
        {
            label oldPointi = pointMap_[newPointi];

            if (oldPointi >= 0)
            {
                renumberedMeshPoints[newPointi] = mesh.points()[oldPointi];
            }
            else
            {
                renumberedMeshPoints[newPointi] = Zero;
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
    List<Map<label>> oldPatchMeshPointMaps;
    labelList oldPatchNMeshPoints;
    labelList oldPatchStarts;
    List<Map<label>> oldFaceZoneMeshPointMaps;

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

        forAll(oldPatches, patchi)
        {
            newBoundary[patchi] = oldPatches[patchi].clone
            (
                newMesh.boundaryMesh(),
                patchi,
                patchSizes[patchi],
                patchStarts[patchi]
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
