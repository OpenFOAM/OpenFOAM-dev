/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    label& nSplit,
    label& nInserted,
    label& nMerge,
    label& nRemove
)
{
    nSplit = 0;
    nInserted = 0;
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
                // Added from another cell
                nSplit++;
            }
        }
        else if (oldCelli == -1)
        {
            // Created from nothing
            nInserted++;
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
    const label patchi
) const
{
    if (nei == -1)
    {
        if (own == -1)
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
        // SortableList<label> nbr(nFaces);
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

        // nbr.sort();
        order.setSize(nFaces);
        sortedOrder(nbr, order);

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

    inplaceReorder(oldToNew, flipFaceFlux_);
    flipFaceFlux_.setCapacity(newSize);
}


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

        renumberKey(localPointMap, oldPoints_);
        renumber(localPointMap, retiredPoints_);

        // Use map to relabel face vertices
        forAll(faces_, facei)
        {
            face& f = faces_[facei];

            // labelList oldF(f);
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
    labelList oldToNew(identityMap(faces_.size()));

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
    List<objectMap>& facesFromFaces,
    List<objectMap>& cellsFromCells,
    List<Map<label>>& oldPatchMeshPointMaps,
    labelList& oldPatchNMeshPoints,
    labelList& oldPatchSizes,
    labelList& oldPatchStarts
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


    // Calculate merging maps
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // These are for the new face(/point/cell) the old faces whose value
    // needs to be
    // averaged/summed to get the new value. These old faces come either from
    // merged old faces (face remove into other face).
    // As an additional complexity will use only internal faces
    // to create new value for internal face and vice versa only patch
    // faces to to create patch face value.

    // For point only point merging
    getMergeSets
    (
        reversePointMap_,
        pointMap_,
        pointsFromPoints
    );

    getMergeSets
    (
        reverseFaceMap_,
        faceMap_,
        facesFromFaces
    );

    getMergeSets
    (
        reverseCellMap_,
        cellMap_,
        cellsFromCells
    );

    const polyBoundaryMesh& boundary = mesh.boundaryMesh();

    // Grab patch mesh point maps
    oldPatchMeshPointMaps.setSize(boundary.size());
    oldPatchNMeshPoints.setSize(boundary.size());
    oldPatchSizes.setSize(boundary.size());
    oldPatchStarts.setSize(boundary.size());

    forAll(boundary, patchi)
    {
        oldPatchMeshPointMaps[patchi] = boundary[patchi].meshPointMap();
        oldPatchNMeshPoints[patchi] = boundary[patchi].meshPoints().size();
        oldPatchSizes[patchi] = boundary[patchi].size();
        oldPatchStarts[patchi] = boundary[patchi].start();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyTopoChange::polyTopoChange(const label nPatches, const bool strict)
:
    strict_(strict),
    nPatches_(nPatches),
    points_(0),
    pointMap_(0),
    reversePointMap_(0),
    retiredPoints_(0),
    oldPoints_(0),
    faces_(0),
    region_(0),
    faceOwner_(0),
    faceNeighbour_(0),
    faceMap_(0),
    reverseFaceMap_(0),
    flipFaceFlux_(0),
    nActiveFaces_(0),
    cellMap_(0),
    reverseCellMap_(0)
{}


Foam::polyTopoChange::polyTopoChange
(
    const polyMesh& mesh,
    const bool strict
)
:
    strict_(strict),
    nPatches_(mesh.boundaryMesh().size()),
    points_(0),
    pointMap_(0),
    reversePointMap_(0),
    retiredPoints_(0),
    oldPoints_(0),
    faces_(0),
    region_(0),
    faceOwner_(0),
    faceNeighbour_(0),
    faceMap_(0),
    reverseFaceMap_(0),
    flipFaceFlux_(0),
    nActiveFaces_(0),
    cellMap_(0),
    reverseCellMap_(0)
{
    // Add points
    {
        const pointField& points = mesh.points();

        // Extend
        points_.setCapacity(points_.size() + points.size());
        pointMap_.setCapacity(pointMap_.size() + points.size());
        reversePointMap_.setCapacity(reversePointMap_.size() + points.size());
        // No need to extend oldPoints_

        // Add points in mesh order
        for (label pointi = 0; pointi < mesh.nPoints(); pointi++)
        {
            addPoint(points[pointi], pointi, true);
        }
    }

    // Add cells
    {
        // Resize

        // Note: polyMesh does not allow retired cells anymore. So allCells
        // always equals nCells
        label nAllCells = mesh.nCells();

        cellMap_.setCapacity(cellMap_.size() + nAllCells);
        reverseCellMap_.setCapacity(reverseCellMap_.size() + nAllCells);

        // Add cells in mesh order
        for (label celli = 0; celli < nAllCells; celli++)
        {
            addCell(celli);
        }
    }

    // Add faces
    {
        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        const faceList& faces = mesh.faces();
        const labelList& faceOwner = mesh.faceOwner();
        const labelList& faceNeighbour = mesh.faceNeighbour();

        // Resize
        label nAllFaces = mesh.faces().size();

        faces_.setCapacity(faces_.size() + nAllFaces);
        region_.setCapacity(region_.size() + nAllFaces);
        faceOwner_.setCapacity(faceOwner_.size() + nAllFaces);
        faceNeighbour_.setCapacity(faceNeighbour_.size() + nAllFaces);
        faceMap_.setCapacity(faceMap_.size() + nAllFaces);
        reverseFaceMap_.setCapacity(reverseFaceMap_.size() + nAllFaces);
        flipFaceFlux_.setCapacity(faces_.size() + nAllFaces);


        // Add faces in order

        // 1. Add internal faces in increasing face order
        for (label facei=0; facei<mesh.nInternalFaces(); facei++)
        {
            addFace
            (
                faces[facei],
                faceOwner[facei],
                faceNeighbour[facei],
                facei,                      // masterFaceID
                false,                      // flipFaceFlux
                -1                          // patchID
            );
        }

        // Find patch order with increasing face order
        SortableList<label> patchStartOrder(patches.size());
        forAll(patches, patchi)
        {
            patchStartOrder[patchi] = patches[patchi].start();
        }
        patchStartOrder.sort();

        // 2. Add patch faces in increasing face order
        forAll(patchStartOrder, i)
        {
            const label patchi = patchStartOrder.indices()[i];
            const polyPatch& pp = patches[patchi];

            forAll(pp, patchFacei)
            {
                const label facei = pp.start() + patchFacei;

                addFace
                (
                    faces[facei],
                    faceOwner[facei],
                    -1,                         // neighbour
                    facei,                      // masterFaceID
                    false,                      // flipFaceFlux
                    patchi                      // patchID
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyTopoChange::clear()
{
    points_.clearStorage();
    pointMap_.clearStorage();
    reversePointMap_.clearStorage();
    retiredPoints_.clearStorage();
    oldPoints_.clearStorage();

    faces_.clearStorage();
    region_.clearStorage();
    faceOwner_.clearStorage();
    faceNeighbour_.clearStorage();
    faceMap_.clearStorage();
    reverseFaceMap_.clearStorage();
    flipFaceFlux_.clearStorage();
    nActiveFaces_ = 0;

    cellMap_.clearStorage();
    reverseCellMap_.clearStorage();
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

    faces_.setCapacity(nFaces);
    region_.setCapacity(nFaces);
    faceOwner_.setCapacity(nFaces);
    faceNeighbour_.setCapacity(nFaces);
    faceMap_.setCapacity(nFaces);
    reverseFaceMap_.setCapacity(nFaces);
    flipFaceFlux_.setCapacity(nFaces);

    cellMap_.setCapacity(nCells);
    reverseCellMap_.setCapacity(nCells);
}


Foam::label Foam::polyTopoChange::addPoint
(
    const point& pt,
    const label masterPointID,
    const bool inCell
)
{
    label pointi = points_.size();

    points_.append(pt);
    pointMap_.append(masterPointID);
    reversePointMap_.append(pointi);

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

    if (inCell)
    {
        retiredPoints_.erase(pointi);
    }
    else
    {
        retiredPoints_.insert(pointi);
    }

    oldPoints_.erase(pointi);
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
    retiredPoints_.erase(pointi);
    oldPoints_.erase(pointi);
}


Foam::label Foam::polyTopoChange::addFace
(
    const face& f,
    const label own,
    const label nei,
    const label masterFaceID,
    const bool flipFaceFlux,
    const label patchID
)
{
    // Check validity
    if (debug)
    {
        checkFace(f, -1, own, nei, patchID);
    }

    label facei = faces_.size();

    faces_.append(f);
    region_.append(patchID);
    faceOwner_.append(own);
    faceNeighbour_.append(nei);

    if (masterFaceID >= 0)
    {
        faceMap_.append(masterFaceID);
    }
    else
    {
        // Allow insert-from-nothing?
        // FatalErrorInFunction
        //    << "Need to specify a master point, edge or face"
        //    << "face:" << f << " own:" << own << " nei:" << nei
        //    << abort(FatalError);
        faceMap_.append(-1);
    }
    reverseFaceMap_.append(facei);

    flipFaceFlux_[facei] = (flipFaceFlux ? 1 : 0);

    return facei;
}


void Foam::polyTopoChange::modifyFace
(
    const face& f,
    const label facei,
    const label own,
    const label nei,
    const bool flipFaceFlux,
    const label patchID
)
{
    // Check validity
    if (debug)
    {
        checkFace(f, facei, own, nei, patchID);
    }

    faces_[facei] = f;
    faceOwner_[facei] = own;
    faceNeighbour_[facei] = nei;
    region_[facei] = patchID;

    flipFaceFlux_[facei] = (flipFaceFlux ? 1 : 0);
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
    flipFaceFlux_[facei] = 0;
}


Foam::label Foam::polyTopoChange::addCell(const label masterCellID)
{
    const label celli = cellMap_.size();
    cellMap_.append(masterCellID);
    reverseCellMap_.append(celli);

    return celli;
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
}


Foam::autoPtr<Foam::polyTopoChangeMap> Foam::polyTopoChange::changeMesh
(
    polyMesh& mesh,
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

    // maps
    List<objectMap> pointsFromPoints;
    List<objectMap> facesFromFaces;
    List<objectMap> cellsFromCells;

    // old mesh info
    List<Map<label>> oldPatchMeshPointMaps;
    labelList oldPatchNMeshPoints;
    labelList oldPatchSizes;
    labelList oldPatchStarts;

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
        facesFromFaces,
        cellsFromCells,
        oldPatchMeshPointMaps,
        oldPatchNMeshPoints,
        oldPatchSizes,
        oldPatchStarts
    );

    const label nOldPoints(mesh.nPoints());
    const label nOldFaces(mesh.nFaces());
    const label nOldCells(mesh.nCells());
    autoPtr<scalarField> oldCellVolumes(new scalarField(mesh.cellVolumes()));


    // Change the mesh
    // ~~~~~~~~~~~~~~~
    // This will invalidate any addressing so better make sure you have
    // all the information you need!!!

    // Set new points.
    mesh.resetPrimitives
    (
        move(newPoints),
        move(faces_),
        move(faceOwner_),
        move(faceNeighbour_),
        patchSizes,
        patchStarts,
        syncParallel
    );

    // Clear out primitives
    {
        retiredPoints_.clearStorage();
        oldPoints_.clearStorage();
        region_.clearStorage();
    }


    if (debug)
    {
        // Some stats on changes
        label nSplit, nInserted, nMerge, nRemove;
        countMap
        (
            pointMap_,
            reversePointMap_,
            nSplit,
            nInserted,
            nMerge,
            nRemove
        );
        Pout<< "Points:"
            << "  added(from point):" << nSplit
            << "  added(from nothing):" << nInserted
            << "  merged(into other point):" << nMerge
            << "  removed:" << nRemove
            << nl;

        countMap(faceMap_, reverseFaceMap_, nSplit, nInserted, nMerge, nRemove);
        Pout<< "Faces:"
            << "  added(from face):" << nSplit
            << "  added(from nothing):" << nInserted
            << "  merged(into other face):" << nMerge
            << "  removed:" << nRemove
            << nl;

        countMap(cellMap_, reverseCellMap_, nSplit, nInserted, nMerge, nRemove);
        Pout<< "Cells:"
            << "  added(from cell):" << nSplit
            << "  added(from nothing):" << nInserted
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

    labelHashSet flipFaceFluxSet(getSetIndices(flipFaceFlux_));

    return autoPtr<polyTopoChangeMap>
    (
        new polyTopoChangeMap
        (
            mesh,
            nOldPoints,
            nOldFaces,
            nOldCells,

            move(pointMap_),
            move(pointsFromPoints),

            move(faceMap_),
            move(facesFromFaces),

            move(cellMap_),
            move(cellsFromCells),

            move(reversePointMap_),
            move(reverseFaceMap_),
            move(reverseCellMap_),

            move(flipFaceFluxSet),

            move(patchPointMap),

            move(oldPatchSizes),
            move(oldPatchStarts),
            move(oldPatchNMeshPoints),

            move(oldCellVolumes)
        )
    );

    // At this point all member DynamicList (pointMap_, cellMap_ etc.) will
    // be invalid.
}


Foam::autoPtr<Foam::polyTopoChangeMap> Foam::polyTopoChange::makeMesh
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
        Pout<< "polyTopoChange::makeMesh"
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
    // maps
    List<objectMap> pointsFromPoints;
    List<objectMap> facesFromFaces;
    List<objectMap> cellsFromCells;

    // old mesh info
    List<Map<label>> oldPatchMeshPointMaps;
    labelList oldPatchNMeshPoints;
    labelList oldPatchSizes;
    labelList oldPatchStarts;

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
        facesFromFaces,
        cellsFromCells,
        oldPatchMeshPointMaps,
        oldPatchNMeshPoints,
        oldPatchSizes,
        oldPatchStarts
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
            move(newPoints),
            move(faces_),
            move(faceOwner_),
            move(faceNeighbour_)
        )
    );
    fvMesh& newMesh = newMeshPtr();

    // Clear out primitives
    {
        retiredPoints_.clearStorage();
        oldPoints_.clearStorage();
        region_.clearStorage();
    }


    if (debug)
    {
        // Some stats on changes
        label nSplit, nInserted, nMerge, nRemove;
        countMap
        (
            pointMap_,
            reversePointMap_,
            nSplit,
            nInserted,
            nMerge,
            nRemove
        );

        Pout<< "Points:"
            << "  added(from point):" << nSplit
            << "  added(from nothing):" << nInserted
            << "  merged(into other point):" << nMerge
            << "  removed:" << nRemove
            << nl;

        countMap(faceMap_, reverseFaceMap_, nSplit, nInserted, nMerge, nRemove);
        Pout<< "Faces:"
            << "  added(from face):" << nSplit
            << "  added(from nothing):" << nInserted
            << "  merged(into other face):" << nMerge
            << "  removed:" << nRemove
            << nl;

        countMap(cellMap_, reverseCellMap_, nSplit, nInserted, nMerge, nRemove);
        Pout<< "Cells:"
            << "  added(from cell):" << nSplit
            << "  added(from nothing):" << nInserted
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

    // Copy pointZone from old mesh
    const pointZoneList& oldPointZones = mesh.pointZones();
    List<pointZone*> pZonePtrs(oldPointZones.size());
    {
        forAll(oldPointZones, i)
        {
            pZonePtrs[i] = oldPointZones[i].clone(newMesh.pointZones()).ptr();
        }
    }

    // Copy faceZone from old mesh
    const faceZoneList& oldFaceZones = mesh.faceZones();
    List<faceZone*> fZonePtrs(oldFaceZones.size());
    {
        forAll(oldFaceZones, i)
        {
            fZonePtrs[i] = oldFaceZones[i].clone(newMesh.faceZones()).ptr();
        }
    }

    // Copy cellZone from old mesh
    const cellZoneList& oldCellZones = mesh.cellZones();
    List<cellZone*> cZonePtrs(oldCellZones.size());
    {
        forAll(oldCellZones, i)
        {
            cZonePtrs[i] = oldCellZones[i].clone(newMesh.cellZones()).ptr();
        }
    }

    newMesh.addZones(pZonePtrs, fZonePtrs, cZonePtrs);

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

    if (debug)
    {
        Pout<< "New mesh:" << nl;
        writeMeshStats(mesh, Pout);
    }

    labelHashSet flipFaceFluxSet(getSetIndices(flipFaceFlux_));

    return autoPtr<polyTopoChangeMap>
    (
        new polyTopoChangeMap
        (
            newMesh,
            nOldPoints,
            nOldFaces,
            nOldCells,

            move(pointMap_),
            move(pointsFromPoints),

            move(faceMap_),
            move(facesFromFaces),

            move(cellMap_),
            move(cellsFromCells),

            move(reversePointMap_),
            move(reverseFaceMap_),
            move(reverseCellMap_),

            move(flipFaceFluxSet),

            move(patchPointMap),

            move(oldPatchSizes),
            move(oldPatchStarts),
            move(oldPatchNMeshPoints),

            move(oldCellVolumes)
        )
    );

    // At this point all member DynamicList (pointMap_, cellMap_ etc.) will
    // be invalid.
}


// ************************************************************************* //
