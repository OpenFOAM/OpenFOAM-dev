/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "meshToMesh0.H"
#include "treeDataCell.H"
#include "treeDataFace.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::meshToMesh0::calcAddressing()
{
    if (debug)
    {
        InfoInFunction
            << "Calculating mesh-to-mesh cell addressing" << endl;
    }

    // Set reference to cells
    const cellList& fromCells = fromMesh_.cells();
    const pointField& fromPoints = fromMesh_.points();

    // In an attempt to preserve the efficiency of linear search and prevent
    // failure, a RESCUE mechanism will be set up. Here, we shall mark all
    // cells next to the solid boundaries. If such a cell is found as the
    // closest, the relationship between the origin and cell will be examined.
    // If the origin is outside the cell, a global n-squared search is
    // triggered.

    // SETTING UP RESCUE

    // Visit all boundaries and mark the cell next to the boundary.

    if (debug)
    {
        InfoInFunction << "Setting up rescue" << endl;
    }

    List<bool> boundaryCell(fromCells.size(), false);

    // Set reference to boundary
    const polyPatchList& patchesFrom = fromMesh_.boundaryMesh();

    forAll(patchesFrom, patchi)
    {
        // Get reference to cells next to the boundary
        const labelUList& bCells = patchesFrom[patchi].faceCells();

        forAll(bCells, facei)
        {
            boundaryCell[bCells[facei]] = true;
        }
    }

    treeBoundBox meshBb(fromPoints);

    scalar typDim = meshBb.avgDim()/(2.0*cbrt(scalar(fromCells.size())));

    treeBoundBox shiftedBb
    (
        meshBb.min(),
        meshBb.max() + vector(typDim, typDim, typDim)
    );

    if (debug)
    {
        Info<< "\nMesh\n"
            << "   bounding box           : " << meshBb << nl
            << "   bounding box (shifted) : " << shiftedBb << nl
            << "   typical dimension      :" << shiftedBb.typDim() << endl;
    }

    indexedOctree<treeDataCell> oc
    (
        treeDataCell(false, fromMesh_, polyMesh::CELL_TETS),
        shiftedBb,      // overall bounding box
        8,              // maxLevel
        10,             // leafsize
        6.0             // duplicity
    );

    if (debug)
    {
        oc.print(Pout, false, 0);
    }

    cellAddresses
    (
        cellAddressing_,
        toMesh_.cellCentres(),
        fromMesh_,
        boundaryCell,
        oc
    );

    forAll(toMesh_.boundaryMesh(), patchi)
    {
        const polyPatch& toPatch = toMesh_.boundaryMesh()[patchi];

        if (cuttingPatches_.found(toPatch.name()))
        {
            boundaryAddressing_[patchi].setSize(toPatch.size());

            cellAddresses
            (
                boundaryAddressing_[patchi],
                toPatch.faceCentres(),
                fromMesh_,
                boundaryCell,
                oc
            );
        }
        else if
        (
            patchMap_.found(toPatch.name())
         && fromMeshPatches_.found(patchMap_.find(toPatch.name())())
        )
        {
            const polyPatch& fromPatch = fromMesh_.boundaryMesh()
            [
                fromMeshPatches_.find(patchMap_.find(toPatch.name())())()
            ];

            if (fromPatch.empty())
            {
                WarningInFunction
                    << "Source patch " << fromPatch.name()
                    << " has no faces. Not performing mapping for it."
                    << endl;
                boundaryAddressing_[patchi].setSize(toPatch.size());
                boundaryAddressing_[patchi] = -1;
            }
            else
            {
                treeBoundBox wallBb(fromPatch.localPoints());
                scalar typDim =
                    wallBb.avgDim()/(2.0*sqrt(scalar(fromPatch.size())));

                treeBoundBox shiftedBb
                (
                    wallBb.min(),
                    wallBb.max() + vector(typDim, typDim, typDim)
                );

                // Note: allow more levels than in meshSearch. Assume patch
                // is not as big as all boundary faces
                indexedOctree<treeDataFace> oc
                (
                    treeDataFace(false, fromPatch),
                    shiftedBb,  // overall search domain
                    12,         // maxLevel
                    10,         // leafsize
                    6.0         // duplicity
                );

                const vectorField::subField centresToBoundary =
                    toPatch.faceCentres();

                boundaryAddressing_[patchi].setSize(toPatch.size());

                const scalar distSqr = sqr(wallBb.mag());

                forAll(toPatch, toi)
                {
                    boundaryAddressing_[patchi][toi] = oc.findNearest
                    (
                        centresToBoundary[toi],
                        distSqr
                    ).index();
                }
            }
        }
    }

    if (debug)
    {
        InfoInFunction
            << "Finished calculating mesh-to-mesh cell addressing" << endl;
    }
}


void Foam::meshToMesh0::cellAddresses
(
    labelList& cellAddressing_,
    const pointField& points,
    const fvMesh& fromMesh,
    const List<bool>& boundaryCell,
    const indexedOctree<treeDataCell>& oc
) const
{
    // The implemented search method is a simple neighbour array search.
    // It starts from a cell zero, searches its neighbours and finds one
    // which is nearer to the target point than the current position.
    // The location of the "current position" is reset to that cell and
    // search through the neighbours continues. The search is finished
    // when all the neighbours of the cell are farther from the target
    // point than the current cell

    // Set curCell label to zero (start)
    label curCell = 0;

    // Set reference to cell to cell addressing
    const vectorField& centresFrom = fromMesh.cellCentres();
    const labelListList& cc = fromMesh.cellCells();

    forAll(points, toi)
    {
        // Pick up target position
        const vector& p = points[toi];

        // Set the sqr-distance
        scalar distSqr = magSqr(p - centresFrom[curCell]);

        bool closer;

        do
        {
            closer = false;

            // Set the current list of neighbouring cells
            const labelList& neighbours = cc[curCell];

            forAll(neighbours, ni)
            {
                const scalar curDistSqr =
                    magSqr(p - centresFrom[neighbours[ni]]);

                // Search through all the neighbours.
                // If the cell is closer, reset current cell and distance
                if (curDistSqr < (1 - small)*distSqr)
                {
                    curCell = neighbours[ni];
                    distSqr = curDistSqr;
                    closer = true;    // A closer neighbour has been found
                }
            }
        } while (closer);

        cellAddressing_[toi] = -1;

        // Check point is actually in the nearest cell
        if (fromMesh.pointInCell(p, curCell))
        {
            cellAddressing_[toi] = curCell;
        }
        else
        {
            // If curCell is a boundary cell then the point maybe either outside
            // the domain or in an other region of the doamin, either way use
            // the octree search to find it.
            if (boundaryCell[curCell])
            {
                cellAddressing_[toi] = oc.findInside(p);
                if (cellAddressing_[toi] != -1)
                {
                    curCell = cellAddressing_[toi];
                }
            }
            else
            {
                // If not on the boundary search the neighbours
                bool found = false;

                // Set the current list of neighbouring cells
                const labelList& neighbours = cc[curCell];

                forAll(neighbours, ni)
                {
                    // Search through all the neighbours.
                    // If point is in neighbour reset current cell
                    if (fromMesh.pointInCell(p, neighbours[ni]))
                    {
                        cellAddressing_[toi] = neighbours[ni];
                        found = true;
                        break;
                    }
                }

                if (!found)
                {
                    // If still not found search the neighbour-neighbours

                    // Set the current list of neighbouring cells
                    const labelList& neighbours = cc[curCell];

                    forAll(neighbours, ni)
                    {
                        // Set the current list of neighbour-neighbouring cells
                        const labelList& nn = cc[neighbours[ni]];

                        forAll(nn, ni)
                        {
                            // Search through all the neighbours.
                            // If point is in neighbour reset current cell
                            if (fromMesh.pointInCell(p, nn[ni]))
                            {
                                cellAddressing_[toi] = nn[ni];
                                found = true;
                                break;
                            }
                        }
                        if (found) break;
                    }
                }

                if (!found)
                {
                    // Still not found so use the octree
                    cellAddressing_[toi] = oc.findInside(p);

                    if (cellAddressing_[toi] != -1)
                    {
                        curCell = cellAddressing_[toi];
                    }
                }
            }
        }
    }
}


// ************************************************************************* //
