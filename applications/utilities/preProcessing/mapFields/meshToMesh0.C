/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
#include "processorFvPatch.H"
#include "demandDrivenData.H"
#include "treeDataCell.H"
#include "treeDataFace.H"
#include "tetOverlapVolume.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(meshToMesh0, 0);
}

const Foam::scalar Foam::meshToMesh0::directHitTol = 1e-5;






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


void Foam::meshToMesh0::calculateInverseDistanceWeights() const
{
    if (debug)
    {
        InfoInFunction
            << "Calculating inverse distance weighting factors" << endl;
    }

    if (inverseDistanceWeightsPtr_)
    {
        FatalErrorInFunction
            << "weighting factors already calculated"
            << exit(FatalError);
    }

    //- Initialise overlap volume to zero
    V_ = 0.0;

    inverseDistanceWeightsPtr_ = new scalarListList(toMesh_.nCells());
    scalarListList& invDistCoeffs = *inverseDistanceWeightsPtr_;

    // get reference to source mesh data
    const labelListList& cc = fromMesh_.cellCells();
    const vectorField& centreFrom = fromMesh_.C();
    const vectorField& centreTo = toMesh_.C();

    forAll(cellAddressing_, celli)
    {
        if (cellAddressing_[celli] != -1)
        {
            const vector& target = centreTo[celli];
            scalar m = mag(target - centreFrom[cellAddressing_[celli]]);

            const labelList& neighbours = cc[cellAddressing_[celli]];

            // if the nearest cell is a boundary cell or there is a direct hit,
            // pick up the value
            label directCelli = -1;
            if (m < directHitTol || neighbours.empty())
            {
                directCelli = celli;
            }
            else
            {
                forAll(neighbours, ni)
                {
                    scalar nm = mag(target - centreFrom[neighbours[ni]]);
                    if (nm < directHitTol)
                    {
                        directCelli = neighbours[ni];
                        break;
                    }
                }
            }


            if (directCelli != -1)
            {
                // Direct hit
                invDistCoeffs[directCelli].setSize(1);
                invDistCoeffs[directCelli][0] = 1.0;
                V_ += fromMesh_.V()[cellAddressing_[directCelli]];
            }
            else
            {
                invDistCoeffs[celli].setSize(neighbours.size() + 1);

                // The first coefficient corresponds to the centre cell.
                // The rest is ordered in the same way as the cellCells list.
                scalar invDist = 1.0/m;
                invDistCoeffs[celli][0] = invDist;
                scalar sumInvDist = invDist;

                // now add the neighbours
                forAll(neighbours, ni)
                {
                    invDist = 1.0/mag(target - centreFrom[neighbours[ni]]);
                    invDistCoeffs[celli][ni + 1] = invDist;
                    sumInvDist += invDist;
                }

                // divide by the total inverse-distance
                forAll(invDistCoeffs[celli], i)
                {
                    invDistCoeffs[celli][i] /= sumInvDist;
                }


                V_ +=
                    invDistCoeffs[celli][0]
                   *fromMesh_.V()[cellAddressing_[celli]];
                for (label i = 1; i < invDistCoeffs[celli].size(); i++)
                {
                    V_ +=
                        invDistCoeffs[celli][i]*fromMesh_.V()[neighbours[i-1]];
                }
            }
        }
    }
}


void Foam::meshToMesh0::calculateInverseVolumeWeights() const
{
    if (debug)
    {
        InfoInFunction
            << "Calculating inverse volume weighting factors" << endl;
    }

    if (inverseVolumeWeightsPtr_)
    {
        FatalErrorInFunction
            << "weighting factors already calculated"
            << exit(FatalError);
    }

    //- Initialise overlap volume to zero
    V_ = 0.0;

    inverseVolumeWeightsPtr_ = new scalarListList(toMesh_.nCells());
    scalarListList& invVolCoeffs = *inverseVolumeWeightsPtr_;

    const labelListList& cellToCell = cellToCellAddressing();

    tetOverlapVolume overlapEngine;

    forAll(cellToCell, celli)
    {
        const labelList& overlapCells = cellToCell[celli];

        if (overlapCells.size() > 0)
        {
            invVolCoeffs[celli].setSize(overlapCells.size());

            forAll(overlapCells, j)
            {
                label cellFrom = overlapCells[j];
                treeBoundBox bbFromMesh
                (
                    pointField
                    (
                        fromMesh_.points(),
                        fromMesh_.cellPoints()[cellFrom]
                    )
                );

                scalar v = overlapEngine.cellCellOverlapVolumeMinDecomp
                (
                    toMesh_,
                    celli,

                    fromMesh_,
                    cellFrom,
                    bbFromMesh
                );
                invVolCoeffs[celli][j] = v/toMesh_.V()[celli];

                V_ += v;
            }
        }
    }
}


void Foam::meshToMesh0::calculateCellToCellAddressing() const
{
    if (debug)
    {
        InfoInFunction
            << "Calculating cell to cell addressing" << endl;
    }

    if (cellToCellAddressingPtr_)
    {
        FatalErrorInFunction
            << "addressing already calculated"
            << exit(FatalError);
    }

    //- Initialise overlap volume to zero
    V_ = 0.0;

    tetOverlapVolume overlapEngine;

    cellToCellAddressingPtr_ = new labelListList(toMesh_.nCells());
    labelListList& cellToCell = *cellToCellAddressingPtr_;


    forAll(cellToCell, iTo)
    {
        const labelList overLapCells =
            overlapEngine.overlappingCells(fromMesh_, toMesh_, iTo);
        if (overLapCells.size() > 0)
        {
            cellToCell[iTo].setSize(overLapCells.size());
            forAll(overLapCells, j)
            {
                cellToCell[iTo][j] = overLapCells[j];
                V_ += fromMesh_.V()[overLapCells[j]];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshToMesh0::meshToMesh0
(
    const fvMesh& meshFrom,
    const fvMesh& meshTo,
    const HashTable<word>& patchMap,
    const wordList& cuttingPatchNames
)
:
    fromMesh_(meshFrom),
    toMesh_(meshTo),
    patchMap_(patchMap),
    cellAddressing_(toMesh_.nCells()),
    boundaryAddressing_(toMesh_.boundaryMesh().size()),
    inverseDistanceWeightsPtr_(nullptr),
    inverseVolumeWeightsPtr_(nullptr),
    cellToCellAddressingPtr_(nullptr),
    V_(0.0)
{
    forAll(fromMesh_.boundaryMesh(), patchi)
    {
        fromMeshPatches_.insert
        (
            fromMesh_.boundaryMesh()[patchi].name(),
            patchi
        );
    }

    forAll(toMesh_.boundaryMesh(), patchi)
    {
        toMeshPatches_.insert
        (
            toMesh_.boundaryMesh()[patchi].name(),
            patchi
        );
    }

    forAll(cuttingPatchNames, i)
    {
        if (toMeshPatches_.found(cuttingPatchNames[i]))
        {
            cuttingPatches_.insert
            (
                cuttingPatchNames[i],
                toMeshPatches_.find(cuttingPatchNames[i])()
            );
        }
        else
        {
            WarningInFunction
                << "Cannot find cutting-patch " << cuttingPatchNames[i]
                << " in destination mesh" << endl;
        }
    }

    forAll(toMesh_.boundaryMesh(), patchi)
    {
        // Add the processor patches in the toMesh to the cuttingPatches list
        if (isA<processorPolyPatch>(toMesh_.boundaryMesh()[patchi]))
        {
            cuttingPatches_.insert
            (
                toMesh_.boundaryMesh()[patchi].name(),
                patchi
            );
        }
    }

    calcAddressing();
}


Foam::meshToMesh0::meshToMesh0
(
    const fvMesh& meshFrom,
    const fvMesh& meshTo
)
:
    fromMesh_(meshFrom),
    toMesh_(meshTo),
    cellAddressing_(toMesh_.nCells()),
    boundaryAddressing_(toMesh_.boundaryMesh().size()),
    inverseDistanceWeightsPtr_(nullptr),
    inverseVolumeWeightsPtr_(nullptr),
    cellToCellAddressingPtr_(nullptr),
    V_(0.0)
{
    // check whether both meshes have got the same number
    // of boundary patches
    if (fromMesh_.boundary().size() != toMesh_.boundary().size())
    {
        FatalErrorInFunction
            << "Incompatible meshes: different number of patches, "
            << "fromMesh = " << fromMesh_.boundary().size()
            << ", toMesh = " << toMesh_.boundary().size()
            << exit(FatalError);
    }

    forAll(fromMesh_.boundaryMesh(), patchi)
    {
        if
        (
            fromMesh_.boundaryMesh()[patchi].name()
         != toMesh_.boundaryMesh()[patchi].name()
        )
        {
            FatalErrorInFunction
                << "Incompatible meshes: different patch names for patch "
                << patchi
                << ", fromMesh = " << fromMesh_.boundary()[patchi].name()
                << ", toMesh = " << toMesh_.boundary()[patchi].name()
                << exit(FatalError);
        }

        if
        (
            fromMesh_.boundaryMesh()[patchi].type()
         != toMesh_.boundaryMesh()[patchi].type()
        )
        {
            FatalErrorInFunction
                << "Incompatible meshes: different patch types for patch "
                << patchi
                << ", fromMesh = " << fromMesh_.boundary()[patchi].type()
                << ", toMesh = " << toMesh_.boundary()[patchi].type()
                << exit(FatalError);
        }

        fromMeshPatches_.insert
        (
            fromMesh_.boundaryMesh()[patchi].name(),
            patchi
        );

        toMeshPatches_.insert
        (
            toMesh_.boundaryMesh()[patchi].name(),
            patchi
        );

        patchMap_.insert
        (
            toMesh_.boundaryMesh()[patchi].name(),
            fromMesh_.boundaryMesh()[patchi].name()
        );
    }

    calcAddressing();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshToMesh0::~meshToMesh0()
{
    deleteDemandDrivenData(inverseDistanceWeightsPtr_);
    deleteDemandDrivenData(inverseVolumeWeightsPtr_);
    deleteDemandDrivenData(cellToCellAddressingPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::scalarListList& Foam::meshToMesh0::inverseDistanceWeights() const
{
    if (!inverseDistanceWeightsPtr_)
    {
        calculateInverseDistanceWeights();
    }

    return *inverseDistanceWeightsPtr_;
}


const Foam::scalarListList& Foam::meshToMesh0::inverseVolumeWeights() const
{
    if (!inverseVolumeWeightsPtr_)
    {
        calculateInverseVolumeWeights();
    }

    return *inverseVolumeWeightsPtr_;
}


const Foam::labelListList& Foam::meshToMesh0::cellToCellAddressing() const
{
    if (!cellToCellAddressingPtr_)
    {
        calculateCellToCellAddressing();
    }

    return *cellToCellAddressingPtr_;
}


// ************************************************************************* //
