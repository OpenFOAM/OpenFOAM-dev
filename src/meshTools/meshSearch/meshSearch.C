/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "meshSearch.H"
#include "polyMesh.H"
#include "indexedOctree.H"
#include "DynamicList.H"
#include "demandDrivenData.H"
#include "treeDataCell.H"
#include "treeDataFace.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(meshSearch, 0);

scalar meshSearch::tol_ = 1e-3;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::meshSearch::findNearer
(
    const point& sample,
    const pointField& points,
    label& nearestI,
    scalar& nearestDistSqr
)
{
    bool nearer = false;

    forAll(points, pointI)
    {
        scalar distSqr = magSqr(points[pointI] - sample);

        if (distSqr < nearestDistSqr)
        {
            nearestDistSqr = distSqr;
            nearestI = pointI;
            nearer = true;
        }
    }

    return nearer;
}


bool Foam::meshSearch::findNearer
(
    const point& sample,
    const pointField& points,
    const labelList& indices,
    label& nearestI,
    scalar& nearestDistSqr
)
{
    bool nearer = false;

    forAll(indices, i)
    {
        label pointI = indices[i];

        scalar distSqr = magSqr(points[pointI] - sample);

        if (distSqr < nearestDistSqr)
        {
            nearestDistSqr = distSqr;
            nearestI = pointI;
            nearer = true;
        }
    }

    return nearer;
}


// tree based searching
Foam::label Foam::meshSearch::findNearestCellTree(const point& location) const
{
    const indexedOctree<treeDataCell>& tree = cellTree();

    pointIndexHit info = tree.findNearest
    (
        location,
        magSqr(tree.bb().max()-tree.bb().min())
    );

    if (!info.hit())
    {
        info = tree.findNearest(location, Foam::sqr(GREAT));
    }
    return info.index();
}


// linear searching
Foam::label Foam::meshSearch::findNearestCellLinear(const point& location) const
{
    const vectorField& centres = mesh_.cellCentres();

    label nearestIndex = 0;
    scalar minProximity = magSqr(centres[nearestIndex] - location);

    findNearer
    (
        location,
        centres,
        nearestIndex,
        minProximity
    );

    return nearestIndex;
}


// walking from seed
Foam::label Foam::meshSearch::findNearestCellWalk
(
    const point& location,
    const label seedCellI
) const
{
    if (seedCellI < 0)
    {
        FatalErrorIn
        (
            "meshSearch::findNearestCellWalk(const point&, const label)"
        )   << "illegal seedCell:" << seedCellI << exit(FatalError);
    }

    // Walk in direction of face that decreases distance

    label curCellI = seedCellI;
    scalar distanceSqr = magSqr(mesh_.cellCentres()[curCellI] - location);

    bool closer;

    do
    {
        // Try neighbours of curCellI
        closer = findNearer
        (
            location,
            mesh_.cellCentres(),
            mesh_.cellCells()[curCellI],
            curCellI,
            distanceSqr
        );
    } while (closer);

    return curCellI;
}


// tree based searching
Foam::label Foam::meshSearch::findNearestFaceTree(const point& location) const
{
    // Search nearest cell centre.
    const indexedOctree<treeDataCell>& tree = cellTree();

    // Search with decent span
    pointIndexHit info = tree.findNearest
    (
        location,
        magSqr(tree.bb().max()-tree.bb().min())
    );

    if (!info.hit())
    {
        // Search with desparate span
        info = tree.findNearest(location, Foam::sqr(GREAT));
    }


    // Now check any of the faces of the nearest cell
    const vectorField& centres = mesh_.faceCentres();
    const cell& ownFaces = mesh_.cells()[info.index()];

    label nearestFaceI = ownFaces[0];
    scalar minProximity = magSqr(centres[nearestFaceI] - location);

    findNearer
    (
        location,
        centres,
        ownFaces,
        nearestFaceI,
        minProximity
    );

    return nearestFaceI;
}


// linear searching
Foam::label Foam::meshSearch::findNearestFaceLinear(const point& location) const
{
    const vectorField& centres = mesh_.faceCentres();

    label nearestFaceI = 0;
    scalar minProximity = magSqr(centres[nearestFaceI] - location);

    findNearer
    (
        location,
        centres,
        nearestFaceI,
        minProximity
    );

    return nearestFaceI;
}


// walking from seed
Foam::label Foam::meshSearch::findNearestFaceWalk
(
    const point& location,
    const label seedFaceI
) const
{
    if (seedFaceI < 0)
    {
        FatalErrorIn
        (
            "meshSearch::findNearestFaceWalk(const point&, const label)"
        )   << "illegal seedFace:" << seedFaceI << exit(FatalError);
    }

    const vectorField& centres = mesh_.faceCentres();


    // Walk in direction of face that decreases distance

    label curFaceI = seedFaceI;
    scalar distanceSqr = magSqr(centres[curFaceI] - location);

    while (true)
    {
        label betterFaceI = curFaceI;

        findNearer
        (
            location,
            centres,
            mesh_.cells()[mesh_.faceOwner()[curFaceI]],
            betterFaceI,
            distanceSqr
        );

        if (mesh_.isInternalFace(curFaceI))
        {
            findNearer
            (
                location,
                centres,
                mesh_.cells()[mesh_.faceNeighbour()[curFaceI]],
                betterFaceI,
                distanceSqr
            );
        }

        if (betterFaceI == curFaceI)
        {
            break;
        }

        curFaceI = betterFaceI;
    }

    return curFaceI;
}


Foam::label Foam::meshSearch::findCellLinear(const point& location) const
{
    bool cellFound = false;
    label n = 0;

    label cellI = -1;

    while ((!cellFound) && (n < mesh_.nCells()))
    {
        if (mesh_.pointInCell(location, n, cellDecompMode_))
        {
            cellFound = true;
            cellI = n;
        }
        else
        {
            n++;
        }
    }
    if (cellFound)
    {
        return cellI;
    }
    else
    {
        return -1;
    }
}


// walking from seed
Foam::label Foam::meshSearch::findCellWalk
(
    const point& location,
    const label seedCellI
) const
{
    if (seedCellI < 0)
    {
        FatalErrorIn
        (
            "meshSearch::findCellWalk(const point&, const label)"
        )   << "illegal seedCell:" << seedCellI << exit(FatalError);
    }

    if (mesh_.pointInCell(location, seedCellI, cellDecompMode_))
    {
        return seedCellI;
    }

    // Walk in direction of face that decreases distance
    label curCellI = seedCellI;
    scalar nearestDistSqr = magSqr(mesh_.cellCentres()[curCellI] - location);

    while(true)
    {
        // Try neighbours of curCellI

        const cell& cFaces = mesh_.cells()[curCellI];

        label nearestCellI = -1;

        forAll(cFaces, i)
        {
            label faceI = cFaces[i];

            if (mesh_.isInternalFace(faceI))
            {
                label cellI = mesh_.faceOwner()[faceI];
                if (cellI == curCellI)
                {
                    cellI = mesh_.faceNeighbour()[faceI];
                }

                // Check if this is the correct cell
                if (mesh_.pointInCell(location, cellI, cellDecompMode_))
                {
                    return cellI;
                }

                // Also calculate the nearest cell
                scalar distSqr = magSqr(mesh_.cellCentres()[cellI] - location);

                if (distSqr < nearestDistSqr)
                {
                    nearestDistSqr = distSqr;
                    nearestCellI = cellI;
                }
            }
        }

        if (nearestCellI == -1)
        {
            return -1;
        }

        // Continue with the nearest cell
        curCellI = nearestCellI;
    }

    return -1;
}


Foam::label Foam::meshSearch::findNearestBoundaryFaceWalk
(
    const point& location,
    const label seedFaceI
) const
{
    if (seedFaceI < 0)
    {
        FatalErrorIn
        (
            "meshSearch::findNearestBoundaryFaceWalk"
            "(const point&, const label)"
        )   << "illegal seedFace:" << seedFaceI << exit(FatalError);
    }

    // Start off from seedFaceI

    label curFaceI = seedFaceI;

    const face& f =  mesh_.faces()[curFaceI];

    scalar minDist = f.nearestPoint
    (
        location,
        mesh_.points()
    ).distance();

    bool closer;

    do
    {
        closer = false;

        // Search through all neighbouring boundary faces by going
        // across edges

        label lastFaceI = curFaceI;

        const labelList& myEdges = mesh_.faceEdges()[curFaceI];

        forAll(myEdges, myEdgeI)
        {
            const labelList& neighbours = mesh_.edgeFaces()[myEdges[myEdgeI]];

            // Check any face which uses edge, is boundary face and
            // is not curFaceI itself.

            forAll(neighbours, nI)
            {
                label faceI = neighbours[nI];

                if
                (
                    (faceI >= mesh_.nInternalFaces())
                 && (faceI != lastFaceI)
                )
                {
                    const face& f =  mesh_.faces()[faceI];

                    pointHit curHit = f.nearestPoint
                    (
                        location,
                        mesh_.points()
                    );

                    // If the face is closer, reset current face and distance
                    if (curHit.distance() < minDist)
                    {
                        minDist = curHit.distance();
                        curFaceI = faceI;
                        closer = true;  // a closer neighbour has been found
                    }
                }
            }
        }
    } while (closer);

    return curFaceI;
}


Foam::vector Foam::meshSearch::offset
(
    const point& bPoint,
    const label bFaceI,
    const vector& dir
) const
{
    // Get the neighbouring cell
    label ownerCellI = mesh_.faceOwner()[bFaceI];

    const point& c = mesh_.cellCentres()[ownerCellI];

    // Typical dimension: distance from point on face to cell centre
    scalar typDim = mag(c - bPoint);

    return tol_*typDim*dir;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::meshSearch::meshSearch
(
    const polyMesh& mesh,
    const polyMesh::cellRepresentation cellDecompMode
)
:
    mesh_(mesh),
    cellDecompMode_(cellDecompMode)
{
    if (cellDecompMode_ == polyMesh::FACEDIAGTETS)
    {
        // Force construction of face diagonals
        (void)mesh.tetBasePtIs();
    }
}


// Construct with a custom bounding box
Foam::meshSearch::meshSearch
(
    const polyMesh& mesh,
    const treeBoundBox& bb,
    const polyMesh::cellRepresentation cellDecompMode
)
:
    mesh_(mesh),
    cellDecompMode_(cellDecompMode)
{
    overallBbPtr_.reset(new treeBoundBox(bb));

    if (cellDecompMode_ == polyMesh::FACEDIAGTETS)
    {
        // Force construction of face diagonals
        (void)mesh.tetBasePtIs();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshSearch::~meshSearch()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::indexedOctree<Foam::treeDataFace>& Foam::meshSearch::boundaryTree()
 const
{
    if (!boundaryTreePtr_.valid())
    {
        //
        // Construct tree
        //

        if (!overallBbPtr_.valid())
        {
            Random rndGen(261782);
            overallBbPtr_.reset
            (
                new treeBoundBox(mesh_.points())
            );

            treeBoundBox& overallBb = overallBbPtr_();
            // Extend slightly and make 3D
            overallBb = overallBb.extend(rndGen, 1e-4);
            overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
            overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        }

        // all boundary faces (not just walls)
        labelList bndFaces(mesh_.nFaces()-mesh_.nInternalFaces());
        forAll(bndFaces, i)
        {
            bndFaces[i] = mesh_.nInternalFaces() + i;
        }

        boundaryTreePtr_.reset
        (
            new indexedOctree<treeDataFace>
            (
                treeDataFace    // all information needed to search faces
                (
                    false,                      // do not cache bb
                    mesh_,
                    bndFaces                    // boundary faces only
                ),
                overallBbPtr_(),                // overall search domain
                8,                              // maxLevel
                10,                             // leafsize
                3.0                             // duplicity
            )
        );
    }

    return boundaryTreePtr_();
}


const Foam::indexedOctree<Foam::treeDataCell>& Foam::meshSearch::cellTree()
const
{
    if (!cellTreePtr_.valid())
    {
        //
        // Construct tree
        //

        if (!overallBbPtr_.valid())
        {
            Random rndGen(261782);
            overallBbPtr_.reset
            (
                new treeBoundBox(mesh_.points())
            );

            treeBoundBox& overallBb = overallBbPtr_();
            // Extend slightly and make 3D
            overallBb = overallBb.extend(rndGen, 1e-4);
            overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
            overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        }

        cellTreePtr_.reset
        (
            new indexedOctree<treeDataCell>
            (
                treeDataCell
                (
                    false,          // not cache bb
                    mesh_,
                    cellDecompMode_ // cell decomposition mode for inside tests
                ),
                overallBbPtr_(),
                8,              // maxLevel
                10,             // leafsize
                6.0             // duplicity
            )
        );
    }

    return cellTreePtr_();
}


//// Is the point in the cell
//// Works by checking if there is a face inbetween the point and the cell
//// centre.
//// Check for internal uses proper face decomposition or just average normal.
//bool Foam::meshSearch::pointInCell(const point& p, label cellI) const
//{
//    if (faceDecomp_)
//    {
//        const point& ctr = mesh_.cellCentres()[cellI];
//
//        vector dir(p - ctr);
//        scalar magDir = mag(dir);
//
//        // Check if any faces are hit by ray from cell centre to p.
//        // If none -> p is in cell.
//        const labelList& cFaces = mesh_.cells()[cellI];
//
//        // Make sure half_ray does not pick up any faces on the wrong
//        // side of the ray.
//        scalar oldTol = intersection::setPlanarTol(0.0);
//
//        forAll(cFaces, i)
//        {
//            label faceI = cFaces[i];
//
//            pointHit inter = mesh_.faces()[faceI].ray
//            (
//                ctr,
//                dir,
//                mesh_.points(),
//                intersection::HALF_RAY,
//                intersection::VECTOR
//            );
//
//            if (inter.hit())
//            {
//                scalar dist = inter.distance();
//
//                if (dist < magDir)
//                {
//                    // Valid hit. Hit face so point is not in cell.
//                    intersection::setPlanarTol(oldTol);
//
//                    return false;
//                }
//            }
//        }
//
//        intersection::setPlanarTol(oldTol);
//
//        // No face inbetween point and cell centre so point is inside.
//        return true;
//    }
//    else
//    {
//        const labelList& f = mesh_.cells()[cellI];
//        const labelList& owner = mesh_.faceOwner();
//        const vectorField& cf = mesh_.faceCentres();
//        const vectorField& Sf = mesh_.faceAreas();
//
//        forAll(f, facei)
//        {
//            label nFace = f[facei];
//            vector proj = p - cf[nFace];
//            vector normal = Sf[nFace];
//            if (owner[nFace] == cellI)
//            {
//                if ((normal & proj) > 0)
//                {
//                    return false;
//                }
//            }
//            else
//            {
//                if ((normal & proj) < 0)
//                {
//                    return false;
//                }
//            }
//        }
//
//        return true;
//    }
//}


Foam::label Foam::meshSearch::findNearestCell
(
    const point& location,
    const label seedCellI,
    const bool useTreeSearch
) const
{
    if (seedCellI == -1)
    {
        if (useTreeSearch)
        {
            return findNearestCellTree(location);
        }
        else
        {
            return findNearestCellLinear(location);
        }
    }
    else
    {
        return findNearestCellWalk(location, seedCellI);
    }
}


Foam::label Foam::meshSearch::findNearestFace
(
    const point& location,
    const label seedFaceI,
    const bool useTreeSearch
) const
{
    if (seedFaceI == -1)
    {
        if (useTreeSearch)
        {
            return findNearestFaceTree(location);
        }
        else
        {
            return findNearestFaceLinear(location);
        }
    }
    else
    {
        return findNearestFaceWalk(location, seedFaceI);
    }
}


Foam::label Foam::meshSearch::findCell
(
    const point& location,
    const label seedCellI,
    const bool useTreeSearch
) const
{
    // Find the nearest cell centre to this location
    if (seedCellI == -1)
    {
        if (useTreeSearch)
        {
            return cellTree().findInside(location);
        }
        else
        {
            return findCellLinear(location);
        }
    }
    else
    {
        return findCellWalk(location, seedCellI);
    }
}


Foam::label Foam::meshSearch::findNearestBoundaryFace
(
    const point& location,
    const label seedFaceI,
    const bool useTreeSearch
) const
{
    if (seedFaceI == -1)
    {
        if (useTreeSearch)
        {
            const indexedOctree<treeDataFace>& tree = boundaryTree();

            pointIndexHit info = boundaryTree().findNearest
            (
                location,
                magSqr(tree.bb().max()-tree.bb().min())
            );

            if (!info.hit())
            {
                info = boundaryTree().findNearest
                (
                    location,
                    Foam::sqr(GREAT)
                );
            }

            return tree.shapes().faceLabels()[info.index()];
        }
        else
        {
            scalar minDist = GREAT;

            label minFaceI = -1;

            for
            (
                label faceI = mesh_.nInternalFaces();
                faceI < mesh_.nFaces();
                faceI++
            )
            {
                const face& f =  mesh_.faces()[faceI];

                pointHit curHit =
                    f.nearestPoint
                    (
                        location,
                        mesh_.points()
                    );

                if (curHit.distance() < minDist)
                {
                    minDist = curHit.distance();
                    minFaceI = faceI;
                }
            }
            return minFaceI;
        }
    }
    else
    {
        return findNearestBoundaryFaceWalk(location, seedFaceI);
    }
}


Foam::pointIndexHit Foam::meshSearch::intersection
(
    const point& pStart,
    const point& pEnd
) const
{
    pointIndexHit curHit = boundaryTree().findLine(pStart, pEnd);

    if (curHit.hit())
    {
        // Change index into octreeData into face label
        curHit.setIndex(boundaryTree().shapes().faceLabels()[curHit.index()]);
    }
    return curHit;
}


Foam::List<Foam::pointIndexHit> Foam::meshSearch::intersections
(
    const point& pStart,
    const point& pEnd
) const
{
    DynamicList<pointIndexHit> hits;

    vector edgeVec = pEnd - pStart;
    edgeVec /= mag(edgeVec);

    point pt = pStart;

    pointIndexHit bHit;
    do
    {
        bHit = intersection(pt, pEnd);

        if (bHit.hit())
        {
            hits.append(bHit);

            const vector& area = mesh_.faceAreas()[bHit.index()];

            scalar typDim = Foam::sqrt(mag(area));

            if ((mag(bHit.hitPoint() - pEnd)/typDim) < SMALL)
            {
                break;
            }

            // Restart from hitPoint shifted a little bit in the direction
            // of the destination

            pt =
                bHit.hitPoint()
              + offset(bHit.hitPoint(), bHit.index(), edgeVec);
        }

    } while (bHit.hit());


    hits.shrink();

    return hits;
}


bool Foam::meshSearch::isInside(const point& p) const
{
    return (boundaryTree().getVolumeType(p) == volumeType::INSIDE);
}


// Delete all storage
void Foam::meshSearch::clearOut()
{
    boundaryTreePtr_.clear();
    cellTreePtr_.clear();
    overallBbPtr_.clear();
}


void Foam::meshSearch::correct()
{
    clearOut();
}


// ************************************************************************* //
