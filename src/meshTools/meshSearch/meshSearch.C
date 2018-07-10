/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

    forAll(points, pointi)
    {
        scalar distSqr = magSqr(points[pointi] - sample);

        if (distSqr < nearestDistSqr)
        {
            nearestDistSqr = distSqr;
            nearestI = pointi;
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
        label pointi = indices[i];

        scalar distSqr = magSqr(points[pointi] - sample);

        if (distSqr < nearestDistSqr)
        {
            nearestDistSqr = distSqr;
            nearestI = pointi;
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
        info = tree.findNearest(location, Foam::sqr(great));
    }
    return info.index();
}


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


Foam::label Foam::meshSearch::findNearestCellWalk
(
    const point& location,
    const label seedCelli
) const
{
    if (seedCelli < 0)
    {
        FatalErrorInFunction
            << "illegal seedCell:" << seedCelli << exit(FatalError);
    }

    // Walk in direction of face that decreases distance

    label curCelli = seedCelli;
    scalar distanceSqr = magSqr(mesh_.cellCentres()[curCelli] - location);

    bool closer;

    do
    {
        // Try neighbours of curCelli
        closer = findNearer
        (
            location,
            mesh_.cellCentres(),
            mesh_.cellCells()[curCelli],
            curCelli,
            distanceSqr
        );
    } while (closer);

    return curCelli;
}


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
        // Search with desperate span
        info = tree.findNearest(location, Foam::sqr(great));
    }


    // Now check any of the faces of the nearest cell
    const vectorField& centres = mesh_.faceCentres();
    const cell& ownFaces = mesh_.cells()[info.index()];

    label nearestFacei = ownFaces[0];
    scalar minProximity = magSqr(centres[nearestFacei] - location);

    findNearer
    (
        location,
        centres,
        ownFaces,
        nearestFacei,
        minProximity
    );

    return nearestFacei;
}


Foam::label Foam::meshSearch::findNearestFaceLinear(const point& location) const
{
    const vectorField& centres = mesh_.faceCentres();

    label nearestFacei = 0;
    scalar minProximity = magSqr(centres[nearestFacei] - location);

    findNearer
    (
        location,
        centres,
        nearestFacei,
        minProximity
    );

    return nearestFacei;
}


Foam::label Foam::meshSearch::findNearestFaceWalk
(
    const point& location,
    const label seedFacei
) const
{
    if (seedFacei < 0)
    {
        FatalErrorInFunction
            << "illegal seedFace:" << seedFacei << exit(FatalError);
    }

    const vectorField& centres = mesh_.faceCentres();


    // Walk in direction of face that decreases distance

    label curFacei = seedFacei;
    scalar distanceSqr = magSqr(centres[curFacei] - location);

    while (true)
    {
        label betterFacei = curFacei;

        findNearer
        (
            location,
            centres,
            mesh_.cells()[mesh_.faceOwner()[curFacei]],
            betterFacei,
            distanceSqr
        );

        if (mesh_.isInternalFace(curFacei))
        {
            findNearer
            (
                location,
                centres,
                mesh_.cells()[mesh_.faceNeighbour()[curFacei]],
                betterFacei,
                distanceSqr
            );
        }

        if (betterFacei == curFacei)
        {
            break;
        }

        curFacei = betterFacei;
    }

    return curFacei;
}


Foam::label Foam::meshSearch::findCellLinear(const point& location) const
{
    bool cellFound = false;
    label n = 0;

    label celli = -1;

    while ((!cellFound) && (n < mesh_.nCells()))
    {
        if (mesh_.pointInCell(location, n, cellDecompMode_))
        {
            cellFound = true;
            celli = n;
        }
        else
        {
            n++;
        }
    }
    if (cellFound)
    {
        return celli;
    }
    else
    {
        return -1;
    }
}


Foam::label Foam::meshSearch::findCellWalk
(
    const point& location,
    const label seedCelli
) const
{
    if (seedCelli < 0)
    {
        FatalErrorInFunction
            << "illegal seedCell:" << seedCelli << exit(FatalError);
    }

    if (mesh_.pointInCell(location, seedCelli, cellDecompMode_))
    {
        return seedCelli;
    }

    // Walk in direction of face that decreases distance
    label curCelli = seedCelli;
    scalar nearestDistSqr = magSqr(mesh_.cellCentres()[curCelli] - location);

    while(true)
    {
        // Try neighbours of curCelli

        const cell& cFaces = mesh_.cells()[curCelli];

        label nearestCelli = -1;

        forAll(cFaces, i)
        {
            label facei = cFaces[i];

            if (mesh_.isInternalFace(facei))
            {
                label celli = mesh_.faceOwner()[facei];
                if (celli == curCelli)
                {
                    celli = mesh_.faceNeighbour()[facei];
                }

                // Check if this is the correct cell
                if (mesh_.pointInCell(location, celli, cellDecompMode_))
                {
                    return celli;
                }

                // Also calculate the nearest cell
                scalar distSqr = magSqr(mesh_.cellCentres()[celli] - location);

                if (distSqr < nearestDistSqr)
                {
                    nearestDistSqr = distSqr;
                    nearestCelli = celli;
                }
            }
        }

        if (nearestCelli == -1)
        {
            return -1;
        }

        // Continue with the nearest cell
        curCelli = nearestCelli;
    }

    return -1;
}


Foam::label Foam::meshSearch::findNearestBoundaryFaceWalk
(
    const point& location,
    const label seedFacei
) const
{
    if (seedFacei < 0)
    {
        FatalErrorInFunction
            << "illegal seedFace:" << seedFacei << exit(FatalError);
    }

    // Start off from seedFacei

    label curFacei = seedFacei;

    const face& f =  mesh_.faces()[curFacei];

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

        label lastFacei = curFacei;

        const labelList& myEdges = mesh_.faceEdges()[curFacei];

        forAll(myEdges, myEdgeI)
        {
            const labelList& neighbours = mesh_.edgeFaces()[myEdges[myEdgeI]];

            // Check any face which uses edge, is boundary face and
            // is not curFacei itself.

            forAll(neighbours, nI)
            {
                label facei = neighbours[nI];

                if
                (
                    (facei >= mesh_.nInternalFaces())
                 && (facei != lastFacei)
                )
                {
                    const face& f =  mesh_.faces()[facei];

                    pointHit curHit = f.nearestPoint
                    (
                        location,
                        mesh_.points()
                    );

                    // If the face is closer, reset current face and distance
                    if (curHit.distance() < minDist)
                    {
                        minDist = curHit.distance();
                        curFacei = facei;
                        closer = true;  // a closer neighbour has been found
                    }
                }
            }
        }
    } while (closer);

    return curFacei;
}


Foam::vector Foam::meshSearch::offset
(
    const point& bPoint,
    const label bFacei,
    const vector& dir
) const
{
    // Get the neighbouring cell
    label ownerCelli = mesh_.faceOwner()[bFacei];

    const point& c = mesh_.cellCentres()[ownerCelli];

    // Typical dimension: distance from point on face to cell centre
    scalar typDim = mag(c - bPoint);

    return tol_*typDim*dir;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshSearch::meshSearch
(
    const polyMesh& mesh,
    const polyMesh::cellDecomposition cellDecompMode
)
:
    mesh_(mesh),
    cellDecompMode_(cellDecompMode)
{
    if
    (
        cellDecompMode_ == polyMesh::FACE_DIAG_TRIS
     || cellDecompMode_ == polyMesh::CELL_TETS)
    {
        // Force construction of face diagonals
        (void)mesh.tetBasePtIs();
    }
}


Foam::meshSearch::meshSearch
(
    const polyMesh& mesh,
    const treeBoundBox& bb,
    const polyMesh::cellDecomposition cellDecompMode
)
:
    mesh_(mesh),
    cellDecompMode_(cellDecompMode)
{
    overallBbPtr_.reset(new treeBoundBox(bb));

    if
    (
        cellDecompMode_ == polyMesh::FACE_DIAG_TRIS
     || cellDecompMode_ == polyMesh::CELL_TETS
    )
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

const Foam::indexedOctree<Foam::treeDataFace>&
Foam::meshSearch::boundaryTree() const
{
    if (!boundaryTreePtr_.valid())
    {
        //
        // Construct tree
        //

        if (!overallBbPtr_.valid())
        {
            overallBbPtr_.reset
            (
                new treeBoundBox(mesh_.points())
            );

            treeBoundBox& overallBb = overallBbPtr_();

            // Extend slightly and make 3D
            overallBb = overallBb.extend(1e-4);
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


const Foam::indexedOctree<Foam::treeDataCell>&
Foam::meshSearch::cellTree() const
{
    if (!cellTreePtr_.valid())
    {
        //
        // Construct tree
        //

        if (!overallBbPtr_.valid())
        {
            overallBbPtr_.reset
            (
                new treeBoundBox(mesh_.points())
            );

            treeBoundBox& overallBb = overallBbPtr_();

            // Extend slightly and make 3D
            overallBb = overallBb.extend(1e-4);
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


Foam::label Foam::meshSearch::findNearestCell
(
    const point& location,
    const label seedCelli,
    const bool useTreeSearch
) const
{
    if (seedCelli == -1)
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
        return findNearestCellWalk(location, seedCelli);
    }
}


Foam::label Foam::meshSearch::findNearestFace
(
    const point& location,
    const label seedFacei,
    const bool useTreeSearch
) const
{
    if (seedFacei == -1)
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
        return findNearestFaceWalk(location, seedFacei);
    }
}


Foam::label Foam::meshSearch::findCell
(
    const point& location,
    const label seedCelli,
    const bool useTreeSearch
) const
{
    // Find the nearest cell centre to this location
    if (seedCelli == -1)
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
        return findCellWalk(location, seedCelli);
    }
}


Foam::label Foam::meshSearch::findNearestBoundaryFace
(
    const point& location,
    const label seedFacei,
    const bool useTreeSearch
) const
{
    if (seedFacei == -1)
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
                    Foam::sqr(great)
                );
            }

            return tree.shapes().faceLabels()[info.index()];
        }
        else
        {
            scalar minDist = great;

            label minFacei = -1;

            for
            (
                label facei = mesh_.nInternalFaces();
                facei < mesh_.nFaces();
                facei++
            )
            {
                const face& f =  mesh_.faces()[facei];

                pointHit curHit =
                    f.nearestPoint
                    (
                        location,
                        mesh_.points()
                    );

                if (curHit.distance() < minDist)
                {
                    minDist = curHit.distance();
                    minFacei = facei;
                }
            }
            return minFacei;
        }
    }
    else
    {
        return findNearestBoundaryFaceWalk(location, seedFacei);
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

            if ((mag(bHit.hitPoint() - pEnd)/typDim) < small)
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
