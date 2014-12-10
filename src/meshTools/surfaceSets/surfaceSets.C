/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "surfaceSets.H"
#include "polyMesh.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "pointSet.H"
#include "cellSet.H"
#include "surfaceToCell.H"
#include "cellToPoint.H"
#include "cellToCell.H"
#include "pointToCell.H"
#include "meshSearch.H"
#include "cellClassification.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//Foam::scalar Foam::surfaceSets::minEdgeLen
//(
//    const primitiveMesh& mesh,
//    const label pointI
//)
//{
//    const edgeList& edges = mesh.edges();
//
//    const pointField& points = mesh.points();
//
//    const labelList& pEdges = mesh.pointEdges()[pointI];
//
//    scalar minLen = GREAT;
//
//    forAll(pEdges, i)
//    {
//        minLen = min(minLen, edges[pEdges[i]].mag(points));
//    }
//    return minLen;
//}
//
//
//// Returns true if cell uses at least one selected point
//bool Foam::surfaceSets::usesPoint
//(
//    const primitiveMesh& mesh,
//    const boolList& selectedPoint,
//    const label cellI
//)
//{
//    const labelList& cFaces = mesh.cells()[cellI];
//
//    forAll(cFaces, cFaceI)
//    {
//        label faceI = cFaces[cFaceI];
//
//        const face& f = mesh.faces()[faceI];
//
//        forAll(f, fp)
//        {
//            if (selectedPoint[f[fp]])
//            {
//                return true;
//            }
//        }
//    }
//    return false;
//}



//// Remove cells in allCells which are connected to other cells in allCells
//// by outside vertices only. Since these outside vertices will be moved onto
//// a surface they might result in flat cells.
//Foam::label Foam::surfaceSets::removeHangingCells
//(
//    const primitiveMesh& mesh,
//    const triSurfaceSearch& querySurf,
//    labelHashSet& internalCells
//)
//{
//    const pointField& points = mesh.points();
//    const cellList& cells = mesh.cells();
//    const faceList& faces = mesh.faces();
//
//    // Determine cells that have all points on the boundary.
//    labelHashSet flatCandidates(getHangingCells(mesh, internalCells));
//
//    // All boundary points will become visible after subsetting and will be
//    // snapped
//    // to surface. Calculate new volume for cells using only these points and
//    // check if it does not become too small.
//
//    // Get points used by flatCandidates.
//    labelHashSet outsidePoints(flatCandidates.size());
//
//    // Snap outside points to surface
//    pointField newPoints(points);
//
//    forAllConstIter(labelHashSet, flatCandidates, iter)
//    {
//        const cell& cFaces = cells[iter.key()];
//
//        forAll(cFaces, cFaceI)
//        {
//            const face& f = faces[cFaces[cFaceI]];
//
//            forAll(f, fp)
//            {
//                label pointI = f[fp];
//
//                if (outsidePoints.insert(pointI))
//                {
//                    // Calculate new position for this outside point
//                    tmp<pointField> tnearest =
//                        querySurf.calcNearest(pointField(1, points[pointI]));
//                    newPoints[pointI] = tnearest()[0];
//                }
//            }
//        }
//    }
//
//
//    // Calculate new volume for mixed cells
//    label nRemoved = 0;
//    forAllConstIter(labelHashSet, flatCandidates, iter)
//    {
//        label cellI = iter.key();
//
//        const cell& cll = cells[cellI];
//
//        scalar newVol = cll.mag(newPoints, faces);
//        scalar oldVol = mesh.cellVolumes()[cellI];
//
//        if (newVol < 0.1 * oldVol)
//        {
//            internalCells.erase(cellI);
//            nRemoved++;
//        }
//    }
//
//    return nRemoved;
//}


//// Select all points out of pointSet where the distance to the surface
//// is less than a factor times a local length scale (minimum length of
//// connected edges)
//void Foam::surfaceSets::getNearPoints
//(
//    const primitiveMesh& mesh,
//    const triSurface&,
//    const triSurfaceSearch& querySurf,
//    const scalar edgeFactor,
//    const pointSet& candidateSet,
//    pointSet& nearPointSet
//)
//{
//    if (edgeFactor <= 0)
//    {
//        return;
//    }
//
//    labelList candidates(candidateSet.toc());
//
//    pointField candidatePoints(candidates.size());
//    forAll(candidates, i)
//    {
//        candidatePoints[i] = mesh.points()[candidates[i]];
//    }
//
//    tmp<pointField> tnearest = querySurf.calcNearest(candidatePoints);
//    const pointField& nearest = tnearest();
//
//    const pointField& points = mesh.points();
//
//    forAll(candidates, i)
//    {
//        label pointI = candidates[i];
//
//        scalar minLen = minEdgeLen(mesh, pointI);
//
//        scalar dist = mag(nearest[i] - points[pointI]);
//
//        if (dist < edgeFactor * minLen)
//        {
//            nearPointSet.insert(pointI);
//        }
//    }
//}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfaceSets::getSurfaceSets
(
    const polyMesh& mesh,
    const fileName&,
    const triSurface&,
    const triSurfaceSearch& querySurf,
    const pointField& outsidePts,

    const label nCutLayers,

    labelHashSet& inside,
    labelHashSet& outside,
    labelHashSet& cut
)
{
    // Construct search engine on mesh
    meshSearch queryMesh(mesh);

    // Cut faces with surface and classify cells
    cellClassification cellType
    (
        mesh,
        queryMesh,
        querySurf,
        outsidePts
    );

    if (nCutLayers > 0)
    {
        // Trim cutCells so they are max nCutLayers away (as seen in point-cell
        // walk) from outside cells.
        cellType.trimCutCells
        (
            nCutLayers,
            cellClassification::OUTSIDE,
            cellClassification::INSIDE
        );
    }

    forAll(cellType, cellI)
    {
        label cType = cellType[cellI];

        if (cType == cellClassification::CUT)
        {
            cut.insert(cellI);
        }
        else if (cType == cellClassification::INSIDE)
        {
            inside.insert(cellI);
        }
        else if (cType == cellClassification::OUTSIDE)
        {
            outside.insert(cellI);
        }
    }
}


Foam::labelHashSet Foam::surfaceSets::getHangingCells
(
    const primitiveMesh& mesh,
    const labelHashSet& internalCells
)
{
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();


    // Divide points into
    // -referenced by internal only
    // -referenced by outside only
    // -mixed

    List<pointStatus> pointSide(mesh.nPoints(), NOTSET);

    for (label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        if (internalCells.found(cellI))
        {
            // Inside cell. Mark all vertices seen from this cell.
            const labelList& cFaces = cells[cellI];

            forAll(cFaces, cFaceI)
            {
                const face& f = faces[cFaces[cFaceI]];

                forAll(f, fp)
                {
                    label pointI = f[fp];

                    if (pointSide[pointI] == NOTSET)
                    {
                        pointSide[pointI] = INSIDE;
                    }
                    else if (pointSide[pointI] == OUTSIDE)
                    {
                        pointSide[pointI] = MIXED;
                    }
                    else
                    {
                        // mixed or inside stay same
                    }
                }
            }
        }
        else
        {
            // Outside cell
            const labelList& cFaces = cells[cellI];

            forAll(cFaces, cFaceI)
            {
                const face& f = faces[cFaces[cFaceI]];

                forAll(f, fp)
                {
                    label pointI = f[fp];

                    if (pointSide[pointI] == NOTSET)
                    {
                        pointSide[pointI] = OUTSIDE;
                    }
                    else if (pointSide[pointI] == INSIDE)
                    {
                        pointSide[pointI] = MIXED;
                    }
                    else
                    {
                        // mixed or outside stay same
                    }
                }
            }
        }
    }


    //OFstream mixedStr("mixed.obj");
    //
    //forAll(pointSide, pointI)
    //{
    //    if (pointSide[pointI] == MIXED)
    //    {
    //        const point& pt = points[pointI];
    //
    //        mixedStr << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z()
    //            << endl;
    //    }
    //}


    // Determine cells using mixed points only

    labelHashSet mixedOnlyCells(internalCells.size());

    forAllConstIter(labelHashSet, internalCells, iter)
    {
        const label cellI = iter.key();
        const cell& cFaces = cells[cellI];

        label usesMixedOnly = true;

        forAll(cFaces, i)
        {
            const face& f = faces[cFaces[i]];

            forAll(f, fp)
            {
                if (pointSide[f[fp]] != MIXED)
                {
                    usesMixedOnly = false;
                    break;
                }
            }

            if (!usesMixedOnly)
            {
                break;
            }
        }
        if (usesMixedOnly)
        {
            mixedOnlyCells.insert(cellI);
        }
    }

    return mixedOnlyCells;
}


//void Foam::surfaceSets::writeSurfaceSets
//(
//    const polyMesh& mesh,
//    const fileName& surfName,
//    const triSurface& surf,
//    const triSurfaceSearch& querySurf,
//    const pointField& outsidePts,
//    const scalar edgeFactor
//)
//{
//    // Cellsets for inside/outside determination
//    cellSet rawInside(mesh, "rawInside", mesh.nCells()/10);
//    cellSet rawOutside(mesh, "rawOutside", mesh.nCells()/10);
//    cellSet cutCells(mesh, "cutCells", mesh.nCells()/10);
//
//    // Get inside/outside/cut cells
//    getSurfaceSets
//    (
//        mesh,
//        surfName,
//        surf,
//        querySurf,
//        outsidePts,
//
//        rawInside,
//        rawOutside,
//        cutCells
//    );
//
//
//    Pout<< "rawInside:" << rawInside.size() << endl;
//
//    label nRemoved;
//    do
//    {
//        nRemoved = removeHangingCells(mesh, querySurf, rawInside);
//
//        Pout<< nl
//            << "Removed " << nRemoved
//            << " rawInside cells that have all their points on the outside"
//            << endl;
//    }
//    while (nRemoved != 0);
//
//    Pout<< "Writing inside cells (" << rawInside.size() << ") to cellSet "
//        << rawInside.instance()/rawInside.local()/rawInside.name()
//        << endl << endl;
//    rawInside.write();
//
//
//
//    // Select outside cells
//    surfaceToCell outsideSource
//    (
//        mesh,
//        surfName,
//        surf,
//        querySurf,
//        outsidePts,
//        false,          // includeCut
//        false,          // includeInside
//        true,           // includeOutside
//        -GREAT,         // nearDist
//        -GREAT          // curvature
//    );
//
//    outsideSource.applyToSet(topoSetSource::NEW, rawOutside);
//
//    Pout<< "rawOutside:" << rawOutside.size() << endl;
//
//    do
//    {
//        nRemoved = removeHangingCells(mesh, querySurf, rawOutside);
//
//        Pout<< nl
//            << "Removed " << nRemoved
//            << " rawOutside cells that have all their points on the outside"
//            << endl;
//    }
//    while (nRemoved != 0);
//
//    Pout<< "Writing outside cells (" << rawOutside.size() << ") to cellSet "
//        << rawOutside.instance()/rawOutside.local()/rawOutside.name()
//        << endl << endl;
//    rawOutside.write();
//
//
//    // Select cut cells by negating inside and outside set.
//    cutCells.invert(mesh.nCells());
//
//    cellToCell deleteInsideSource(mesh, rawInside.name());
//
//    deleteInsideSource.applyToSet(topoSetSource::DELETE, cutCells);
//    Pout<< "Writing cut cells (" << cutCells.size() << ") to cellSet "
//        << cutCells.instance()/cutCells.local()/cutCells.name()
//        << endl << endl;
//    cutCells.write();
//
//
//    //
//    // Remove cells with points too close to surface.
//    //
//
//
//    // Get all points in cutCells.
//    pointSet cutPoints(mesh, "cutPoints", 4*cutCells.size());
//    cellToPoint cutSource(mesh, "cutCells", cellToPoint::ALL);
//    cutSource.applyToSet(topoSetSource::NEW, cutPoints);
//
//    // Get all points that are too close to surface.
//    pointSet nearPoints(mesh, "nearPoints", cutPoints.size());
//
//    getNearPoints
//    (
//        mesh,
//        surf,
//        querySurf,
//        edgeFactor,
//        cutPoints,
//        nearPoints
//    );
//
//    Pout<< nl
//        << "Selected " << nearPoints.size()
//        << " points that are closer than " << edgeFactor
//        << " times the local minimum lengthscale to the surface"
//        << nl << endl;
//
//
//    // Remove cells that use any of the points in nearPoints
//    // from the inside and outsideCells.
//    nearPoints.write();
//    pointToCell pToCell(mesh, nearPoints.name(), pointToCell::ANY);
//
//
//
//    // Start off from copy of rawInside, rawOutside
//    cellSet inside(mesh, "inside", rawInside);
//    cellSet outside(mesh, "outside", rawOutside);
//
//    pToCell.applyToSet(topoSetSource::DELETE, inside);
//    pToCell.applyToSet(topoSetSource::DELETE, outside);
//
//    Pout<< nl
//        << "Removed " << rawInside.size() - inside.size()
//        << " inside cells that are too close to the surface" << endl;
//
//    Pout<< nl
//        << "Removed " << rawOutside.size() - outside.size()
//        << " inside cells that are too close to the surface" << nl << endl;
//
//
//
//    //
//    // Remove cells with one or no faces on rest of cellSet. Note: Problem is
//    // not these cells an sich but rather that all of their vertices will be
//    // outside vertices and thus projected onto surface. Which might (if they
//    // project onto same surface) result in flattened cells.
//    //
//
//    do
//    {
//        nRemoved = removeHangingCells(mesh, querySurf, inside);
//
//        Pout<< nl
//            << "Removed " << nRemoved
//            << " inside cells that have all their points on the outside"
//            << endl;
//    }
//    while (nRemoved != 0);
//    do
//    {
//        nRemoved = removeHangingCells(mesh, querySurf, outside);
//
//        Pout<< nl
//            << "Removed " << nRemoved
//            << " outside cells that have all their points on the inside"
//            << endl;
//    }
//    while (nRemoved != 0);
//
//
//    //
//    // Write
//    //
//
//
//    Pout<< "Writing inside cells (" << inside.size() << ") to cellSet "
//        << inside.instance()/inside.local()/inside.name()
//        << endl << endl;
//
//    inside.write();
//
//    Pout<< "Writing outside cells (" << outside.size() << ") to cellSet "
//        << outside.instance()/outside.local()/outside.name()
//        << endl << endl;
//
//    outside.write();
//}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
