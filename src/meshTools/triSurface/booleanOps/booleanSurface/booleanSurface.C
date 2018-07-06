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

#include "booleanSurface.H"
#include "intersectedSurface.H"
#include "orientedSurface.H"
#include "triSurfaceSearch.H"
#include "OFstream.H"
#include "treeBoundBox.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(booleanSurface, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Check whether at least one of faces connected to the intersection has been
// marked.
void Foam::booleanSurface::checkIncluded
(
    const intersectedSurface& surf,
    const labelList& faceZone,
    const label includedFace
)
{
    forAll(surf.intersectionEdges(), intEdgeI)
    {
        label edgeI = surf.intersectionEdges()[intEdgeI];

        const labelList& myFaces = surf.edgeFaces()[edgeI];

        bool usesIncluded = false;

        forAll(myFaces, myFacei)
        {
            if (faceZone[myFaces[myFacei]] == faceZone[includedFace])
            {
                usesIncluded = true;

                break;
            }
        }

        if (!usesIncluded)
        {
            FatalErrorInFunction
                << "None of the faces reachable from face " << includedFace
                << " connects to the intersection."
                << exit(FatalError);
        }
    }
}


// Linear lookup
Foam::label Foam::booleanSurface::index
(
    const labelList& elems,
    const label elem
)
{
    forAll(elems, elemI)
    {
        if (elems[elemI] == elem)
        {
            return elemI;
        }
    }
    return -1;
}


Foam::label Foam::booleanSurface::findEdge
(
    const edgeList& edges,
    const labelList& edgeLabels,
    const edge& e
)
{
    forAll(edgeLabels, edgeLabelI)
    {
        if (edges[edgeLabels[edgeLabelI]] == e)
        {
            return edgeLabels[edgeLabelI];
        }
    }
    FatalErrorInFunction
        << "Cannot find edge " << e << " in edges " << edgeLabels
        << abort(FatalError);

    return -1;
}


// Generate combined patchList (returned). Sets patchMap to map from surf2
// region numbers into combined patchList
Foam::geometricSurfacePatchList Foam::booleanSurface::mergePatches
(
    const triSurface& surf1,
    const triSurface& surf2,
    labelList& patchMap2
)
{
    // Size too big.
    geometricSurfacePatchList combinedPatches
    (
        surf1.patches().size()
      + surf2.patches().size()
    );

    // Copy all patches of surf1
    label combinedPatchi = 0;
    forAll(surf1.patches(), patchi)
    {
        combinedPatches[combinedPatchi++] = surf1.patches()[patchi];
    }

    // (inefficiently) add unique patches from surf2
    patchMap2.setSize(surf2.patches().size());

    forAll(surf2.patches(), patch2I)
    {
        label index = -1;

        forAll(surf1.patches(), patch1I)
        {
            if (surf1.patches()[patch1I] == surf2.patches()[patch2I])
            {
                index = patch1I;

                break;
            }
        }

        if (index == -1)
        {
            combinedPatches[combinedPatchi] = surf2.patches()[patch2I];
            patchMap2[patch2I] = combinedPatchi;
            combinedPatchi++;
        }
        else
        {
            patchMap2[patch2I] = index;
        }
    }

    combinedPatches.setSize(combinedPatchi);

    return combinedPatches;
}


void Foam::booleanSurface::propagateEdgeSide
(
    const triSurface& surf,
    const label prevVert0,
    const label prevFacei,
    const label prevState,
    const label edgeI,
    labelList& side
)
{
    const labelList& eFaces = surf.sortedEdgeFaces()[edgeI];

    // Simple case. Propagate side.
    if (eFaces.size() == 2)
    {
        forAll(eFaces, eFacei)
        {
            propagateSide
            (
                surf,
                prevState,
                eFaces[eFacei],
                side
            );
        }
    }


    if (((eFaces.size() % 2) == 1) && (eFaces.size() != 1))
    {
        FatalErrorInFunction
            << "Don't know how to handle edges with odd number of faces"
            << endl
            << "edge:" << edgeI << " vertices:" << surf.edges()[edgeI]
            << " coming from face:" << prevFacei
            << " edgeFaces:" << eFaces << abort(FatalError);
    }


    // Get position of face in edgeFaces
    label ind = index(eFaces, prevFacei);

    // Determine orientation of faces around edge prevVert0
    // (might be opposite of edge)
    const edge& e = surf.edges()[edgeI];

    // Get next face to include
    label nextInd;
    label prevInd;

    if (e.start() == prevVert0)
    {
        // Edge (and hence eFaces) in same order as prevVert0.
        // Take next face from sorted list
        nextInd = eFaces.fcIndex(ind);
        prevInd = eFaces.rcIndex(ind);
    }
    else
    {
        // Take previous face from sorted neighbours
        nextInd = eFaces.rcIndex(ind);
        prevInd = eFaces.fcIndex(ind);
    }


    if (prevState == OUTSIDE)
    {
        // Coming from outside. nextInd is outside, rest is inside.

        forAll(eFaces, eFacei)
        {
            if (eFacei != ind)
            {
                label nextState;

                if (eFacei == nextInd)
                {
                    nextState = OUTSIDE;
                }
                else
                {
                    nextState = INSIDE;
                }

                propagateSide
                (
                    surf,
                    nextState,
                    eFaces[eFacei],
                    side
                );
            }
        }
    }
    else
    {
        // Coming from inside. prevInd is inside as well, rest is outside.

        forAll(eFaces, eFacei)
        {
            if (eFacei != ind)
            {
                label nextState;

                if (eFacei == prevInd)
                {
                    nextState = INSIDE;
                }
                else
                {
                    nextState = OUTSIDE;
                }

                propagateSide
                (
                    surf,
                    nextState,
                    eFaces[eFacei],
                    side
                );
            }
        }
    }
}


// Face-edge walk. Determines inside/outside for all faces connected to an edge.
void Foam::booleanSurface::propagateSide
(
    const triSurface& surf,
    const label prevState,
    const label facei,
    labelList& side
)
{
    if (side[facei] == UNVISITED)
    {
        side[facei] = prevState;

        const labelledTri& tri = surf.localFaces()[facei];

        // Get copy of face labels
        label a = tri[0];
        label b = tri[1];
        label c = tri[2];

        // Go and visit my edges' face-neighbours.

        const labelList& myEdges = surf.faceEdges()[facei];

        label edgeAB = findEdge(surf.edges(), myEdges, edge(a, b));

        propagateEdgeSide
        (
            surf,
            a,
            facei,
            prevState,
            edgeAB,
            side
        );

        label edgeBC = findEdge(surf.edges(), myEdges, edge(b, c));

        propagateEdgeSide
        (
            surf,
            b,
            facei,
            prevState,
            edgeBC,
            side
        );

        label edgeCA = findEdge(surf.edges(), myEdges, edge(c, a));

        propagateEdgeSide
        (
            surf,
            c,
            facei,
            prevState,
            edgeCA,
            side
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
Foam::booleanSurface::booleanSurface()
:
    triSurface()
{}


// Construct from surfaces and face to include for every surface
Foam::booleanSurface::booleanSurface
(
    const triSurface& surf1,
    const triSurface& surf2,
    const surfaceIntersection& inter,
    const label includeFace1,
    const label includeFace2
)
:
    triSurface(),
    faceMap_()
{
    if (debug)
    {
        Pout<< "booleanSurface : Generating intersected surface for surf1"
            << endl;
    }

    // Add intersection to surface1 (retriangulates cut faces)
    intersectedSurface cutSurf1(surf1, true, inter);


    if (debug)
    {
        Pout<< "booleanSurface : Generated cutSurf1: " << endl;
        cutSurf1.writeStats(Pout);

        Pout<< "Writing to file cutSurf1.obj" << endl;
        cutSurf1.write("cutSurf1.obj");
    }

    if (debug)
    {
        Pout<< "booleanSurface : Generating intersected surface for surf2"
            << endl;
    }

    // Add intersection to surface2
    intersectedSurface cutSurf2(surf2, false, inter);

    if (debug)
    {
        Pout<< "booleanSurface : Generated cutSurf2: " << endl;
        cutSurf2.writeStats(Pout);

        Pout<< "Writing to file cutSurf2.obj" << endl;
        cutSurf2.write("cutSurf2.obj");
    }


    // Find (first) face of cutSurf1 that originates from includeFace1
    label cutSurf1Facei = index(cutSurf1.faceMap(), includeFace1);

    if (debug)
    {
        Pout<< "cutSurf1 : starting to fill from face:" << cutSurf1Facei
            << endl;
    }

    if (cutSurf1Facei == -1)
    {
        FatalErrorInFunction
            << "Did not find face with label " << includeFace1
            << " in intersectedSurface."
            << exit(FatalError);
    }

    // Find (first) face of cutSurf2 that originates from includeFace1
    label cutSurf2Facei = index(cutSurf2.faceMap(), includeFace2);

    if (debug)
    {
        Pout<< "cutSurf2 : starting to fill from face:" << cutSurf2Facei
            << endl;
    }
    if (cutSurf2Facei == -1)
    {
        FatalErrorInFunction
            << "Did not find face with label " << includeFace2
            << " in intersectedSurface."
            << exit(FatalError);
    }


    //
    // Mark faces of cutSurf1 that need to be kept by walking from includeFace1
    // without crossing any edges of the intersection.
    //

    // Mark edges on intersection
    const labelList& int1Edges = cutSurf1.intersectionEdges();

    boolList isIntersectionEdge1(cutSurf1.nEdges(), false);
    forAll(int1Edges, intEdgeI)
    {
        label edgeI = int1Edges[intEdgeI];
        isIntersectionEdge1[edgeI] = true;
    }

    labelList faceZone1;
    cutSurf1.markZones(isIntersectionEdge1, faceZone1);


    // Check whether at least one of sides of intersection has been marked.
    checkIncluded(cutSurf1, faceZone1, cutSurf1Facei);

    // Subset zone which includes cutSurf2Facei
    boolList includedFaces1(cutSurf1.size(), false);

    forAll(faceZone1, facei)
    {
        if (faceZone1[facei] == faceZone1[cutSurf1Facei])
        {
            includedFaces1[facei] = true;
        }
    }

    // Subset to include only interesting part
    labelList pointMap1;
    labelList faceMap1;

    triSurface subSurf1
    (
        cutSurf1.subsetMesh
        (
            includedFaces1,
            pointMap1,
            faceMap1
        )
    );


    //
    // Mark faces of cutSurf2 that need to be kept by walking from includeFace2
    // without crossing any edges of the intersection.
    //

    // Mark edges and points on intersection
    const labelList& int2Edges = cutSurf2.intersectionEdges();

    boolList isIntersectionEdge2(cutSurf2.nEdges(), false);
    forAll(int2Edges, intEdgeI)
    {
        label edgeI = int2Edges[intEdgeI];
        isIntersectionEdge2[edgeI] = true;
    }

    labelList faceZone2;
    cutSurf2.markZones(isIntersectionEdge2, faceZone2);


    // Check whether at least one of sides of intersection has been marked.
    checkIncluded(cutSurf2, faceZone2, cutSurf2Facei);

    // Subset zone which includes cutSurf2Facei
    boolList includedFaces2(cutSurf2.size(), false);

    forAll(faceZone2, facei)
    {
        if (faceZone2[facei] == faceZone2[cutSurf2Facei])
        {
            includedFaces2[facei] = true;
        }
    }

    labelList pointMap2;
    labelList faceMap2;

    triSurface subSurf2
    (
        cutSurf2.subsetMesh
        (
            includedFaces2,
            pointMap2,
            faceMap2
        )
    );


    //
    // Now match up the corresponding points on the intersection. The
    // intersectedSurfaces will have the points resulting from the
    // intersection last in their points and in the same
    // order so we can use the pointMaps from the subsets to find them.
    //
    // We keep the vertices on the first surface and renumber those on the
    // second one.


    //
    // points
    //
    pointField combinedPoints
    (
        subSurf1.nPoints()
      + subSurf2.nPoints()
      - (cutSurf2.nPoints() - cutSurf2.nSurfacePoints())
    );

    // Copy points from subSurf1 and remember the labels of the ones in
    // the intersection
    labelList intersectionLabels
    (
        cutSurf1.nPoints() - cutSurf1.nSurfacePoints()
    );

    label combinedPointi = 0;

    forAll(subSurf1.points(), pointi)
    {
        // Label in cutSurf
        label cutSurfPointi = pointMap1[pointi];

        if (!cutSurf1.isSurfacePoint(cutSurfPointi))
        {
            // Label in original intersection is equal to the cutSurfPointi

            // Remember label in combinedPoints for intersection point.
            intersectionLabels[cutSurfPointi] = combinedPointi;
        }

        // Copy point
        combinedPoints[combinedPointi++] = subSurf1.points()[pointi];
    }

    // Append points from subSurf2 (if they are not intersection points)
    // and construct mapping
    labelList pointMap(subSurf2.nPoints());

    forAll(subSurf2.points(), pointi)
    {
        // Label in cutSurf
        label cutSurfPointi = pointMap2[pointi];

        if (!cutSurf2.isSurfacePoint(cutSurfPointi))
        {
            // Lookup its label in combined point list.
            pointMap[pointi] = intersectionLabels[cutSurfPointi];
        }
        else
        {
            pointMap[pointi] = combinedPointi;

            combinedPoints[combinedPointi++] = subSurf2.points()[pointi];
        }
    }


    //
    // patches
    //

    labelList patchMap2;

    geometricSurfacePatchList combinedPatches
    (
        mergePatches
        (
            surf1,
            surf2,
            patchMap2
        )
    );


    //
    // faces
    //

    List<labelledTri> combinedFaces(subSurf1.size() + subSurf2.size());

    faceMap_.setSize(combinedFaces.size());

    // Copy faces from subSurf1. No need for renumbering.
    label combinedFacei = 0;
    forAll(subSurf1, facei)
    {
        faceMap_[combinedFacei] = faceMap1[facei];
        combinedFaces[combinedFacei++] = subSurf1[facei];
    }

    // Copy and renumber faces from subSurf2.
    forAll(subSurf2, facei)
    {
        const labelledTri& f = subSurf2[facei];

        faceMap_[combinedFacei] = -faceMap2[facei]-1;

        combinedFaces[combinedFacei++] =
            labelledTri
            (
                pointMap[f[0]],
                pointMap[f[1]],
                pointMap[f[2]],
                patchMap2[f.region()]
            );
    }

    triSurface::operator=
    (
        triSurface
        (
            combinedFaces,
            combinedPatches,
            combinedPoints
        )
    );
}


// Construct from surfaces and boolean operation
Foam::booleanSurface::booleanSurface
(
    const triSurface& surf1,
    const triSurface& surf2,
    const surfaceIntersection& inter,
    const label booleanOp
)
:
    triSurface(),
    faceMap_()
{
    if (debug)
    {
        Pout<< "booleanSurface : Testing surf1 and surf2" << endl;

        {
            const labelListList& edgeFaces = surf1.edgeFaces();

            forAll(edgeFaces, edgeI)
            {
                const labelList& eFaces = edgeFaces[edgeI];

                if (eFaces.size() == 1)
                {
                    WarningInFunction
                        << "surf1 is open surface at edge " << edgeI
                        << " verts:" << surf1.edges()[edgeI]
                        << " connected to faces " << eFaces << endl;
                }
            }
        }
        {
            const labelListList& edgeFaces = surf2.edgeFaces();

            forAll(edgeFaces, edgeI)
            {
                const labelList& eFaces = edgeFaces[edgeI];

                if (eFaces.size() == 1)
                {
                    WarningInFunction
                        << "surf2 is open surface at edge " << edgeI
                        << " verts:" << surf2.edges()[edgeI]
                        << " connected to faces " << eFaces << endl;
                }
            }
        }
    }


    //
    // Surface 1
    //

    if (debug)
    {
        Pout<< "booleanSurface : Generating intersected surface for surf1"
            << endl;
    }

    // Add intersection to surface1 (retriangulates cut faces)
    intersectedSurface cutSurf1(surf1, true, inter);

    if (debug)
    {
        Pout<< "booleanSurface : Generated cutSurf1: " << endl;
        cutSurf1.writeStats(Pout);

        Pout<< "Writing to file cutSurf1.obj" << endl;
        cutSurf1.write("cutSurf1.obj");
    }


    //
    // Surface 2
    //

    if (debug)
    {
        Pout<< "booleanSurface : Generating intersected surface for surf2"
            << endl;
    }

    // Add intersection to surface2
    intersectedSurface cutSurf2(surf2, false, inter);

    if (debug)
    {
        Pout<< "booleanSurface : Generated cutSurf2: " << endl;
        cutSurf2.writeStats(Pout);

        Pout<< "Writing to file cutSurf2.obj" << endl;
        cutSurf2.write("cutSurf2.obj");
    }


    //
    // patches
    //

    labelList patchMap2;

    geometricSurfacePatchList combinedPatches
    (
        mergePatches
        (
            surf1,
            surf2,
            patchMap2
        )
    );


    //
    // Now match up the corresponding points on the intersection. The
    // intersectedSurfaces will have the points resulting from the
    // intersection first in their points and in the same
    // order
    //
    // We keep the vertices on the first surface and renumber those on the
    // second one.

    pointField combinedPoints(cutSurf1.nPoints() + cutSurf2.nSurfacePoints());

    // Copy all points from 1 and non-intersection ones from 2.

    label combinedPointi = 0;

    forAll(cutSurf1.points(), pointi)
    {
        combinedPoints[combinedPointi++] = cutSurf1.points()[pointi];
    }

    for
    (
        label pointi = 0;
        pointi < cutSurf2.nSurfacePoints();
        pointi++
    )
    {
        combinedPoints[combinedPointi++] = cutSurf2.points()[pointi];
    }

    // Point order is now
    // - 0.. cutSurf1.nSurfacePoints : original surf1 points
    // -  .. cutSurf1.nPoints        : intersection points
    // -  .. cutSurf2.nSurfacePoints : original surf2 points

    if (debug)
    {
        Pout<< "booleanSurface : generated points:" << nl
            << "    0 .. " << cutSurf1.nSurfacePoints()-1
            << " : original surface1"
            << nl
            << "    " << cutSurf1.nSurfacePoints()
            << " .. " << cutSurf1.nPoints()-1
            << " : intersection points"
            << nl
            << "    " << cutSurf1.nPoints() << " .. "
            << cutSurf2.nSurfacePoints()-1
            << " : surface2 points"
            << nl
            << endl;
    }

    // Copy faces. Faces from surface 1 keep vertex numbering and region info.
    // Faces from 2 get vertices and region renumbered.
    List<labelledTri> combinedFaces(cutSurf1.size() + cutSurf2.size());

    label combinedFacei = 0;

    forAll(cutSurf1, facei)
    {
        combinedFaces[combinedFacei++] = cutSurf1[facei];
    }

    forAll(cutSurf2, facei)
    {
        labelledTri& combinedTri = combinedFaces[combinedFacei++];

        const labelledTri& tri = cutSurf2[facei];

        forAll(tri, fp)
        {
            if (cutSurf2.isSurfacePoint(tri[fp]))
            {
                // Renumber. Surface2 points are after ones from surf 1.
                combinedTri[fp] = tri[fp] + cutSurf1.nPoints();
            }
            else
            {
                // Is intersection.
                combinedTri[fp] =
                    tri[fp]
                  - cutSurf2.nSurfacePoints()
                  + cutSurf1.nSurfacePoints();
            }
        }
        combinedTri.region() = patchMap2[tri.region()];
    }


    // Now we have surface in combinedFaces and combinedPoints. Use
    // booleanOp to determine which part of what to keep.

    // Construct addressing for whole part.
    triSurface combinedSurf
    (
        combinedFaces,
        combinedPatches,
        combinedPoints
    );

    if (debug)
    {
        Pout<< "booleanSurface : Generated combinedSurf: " << endl;
        combinedSurf.writeStats(Pout);

        Pout<< "Writing to file combinedSurf.obj" << endl;
        combinedSurf.write("combinedSurf.obj");
    }


    if (booleanOp == booleanSurface::ALL)
    {
        // Special case: leave surface multiply connected

        faceMap_.setSize(combinedSurf.size());

        label combinedFacei = 0;

        forAll(cutSurf1, facei)
        {
            faceMap_[combinedFacei++] = cutSurf1.faceMap()[facei];
        }
        forAll(cutSurf2, facei)
        {
            faceMap_[combinedFacei++] = -cutSurf2.faceMap()[facei] - 1;
        }

        triSurface::operator=(combinedSurf);

        return;
    }


    // Get outside point.
    point outsidePoint = 2 * treeBoundBox(combinedSurf.localPoints()).span();

    //
    // Linear search for nearest point on surface.
    //

    const pointField& pts = combinedSurf.points();

    label minFacei = -1;
    pointHit minHit(false, Zero, great, true);

    forAll(combinedSurf, facei)
    {
        pointHit curHit = combinedSurf[facei].nearestPoint(outsidePoint, pts);

        if (curHit.distance() < minHit.distance())
        {
            minHit = curHit;
            minFacei = facei;
        }
    }

    if (debug)
    {
        Pout<< "booleanSurface : found for point:" << outsidePoint
            << "  nearest face:" << minFacei
            << "  nearest point:" << minHit.rawPoint()
            << endl;
    }

    // Visibility/side of face:
    //      UNVISITED: unvisited
    //      OUTSIDE: visible from outside
    //      INSIDE: invisible from outside
    labelList side(combinedSurf.size(), UNVISITED);

    // Walk face-edge-face and propagate inside/outside status.
    propagateSide(combinedSurf, OUTSIDE, minFacei, side);


    // Depending on operation include certain faces.
    //  INTERSECTION: faces on inside of 1 and of 2
    //  UNION: faces on outside of 1 and of 2
    //  DIFFERENCE: faces on outside of 1 and inside of 2

    boolList include(combinedSurf.size(), false);

    forAll(side, facei)
    {
        if (side[facei] == UNVISITED)
        {
            FatalErrorInFunction
                << "Face " << facei << " has not been reached by walking from"
                << " nearest point " << minHit.rawPoint()
                << " nearest face " << minFacei << exit(FatalError);
        }
        else if (side[facei] == OUTSIDE)
        {
            if (booleanOp == booleanSurface::UNION)
            {
                include[facei] = true;
            }
            else if (booleanOp == booleanSurface::INTERSECTION)
            {
                include[facei] = false;
            }
            else    // difference
            {
                include[facei] = (facei < cutSurf1.size()); // face from surf1
            }
        }
        else    // inside
        {
            if (booleanOp == booleanSurface::UNION)
            {
                include[facei] = false;
            }
            else if (booleanOp == booleanSurface::INTERSECTION)
            {
                include[facei] = true;
            }
            else    // difference
            {
                include[facei] = (facei >= cutSurf1.size()); // face from surf2
            }
        }
    }

    // Create subsetted surface
    labelList subToCombinedPoint;
    labelList subToCombinedFace;
    triSurface subSurf
    (
        combinedSurf.subsetMesh
        (
            include,
            subToCombinedPoint,
            subToCombinedFace
        )
    );

    // Create face map
    faceMap_.setSize(subSurf.size());

    forAll(subToCombinedFace, facei)
    {
        // Get label in combinedSurf
        label combinedFacei = subToCombinedFace[facei];

        // First faces in combinedSurf come from cutSurf1.

        if (combinedFacei < cutSurf1.size())
        {
            label cutSurf1Face = combinedFacei;

            faceMap_[facei] = cutSurf1.faceMap()[cutSurf1Face];
        }
        else
        {
            label cutSurf2Face = combinedFacei - cutSurf1.size();

            faceMap_[facei] = - cutSurf2.faceMap()[cutSurf2Face] - 1;
        }
    }

    // Orient outwards
    orientedSurface outSurf(subSurf);

    // Assign.
    triSurface::operator=(outSurf);
}


// ************************************************************************* //
