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

Application
    surfaceSplitNonManifolds

Description
    Takes multiply connected surface and tries to split surface at
    multiply connected edges by duplicating points. Introduces concept of
    - borderEdge. Edge with 4 faces connected to it.
    - borderPoint. Point connected to exactly 2 borderEdges.
    - borderLine. Connected list of borderEdges.

    By duplicating borderPoints this will split 'borderLines'. As a
    preprocessing step it can detect borderEdges without any borderPoints
    and explicitly split these triangles.

    The problems in this algorithm are:
    - determining which two (of the four) faces form a surface. Done by walking
      face-edge-face while keeping and edge or point on the borderEdge
      borderPoint.
    - determining the outwards pointing normal to be used to slightly offset the
      duplicated point.

    Uses sortedEdgeFaces quite a bit.

    Is tested on simple borderLines resulting from extracting a surface
    from a hex mesh. Will quite possibly go wrong on more complicated border
    lines (i.e. ones forming a loop).

    Dumps surface every so often since might take a long time to complete.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "triSurface.H"
#include "OFstream.H"
#include "ListOps.H"
#include "triSurfaceTools.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writeOBJ(Ostream& os, const pointField& pts)
{
    forAll(pts, i)
    {
        const point& pt = pts[i];

        os  << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
    }
}


void dumpPoints(const triSurface& surf, const labelList& borderPoint)
{
    fileName fName("borderPoints.obj");

    Info<< "Dumping borderPoints as Lightwave .obj file to " << fName
        << "\nThis can be visualized with e.g. javaview (www.javaview.de)\n\n";

    OFstream os(fName);

    forAll(borderPoint, pointI)
    {
        if (borderPoint[pointI] != -1)
        {
            const point& pt = surf.localPoints()[pointI];

            os  << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
        }
    }
}


void dumpEdges(const triSurface& surf, const boolList& borderEdge)
{
    fileName fName("borderEdges.obj");

    Info<< "Dumping borderEdges as Lightwave .obj file to " << fName
        << "\nThis can be visualized with e.g. javaview (www.javaview.de)\n\n";

    OFstream os(fName);

    writeOBJ(os, surf.localPoints());

    forAll(borderEdge, edgeI)
    {
        if (borderEdge[edgeI])
        {
            const edge& e = surf.edges()[edgeI];

            os  << "l " << e.start()+1 << ' ' << e.end()+1 << endl;
        }
    }
}


void dumpFaces
(
    const fileName& fName,
    const triSurface& surf,
    const Map<label>& connectedFaces
)
{
    Info<< "Dumping connectedFaces as Lightwave .obj file to " << fName
        << "\nThis can be visualized with e.g. javaview (www.javaview.de)\n\n";

    OFstream os(fName);

    forAllConstIter(Map<label>, connectedFaces, iter)
    {
        point ctr = surf.localFaces()[iter.key()].centre(surf.localPoints());

        os  << "v " << ctr.x() << ' ' << ctr.y() << ' ' << ctr.z() << endl;
    }
}


void testSortedEdgeFaces(const triSurface& surf)
{
    const labelListList& edgeFaces = surf.edgeFaces();
    const labelListList& sortedEdgeFaces = surf.sortedEdgeFaces();

    forAll(edgeFaces, edgeI)
    {
        const labelList& myFaces = edgeFaces[edgeI];
        const labelList& sortMyFaces = sortedEdgeFaces[edgeI];

        forAll(myFaces, i)
        {
            if (findIndex(sortMyFaces, myFaces[i]) == -1)
            {
                FatalErrorIn("testSortedEdgeFaces(..)") << abort(FatalError);
            }
        }
        forAll(sortMyFaces, i)
        {
            if (findIndex(myFaces, sortMyFaces[i]) == -1)
            {
                FatalErrorIn("testSortedEdgeFaces(..)") << abort(FatalError);
            }
        }
    }
}


// Mark off all border edges. Return number of border edges.
label markBorderEdges
(
    const bool debug,
    const triSurface& surf,
    boolList& borderEdge
)
{
    label nBorderEdges = 0;

    const labelListList& edgeFaces = surf.edgeFaces();

    forAll(edgeFaces, edgeI)
    {
        if (edgeFaces[edgeI].size() == 4)
        {
            borderEdge[edgeI] = true;

            nBorderEdges++;
        }
    }

    if (debug)
    {
        dumpEdges(surf, borderEdge);
    }

    return nBorderEdges;
}


// Mark off all border points. Return number of border points. Border points
// marked by setting value to newly introduced point.
label markBorderPoints
(
    const bool debug,
    const triSurface& surf,
    const boolList& borderEdge,
    labelList& borderPoint
)
{
    label nPoints = surf.nPoints();

    const labelListList& pointEdges = surf.pointEdges();

    forAll(pointEdges, pointI)
    {
        const labelList& pEdges = pointEdges[pointI];

        label nBorderEdges = 0;

        forAll(pEdges, i)
        {
            if (borderEdge[pEdges[i]])
            {
                nBorderEdges++;
            }
        }

        if (nBorderEdges == 2 && borderPoint[pointI] == -1)
        {
            borderPoint[pointI] = nPoints++;
        }
    }

    label nBorderPoints = nPoints - surf.nPoints();

    if (debug)
    {
        dumpPoints(surf, borderPoint);
    }

    return nBorderPoints;
}


// Get minumum length of edges connected to pointI
// Serves to get some local length scale.
scalar minEdgeLen(const triSurface& surf, const label pointI)
{
    const pointField& points = surf.localPoints();

    const labelList& pEdges = surf.pointEdges()[pointI];

    scalar minLen = GREAT;

    forAll(pEdges, i)
    {
        label edgeI = pEdges[i];

        scalar len = surf.edges()[edgeI].mag(points);

        if (len < minLen)
        {
            minLen = len;
        }
    }
    return minLen;
}


// Find edge among edgeLabels with endpoints v0,v1
label findEdge
(
    const triSurface& surf,
    const labelList& edgeLabels,
    const label v0,
    const label v1
)
{
    forAll(edgeLabels, i)
    {
        label edgeI = edgeLabels[i];

        const edge& e = surf.edges()[edgeI];

        if
        (
            (
                e.start() == v0
             && e.end() == v1
            )
         || (
                e.start() == v1
             && e.end() == v0
            )
        )
        {
            return edgeI;
        }
    }


    FatalErrorIn("findEdge(..)") << "Cannot find edge with labels " << v0
        << ' ' << v1 << " in candidates " << edgeLabels
        << " with vertices:" << UIndirectList<edge>(surf.edges(), edgeLabels)()
        << abort(FatalError);

    return -1;
}


// Get the other edge connected to pointI on faceI.
label otherEdge
(
    const triSurface& surf,
    const label faceI,
    const label otherEdgeI,
    const label pointI
)
{
    const labelList& fEdges = surf.faceEdges()[faceI];

    forAll(fEdges, i)
    {
        label edgeI = fEdges[i];

        const edge& e = surf.edges()[edgeI];

        if
        (
            edgeI != otherEdgeI
         && (
                e.start() == pointI
             || e.end() == pointI
            )
        )
        {
            return edgeI;
        }
    }

    FatalErrorIn("otherEdge(..)") << "Cannot find other edge on face " << faceI
        << " verts:" << surf.localPoints()[faceI]
        << " connected to point " << pointI
        << " faceEdges:" << UIndirectList<edge>(surf.edges(), fEdges)()
        << abort(FatalError);

    return -1;
}


// Starting from startPoint on startEdge on startFace walk along border
// and insert faces along the way. Walk keeps always one point or one edge
// on the border line.
void walkSplitLine
(
    const triSurface& surf,
    const boolList& borderEdge,
    const labelList& borderPoint,

    const label startFaceI,
    const label startEdgeI,     // is border edge
    const label startPointI,    // is border point

    Map<label>& faceToEdge,
    Map<label>& faceToPoint
)
{
    label faceI = startFaceI;
    label edgeI = startEdgeI;
    label pointI = startPointI;

    do
    {
        //
        // Stick to pointI and walk face-edge-face until back on border edge.
        //

        do
        {
            // Cross face to next edge.
            edgeI = otherEdge(surf, faceI, edgeI, pointI);

            if (borderEdge[edgeI])
            {
                if (!faceToEdge.insert(faceI, edgeI))
                {
                    // Was already visited.
                    return;
                }
                else
                {
                    // First visit to this borderEdge. We're back on borderline.
                    break;
                }
            }
            else if (!faceToPoint.insert(faceI, pointI))
            {
                // Was already visited.
                return;
            }

            // Cross edge to other face
            const labelList& eFaces = surf.edgeFaces()[edgeI];

            if (eFaces.size() != 2)
            {
                FatalErrorIn("walkSplitPoint(..)")
                    << "Can only handle edges with 2 or 4 edges for now."
                    << abort(FatalError);
            }

            if (eFaces[0] == faceI)
            {
                faceI = eFaces[1];
            }
            else if (eFaces[1] == faceI)
            {
                faceI = eFaces[0];
            }
            else
            {
                FatalErrorIn("walkSplitPoint(..)") << abort(FatalError);
            }
        }
        while (true);


        //
        // Back on border edge. Cross to other point on edge.
        //

        pointI = surf.edges()[edgeI].otherVertex(pointI);

        if (borderPoint[pointI] == -1)
        {
            return;
        }

    }
    while (true);
}


// Find second face which is from same surface i.e. has points on the
// shared edge in reverse order.
label sharedFace
(
    const triSurface& surf,
    const label firstFaceI,
    const label sharedEdgeI
)
{
    // Find ordering of face points in edge.

    const edge& e = surf.edges()[sharedEdgeI];

    const triSurface::FaceType& f = surf.localFaces()[firstFaceI];

    label startIndex = findIndex(f, e.start());

    // points in face in same order as edge
    bool edgeOrder = (f[f.fcIndex(startIndex)] == e.end());

    // Get faces using edge in sorted order. (sorted such that walking
    // around them in anti-clockwise order corresponds to edge vector
    // acc. to right-hand rule)
    const labelList& eFaces = surf.sortedEdgeFaces()[sharedEdgeI];

    // Get position of face in sorted edge faces
    label faceIndex = findIndex(eFaces, firstFaceI);

    if (edgeOrder)
    {
        // Get face before firstFaceI
        return eFaces[eFaces.rcIndex(faceIndex)];
    }
    else
    {
        // Get face after firstFaceI
        return eFaces[eFaces.fcIndex(faceIndex)];
    }
}


// Calculate (inward pointing) normals on edges shared by faces in
// faceToEdge and averages them to pointNormals.
void calcPointVecs
(
    const triSurface& surf,
    const Map<label>& faceToEdge,
    const Map<label>& faceToPoint,
    vectorField& borderPointVec
)
{
    const labelListList& sortedEdgeFaces = surf.sortedEdgeFaces();
    const edgeList& edges = surf.edges();
    const pointField& points = surf.localPoints();

    boolList edgeDone(surf.nEdges(), false);

    forAllConstIter(Map<label>, faceToEdge, iter)
    {
        const label edgeI = iter();

        if (!edgeDone[edgeI])
        {
            edgeDone[edgeI] = true;

            // Get the two connected faces in sorted order
            // Note: should have stored this when walking ...

            label face0I = -1;
            label face1I = -1;

            const labelList& eFaces = sortedEdgeFaces[edgeI];

            forAll(eFaces, i)
            {
                label faceI = eFaces[i];

                if (faceToEdge.found(faceI))
                {
                    if (face0I == -1)
                    {
                        face0I = faceI;
                    }
                    else if (face1I == -1)
                    {
                        face1I = faceI;

                        break;
                    }
                }
            }

            if (face0I == -1 && face1I == -1)
            {
                Info<< "Writing surface to errorSurf.obj" << endl;

                surf.write("errorSurf.obj");

                FatalErrorIn("calcPointVecs(..)")
                    << "Cannot find two faces using border edge " << edgeI
                    << " verts:" << edges[edgeI]
                    << " eFaces:" << eFaces << endl
                    << "face0I:" << face0I
                    << " face1I:" << face1I << nl
                    << "faceToEdge:" << faceToEdge << nl
                    << "faceToPoint:" << faceToPoint
                    << "Written surface to errorSurf.obj"
                    << abort(FatalError);
            }

            // Now we have edge and the two faces in counter-clockwise order
            // as seen from edge vector. Calculate normal.
            const edge& e = edges[edgeI];
            vector eVec = e.vec(points);

            // Determine vector as average of the vectors in the two faces.
            // If there is only one face available use only one vector.
            vector midVec(vector::zero);

            if (face0I != -1)
            {
                label v0 = triSurfaceTools::oppositeVertex(surf, face0I, edgeI);
                vector e0 = (points[v0] - points[e.start()]) ^ eVec;
                e0 /= mag(e0);
                midVec = e0;
            }

            if (face1I != -1)
            {
                label v1 = triSurfaceTools::oppositeVertex(surf, face1I, edgeI);
                vector e1 = (points[e.start()] - points[v1]) ^ eVec;
                e1 /= mag(e1);
                midVec += e1;
            }

            scalar magMidVec = mag(midVec);

            if (magMidVec > SMALL)
            {
                midVec /= magMidVec;

                // Average to pointVec
                borderPointVec[e.start()] += midVec;
                borderPointVec[e.end()] += midVec;
            }
        }
    }
}


// Renumbers vertices (of triangles in faceToEdge) of which the pointMap is
// not -1.
void renumberFaces
(
    const triSurface& surf,
    const labelList& pointMap,
    const Map<label>& faceToEdge,
    List<triSurface::FaceType>& newTris
)
{
    forAllConstIter(Map<label>, faceToEdge, iter)
    {
        const label faceI = iter.key();
        const triSurface::FaceType& f = surf.localFaces()[faceI];

        forAll(f, fp)
        {
            if (pointMap[f[fp]] != -1)
            {
                newTris[faceI][fp] = pointMap[f[fp]];
            }
        }
    }
}


// Split all borderEdges that don't have borderPoint. Return true if split
// anything.
bool splitBorderEdges
(
    triSurface& surf,
    const boolList& borderEdge,
    const labelList& borderPoint
)
{
    labelList edgesToBeSplit(surf.nEdges());
    label nSplit = 0;

    forAll(borderEdge, edgeI)
    {
        if (borderEdge[edgeI])
        {
            const edge& e = surf.edges()[edgeI];

            if (borderPoint[e.start()] == -1 && borderPoint[e.end()] == -1)
            {
                // None of the points of the edge is borderPoint. Split edge
                // to introduce border point.
                edgesToBeSplit[nSplit++] = edgeI;
            }
        }
    }
    edgesToBeSplit.setSize(nSplit);

    if (nSplit > 0)
    {
        Info<< "Splitting surface along " << nSplit << " borderEdges that don't"
            << " neighbour other borderEdges" << nl << endl;

        surf = triSurfaceTools::greenRefine(surf, edgesToBeSplit);

        return true;
    }
    else
    {
        Info<< "No edges to be split" <<nl << endl;

        return false;
    }
}



int main(int argc, char *argv[])
{
    argList::addNote
    (
        "split multiply connected surface edges by duplicating points"
    );
    argList::noParallel();
    argList::validArgs.append("surfaceFile");
    argList::validArgs.append("output surfaceFile");
    argList::addBoolOption
    (
        "debug",
        "add debugging output"
    );

    argList args(argc, argv);

    const fileName inSurfName  = args[1];
    const fileName outSurfName = args[2];
    const bool debug = args.optionFound("debug");

    Info<< "Reading surface from " << inSurfName << endl;

    triSurface surf(inSurfName);

    // Make sure sortedEdgeFaces is calculated correctly
    testSortedEdgeFaces(surf);

    // Get all quad connected edges. These are seen as borders when walking.
    boolList borderEdge(surf.nEdges(), false);
    markBorderEdges(debug, surf, borderEdge);

    // Points on two sides connected to borderEdges are called
    // borderPoints and will be duplicated. borderPoint contains label
    // of newly introduced vertex.
    labelList borderPoint(surf.nPoints(), -1);
    markBorderPoints(debug, surf, borderEdge, borderPoint);

    // Split edges where there would be no borderPoint to duplicate.
    splitBorderEdges(surf, borderEdge, borderPoint);


    Info<< "Writing split surface to " << outSurfName << nl << endl;
    surf.write(outSurfName);
    Info<< "Finished writing surface to " << outSurfName << nl << endl;


    // Last iteration values.
    label nOldBorderEdges = -1;
    label nOldBorderPoints = -1;

    label iteration = 0;

    do
    {
        // Redo borderEdge/borderPoint calculation.
        boolList borderEdge(surf.nEdges(), false);

        label nBorderEdges = markBorderEdges(debug, surf, borderEdge);

        if (nBorderEdges == 0)
        {
            Info<< "Found no border edges. Exiting." << nl << nl;

            break;
        }

        // Label of newly introduced duplicate.
        labelList borderPoint(surf.nPoints(), -1);

        label nBorderPoints =
            markBorderPoints
            (
                debug,
                surf,
                borderEdge,
                borderPoint
            );

        if (nBorderPoints == 0)
        {
            Info<< "Found no border points. Exiting." << nl << nl;

            break;
        }

        Info<< "Found:\n"
            << "    border edges :" << nBorderEdges << nl
            << "    border points:" << nBorderPoints << nl
            << endl;

        if
        (
            nBorderPoints == nOldBorderPoints
         && nBorderEdges == nOldBorderEdges
        )
        {
            Info<< "Stopping since number of border edges and point is same"
                << " as in previous iteration" << nl << endl;

            break;
        }

        //
        // Define splitLine as a series of connected borderEdges. Find start
        // of one (as edge and point on edge)
        //

        label startEdgeI = -1;
        label startPointI = -1;

        forAll(borderEdge, edgeI)
        {
            if (borderEdge[edgeI])
            {
                const edge& e = surf.edges()[edgeI];

                if ((borderPoint[e[0]] != -1) && (borderPoint[e[1]] == -1))
                {
                    startEdgeI = edgeI;
                    startPointI = e[0];

                    break;
                }
                else if ((borderPoint[e[0]] == -1) && (borderPoint[e[1]] != -1))
                {
                    startEdgeI = edgeI;
                    startPointI = e[1];

                    break;
                }
            }
        }

        if (startEdgeI == -1)
        {
            Info<< "Cannot find starting point of splitLine\n" << endl;

            break;
        }

        // Pick any face using edge to start from.
        const labelList& eFaces = surf.edgeFaces()[startEdgeI];

        label firstFaceI = eFaces[0];

        // Find second face which is from same surface i.e. has outwards
        // pointing normal as well (actually bit more complex than this)
        label secondFaceI = sharedFace(surf, firstFaceI, startEdgeI);

        Info<< "Starting local walk from:" << endl
            << "    edge :" << startEdgeI << endl
            << "    point:" << startPointI << endl
            << "    face0:" << firstFaceI << endl
            << "    face1:" << secondFaceI << endl
            << endl;

        // From face on border edge to edge.
        Map<label> faceToEdge(2*nBorderEdges);
        // From face connected to border point (but not border edge) to point.
        Map<label> faceToPoint(2*nBorderPoints);

        faceToEdge.insert(firstFaceI, startEdgeI);

        walkSplitLine
        (
            surf,
            borderEdge,
            borderPoint,

            firstFaceI,
            startEdgeI,
            startPointI,

            faceToEdge,
            faceToPoint
        );

        faceToEdge.insert(secondFaceI, startEdgeI);

        walkSplitLine
        (
            surf,
            borderEdge,
            borderPoint,

            secondFaceI,
            startEdgeI,
            startPointI,

            faceToEdge,
            faceToPoint
        );

        Info<< "Finished local walk and visited" << nl
            << "    border edges                :" << faceToEdge.size() << nl
            << "    border points(but not edges):" << faceToPoint.size() << nl
            << endl;

        if (debug)
        {
            dumpFaces("faceToEdge.obj", surf, faceToEdge);
            dumpFaces("faceToPoint.obj", surf, faceToPoint);
        }

        //
        // Create coordinates for borderPoints by duplicating the existing
        // point and then slightly shifting it inwards. To determine the
        // inwards direction get the average normal of both connectedFaces on
        // the edge and then interpolate this to the (border)point.
        //

        vectorField borderPointVec(surf.nPoints(), vector(GREAT, GREAT, GREAT));

        calcPointVecs(surf, faceToEdge, faceToPoint, borderPointVec);


        // New position. Start off from copy of old points.
        pointField newPoints(surf.localPoints());
        newPoints.setSize(newPoints.size() + nBorderPoints);

        forAll(borderPoint, pointI)
        {
            label newPointI = borderPoint[pointI];

            if (newPointI != -1)
            {
                scalar minLen = minEdgeLen(surf, pointI);

                vector n = borderPointVec[pointI];
                n /= mag(n);

                newPoints[newPointI] = newPoints[pointI] + 0.1 * minLen * n;
            }
        }


        //
        // Renumber all faces in connectedFaces
        //

        // Start off from copy of faces.
        List<labelledTri> newTris(surf.size());

        forAll(surf, faceI)
        {
            newTris[faceI] = surf.localFaces()[faceI];
            newTris[faceI].region() = surf[faceI].region();
        }

        // Renumber all faces in faceToEdge
        renumberFaces(surf, borderPoint, faceToEdge, newTris);
        // Renumber all faces in faceToPoint
        renumberFaces(surf, borderPoint, faceToPoint, newTris);


        // Check if faces use unmoved points.
        forAll(newTris, faceI)
        {
            const triSurface::FaceType& f = newTris[faceI];

            forAll(f, fp)
            {
                const point& pt = newPoints[f[fp]];

                if (mag(pt) >= GREAT/2)
                {
                    Info<< "newTri:" << faceI << " verts:" << f
                        << " vert:" << f[fp] << " point:" << pt << endl;
                }
            }
        }


        surf = triSurface(newTris, surf.patches(), newPoints);

        if (debug || (iteration != 0 && (iteration % 20) == 0))
        {
            Info<< "Writing surface to " << outSurfName << nl << endl;

            surf.write(outSurfName);

            Info<< "Finished writing surface to " << outSurfName << nl << endl;
        }

        // Save prev iteration values.
        nOldBorderEdges = nBorderEdges;
        nOldBorderPoints = nBorderPoints;

        iteration++;
    }
    while (true);

    Info<< "Writing final surface to " << outSurfName << nl << endl;

    surf.write(outSurfName);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
