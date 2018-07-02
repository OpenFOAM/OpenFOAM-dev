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

    forAll(borderPoint, pointi)
    {
        if (borderPoint[pointi] != -1)
        {
            const point& pt = surf.localPoints()[pointi];

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
                FatalErrorInFunction << abort(FatalError);
            }
        }
        forAll(sortMyFaces, i)
        {
            if (findIndex(myFaces, sortMyFaces[i]) == -1)
            {
                FatalErrorInFunction << abort(FatalError);
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

    forAll(pointEdges, pointi)
    {
        const labelList& pEdges = pointEdges[pointi];

        label nBorderEdges = 0;

        forAll(pEdges, i)
        {
            if (borderEdge[pEdges[i]])
            {
                nBorderEdges++;
            }
        }

        if (nBorderEdges == 2 && borderPoint[pointi] == -1)
        {
            borderPoint[pointi] = nPoints++;
        }
    }

    label nBorderPoints = nPoints - surf.nPoints();

    if (debug)
    {
        dumpPoints(surf, borderPoint);
    }

    return nBorderPoints;
}


// Get minimum length of edges connected to pointi
// Serves to get some local length scale.
scalar minEdgeLen(const triSurface& surf, const label pointi)
{
    const pointField& points = surf.localPoints();

    const labelList& pEdges = surf.pointEdges()[pointi];

    scalar minLen = great;

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


    FatalErrorInFunction
        << ' ' << v1 << " in candidates " << edgeLabels
        << " with vertices:" << UIndirectList<edge>(surf.edges(), edgeLabels)()
        << abort(FatalError);

    return -1;
}


// Get the other edge connected to pointi on facei.
label otherEdge
(
    const triSurface& surf,
    const label facei,
    const label otherEdgeI,
    const label pointi
)
{
    const labelList& fEdges = surf.faceEdges()[facei];

    forAll(fEdges, i)
    {
        label edgeI = fEdges[i];

        const edge& e = surf.edges()[edgeI];

        if
        (
            edgeI != otherEdgeI
         && (
                e.start() == pointi
             || e.end() == pointi
            )
        )
        {
            return edgeI;
        }
    }

    FatalErrorInFunction
        << " verts:" << surf.localPoints()[facei]
        << " connected to point " << pointi
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

    const label startFacei,
    const label startEdgeI,     // is border edge
    const label startPointi,    // is border point

    Map<label>& faceToEdge,
    Map<label>& faceToPoint
)
{
    label facei = startFacei;
    label edgeI = startEdgeI;
    label pointi = startPointi;

    do
    {
        //
        // Stick to pointi and walk face-edge-face until back on border edge.
        //

        do
        {
            // Cross face to next edge.
            edgeI = otherEdge(surf, facei, edgeI, pointi);

            if (borderEdge[edgeI])
            {
                if (!faceToEdge.insert(facei, edgeI))
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
            else if (!faceToPoint.insert(facei, pointi))
            {
                // Was already visited.
                return;
            }

            // Cross edge to other face
            const labelList& eFaces = surf.edgeFaces()[edgeI];

            if (eFaces.size() != 2)
            {
                FatalErrorInFunction
                    << "Can only handle edges with 2 or 4 edges for now."
                    << abort(FatalError);
            }

            if (eFaces[0] == facei)
            {
                facei = eFaces[1];
            }
            else if (eFaces[1] == facei)
            {
                facei = eFaces[0];
            }
            else
            {
                FatalErrorInFunction << abort(FatalError);
            }
        }
        while (true);


        //
        // Back on border edge. Cross to other point on edge.
        //

        pointi = surf.edges()[edgeI].otherVertex(pointi);

        if (borderPoint[pointi] == -1)
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
    const label firstFacei,
    const label sharedEdgeI
)
{
    // Find ordering of face points in edge.

    const edge& e = surf.edges()[sharedEdgeI];

    const triSurface::FaceType& f = surf.localFaces()[firstFacei];

    label startIndex = findIndex(f, e.start());

    // points in face in same order as edge
    bool edgeOrder = (f[f.fcIndex(startIndex)] == e.end());

    // Get faces using edge in sorted order. (sorted such that walking
    // around them in anti-clockwise order corresponds to edge vector
    // acc. to right-hand rule)
    const labelList& eFaces = surf.sortedEdgeFaces()[sharedEdgeI];

    // Get position of face in sorted edge faces
    label faceIndex = findIndex(eFaces, firstFacei);

    if (edgeOrder)
    {
        // Get face before firstFacei
        return eFaces[eFaces.rcIndex(faceIndex)];
    }
    else
    {
        // Get face after firstFacei
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
                label facei = eFaces[i];

                if (faceToEdge.found(facei))
                {
                    if (face0I == -1)
                    {
                        face0I = facei;
                    }
                    else if (face1I == -1)
                    {
                        face1I = facei;

                        break;
                    }
                }
            }

            if (face0I == -1 && face1I == -1)
            {
                Info<< "Writing surface to errorSurf.obj" << endl;

                surf.write("errorSurf.obj");

                FatalErrorInFunction
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
            vector midVec(Zero);

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

            if (magMidVec > small)
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
        const label facei = iter.key();
        const triSurface::FaceType& f = surf.localFaces()[facei];

        forAll(f, fp)
        {
            if (pointMap[f[fp]] != -1)
            {
                newTris[facei][fp] = pointMap[f[fp]];
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
    #include "removeCaseOptions.H"

    argList::addNote
    (
        "split multiply connected surface edges by duplicating points"
    );
    argList::validArgs.append("surface file");
    argList::validArgs.append("output surface file");
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
        label startPointi = -1;

        forAll(borderEdge, edgeI)
        {
            if (borderEdge[edgeI])
            {
                const edge& e = surf.edges()[edgeI];

                if ((borderPoint[e[0]] != -1) && (borderPoint[e[1]] == -1))
                {
                    startEdgeI = edgeI;
                    startPointi = e[0];

                    break;
                }
                else if ((borderPoint[e[0]] == -1) && (borderPoint[e[1]] != -1))
                {
                    startEdgeI = edgeI;
                    startPointi = e[1];

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

        label firstFacei = eFaces[0];

        // Find second face which is from same surface i.e. has outwards
        // pointing normal as well (actually bit more complex than this)
        label secondFacei = sharedFace(surf, firstFacei, startEdgeI);

        Info<< "Starting local walk from:" << endl
            << "    edge :" << startEdgeI << endl
            << "    point:" << startPointi << endl
            << "    face0:" << firstFacei << endl
            << "    face1:" << secondFacei << endl
            << endl;

        // From face on border edge to edge.
        Map<label> faceToEdge(2*nBorderEdges);
        // From face connected to border point (but not border edge) to point.
        Map<label> faceToPoint(2*nBorderPoints);

        faceToEdge.insert(firstFacei, startEdgeI);

        walkSplitLine
        (
            surf,
            borderEdge,
            borderPoint,

            firstFacei,
            startEdgeI,
            startPointi,

            faceToEdge,
            faceToPoint
        );

        faceToEdge.insert(secondFacei, startEdgeI);

        walkSplitLine
        (
            surf,
            borderEdge,
            borderPoint,

            secondFacei,
            startEdgeI,
            startPointi,

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

        vectorField borderPointVec(surf.nPoints(), vector(great, great, great));

        calcPointVecs(surf, faceToEdge, faceToPoint, borderPointVec);


        // New position. Start off from copy of old points.
        pointField newPoints(surf.localPoints());
        newPoints.setSize(newPoints.size() + nBorderPoints);

        forAll(borderPoint, pointi)
        {
            label newPointi = borderPoint[pointi];

            if (newPointi != -1)
            {
                scalar minLen = minEdgeLen(surf, pointi);

                vector n = borderPointVec[pointi];
                n /= mag(n);

                newPoints[newPointi] = newPoints[pointi] + 0.1 * minLen * n;
            }
        }


        //
        // Renumber all faces in connectedFaces
        //

        // Start off from copy of faces.
        List<labelledTri> newTris(surf.size());

        forAll(surf, facei)
        {
            newTris[facei] = surf.localFaces()[facei];
            newTris[facei].region() = surf[facei].region();
        }

        // Renumber all faces in faceToEdge
        renumberFaces(surf, borderPoint, faceToEdge, newTris);
        // Renumber all faces in faceToPoint
        renumberFaces(surf, borderPoint, faceToPoint, newTris);


        // Check if faces use unmoved points.
        forAll(newTris, facei)
        {
            const triSurface::FaceType& f = newTris[facei];

            forAll(f, fp)
            {
                const point& pt = newPoints[f[fp]];

                if (mag(pt) >= great/2)
                {
                    Info<< "newTri:" << facei << " verts:" << f
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
