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

#include "triSurfaceTools.H"

#include "triSurface.H"
#include "OFstream.H"
#include "mergePoints.H"
#include "polyMesh.H"
#include "plane.H"
#include "geompack.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::label Foam::triSurfaceTools::ANYEDGE = -1;
const Foam::label Foam::triSurfaceTools::NOEDGE = -2;
const Foam::label Foam::triSurfaceTools::COLLAPSED = -3;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

/*
    Refine by splitting all three edges of triangle ('red' refinement).
    Neighbouring triangles (which are not marked for refinement get split
    in half ('green') refinement. (R. Verfuerth, "A review of a posteriori
    error estimation and adaptive mesh refinement techniques",
    Wiley-Teubner, 1996)
*/

// FaceI gets refined ('red'). Propagate edge refinement.
void Foam::triSurfaceTools::calcRefineStatus
(
    const triSurface& surf,
    const label faceI,
    List<refineType>& refine
)
{
    if (refine[faceI] == RED)
    {
        // Already marked for refinement. Do nothing.
    }
    else
    {
        // Not marked or marked for 'green' refinement. Refine.
        refine[faceI] = RED;

        const labelList& myNeighbours = surf.faceFaces()[faceI];

        forAll(myNeighbours, myNeighbourI)
        {
            label neighbourFaceI = myNeighbours[myNeighbourI];

            if (refine[neighbourFaceI] == GREEN)
            {
                // Change to red refinement and propagate
                calcRefineStatus(surf, neighbourFaceI, refine);
            }
            else if (refine[neighbourFaceI] == NONE)
            {
                refine[neighbourFaceI] = GREEN;
            }
        }
    }
}


// Split faceI along edgeI at position newPointI
void Foam::triSurfaceTools::greenRefine
(
    const triSurface& surf,
    const label faceI,
    const label edgeI,
    const label newPointI,
    DynamicList<labelledTri>& newFaces
)
{
    const labelledTri& f = surf.localFaces()[faceI];
    const edge& e = surf.edges()[edgeI];

    // Find index of edge in face.

    label fp0 = findIndex(f, e[0]);
    label fp1 = f.fcIndex(fp0);
    label fp2 = f.fcIndex(fp1);

    if (f[fp1] == e[1])
    {
        // Edge oriented like face
        newFaces.append
        (
            labelledTri
            (
                f[fp0],
                newPointI,
                f[fp2],
                f.region()
            )
        );
        newFaces.append
        (
            labelledTri
            (
                newPointI,
                f[fp1],
                f[fp2],
                f.region()
            )
        );
    }
    else
    {
        newFaces.append
        (
            labelledTri
            (
                f[fp2],
                newPointI,
                f[fp1],
                f.region()
            )
        );
        newFaces.append
        (
            labelledTri
            (
                newPointI,
                f[fp0],
                f[fp1],
                f.region()
            )
        );
    }
}


// Refine all triangles marked for refinement.
Foam::triSurface Foam::triSurfaceTools::doRefine
(
    const triSurface& surf,
    const List<refineType>& refineStatus
)
{
    // Storage for new points. (start after old points)
    DynamicList<point> newPoints(surf.nPoints());
    forAll(surf.localPoints(), pointI)
    {
        newPoints.append(surf.localPoints()[pointI]);
    }
    label newVertI = surf.nPoints();

    // Storage for new faces
    DynamicList<labelledTri> newFaces(surf.size());


    // Point index for midpoint on edge
    labelList edgeMid(surf.nEdges(), -1);

    forAll(refineStatus, faceI)
    {
        if (refineStatus[faceI] == RED)
        {
            // Create new vertices on all edges to be refined.
            const labelList& fEdges = surf.faceEdges()[faceI];

            forAll(fEdges, i)
            {
                label edgeI = fEdges[i];

                if (edgeMid[edgeI] == -1)
                {
                    const edge& e = surf.edges()[edgeI];

                    // Create new point on mid of edge
                    newPoints.append
                    (
                        0.5
                      * (
                            surf.localPoints()[e.start()]
                          + surf.localPoints()[e.end()]
                        )
                    );
                    edgeMid[edgeI] = newVertI++;
                }
            }

            // Now we have new mid edge vertices for all edges on face
            // so create triangles for RED rerfinement.

            const edgeList& edges = surf.edges();

            // Corner triangles
            newFaces.append
            (
                labelledTri
                (
                    edgeMid[fEdges[0]],
                    edges[fEdges[0]].commonVertex(edges[fEdges[1]]),
                    edgeMid[fEdges[1]],
                    surf[faceI].region()
                )
            );

            newFaces.append
            (
                labelledTri
                (
                    edgeMid[fEdges[1]],
                    edges[fEdges[1]].commonVertex(edges[fEdges[2]]),
                    edgeMid[fEdges[2]],
                    surf[faceI].region()
                )
            );

            newFaces.append
            (
                labelledTri
                (
                    edgeMid[fEdges[2]],
                    edges[fEdges[2]].commonVertex(edges[fEdges[0]]),
                    edgeMid[fEdges[0]],
                    surf[faceI].region()
                )
            );

            // Inner triangle
            newFaces.append
            (
                labelledTri
                (
                    edgeMid[fEdges[0]],
                    edgeMid[fEdges[1]],
                    edgeMid[fEdges[2]],
                    surf[faceI].region()
                )
            );


            // Create triangles for GREEN refinement.
            forAll(fEdges, i)
            {
                const label edgeI = fEdges[i];

                label otherFaceI = otherFace(surf, faceI, edgeI);

                if ((otherFaceI != -1) && (refineStatus[otherFaceI] == GREEN))
                {
                    greenRefine
                    (
                        surf,
                        otherFaceI,
                        edgeI,
                        edgeMid[edgeI],
                        newFaces
                    );
                }
            }
        }
    }

    // Copy unmarked triangles since keep original vertex numbering.
    forAll(refineStatus, faceI)
    {
        if (refineStatus[faceI] == NONE)
        {
            newFaces.append(surf.localFaces()[faceI]);
        }
    }

    newFaces.shrink();
    newPoints.shrink();


    // Transfer DynamicLists to straight ones.
    pointField allPoints;
    allPoints.transfer(newPoints);
    newPoints.clear();

    return triSurface(newFaces, surf.patches(), allPoints, true);
}


// Angle between two neighbouring triangles,
// angle between shared-edge end points and left and right face end points
Foam::scalar Foam::triSurfaceTools::faceCosAngle
(
    const point& pStart,
    const point& pEnd,
    const point& pLeft,
    const point& pRight
)
{
    const vector common(pEnd - pStart);
    const vector base0(pLeft - pStart);
    const vector base1(pRight - pStart);

    vector n0(common ^ base0);
    n0 /= Foam::mag(n0);

    vector n1(base1 ^ common);
    n1 /= Foam::mag(n1);

    return n0 & n1;
}


// Protect edges around vertex from collapsing.
// Moving vertI to a different position
// - affects obviously the cost of the faces using it
// - but also their neighbours since the edge between the faces changes
void Foam::triSurfaceTools::protectNeighbours
(
    const triSurface& surf,
    const label vertI,
    labelList& faceStatus
)
{
//    const labelList& myFaces = surf.pointFaces()[vertI];
//    forAll(myFaces, i)
//    {
//        label faceI = myFaces[i];
//
//        if ((faceStatus[faceI] == ANYEDGE) || (faceStatus[faceI] >= 0))
//        {
//            faceStatus[faceI] = NOEDGE;
//        }
//    }

    const labelList& myEdges = surf.pointEdges()[vertI];
    forAll(myEdges, i)
    {
        const labelList& myFaces = surf.edgeFaces()[myEdges[i]];

        forAll(myFaces, myFaceI)
        {
            label faceI = myFaces[myFaceI];

            if ((faceStatus[faceI] == ANYEDGE) || (faceStatus[faceI] >= 0))
            {
                faceStatus[faceI] = NOEDGE;
            }
       }
    }
}


//
// Edge collapse helper functions
//

// Get all faces that will get collapsed if edgeI collapses.
Foam::labelHashSet Foam::triSurfaceTools::getCollapsedFaces
(
    const triSurface& surf,
    label edgeI
)
{
    const edge& e = surf.edges()[edgeI];
    label v1 = e.start();
    label v2 = e.end();

    // Faces using edge will certainly get collapsed.
    const labelList& myFaces = surf.edgeFaces()[edgeI];

    labelHashSet facesToBeCollapsed(2*myFaces.size());

    forAll(myFaces, myFaceI)
    {
        facesToBeCollapsed.insert(myFaces[myFaceI]);
    }

    // From faces using v1 check if they share an edge with faces
    // using v2.
    //  - share edge: are part of 'splay' tree and will collapse if edge
    //    collapses
    const labelList& v1Faces = surf.pointFaces()[v1];

    forAll(v1Faces, v1FaceI)
    {
        label face1I = v1Faces[v1FaceI];

        label otherEdgeI = oppositeEdge(surf, face1I, v1);

        // Step across edge to other face
        label face2I = otherFace(surf, face1I, otherEdgeI);

        if (face2I != -1)
        {
            // found face on other side of edge. Now check if uses v2.
            if (oppositeVertex(surf, face2I, otherEdgeI) == v2)
            {
                // triangles face1I and face2I form splay tree and will
                // collapse.
                facesToBeCollapsed.insert(face1I);
                facesToBeCollapsed.insert(face2I);
            }
        }
    }

    return facesToBeCollapsed;
}


// Return value of faceUsed for faces using vertI
Foam::label Foam::triSurfaceTools::vertexUsesFace
(
    const triSurface& surf,
    const labelHashSet& faceUsed,
    const label vertI
)
{
    const labelList& myFaces = surf.pointFaces()[vertI];

    forAll(myFaces, myFaceI)
    {
        label face1I = myFaces[myFaceI];

        if (faceUsed.found(face1I))
        {
            return face1I;
        }
    }
    return -1;
}


// Calculate new edge-face addressing resulting from edge collapse
void Foam::triSurfaceTools::getMergedEdges
(
    const triSurface& surf,
    const label edgeI,
    const labelHashSet& collapsedFaces,
    HashTable<label, label, Hash<label> >& edgeToEdge,
    HashTable<label, label, Hash<label> >& edgeToFace
)
{
    const edge& e = surf.edges()[edgeI];
    label v1 = e.start();
    label v2 = e.end();

    const labelList& v1Faces = surf.pointFaces()[v1];
    const labelList& v2Faces = surf.pointFaces()[v2];

    // Mark all (non collapsed) faces using v2
    labelHashSet v2FacesHash(v2Faces.size());

    forAll(v2Faces, v2FaceI)
    {
        if (!collapsedFaces.found(v2Faces[v2FaceI]))
        {
            v2FacesHash.insert(v2Faces[v2FaceI]);
        }
    }


    forAll(v1Faces, v1FaceI)
    {
        label face1I = v1Faces[v1FaceI];

        if (collapsedFaces.found(face1I))
        {
            continue;
        }

        //
        // Check if face1I uses a vertex connected to a face using v2
        //

        label vert1I = -1;
        label vert2I = -1;
        otherVertices
        (
            surf,
            face1I,
            v1,
            vert1I,
            vert2I
        );
        //Pout<< "Face:" << surf.localFaces()[face1I] << " other vertices:"
        //    << vert1I << ' ' << vert2I << endl;

        // Check vert1, vert2 for usage by v2Face.

        label commonVert = vert1I;
        label face2I = vertexUsesFace(surf, v2FacesHash, commonVert);
        if (face2I == -1)
        {
            commonVert = vert2I;
            face2I = vertexUsesFace(surf, v2FacesHash, commonVert);
        }

        if (face2I != -1)
        {
            // Found one: commonVert is used by both face1I and face2I
            label edge1I = getEdge(surf, v1, commonVert);
            label edge2I = getEdge(surf, v2, commonVert);

            edgeToEdge.insert(edge1I, edge2I);
            edgeToEdge.insert(edge2I, edge1I);

            edgeToFace.insert(edge1I, face2I);
            edgeToFace.insert(edge2I, face1I);
        }
    }
}


// Calculates (cos of) angle across edgeI of faceI,
// taking into account updated addressing (resulting from edge collapse)
Foam::scalar Foam::triSurfaceTools::edgeCosAngle
(
    const triSurface& surf,
    const label v1,
    const point& pt,
    const labelHashSet& collapsedFaces,
    const HashTable<label, label, Hash<label> >& edgeToEdge,
    const HashTable<label, label, Hash<label> >& edgeToFace,
    const label faceI,
    const label edgeI
)
{
    const pointField& localPoints = surf.localPoints();

    label A = surf.edges()[edgeI].start();
    label B = surf.edges()[edgeI].end();
    label C = oppositeVertex(surf, faceI, edgeI);

    label D = -1;

    label face2I = -1;

    if (edgeToEdge.found(edgeI))
    {
        // Use merged addressing
        label edge2I = edgeToEdge[edgeI];
        face2I = edgeToFace[edgeI];

        D = oppositeVertex(surf, face2I, edge2I);
    }
    else
    {
        // Use normal edge-face addressing
        face2I = otherFace(surf, faceI, edgeI);

        if ((face2I != -1) && !collapsedFaces.found(face2I))
        {
            D = oppositeVertex(surf, face2I, edgeI);
        }
    }

    scalar cosAngle = 1;

    if (D != -1)
    {
        if (A == v1)
        {
            cosAngle = faceCosAngle
            (
                pt,
                localPoints[B],
                localPoints[C],
                localPoints[D]
            );
        }
        else if (B == v1)
        {
            cosAngle = faceCosAngle
            (
                localPoints[A],
                pt,
                localPoints[C],
                localPoints[D]
            );
        }
        else if (C == v1)
        {
            cosAngle = faceCosAngle
            (
                localPoints[A],
                localPoints[B],
                pt,
                localPoints[D]
            );
        }
        else if (D == v1)
        {
            cosAngle = faceCosAngle
            (
                localPoints[A],
                localPoints[B],
                localPoints[C],
                pt
            );
        }
        else
        {
            FatalErrorIn("edgeCosAngle")
                << "face " << faceI << " does not use vertex "
                << v1 << " of collapsed edge" << abort(FatalError);
        }
    }
    return cosAngle;
}


//- Calculate minimum (cos of) edge angle using addressing from collapsing
//  edge to v1 at pt.
Foam::scalar Foam::triSurfaceTools::collapseMinCosAngle
(
    const triSurface& surf,
    const label v1,
    const point& pt,
    const labelHashSet& collapsedFaces,
    const HashTable<label, label, Hash<label> >& edgeToEdge,
    const HashTable<label, label, Hash<label> >& edgeToFace
)
{
    const labelList& v1Faces = surf.pointFaces()[v1];

    scalar minCos = 1;

    forAll(v1Faces, v1FaceI)
    {
        label faceI = v1Faces[v1FaceI];

        if (collapsedFaces.found(faceI))
        {
            continue;
        }

        const labelList& myEdges = surf.faceEdges()[faceI];

        forAll(myEdges, myEdgeI)
        {
            label edgeI = myEdges[myEdgeI];

            minCos =
                min
                (
                    minCos,
                    edgeCosAngle
                    (
                        surf,
                        v1,
                        pt,
                        collapsedFaces,
                        edgeToEdge,
                        edgeToFace,
                        faceI,
                        edgeI
                    )
                );
        }
    }

    return minCos;
}


// Calculate max dihedral angle after collapsing edge to v1 (at pt).
// Note that all edges of all faces using v1 are affected.
bool Foam::triSurfaceTools::collapseCreatesFold
(
    const triSurface& surf,
    const label v1,
    const point& pt,
    const labelHashSet& collapsedFaces,
    const HashTable<label, label, Hash<label> >& edgeToEdge,
    const HashTable<label, label, Hash<label> >& edgeToFace,
    const scalar minCos
)
{
    const labelList& v1Faces = surf.pointFaces()[v1];

    forAll(v1Faces, v1FaceI)
    {
        label faceI = v1Faces[v1FaceI];

        if (collapsedFaces.found(faceI))
        {
            continue;
        }

        const labelList& myEdges = surf.faceEdges()[faceI];

        forAll(myEdges, myEdgeI)
        {
            label edgeI = myEdges[myEdgeI];

            if
            (
                edgeCosAngle
                (
                    surf,
                    v1,
                    pt,
                    collapsedFaces,
                    edgeToEdge,
                    edgeToFace,
                    faceI,
                    edgeI
                )
              < minCos
            )
            {
                return true;
            }
        }
    }

    return false;
}


//// Return true if collapsing edgeI creates triangles on top of each other.
//// Is when the triangles neighbouring the collapsed one already share
// a vertex.
//bool Foam::triSurfaceTools::collapseCreatesDuplicates
//(
//    const triSurface& surf,
//    const label edgeI,
//    const labelHashSet& collapsedFaces
//)
//{
//Pout<< "duplicate : edgeI:" << edgeI << surf.edges()[edgeI]
//    << " collapsedFaces:" << collapsedFaces.toc() << endl;
//
//    // Mark neighbours of faces to be collapsed, i.e. get the first layer
//    // of triangles outside collapsedFaces.
//    // neighbours actually contains the
//    // edge with which triangle connects to collapsedFaces.
//
//    HashTable<label, label, Hash<label> > neighbours;
//
//    labelList collapsed = collapsedFaces.toc();
//
//    forAll(collapsed, collapseI)
//    {
//        const label faceI = collapsed[collapseI];
//
//        const labelList& myEdges = surf.faceEdges()[faceI];
//
//        Pout<< "collapsing faceI:" << faceI << " uses edges:" << myEdges
//            << endl;
//
//        forAll(myEdges, myEdgeI)
//        {
//            const labelList& myFaces = surf.edgeFaces()[myEdges[myEdgeI]];
//
//            Pout<< "Edge " << myEdges[myEdgeI] << " is used by faces "
//                << myFaces << endl;
//
//            if ((myEdges[myEdgeI] != edgeI) && (myFaces.size() == 2))
//            {
//                // Get the other face
//                label neighbourFaceI = myFaces[0];
//
//                if (neighbourFaceI == faceI)
//                {
//                    neighbourFaceI = myFaces[1];
//                }
//
//                // Is 'outside' face if not in collapsedFaces itself
//                if (!collapsedFaces.found(neighbourFaceI))
//                {
//                    neighbours.insert(neighbourFaceI, myEdges[myEdgeI]);
//                }
//            }
//        }
//    }
//
//    // Now neighbours contains first layer of triangles outside of
//    // collapseFaces
//    // There should be
//    // -two if edgeI is a boundary edge
//    // since the outside 'edge' of collapseFaces should
//    // form a triangle and the face connected to edgeI is not inserted.
//    // -four if edgeI is not a boundary edge since then the outside edge forms
//    // a diamond.
//
//    // Check if any of the faces in neighbours share an edge. (n^2)
//
//    labelList neighbourList = neighbours.toc();
//
//Pout<< "edgeI:" << edgeI << "  neighbourList:" << neighbourList << endl;
//
//
//    forAll(neighbourList, i)
//    {
//        const labelList& faceIEdges = surf.faceEdges()[neighbourList[i]];
//
//        for (label j = i+1; j < neighbourList.size(); i++)
//        {
//            const labelList& faceJEdges = surf.faceEdges()[neighbourList[j]];
//
//            // Check if faceI and faceJ share an edge
//            forAll(faceIEdges, fI)
//            {
//                forAll(faceJEdges, fJ)
//                {
//                    Pout<< " comparing " << faceIEdges[fI] << " to "
//                        << faceJEdges[fJ] << endl;
//
//                    if (faceIEdges[fI] == faceJEdges[fJ])
//                    {
//                        return true;
//                    }
//                }
//            }
//        }
//    }
//    Pout<< "Found no match. Returning false" << endl;
//    return false;
//}


// Finds the triangle edge cut by the plane between a point inside / on edge
// of a triangle and a point outside. Returns:
//  - cut through edge/point
//  - miss()
Foam::surfaceLocation Foam::triSurfaceTools::cutEdge
(
    const triSurface& s,
    const label triI,
    const label excludeEdgeI,
    const label excludePointI,

    const point& triPoint,
    const plane& cutPlane,
    const point& toPoint
)
{
    const pointField& points = s.points();
    const labelledTri& f = s[triI];
    const labelList& fEdges = s.faceEdges()[triI];

    // Get normal distance to planeN
    FixedList<scalar, 3> d;

    scalar norm = 0;
    forAll(d, fp)
    {
        d[fp] = (points[f[fp]]-cutPlane.refPoint()) & cutPlane.normal();
        norm += mag(d[fp]);
    }

    // Normalise and truncate
    forAll(d, i)
    {
        d[i] /= norm;

        if (mag(d[i]) < 1e-6)
        {
            d[i] = 0.0;
        }
    }

    // Return information
    surfaceLocation cut;

    if (excludePointI != -1)
    {
        // Excluded point. Test only opposite edge.

        label fp0 = findIndex(s.localFaces()[triI], excludePointI);

        if (fp0 == -1)
        {
            FatalErrorIn("cutEdge(..)") << "excludePointI:" << excludePointI
                << " localF:" << s.localFaces()[triI] << abort(FatalError);
        }

        label fp1 = f.fcIndex(fp0);
        label fp2 = f.fcIndex(fp1);


        if (d[fp1] == 0.0)
        {
            cut.setHit();
            cut.setPoint(points[f[fp1]]);
            cut.elementType() = triPointRef::POINT;
            cut.setIndex(s.localFaces()[triI][fp1]);
        }
        else if (d[fp2] == 0.0)
        {
            cut.setHit();
            cut.setPoint(points[f[fp2]]);
            cut.elementType() = triPointRef::POINT;
            cut.setIndex(s.localFaces()[triI][fp2]);
        }
        else if
        (
            (d[fp1] < 0 && d[fp2] < 0)
         || (d[fp1] > 0 && d[fp2] > 0)
        )
        {
            // Both same sign. Not crossing edge at all.
            // cut already set to miss().
        }
        else
        {
            cut.setHit();
            cut.setPoint
            (
                (d[fp2]*points[f[fp1]] - d[fp1]*points[f[fp2]])
              / (d[fp2] - d[fp1])
            );
            cut.elementType() = triPointRef::EDGE;
            cut.setIndex(fEdges[fp1]);
        }
    }
    else
    {
        // Find the two intersections
        FixedList<surfaceLocation, 2> inters;
        label interI = 0;

        forAll(f, fp0)
        {
            label fp1 = f.fcIndex(fp0);

            if (d[fp0] == 0)
            {
                if (interI >= 2)
                {
                    FatalErrorIn("cutEdge(..)")
                        << "problem : triangle has three intersections." << nl
                        << "triangle:" << f.tri(points)
                        << " d:" << d << abort(FatalError);
                }
                inters[interI].setHit();
                inters[interI].setPoint(points[f[fp0]]);
                inters[interI].elementType() = triPointRef::POINT;
                inters[interI].setIndex(s.localFaces()[triI][fp0]);
                interI++;
            }
            else if
            (
                (d[fp0] < 0 && d[fp1] > 0)
             || (d[fp0] > 0 && d[fp1] < 0)
            )
            {
                if (interI >= 2)
                {
                    FatalErrorIn("cutEdge(..)")
                        << "problem : triangle has three intersections." << nl
                        << "triangle:" << f.tri(points)
                        << " d:" << d << abort(FatalError);
                }
                inters[interI].setHit();
                inters[interI].setPoint
                (
                    (d[fp0]*points[f[fp1]] - d[fp1]*points[f[fp0]])
                  / (d[fp0] - d[fp1])
                );
                inters[interI].elementType() = triPointRef::EDGE;
                inters[interI].setIndex(fEdges[fp0]);
                interI++;
            }
        }


        if (interI == 0)
        {
            // Return miss
        }
        else if (interI == 1)
        {
            // Only one intersection. Should not happen!
            cut = inters[0];
        }
        else if (interI == 2)
        {
            // Handle excludeEdgeI
            if
            (
                inters[0].elementType() == triPointRef::EDGE
             && inters[0].index() == excludeEdgeI
            )
            {
                cut = inters[1];
            }
            else if
            (
                inters[1].elementType() == triPointRef::EDGE
             && inters[1].index() == excludeEdgeI
            )
            {
                cut = inters[0];
            }
            else
            {
                // Two cuts. Find nearest.
                if
                (
                    magSqr(inters[0].rawPoint() - toPoint)
                  < magSqr(inters[1].rawPoint() - toPoint)
                )
                {
                    cut = inters[0];
                }
                else
                {
                    cut = inters[1];
                }
            }
        }
    }
    return cut;
}


// 'Snap' point on to endPoint.
void Foam::triSurfaceTools::snapToEnd
(
    const triSurface& s,
    const surfaceLocation& end,
    surfaceLocation& current
)
{
    if (end.elementType() == triPointRef::NONE)
    {
        if (current.elementType() == triPointRef::NONE)
        {
            // endpoint on triangle; current on triangle
            if (current.index() == end.index())
            {
                //if (debug)
                //{
                //    Pout<< "snapToEnd : snapping:" << current << " onto:"
                //        << end << endl;
                //}
                current = end;
                current.setHit();
            }
        }
        // No need to handle current on edge/point since tracking handles this.
    }
    else if (end.elementType() == triPointRef::EDGE)
    {
        if (current.elementType() == triPointRef::NONE)
        {
            // endpoint on edge; current on triangle
            const labelList& fEdges = s.faceEdges()[current.index()];

            if (findIndex(fEdges, end.index()) != -1)
            {
                //if (debug)
                //{
                //    Pout<< "snapToEnd : snapping:" << current << " onto:"
                //        << end << endl;
                //}
                current = end;
                current.setHit();
            }
        }
        else if (current.elementType() == triPointRef::EDGE)
        {
            // endpoint on edge; current on edge
            if (current.index() == end.index())
            {
                //if (debug)
                //{
                //    Pout<< "snapToEnd : snapping:" << current << " onto:"
                //        << end << endl;
                //}
                current = end;
                current.setHit();
            }
        }
        else
        {
            // endpoint on edge; current on point
            const edge& e = s.edges()[end.index()];

            if (current.index() == e[0] || current.index() == e[1])
            {
                //if (debug)
                //{
                //    Pout<< "snapToEnd : snapping:" << current << " onto:"
                //        << end << endl;
                //}
                current = end;
                current.setHit();
            }
        }
    }
    else    // end.elementType() == POINT
    {
        if (current.elementType() == triPointRef::NONE)
        {
            // endpoint on point; current on triangle
            const triSurface::FaceType& f = s.localFaces()[current.index()];

            if (findIndex(f, end.index()) != -1)
            {
                //if (debug)
                //{
                //    Pout<< "snapToEnd : snapping:" << current << " onto:"
                //        << end << endl;
                //}
                current = end;
                current.setHit();
            }
        }
        else if (current.elementType() == triPointRef::EDGE)
        {
            // endpoint on point; current on edge
            const edge& e = s.edges()[current.index()];

            if (end.index() == e[0] || end.index() == e[1])
            {
                //if (debug)
                //{
                //    Pout<< "snapToEnd : snapping:" << current << " onto:"
                //        << end << endl;
                //}
                current = end;
                current.setHit();
            }
        }
        else
        {
            // endpoint on point; current on point
            if (current.index() == end.index())
            {
                //if (debug)
                //{
                //    Pout<< "snapToEnd : snapping:" << current << " onto:"
                //        << end << endl;
                //}
                current = end;
                current.setHit();
            }
        }
    }
}


// Start:
//  - location
//  - element type (triangle/edge/point) and index
//  - triangle to exclude
Foam::surfaceLocation Foam::triSurfaceTools::visitFaces
(
    const triSurface& s,
    const labelList& eFaces,
    const surfaceLocation& start,
    const label excludeEdgeI,
    const label excludePointI,
    const surfaceLocation& end,
    const plane& cutPlane
)
{
    surfaceLocation nearest;

    scalar minDistSqr = Foam::sqr(GREAT);

    forAll(eFaces, i)
    {
        label triI = eFaces[i];

        // Make sure we don't revisit previous face
        if (triI != start.triangle())
        {
            if (end.elementType() == triPointRef::NONE && end.index() == triI)
            {
                // Endpoint is in this triangle. Jump there.
                nearest = end;
                nearest.setHit();
                nearest.triangle() = triI;
                break;
            }
            else
            {
               // Which edge is cut.

                surfaceLocation cutInfo = cutEdge
                (
                    s,
                    triI,
                    excludeEdgeI,       // excludeEdgeI
                    excludePointI,      // excludePointI
                    start.rawPoint(),
                    cutPlane,
                    end.rawPoint()
                );

                // If crossing an edge we expect next edge to be cut.
                if (excludeEdgeI != -1 && !cutInfo.hit())
                {
                    FatalErrorIn("triSurfaceTools::visitFaces(..)")
                        << "Triangle:" << triI
                        << " excludeEdge:" << excludeEdgeI
                        << " point:" << start.rawPoint()
                        << " plane:" << cutPlane
                        << " . No intersection!" << abort(FatalError);
                }

                if (cutInfo.hit())
                {
                    scalar distSqr = magSqr(cutInfo.rawPoint()-end.rawPoint());

                    if (distSqr < minDistSqr)
                    {
                        minDistSqr = distSqr;
                        nearest = cutInfo;
                        nearest.triangle() = triI;
                        nearest.setMiss();
                    }
                }
            }
        }
    }

    if (nearest.triangle() == -1)
    {
        // Did not move from edge. Give warning? Return something special?
        // For now responsability of caller to make sure that nothing has
        // moved.
    }

    return nearest;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// Write pointField to file
void Foam::triSurfaceTools::writeOBJ
(
    const fileName& fName,
    const pointField& pts
)
{
    OFstream outFile(fName);

    forAll(pts, pointI)
    {
        const point& pt = pts[pointI];

        outFile<< "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
    }
    Pout<< "Written " << pts.size() << " vertices to file " << fName << endl;
}


// Write vertex subset to OBJ format file
void Foam::triSurfaceTools::writeOBJ
(
    const triSurface& surf,
    const fileName& fName,
    const boolList& markedVerts
)
{
    OFstream outFile(fName);

    label nVerts = 0;
    forAll(markedVerts, vertI)
    {
        if (markedVerts[vertI])
        {
            const point& pt = surf.localPoints()[vertI];

            outFile<< "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;

            nVerts++;
        }
    }
    Pout<< "Written " << nVerts << " vertices to file " << fName << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Addressing helper functions:


// Get all triangles using vertices of edge
void Foam::triSurfaceTools::getVertexTriangles
(
    const triSurface& surf,
    const label edgeI,
    labelList& edgeTris
)
{
    const edge& e = surf.edges()[edgeI];
    const labelList& myFaces = surf.edgeFaces()[edgeI];

    label face1I = myFaces[0];
    label face2I = -1;
    if (myFaces.size() == 2)
    {
        face2I = myFaces[1];
    }

    const labelList& startFaces = surf.pointFaces()[e.start()];
    const labelList& endFaces = surf.pointFaces()[e.end()];

    // Number of triangles is sum of pointfaces - common faces
    // (= faces using edge)
    edgeTris.setSize(startFaces.size() + endFaces.size() - myFaces.size());

    label nTris = 0;
    forAll(startFaces, startFaceI)
    {
        edgeTris[nTris++] = startFaces[startFaceI];
    }

    forAll(endFaces, endFaceI)
    {
        label faceI = endFaces[endFaceI];

        if ((faceI != face1I) && (faceI != face2I))
        {
            edgeTris[nTris++] = faceI;
        }
    }
}


// Get all vertices connected to vertices of edge
Foam::labelList Foam::triSurfaceTools::getVertexVertices
(
    const triSurface& surf,
    const edge& e
)
{
    const edgeList& edges = surf.edges();

    label v1 = e.start();
    label v2 = e.end();

    // Get all vertices connected to v1 or v2 through an edge
    labelHashSet vertexNeighbours;

    const labelList& v1Edges = surf.pointEdges()[v1];

    forAll(v1Edges, v1EdgeI)
    {
        const edge& e = edges[v1Edges[v1EdgeI]];
        vertexNeighbours.insert(e.otherVertex(v1));
    }

    const labelList& v2Edges = surf.pointEdges()[v2];

    forAll(v2Edges, v2EdgeI)
    {
        const edge& e = edges[v2Edges[v2EdgeI]];

        label vertI = e.otherVertex(v2);

        vertexNeighbours.insert(vertI);
    }
    return vertexNeighbours.toc();
}


//// Order vertices consistent with face
//void Foam::triSurfaceTools::orderVertices
//(
//    const labelledTri& f,
//    const label v1,
//    const label v2,
//    label& vA,
//    label& vB
//)
//{
//    // Order v1, v2 in anticlockwise order.
//    bool reverse = false;
//
//    if (f[0] == v1)
//    {
//        if (f[1] != v2)
//        {
//            reverse = true;
//        }
//    }
//    else if (f[1] == v1)
//    {
//        if (f[2] != v2)
//        {
//            reverse = true;
//        }
//    }
//    else
//    {
//        if (f[0] != v2)
//        {
//            reverse = true;
//        }
//    }
//
//    if (reverse)
//    {
//        vA = v2;
//        vB = v1;
//    }
//    else
//    {
//        vA = v1;
//        vB = v2;
//    }
//}


// Get the other face using edgeI
Foam::label Foam::triSurfaceTools::otherFace
(
    const triSurface& surf,
    const label faceI,
    const label edgeI
)
{
    const labelList& myFaces = surf.edgeFaces()[edgeI];

    if (myFaces.size() != 2)
    {
        return -1;
    }
    else
    {
        if (faceI == myFaces[0])
        {
            return myFaces[1];
        }
        else
        {
            return myFaces[0];
        }
    }
}


// Get the two edges on faceI counterclockwise after edgeI
void Foam::triSurfaceTools::otherEdges
(
    const triSurface& surf,
    const label faceI,
    const label edgeI,
    label& e1,
    label& e2
)
{
    const labelList& eFaces = surf.faceEdges()[faceI];

    label i0 = findIndex(eFaces, edgeI);

    if (i0 == -1)
    {
        FatalErrorIn
        (
            "otherEdges"
            "(const triSurface&, const label, const label,"
            " label&, label&)"
        )   << "Edge " << surf.edges()[edgeI] << " not in face "
            << surf.localFaces()[faceI] << abort(FatalError);
    }

    label i1 = eFaces.fcIndex(i0);
    label i2 = eFaces.fcIndex(i1);

    e1 = eFaces[i1];
    e2 = eFaces[i2];
}


// Get the two vertices on faceI counterclockwise vertI
void Foam::triSurfaceTools::otherVertices
(
    const triSurface& surf,
    const label faceI,
    const label vertI,
    label& vert1I,
    label& vert2I
)
{
    const labelledTri& f = surf.localFaces()[faceI];

    if (vertI == f[0])
    {
        vert1I = f[1];
        vert2I = f[2];
    }
    else if (vertI == f[1])
    {
        vert1I = f[2];
        vert2I = f[0];
    }
    else if (vertI == f[2])
    {
        vert1I = f[0];
        vert2I = f[1];
    }
    else
    {
        FatalErrorIn
        (
            "otherVertices"
            "(const triSurface&, const label, const label,"
            " label&, label&)"
        )   << "Vertex " << vertI << " not in face " << f << abort(FatalError);
    }
}


// Get edge opposite vertex
Foam::label Foam::triSurfaceTools::oppositeEdge
(
    const triSurface& surf,
    const label faceI,
    const label vertI
)
{
    const labelList& myEdges = surf.faceEdges()[faceI];

    forAll(myEdges, myEdgeI)
    {
        label edgeI = myEdges[myEdgeI];

        const edge& e = surf.edges()[edgeI];

        if ((e.start() != vertI) && (e.end() != vertI))
        {
            return edgeI;
        }
    }

    FatalErrorIn
    (
        "oppositeEdge"
        "(const triSurface&, const label, const label)"
    )   << "Cannot find vertex " << vertI << " in edges of face " << faceI
        << abort(FatalError);

    return -1;
}


// Get vertex opposite edge
Foam::label Foam::triSurfaceTools::oppositeVertex
(
    const triSurface& surf,
    const label faceI,
    const label edgeI
)
{
    const triSurface::FaceType& f = surf.localFaces()[faceI];
    const edge& e = surf.edges()[edgeI];

    forAll(f, fp)
    {
        label vertI = f[fp];

        if (vertI != e.start() && vertI != e.end())
        {
            return vertI;
        }
    }

    FatalErrorIn("triSurfaceTools::oppositeVertex")
        << "Cannot find vertex opposite edge " << edgeI << " vertices " << e
        << " in face " << faceI << " vertices " << f << abort(FatalError);

    return -1;
}


// Returns edge label connecting v1, v2
Foam::label Foam::triSurfaceTools::getEdge
(
    const triSurface& surf,
    const label v1,
    const label v2
)
{
    const labelList& v1Edges = surf.pointEdges()[v1];

    forAll(v1Edges, v1EdgeI)
    {
        label edgeI = v1Edges[v1EdgeI];
        const edge& e = surf.edges()[edgeI];

        if ((e.start() == v2) || (e.end() == v2))
        {
            return edgeI;
        }
    }
    return -1;
}


// Return index of triangle (or -1) using all three edges
Foam::label Foam::triSurfaceTools::getTriangle
(
    const triSurface& surf,
    const label e0I,
    const label e1I,
    const label e2I
)
{
    if ((e0I == e1I) || (e0I == e2I) || (e1I == e2I))
    {
        FatalErrorIn
        (
            "getTriangle"
            "(const triSurface&, const label, const label,"
            " const label)"
        )   << "Duplicate edge labels : e0:" << e0I << " e1:" << e1I
            << " e2:" << e2I
            << abort(FatalError);
    }

    const labelList& eFaces = surf.edgeFaces()[e0I];

    forAll(eFaces, eFaceI)
    {
        label faceI = eFaces[eFaceI];

        const labelList& myEdges = surf.faceEdges()[faceI];

        if
        (
            (myEdges[0] == e1I)
         || (myEdges[1] == e1I)
         || (myEdges[2] == e1I)
        )
        {
            if
            (
                (myEdges[0] == e2I)
             || (myEdges[1] == e2I)
             || (myEdges[2] == e2I)
            )
            {
                return faceI;
            }
        }
    }
    return -1;
}


// Collapse indicated edges. Return new tri surface.
Foam::triSurface Foam::triSurfaceTools::collapseEdges
(
    const triSurface& surf,
    const labelList& collapsableEdges
)
{
    pointField edgeMids(surf.nEdges());

    forAll(edgeMids, edgeI)
    {
        const edge& e = surf.edges()[edgeI];

        edgeMids[edgeI] =
            0.5
          * (
                surf.localPoints()[e.start()]
              + surf.localPoints()[e.end()]
            );
    }


    labelList faceStatus(surf.size(), ANYEDGE);

    //// Protect triangles which are on the border of different regions
    //forAll(edges, edgeI)
    //{
    //    const labelList& neighbours = edgeFaces[edgeI];
    //
    //    if ((neighbours.size() != 2) && (neighbours.size() != 1))
    //    {
    //        FatalErrorIn("collapseEdges")
    //            << abort(FatalError);
    //    }
    //
    //    if (neighbours.size() == 2)
    //    {
    //        if (surf[neighbours[0]].region() != surf[neighbours[1]].region())
    //        {
    //            // Neighbours on different regions. For now dont allow
    //            // any collapse.
    //            //Pout<< "protecting face " << neighbours[0]
    //            //    << ' ' << neighbours[1] << endl;
    //            faceStatus[neighbours[0]] = NOEDGE;
    //            faceStatus[neighbours[1]] = NOEDGE;
    //        }
    //    }
    //}

    return collapseEdges(surf, collapsableEdges, edgeMids, faceStatus);
}


// Collapse indicated edges. Return new tri surface.
Foam::triSurface Foam::triSurfaceTools::collapseEdges
(
    const triSurface& surf,
    const labelList& collapseEdgeLabels,
    const pointField& edgeMids,
    labelList& faceStatus
)
{
    const labelListList& edgeFaces = surf.edgeFaces();
    const pointField& localPoints = surf.localPoints();
    const edgeList& edges = surf.edges();

    // Storage for new points
    pointField newPoints(localPoints);

    // Map for old to new points
    labelList pointMap(localPoints.size());
    forAll(localPoints, pointI)
    {
        pointMap[pointI] = pointI;
    }


    // Do actual 'collapsing' of edges

    forAll(collapseEdgeLabels, collapseEdgeI)
    {
        const label edgeI = collapseEdgeLabels[collapseEdgeI];

        if ((edgeI < 0) || (edgeI >= surf.nEdges()))
        {
            FatalErrorIn("collapseEdges")
                << "Edge label outside valid range." << endl
                << "edge label:" << edgeI << endl
                << "total number of edges:" << surf.nEdges() << endl
                << abort(FatalError);
        }

        const labelList& neighbours = edgeFaces[edgeI];

        if (neighbours.size() == 2)
        {
            const label stat0 = faceStatus[neighbours[0]];
            const label stat1 = faceStatus[neighbours[1]];

            // Check faceStatus to make sure this one can be collapsed
            if
            (
                ((stat0 == ANYEDGE) || (stat0 == edgeI))
             && ((stat1 == ANYEDGE) || (stat1 == edgeI))
            )
            {
                const edge& e = edges[edgeI];

                // Set up mapping to 'collapse' points of edge
                if
                (
                    (pointMap[e.start()] != e.start())
                 || (pointMap[e.end()] != e.end())
                )
                {
                    FatalErrorIn("collapseEdges")
                        << "points already mapped. Double collapse." << endl
                        << "edgeI:" << edgeI
                        << "  start:" << e.start()
                        << "  end:" << e.end()
                        << "  pointMap[start]:" << pointMap[e.start()]
                        << "  pointMap[end]:" << pointMap[e.end()]
                        << abort(FatalError);
                }

                const label minVert = min(e.start(), e.end());
                pointMap[e.start()] = minVert;
                pointMap[e.end()] = minVert;

                // Move shared vertex to mid of edge
                newPoints[minVert] = edgeMids[edgeI];

                // Protect neighbouring faces
                protectNeighbours(surf, e.start(), faceStatus);
                protectNeighbours(surf, e.end(), faceStatus);
                protectNeighbours
                (
                    surf,
                    oppositeVertex(surf, neighbours[0], edgeI),
                    faceStatus
                );
                protectNeighbours
                (
                    surf,
                    oppositeVertex(surf, neighbours[1], edgeI),
                    faceStatus
                );

                // Mark all collapsing faces
                labelList collapseFaces =
                    getCollapsedFaces
                    (
                        surf,
                        edgeI
                    ).toc();

                forAll(collapseFaces, collapseI)
                {
                     faceStatus[collapseFaces[collapseI]] = COLLAPSED;
                }
            }
        }
    }


    // Storage for new triangles
    List<labelledTri> newTris(surf.size());
    label newTriI = 0;

    const List<labelledTri>& localFaces = surf.localFaces();


    // Get only non-collapsed triangles and renumber vertex labels.
    forAll(localFaces, faceI)
    {
        const labelledTri& f = localFaces[faceI];

        const label a = pointMap[f[0]];
        const label b = pointMap[f[1]];
        const label c = pointMap[f[2]];

        if
        (
            (a != b) && (a != c) && (b != c)
         && (faceStatus[faceI] != COLLAPSED)
        )
        {
            // uncollapsed triangle
            newTris[newTriI++] = labelledTri(a, b, c, f.region());
        }
        else
        {
            //Pout<< "Collapsed triangle " << faceI
            //    << " vertices:" << f << endl;
        }
    }
    newTris.setSize(newTriI);



    // Pack faces

    triSurface tempSurf(newTris, surf.patches(), newPoints);

    return
        triSurface
        (
            tempSurf.localFaces(),
            tempSurf.patches(),
            tempSurf.localPoints()
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::triSurface Foam::triSurfaceTools::redGreenRefine
(
    const triSurface& surf,
    const labelList& refineFaces
)
{
    List<refineType> refineStatus(surf.size(), NONE);

    // Mark & propagate refinement
    forAll(refineFaces, refineFaceI)
    {
        calcRefineStatus(surf, refineFaces[refineFaceI], refineStatus);
    }

    // Do actual refinement
    return doRefine(surf, refineStatus);
}


Foam::triSurface Foam::triSurfaceTools::greenRefine
(
    const triSurface& surf,
    const labelList& refineEdges
)
{
    // Storage for marking faces
    List<refineType> refineStatus(surf.size(), NONE);

    // Storage for new faces
    DynamicList<labelledTri> newFaces(0);

    pointField newPoints(surf.localPoints());
    newPoints.setSize(surf.nPoints() + surf.nEdges());
    label newPointI = surf.nPoints();


    // Refine edges
    forAll(refineEdges, refineEdgeI)
    {
        label edgeI = refineEdges[refineEdgeI];

        const labelList& myFaces = surf.edgeFaces()[edgeI];

        bool neighbourIsRefined= false;

        forAll(myFaces, myFaceI)
        {
            if (refineStatus[myFaces[myFaceI]] != NONE)
            {
                neighbourIsRefined =  true;
            }
        }

        // Only refine if none of the faces is refined
        if (!neighbourIsRefined)
        {
            // Refine edge
            const edge& e = surf.edges()[edgeI];

            point mid =
                0.5
              * (
                    surf.localPoints()[e.start()]
                  + surf.localPoints()[e.end()]
                );

            newPoints[newPointI] = mid;

            // Refine faces using edge
            forAll(myFaces, myFaceI)
            {
                // Add faces to newFaces
                greenRefine
                (
                    surf,
                    myFaces[myFaceI],
                    edgeI,
                    newPointI,
                    newFaces
                );

                // Mark as refined
                refineStatus[myFaces[myFaceI]] = GREEN;
            }

            newPointI++;
        }
    }

    // Add unrefined faces
    forAll(surf.localFaces(), faceI)
    {
        if (refineStatus[faceI] == NONE)
        {
            newFaces.append(surf.localFaces()[faceI]);
        }
    }

    newFaces.shrink();
    newPoints.setSize(newPointI);

    return triSurface(newFaces, surf.patches(), newPoints, true);
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Geometric helper functions:


// Returns element in edgeIndices with minimum length
Foam::label Foam::triSurfaceTools::minEdge
(
    const triSurface& surf,
    const labelList& edgeIndices
)
{
    scalar minLength = GREAT;
    label minIndex = -1;
    forAll(edgeIndices, i)
    {
        const edge& e = surf.edges()[edgeIndices[i]];

        scalar length =
            mag
            (
                surf.localPoints()[e.end()]
              - surf.localPoints()[e.start()]
            );

        if (length < minLength)
        {
            minLength = length;
            minIndex = i;
        }
    }
    return edgeIndices[minIndex];
}


// Returns element in edgeIndices with maximum length
Foam::label Foam::triSurfaceTools::maxEdge
(
    const triSurface& surf,
    const labelList& edgeIndices
)
{
    scalar maxLength = -GREAT;
    label maxIndex = -1;
    forAll(edgeIndices, i)
    {
        const edge& e = surf.edges()[edgeIndices[i]];

        scalar length =
            mag
            (
                surf.localPoints()[e.end()]
              - surf.localPoints()[e.start()]
            );

        if (length > maxLength)
        {
            maxLength = length;
            maxIndex = i;
        }
    }
    return edgeIndices[maxIndex];
}


// Merge points and reconstruct surface
Foam::triSurface Foam::triSurfaceTools::mergePoints
(
    const triSurface& surf,
    const scalar mergeTol
)
{
    pointField newPoints(surf.nPoints());

    labelList pointMap(surf.nPoints());

    bool hasMerged = Foam::mergePoints
    (
        surf.localPoints(),
        mergeTol,
        false,
        pointMap,
        newPoints
    );

    if (hasMerged)
    {
        // Pack the triangles

        // Storage for new triangles
        List<labelledTri> newTriangles(surf.size());
        label newTriangleI = 0;

        forAll(surf, faceI)
        {
            const labelledTri& f = surf.localFaces()[faceI];

            label newA = pointMap[f[0]];
            label newB = pointMap[f[1]];
            label newC = pointMap[f[2]];

            if ((newA != newB) && (newA != newC) && (newB != newC))
            {
                newTriangles[newTriangleI++] =
                    labelledTri(newA, newB, newC, f.region());
            }
        }
        newTriangles.setSize(newTriangleI);

        return triSurface
        (
            newTriangles,
            surf.patches(),
            newPoints,
            true                //reuse storage
        );
    }
    else
    {
        return surf;
    }
}


// Calculate normal on triangle
Foam::vector Foam::triSurfaceTools::surfaceNormal
(
    const triSurface& surf,
    const label nearestFaceI,
    const point& nearestPt
)
{
    const triSurface::FaceType& f = surf[nearestFaceI];
    const pointField& points = surf.points();

    label nearType, nearLabel;

    f.nearestPointClassify(nearestPt, points, nearType, nearLabel);

    if (nearType == triPointRef::NONE)
    {
        // Nearest to face
        return surf.faceNormals()[nearestFaceI];
    }
    else if (nearType == triPointRef::EDGE)
    {
        // Nearest to edge. Assume order of faceEdges same as face vertices.
        label edgeI = surf.faceEdges()[nearestFaceI][nearLabel];

        // Calculate edge normal by averaging face normals
        const labelList& eFaces = surf.edgeFaces()[edgeI];

        vector edgeNormal(vector::zero);

        forAll(eFaces, i)
        {
            edgeNormal += surf.faceNormals()[eFaces[i]];
        }
        return edgeNormal/(mag(edgeNormal) + VSMALL);
    }
    else
    {
        // Nearest to point
        const triSurface::FaceType& localF = surf.localFaces()[nearestFaceI];
        return surf.pointNormals()[localF[nearLabel]];
    }
}


Foam::triSurfaceTools::sideType Foam::triSurfaceTools::edgeSide
(
    const triSurface& surf,
    const point& sample,
    const point& nearestPoint,
    const label edgeI
)
{
    const labelList& eFaces = surf.edgeFaces()[edgeI];

    if (eFaces.size() != 2)
    {
        // Surface not closed.
        return UNKNOWN;
    }
    else
    {
        const vectorField& faceNormals = surf.faceNormals();

        // Compare to bisector. This is actually correct since edge is
        // nearest so there is a knife-edge.

        vector n = 0.5*(faceNormals[eFaces[0]] + faceNormals[eFaces[1]]);

        if (((sample - nearestPoint) & n) > 0)
        {
            return OUTSIDE;
        }
        else
        {
            return INSIDE;
        }
    }
}


// Calculate normal on triangle
Foam::triSurfaceTools::sideType Foam::triSurfaceTools::surfaceSide
(
    const triSurface& surf,
    const point& sample,
    const label nearestFaceI
)
{
    const triSurface::FaceType& f = surf[nearestFaceI];
    const pointField& points = surf.points();

    // Find where point is on face
    label nearType, nearLabel;

    pointHit pHit = f.nearestPointClassify(sample, points, nearType, nearLabel);

    const point& nearestPoint(pHit.rawPoint());

    if (nearType == triPointRef::NONE)
    {
        vector sampleNearestVec = (sample - nearestPoint);

        // Nearest to face interior. Use faceNormal to determine side
        scalar c = sampleNearestVec & surf.faceNormals()[nearestFaceI];

        // // If the sample is essentially on the face, do not check for
        // // it being perpendicular.

        // scalar magSampleNearestVec = mag(sampleNearestVec);

        // if (magSampleNearestVec > SMALL)
        // {
        //     c /= magSampleNearestVec*mag(surf.faceNormals()[nearestFaceI]);

        //     if (mag(c) < 0.99)
        //     {
        //         FatalErrorIn("triSurfaceTools::surfaceSide")
        //             << "nearestPoint identified as being on triangle face "
        //             << "but vector from nearestPoint to sample is not "
        //             << "perpendicular to the normal." << nl
        //             << "sample: " << sample << nl
        //             << "nearestPoint: " << nearestPoint << nl
        //             << "sample - nearestPoint: "
        //             << sample - nearestPoint << nl
        //             << "normal: " << surf.faceNormals()[nearestFaceI] << nl
        //             << "mag(sample - nearestPoint): "
        //             << mag(sample - nearestPoint) << nl
        //             << "normalised dot product: " << c << nl
        //             << "triangle vertices: " << nl
        //             << "    " << points[f[0]] << nl
        //             << "    " << points[f[1]] << nl
        //             << "    " << points[f[2]] << nl
        //             << abort(FatalError);
        //     }
        // }

        if (c > 0)
        {
            return OUTSIDE;
        }
        else
        {
            return INSIDE;
        }
    }
    else if (nearType == triPointRef::EDGE)
    {
        // Nearest to edge nearLabel. Note that this can only be a knife-edge
        // situation since otherwise the nearest point could never be the edge.

        // Get the edge. Assume order of faceEdges same as face vertices.
        label edgeI = surf.faceEdges()[nearestFaceI][nearLabel];

        // if (debug)
        // {
        //    // Check order of faceEdges same as face vertices.
        //    const edge& e = surf.edges()[edgeI];
        //    const labelList& meshPoints = surf.meshPoints();
        //    const edge meshEdge(meshPoints[e[0]], meshPoints[e[1]]);

        //    if
        //    (
        //        meshEdge
        //     != edge(f[nearLabel], f[f.fcIndex(nearLabel)])
        //    )
        //    {
        //        FatalErrorIn("triSurfaceTools::surfaceSide")
        //            << "Edge:" << edgeI << " local vertices:" << e
        //            << " mesh vertices:" << meshEdge
        //            << " not at position " << nearLabel
        //            << " in face " << f
        //            << abort(FatalError);
        //    }
        // }

        return edgeSide(surf, sample, nearestPoint, edgeI);
    }
    else
    {
        // Nearest to point. Could use pointNormal here but is not correct.
        // Instead determine which edge using point is nearest and use test
        // above (nearType == triPointRef::EDGE).


        const triSurface::FaceType& localF = surf.localFaces()[nearestFaceI];
        label nearPointI = localF[nearLabel];

        const edgeList& edges = surf.edges();
        const pointField& localPoints = surf.localPoints();
        const point& base = localPoints[nearPointI];

        const labelList& pEdges = surf.pointEdges()[nearPointI];

        scalar minDistSqr = Foam::sqr(GREAT);
        label minEdgeI = -1;

        forAll(pEdges, i)
        {
            label edgeI = pEdges[i];

            const edge& e = edges[edgeI];

            label otherPointI = e.otherVertex(nearPointI);

            // Get edge normal.
            vector eVec(localPoints[otherPointI] - base);
            scalar magEVec = mag(eVec);

            if (magEVec > VSMALL)
            {
                eVec /= magEVec;

                // Get point along vector and determine closest.
                const point perturbPoint = base + eVec;

                scalar distSqr = Foam::magSqr(sample - perturbPoint);

                if (distSqr < minDistSqr)
                {
                    minDistSqr = distSqr;
                    minEdgeI = edgeI;
                }
            }
        }

        if (minEdgeI == -1)
        {
            FatalErrorIn("treeDataTriSurface::getSide")
                << "Problem: did not find edge closer than " << minDistSqr
                << abort(FatalError);
        }

        return edgeSide(surf, sample, nearestPoint, minEdgeI);
    }
}


// triangulation of boundaryMesh
Foam::triSurface Foam::triSurfaceTools::triangulate
(
    const polyBoundaryMesh& bMesh,
    const labelHashSet& includePatches,
    const bool verbose
)
{
    const polyMesh& mesh = bMesh.mesh();

    // Storage for surfaceMesh. Size estimate.
    DynamicList<labelledTri> triangles
    (
        mesh.nFaces() - mesh.nInternalFaces()
    );

    label newPatchI = 0;

    forAllConstIter(labelHashSet, includePatches, iter)
    {
        const label patchI = iter.key();
        const polyPatch& patch = bMesh[patchI];
        const pointField& points = patch.points();

        label nTriTotal = 0;

        forAll(patch, patchFaceI)
        {
            const face& f = patch[patchFaceI];

            faceList triFaces(f.nTriangles(points));

            label nTri = 0;

            f.triangles(points, nTri, triFaces);

            forAll(triFaces, triFaceI)
            {
                const face& f = triFaces[triFaceI];

                triangles.append(labelledTri(f[0], f[1], f[2], newPatchI));

                nTriTotal++;
            }
        }

        if (verbose)
        {
            Pout<< patch.name() << " : generated " << nTriTotal
                << " triangles from " << patch.size() << " faces with"
                << " new patchid " << newPatchI << endl;
        }

        newPatchI++;
    }
    triangles.shrink();

    // Create globally numbered tri surface
    triSurface rawSurface(triangles, mesh.points());

    // Create locally numbered tri surface
    triSurface surface
    (
        rawSurface.localFaces(),
        rawSurface.localPoints()
    );

    // Add patch names to surface
    surface.patches().setSize(newPatchI);

    newPatchI = 0;

    forAllConstIter(labelHashSet, includePatches, iter)
    {
        const label patchI = iter.key();
        const polyPatch& patch = bMesh[patchI];

        surface.patches()[newPatchI].name() = patch.name();
        surface.patches()[newPatchI].geometricType() = patch.type();

        newPatchI++;
    }

    return surface;
}


// triangulation of boundaryMesh
Foam::triSurface Foam::triSurfaceTools::triangulateFaceCentre
(
    const polyBoundaryMesh& bMesh,
    const labelHashSet& includePatches,
    const bool verbose
)
{
    const polyMesh& mesh = bMesh.mesh();

    // Storage for new points = meshpoints + face centres.
    const pointField& points = mesh.points();
    const pointField& faceCentres = mesh.faceCentres();

    pointField newPoints(points.size() + faceCentres.size());

    label newPointI = 0;

    forAll(points, pointI)
    {
        newPoints[newPointI++] = points[pointI];
    }
    forAll(faceCentres, faceI)
    {
        newPoints[newPointI++] = faceCentres[faceI];
    }


    // Count number of faces.
    DynamicList<labelledTri> triangles
    (
        mesh.nFaces() - mesh.nInternalFaces()
    );

    label newPatchI = 0;

    forAllConstIter(labelHashSet, includePatches, iter)
    {
        const label patchI = iter.key();
        const polyPatch& patch = bMesh[patchI];

        label nTriTotal = 0;

        forAll(patch, patchFaceI)
        {
            // Face in global coords.
            const face& f = patch[patchFaceI];

            // Index in newPointI of face centre.
            label fc = points.size() + patchFaceI + patch.start();

            forAll(f, fp)
            {
                label fp1 = f.fcIndex(fp);

                triangles.append(labelledTri(f[fp], f[fp1], fc, newPatchI));

                nTriTotal++;
            }
        }

        if (verbose)
        {
            Pout<< patch.name() << " : generated " << nTriTotal
                << " triangles from " << patch.size() << " faces with"
                << " new patchid " << newPatchI << endl;
        }

        newPatchI++;
    }
    triangles.shrink();


    // Create globally numbered tri surface
    triSurface rawSurface(triangles, newPoints);

    // Create locally numbered tri surface
    triSurface surface
    (
        rawSurface.localFaces(),
        rawSurface.localPoints()
    );

    // Add patch names to surface
    surface.patches().setSize(newPatchI);

    newPatchI = 0;

    forAllConstIter(labelHashSet, includePatches, iter)
    {
        const label patchI = iter.key();
        const polyPatch& patch = bMesh[patchI];

        surface.patches()[newPatchI].name() = patch.name();
        surface.patches()[newPatchI].geometricType() = patch.type();

        newPatchI++;
    }

    return surface;
}


Foam::triSurface Foam::triSurfaceTools::delaunay2D(const List<vector2D>& pts)
{
    // Vertices in geompack notation. Note that could probably just use
    // pts.begin() if double precision.
    List<doubleScalar> geompackVertices(2*pts.size());
    label doubleI = 0;
    forAll(pts, i)
    {
        geompackVertices[doubleI++] = pts[i][0];
        geompackVertices[doubleI++] = pts[i][1];
    }

    // Storage for triangles
    int m2 = 3;
    List<int> triangle_node(m2*3*pts.size());
    List<int> triangle_neighbor(m2*3*pts.size());

    // Triangulate
    int nTris = 0;
    int err = dtris2
    (
        pts.size(),
        geompackVertices.begin(),
        &nTris,
        triangle_node.begin(),
        triangle_neighbor.begin()
    );

    if (err != 0)
    {
        FatalErrorIn("triSurfaceTools::delaunay2D(const List<vector2D>&)")
            << "Failed dtris2 with vertices:" << pts.size()
            << abort(FatalError);
    }

    // Trim
    triangle_node.setSize(3*nTris);
    triangle_neighbor.setSize(3*nTris);

    // Convert to triSurface.
    List<labelledTri> faces(nTris);

    forAll(faces, i)
    {
        faces[i] = labelledTri
        (
            triangle_node[3*i]-1,
            triangle_node[3*i+1]-1,
            triangle_node[3*i+2]-1,
            0
        );
    }

    pointField points(pts.size());
    forAll(pts, i)
    {
        points[i][0] = pts[i][0];
        points[i][1] = pts[i][1];
        points[i][2] = 0.0;
    }

    return triSurface(faces, points);
}


void Foam::triSurfaceTools::calcInterpolationWeights
(
    const triPointRef& tri,
    const point& p,
    FixedList<scalar, 3>& weights
)
{
    // calculate triangle edge vectors and triangle face normal
    // the 'i':th edge is opposite node i
    FixedList<vector, 3> edge;
    edge[0] = tri.c()-tri.b();
    edge[1] = tri.a()-tri.c();
    edge[2] = tri.b()-tri.a();

    vector triangleFaceNormal = edge[1] ^ edge[2];

    // calculate edge normal (pointing inwards)
    FixedList<vector, 3> normal;
    for (label i=0; i<3; i++)
    {
        normal[i] = triangleFaceNormal ^ edge[i];
        normal[i] /= mag(normal[i]) + VSMALL;
    }

    weights[0] = ((p-tri.b()) & normal[0]) / max(VSMALL, normal[0] & edge[1]);
    weights[1] = ((p-tri.c()) & normal[1]) / max(VSMALL, normal[1] & edge[2]);
    weights[2] = ((p-tri.a()) & normal[2]) / max(VSMALL, normal[2] & edge[0]);
}


// Calculate weighting factors from samplePts to triangle it is in.
// Uses linear search.
void Foam::triSurfaceTools::calcInterpolationWeights
(
    const triSurface& s,
    const pointField& samplePts,
    List<FixedList<label, 3> >& allVerts,
    List<FixedList<scalar, 3> >& allWeights
)
{
    allVerts.setSize(samplePts.size());
    allWeights.setSize(samplePts.size());

    const pointField& points = s.points();

    forAll(samplePts, i)
    {
        const point& samplePt = samplePts[i];


        FixedList<label, 3>& verts = allVerts[i];
        FixedList<scalar, 3>& weights = allWeights[i];

        scalar minDistance = GREAT;

        forAll(s, faceI)
        {
            const labelledTri& f = s[faceI];

            triPointRef tri(f.tri(points));

            label nearType, nearLabel;

            pointHit nearest = tri.nearestPointClassify
            (
                samplePt,
                nearType,
                nearLabel
            );

            if (nearest.hit())
            {
                // samplePt inside triangle
                verts[0] = f[0];
                verts[1] = f[1];
                verts[2] = f[2];

                calcInterpolationWeights(tri, nearest.rawPoint(), weights);

                //Pout<< "calcScalingFactors : samplePt:" << samplePt
                //    << " inside triangle:" << faceI
                //    << " verts:" << verts
                //    << " weights:" << weights
                //    << endl;

                break;
            }
            else if (nearest.distance() < minDistance)
            {
                minDistance = nearest.distance();

                // Outside triangle. Store nearest.

                if (nearType == triPointRef::POINT)
                {
                    verts[0] = f[nearLabel];
                    weights[0] = 1;
                    verts[1] = -1;
                    weights[1] = -GREAT;
                    verts[2] = -1;
                    weights[2] = -GREAT;

                    //Pout<< "calcScalingFactors : samplePt:" << samplePt
                    //    << " distance:" << nearest.distance()
                    //    << " from point:" << points[f[nearLabel]]
                    //    << endl;
                }
                else if (nearType == triPointRef::EDGE)
                {
                    verts[0] = f[nearLabel];
                    verts[1] = f[f.fcIndex(nearLabel)];
                    verts[2] = -1;

                    const point& p0 = points[verts[0]];
                    const point& p1 = points[verts[1]];

                    scalar s = min
                    (
                        1,
                        max
                        (
                            0,
                            mag(nearest.rawPoint() - p0)/mag(p1 - p0)
                        )
                    );

                    // Interpolate
                    weights[0] = 1 - s;
                    weights[1] = s;
                    weights[2] = -GREAT;

                    //Pout<< "calcScalingFactors : samplePt:" << samplePt
                    //    << " distance:" << nearest.distance()
                    //    << " from edge:" << p0 << p1 << " s:" << s
                    //    << endl;
                }
                else
                {
                    // triangle. Can only happen because of truncation errors.
                    verts[0] = f[0];
                    verts[1] = f[1];
                    verts[2] = f[2];

                    calcInterpolationWeights(tri, nearest.rawPoint(), weights);

                    //Pout<< "calcScalingFactors : samplePt:" << samplePt
                    //    << " distance:" << nearest.distance()
                    //    << " to verts:" << verts
                    //    << " weights:" << weights
                    //    << endl;

                    break;
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Tracking:


// Test point on surface to see if is on face,edge or point.
Foam::surfaceLocation Foam::triSurfaceTools::classify
(
    const triSurface& s,
    const label triI,
    const point& trianglePoint
)
{
    surfaceLocation nearest;

    // Nearest point could be on point or edge. Retest.
    label index, elemType;
    //bool inside =
    triPointRef(s[triI].tri(s.points())).classify
    (
        trianglePoint,
        elemType,
        index
    );

    nearest.setPoint(trianglePoint);

    if (elemType == triPointRef::NONE)
    {
        nearest.setHit();
        nearest.setIndex(triI);
        nearest.elementType() = triPointRef::NONE;
    }
    else if (elemType == triPointRef::EDGE)
    {
        nearest.setMiss();
        nearest.setIndex(s.faceEdges()[triI][index]);
        nearest.elementType() = triPointRef::EDGE;
    }
    else //if (elemType == triPointRef::POINT)
    {
        nearest.setMiss();
        nearest.setIndex(s.localFaces()[triI][index]);
        nearest.elementType() = triPointRef::POINT;
    }

    return nearest;
}


Foam::surfaceLocation Foam::triSurfaceTools::trackToEdge
(
    const triSurface& s,
    const surfaceLocation& start,
    const surfaceLocation& end,
    const plane& cutPlane
)
{
    // Start off from starting point
    surfaceLocation nearest = start;
    nearest.setMiss();

    // See if in same triangle as endpoint. If so snap.
    snapToEnd(s, end, nearest);

    if (!nearest.hit())
    {
        // Not yet at end point

        if (start.elementType() == triPointRef::NONE)
        {
            // Start point is inside triangle. Trivial cases already handled
            // above.

            // end point is on edge or point so cross currrent triangle to
            // see which edge is cut.

            nearest = cutEdge
            (
                s,
                start.index(),          // triangle
                -1,                     // excludeEdge
                -1,                     // excludePoint
                start.rawPoint(),
                cutPlane,
                end.rawPoint()
            );
            nearest.elementType() = triPointRef::EDGE;
            nearest.triangle() = start.index();
            nearest.setMiss();
        }
        else if (start.elementType() == triPointRef::EDGE)
        {
            // Pick connected triangle that is most in direction.
            const labelList& eFaces = s.edgeFaces()[start.index()];

            nearest = visitFaces
            (
                s,
                eFaces,
                start,
                start.index(),      // excludeEdgeI
                -1,                 // excludePointI
                end,
                cutPlane
            );
        }
        else    // start.elementType() == triPointRef::POINT
        {
            const labelList& pFaces = s.pointFaces()[start.index()];

            nearest = visitFaces
            (
                s,
                pFaces,
                start,
                -1,                 // excludeEdgeI
                start.index(),      // excludePointI
                end,
                cutPlane
            );
        }
        snapToEnd(s, end, nearest);
    }
    return nearest;
}


void Foam::triSurfaceTools::track
(
    const triSurface& s,
    const surfaceLocation& endInfo,
    const plane& cutPlane,
    surfaceLocation& hitInfo
)
{
    //OFstream str("track.obj");
    //label vertI = 0;
    //meshTools::writeOBJ(str, hitInfo.rawPoint());
    //vertI++;

    // Track across surface.
    while (true)
    {
        //Pout<< "Tracking from:" << nl
        //    << "    " << hitInfo.info()
        //    << endl;

        hitInfo = trackToEdge
        (
            s,
            hitInfo,
            endInfo,
            cutPlane
        );

        //meshTools::writeOBJ(str, hitInfo.rawPoint());
        //vertI++;
        //str<< "l " << vertI-1 << ' ' << vertI << nl;

        //Pout<< "Tracked to:" << nl
        //    << "    " << hitInfo.info() << endl;

        if (hitInfo.hit() || hitInfo.triangle() == -1)
        {
            break;
        }
    }
}


// ************************************************************************* //
