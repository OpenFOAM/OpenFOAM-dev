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

\*---------------------------------------------------------------------------*/

#include "collapseBase.H"
#include "triSurfaceTools.H"
#include "argList.H"
#include "OFstream.H"
#include "SubList.H"
#include "labelPair.H"
#include "meshTools.H"
#include "OSspecific.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//// Dump collapse region to .obj file
//static void writeRegionOBJ
//(
//    const triSurface& surf,
//    const label regionI,
//    const labelList& collapseRegion,
//    const labelList& outsideVerts
//)
//{
//    fileName dir("regions");
//
//    mkDir(dir);
//    fileName regionName(dir / "region_" + name(regionI) + ".obj");
//
//    Pout<< "Dumping region " << regionI << " to file " << regionName << endl;
//
//    boolList include(surf.size(), false);
//
//    forAll(collapseRegion, facei)
//    {
//        if (collapseRegion[facei] == regionI)
//        {
//            include[facei] = true;
//        }
//    }
//
//    labelList pointMap, faceMap;
//
//    triSurface regionSurf(surf.subsetMesh(include, pointMap, faceMap));
//
//    Pout<< "Region " << regionI << " surface:" << nl;
//    regionSurf.writeStats(Pout);
//
//    regionSurf.write(regionName);
//
//
//    // Dump corresponding outside vertices.
//    fileName pointsName(dir / "regionPoints_" + name(regionI) + ".obj");
//
//    Pout<< "Dumping region " << regionI << " points to file " << pointsName
//        << endl;
//
//    OFstream str(pointsName);
//
//    forAll(outsideVerts, i)
//    {
//        meshTools::writeOBJ(str, surf.localPoints()[outsideVerts[i]]);
//    }
//}


// Split triangle into multiple triangles because edge e being split
// into multiple edges.
static void splitTri
(
    const labelledTri& f,
    const edge& e,
    const labelList& splitPoints,
    DynamicList<labelledTri>& tris
)
{
    //label oldNTris = tris.size();

    label fp = findIndex(f, e[0]);
    label fp1 = f.fcIndex(fp);
    label fp2 = f.fcIndex(fp1);

    if (f[fp1] == e[1])
    {
        // Split triangle along fp to fp1
        tris.append(labelledTri(f[fp2], f[fp], splitPoints[0], f.region()));

        for (label i = 1; i < splitPoints.size(); i++)
        {
            tris.append
            (
                labelledTri
                (
                    f[fp2],
                    splitPoints[i-1],
                    splitPoints[i],
                    f.region()
                )
            );
        }

        tris.append
        (
            labelledTri
            (
                f[fp2],
                splitPoints.last(),
                f[fp1],
                f.region()
            )
        );
    }
    else if (f[fp2] == e[1])
    {
        // Split triangle along fp2 to fp. (Reverse order of splitPoints)

        tris.append
        (
            labelledTri
            (
                f[fp1],
                f[fp2],
                splitPoints.last(),
                f.region()
            )
        );

        for (label i = splitPoints.size()-1; i > 0; --i)
        {
            tris.append
            (
                labelledTri
                (
                    f[fp1],
                    splitPoints[i],
                    splitPoints[i-1],
                    f.region()
                )
            );
        }

        tris.append
        (
            labelledTri
            (
                f[fp1],
                splitPoints[0],
                f[fp],
                f.region()
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "Edge " << e << " not part of triangle " << f
            << " fp:" << fp
            << " fp1:" << fp1
            << " fp2:" << fp2
            << abort(FatalError);
    }

    //Pout<< "Split face " << f << " along edge " << e
    //    << " into triangles:" << endl;
    //
    //for (label i = oldNTris; i < tris.size(); i++)
    //{
    //    Pout<< "   " << tris[i] << nl;
    //}
}


// Insert scalar into sortedVerts/sortedWeights so the weights are in
// incrementing order.
static bool insertSorted
(
    const label vertI,
    const scalar weight,

    labelList& sortedVerts,
    scalarField& sortedWeights
)
{
    if (findIndex(sortedVerts, vertI) != -1)
    {
        FatalErrorInFunction
            << " which is already in list of sorted vertices "
            << sortedVerts << abort(FatalError);
    }

    if (weight <= 0 || weight >= 1)
    {
        FatalErrorInFunction
            << " with illegal weight " << weight
            << " into list of sorted vertices "
            << sortedVerts << abort(FatalError);
    }


    label insertI = sortedVerts.size();

    forAll(sortedVerts, sortedI)
    {
        scalar w = sortedWeights[sortedI];

        if (mag(w - weight) < small)
        {
            WarningInFunction
                << "Trying to insert weight " << weight << " which is close to"
                << " existing weight " << w << " in " << sortedWeights
                << endl;
        }

        if (w > weight)
        {
            // Start inserting before sortedI.
            insertI = sortedI;
            break;
        }
    }


    label sz = sortedWeights.size();

    sortedWeights.setSize(sz + 1);
    sortedVerts.setSize(sz + 1);

    // Leave everything up to (not including) insertI intact.

    // Make space by copying from insertI up.
    for (label i = sz-1; i >= insertI; --i)
    {
        sortedWeights[i+1] = sortedWeights[i];
        sortedVerts[i+1] = sortedVerts[i];
    }
    sortedWeights[insertI] = weight;
    sortedVerts[insertI] = vertI;

    return true;
}


// Is triangle candidate for collapse? Small height or small quality
bool isSliver
(
    const triSurface& surf,
    const scalar minLen,
    const scalar minQuality,
    const label facei,
    const label edgeI
)
{
    const pointField& localPoints = surf.localPoints();

    // Check
    // - opposite vertex projects onto base edge
    // - normal distance is small
    // - or triangle quality is small

    label opposite0 =
        triSurfaceTools::oppositeVertex
        (
            surf,
            facei,
            edgeI
        );

    const edge& e = surf.edges()[edgeI];
    const labelledTri& f = surf[facei];

    pointHit pHit =
        e.line(localPoints).nearestDist
        (
            localPoints[opposite0]
        );

    if
    (
        pHit.hit()
     && (
            pHit.distance() < minLen
         || f.tri(surf.points()).quality() < minQuality
        )
    )
    {
        // Remove facei and split all other faces using this
        // edge. This is done by 'replacing' the edgeI with the
        // opposite0 vertex
        //Pout<< "Splitting face " << facei << " since distance "
        //    << pHit.distance()
        //    << " from vertex " << opposite0
        //    << " to edge " << edgeI
        //    << "  points "
        //    << localPoints[e[0]]
        //    << localPoints[e[1]]
        //    << " is too small or triangle quality "
        //    << f.tri(surf.points()).quality()
        //    << " too small." << endl;

        return true;
    }
    else
    {
        return false;
    }
}


// Mark all faces that are going to be collapsed.
// faceToEdge: per face -1 or the base edge of the face.
static void markCollapsedFaces
(
    const triSurface& surf,
    const scalar minLen,
    const scalar minQuality,
    labelList& faceToEdge
)
{
    faceToEdge.setSize(surf.size());
    faceToEdge = -1;

    const labelListList& edgeFaces = surf.edgeFaces();

    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = surf.edgeFaces()[edgeI];

        forAll(eFaces, i)
        {
            label facei = eFaces[i];

            bool isCandidate = isSliver(surf, minLen, minQuality, facei, edgeI);

            if (isCandidate)
            {
                // Mark face as being collapsed
                if (faceToEdge[facei] != -1)
                {
                    FatalErrorInFunction
                        << "Cannot collapse face " << facei << " since "
                        << " is marked to be collapsed both to edge "
                        << faceToEdge[facei] << " and " << edgeI
                        << abort(FatalError);
                }

                faceToEdge[facei] = edgeI;
            }
        }
    }
}


// Recurse through collapsed faces marking all of them with regionI (in
// collapseRegion)
static void markRegion
(
    const triSurface& surf,
    const labelList& faceToEdge,
    const label regionI,
    const label facei,
    labelList& collapseRegion
)
{
    if (faceToEdge[facei] == -1 || collapseRegion[facei] != -1)
    {
        FatalErrorInFunction
            << "Problem : crossed into uncollapsed/regionized face"
            << abort(FatalError);
    }

    collapseRegion[facei] = regionI;

    // Recurse across edges to collapsed neighbours

    const labelList& fEdges = surf.faceEdges()[facei];

    forAll(fEdges, fEdgeI)
    {
        label edgeI = fEdges[fEdgeI];

        const labelList& eFaces = surf.edgeFaces()[edgeI];

        forAll(eFaces, i)
        {
            label nbrFacei = eFaces[i];

            if (faceToEdge[nbrFacei] != -1)
            {
                if (collapseRegion[nbrFacei] == -1)
                {
                    markRegion
                    (
                        surf,
                        faceToEdge,
                        regionI,
                        nbrFacei,
                        collapseRegion
                    );
                }
                else if (collapseRegion[nbrFacei] != regionI)
                {
                    FatalErrorInFunction
                        << "Edge:" << edgeI << " between face " << facei
                        << " with region " << regionI
                        << " and face " << nbrFacei
                        << " with region " << collapseRegion[nbrFacei]
                        << endl;
                }
            }
        }
    }
}


// Mark every face with region (in collapseRegion) (or -1).
// Return number of regions.
static label markRegions
(
    const triSurface& surf,
    const labelList& faceToEdge,
    labelList& collapseRegion
)
{
    label regionI = 0;

    forAll(faceToEdge, facei)
    {
        if (collapseRegion[facei] == -1 && faceToEdge[facei] != -1)
        {
            //Pout<< "markRegions : Marking region:" << regionI
            //    << " starting from face " << facei << endl;

            // Collapsed face. Mark connected region with current region number
            markRegion(surf, faceToEdge, regionI++, facei, collapseRegion);
        }
    }
    return regionI;
}


// Type of region.
// -1  : edge in between uncollapsed faces.
// -2  : edge in between collapsed faces
// >=0 : edge in between uncollapsed and collapsed region. Returns region.
static label edgeType
(
    const triSurface& surf,
    const labelList& collapseRegion,
    const label edgeI
)
{
    const labelList& eFaces = surf.edgeFaces()[edgeI];

    // Detect if edge is in between collapseRegion and non-collapse face
    bool usesUncollapsed = false;
    label usesRegion = -1;

    forAll(eFaces, i)
    {
        label facei = eFaces[i];

        label region = collapseRegion[facei];

        if (region == -1)
        {
            usesUncollapsed = true;
        }
        else if (usesRegion == -1)
        {
            usesRegion = region;
        }
        else if (usesRegion != region)
        {
            FatalErrorInFunction << abort(FatalError);
        }
        else
        {
            // Equal regions.
        }
    }

    if (usesUncollapsed)
    {
        if (usesRegion == -1)
        {
            // uncollapsed faces only.
            return -1;
        }
        else
        {
            // between uncollapsed and collapsed.
            return usesRegion;
        }
    }
    else
    {
        if (usesRegion == -1)
        {
            FatalErrorInFunction << abort(FatalError);
            return -2;
        }
        else
        {
            return -2;
        }
    }
}


// Get points on outside edge of region (= outside points)
static labelListList getOutsideVerts
(
    const triSurface& surf,
    const labelList& collapseRegion,
    const label nRegions
)
{
    const labelListList& edgeFaces = surf.edgeFaces();

    // Per region all the outside vertices.
    labelListList outsideVerts(nRegions);

    forAll(edgeFaces, edgeI)
    {
        // Detect if edge is in between collapseRegion and non-collapse face
        label regionI = edgeType(surf, collapseRegion, edgeI);

        if (regionI >= 0)
        {
            // Edge borders both uncollapsed face and collapsed face on region
            // usesRegion.

            const edge& e = surf.edges()[edgeI];

            labelList& regionVerts = outsideVerts[regionI];

            // Add both edge points to regionVerts.
            forAll(e, eI)
            {
                label v = e[eI];

                if (findIndex(regionVerts, v) == -1)
                {
                    label sz = regionVerts.size();
                    regionVerts.setSize(sz+1);
                    regionVerts[sz] = v;
                }
            }
        }
    }

    return outsideVerts;
}


// n^2 search for furthest removed point pair.
static labelPair getSpanPoints
(
    const triSurface& surf,
    const labelList& outsideVerts
)
{
    const pointField& localPoints = surf.localPoints();

    scalar maxDist = -great;
    labelPair maxPair;

    forAll(outsideVerts, i)
    {
        label v0 = outsideVerts[i];

        for (label j = i+1; j < outsideVerts.size(); j++)
        {
            label v1 = outsideVerts[j];

            scalar d = mag(localPoints[v0] - localPoints[v1]);

            if (d > maxDist)
            {
                maxDist = d;
                maxPair[0] = v0;
                maxPair[1] = v1;
            }
        }
    }

    return maxPair;
}


// Project all non-span points onto the span edge.
static void projectNonSpanPoints
(
    const triSurface& surf,
    const labelList& outsideVerts,
    const labelPair& spanPair,
    labelList& sortedVertices,
    scalarField& sortedWeights
)
{
    const point& p0 = surf.localPoints()[spanPair[0]];
    const point& p1 = surf.localPoints()[spanPair[1]];

    forAll(outsideVerts, i)
    {
        label v = outsideVerts[i];

        if (v != spanPair[0] && v != spanPair[1])
        {
            // Is a non-span point. Project onto spanning edge.

            pointHit pHit =
                linePointRef(p0, p1).nearestDist
                (
                    surf.localPoints()[v]
                );

            if (!pHit.hit())
            {
                FatalErrorInFunction
                    << abort(FatalError);
            }

            scalar w = mag(pHit.hitPoint() - p0) / mag(p1 - p0);

            insertSorted(v, w, sortedVertices, sortedWeights);
        }
    }
}


// Slice part of the orderVertices (and optionally reverse) for this edge.
static void getSplitVerts
(
    const triSurface& surf,
    const label regionI,
    const labelPair& spanPoints,
    const labelList& orderedVerts,
    const scalarField& orderedWeights,
    const label edgeI,

    labelList& splitVerts,
    scalarField& splitWeights
)
{
    const edge& e = surf.edges()[edgeI];
    const label sz = orderedVerts.size();

    if (e[0] == spanPoints[0])
    {
        // Edge in same order as spanPoints&orderedVerts. Keep order.

        if (e[1] == spanPoints[1])
        {
            // Copy all.
            splitVerts = orderedVerts;
            splitWeights = orderedWeights;
        }
        else
        {
            // Copy up to (but not including) e[1]
            label i1 = findIndex(orderedVerts, e[1]);
            splitVerts = SubList<label>(orderedVerts, i1, 0);
            splitWeights = SubList<scalar>(orderedWeights, i1, 0);
        }
    }
    else if (e[0] == spanPoints[1])
    {
        // Reverse.

        if (e[1] == spanPoints[0])
        {
            // Copy all.
            splitVerts = orderedVerts;
            reverse(splitVerts);
            splitWeights = orderedWeights;
            reverse(splitWeights);
        }
        else
        {
            // Copy downto (but not including) e[1]

            label i1 = findIndex(orderedVerts, e[1]);
            splitVerts = SubList<label>(orderedVerts, sz-(i1+1), i1+1);
            reverse(splitVerts);
            splitWeights = SubList<scalar>(orderedWeights, sz-(i1+1), i1+1);
            reverse(splitWeights);
        }
    }
    else if (e[1] == spanPoints[0])
    {
        // Reverse.

        // Copy up to (but not including) e[0]

        label i0 = findIndex(orderedVerts, e[0]);
        splitVerts = SubList<label>(orderedVerts, i0, 0);
        reverse(splitVerts);
        splitWeights = SubList<scalar>(orderedWeights, i0, 0);
        reverse(splitWeights);
    }
    else if (e[1] == spanPoints[1])
    {
        // Copy from (but not including) e[0] to end

        label i0 = findIndex(orderedVerts, e[0]);
        splitVerts = SubList<label>(orderedVerts, sz-(i0+1), i0+1);
        splitWeights = SubList<scalar>(orderedWeights, sz-(i0+1), i0+1);
    }
    else
    {
        label i0 = findIndex(orderedVerts, e[0]);
        label i1 = findIndex(orderedVerts, e[1]);

        if (i0 == -1 || i1 == -1)
        {
            FatalErrorInFunction
                << "Did not find edge in projected vertices." << nl
                << "region:" << regionI << nl
                << "spanPoints:" << spanPoints
                << "  coords:" << surf.localPoints()[spanPoints[0]]
                << surf.localPoints()[spanPoints[1]] << nl
                << "edge:" << edgeI
                << "  verts:" << e
                << "  coords:" << surf.localPoints()[e[0]]
                << surf.localPoints()[e[1]] << nl
                << "orderedVerts:" << orderedVerts << nl
                << abort(FatalError);
        }

        if (i0 < i1)
        {
            splitVerts = SubList<label>(orderedVerts, i1-i0-1, i0+1);
            splitWeights = SubList<scalar>(orderedWeights, i1-i0-1, i0+1);
        }
        else
        {
            splitVerts = SubList<label>(orderedVerts, i0-i1-1, i1+1);
            reverse(splitVerts);
            splitWeights = SubList<scalar>(orderedWeights, i0-i1-1, i1+1);
            reverse(splitWeights);
        }
    }
}


label collapseBase
(
    triSurface& surf,
    const scalar minLen,
    const scalar minQuality
)
{
    label nTotalSplit = 0;

    label iter = 0;

    while (true)
    {
        // Detect faces to collapse
        // ~~~~~~~~~~~~~~~~~~~~~~~~

        // -1 or edge the face is collapsed onto.
        labelList faceToEdge(surf.size(), -1);

        // Calculate faceToEdge (face collapses)
        markCollapsedFaces(surf, minLen, minQuality, faceToEdge);


        // Find regions of connected collapsed faces
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // per face -1 or region
        labelList collapseRegion(surf.size(), -1);

        label nRegions = markRegions(surf, faceToEdge, collapseRegion);

        //Pout<< "Detected " << nRegions << " regions of faces to be collapsed"
        //    << nl << endl;

        // Pick up all vertices on outside of region
        labelListList outsideVerts
        (
            getOutsideVerts(surf, collapseRegion, nRegions)
        );

        // For all regions determine maximum distance between points
        List<labelPair> spanPoints(nRegions);
        labelListList orderedVertices(nRegions);
        List<scalarField> orderedWeights(nRegions);

        forAll(spanPoints, regionI)
        {
            spanPoints[regionI] = getSpanPoints(surf, outsideVerts[regionI]);

            //Pout<< "For region " << regionI << " found extrema at points "
            //    << surf.localPoints()[spanPoints[regionI][0]]
            //    << surf.localPoints()[spanPoints[regionI][1]]
            //    << endl;

            // Project all non-span points onto the span edge.
            projectNonSpanPoints
            (
                surf,
                outsideVerts[regionI],
                spanPoints[regionI],
                orderedVertices[regionI],
                orderedWeights[regionI]
            );

            //Pout<< "For region:" << regionI
            //    << " span:" << spanPoints[regionI]
            //    << " orderedVerts:" << orderedVertices[regionI]
            //    << " orderedWeights:" << orderedWeights[regionI]
            //    << endl;

            //writeRegionOBJ
            //(
            //    surf,
            //    regionI,
            //    collapseRegion,
            //    outsideVerts[regionI]
            //);

            //Pout<< endl;
        }



        // Actually split the edges
        // ~~~~~~~~~~~~~~~~~~~~~~~~


        const List<labelledTri>& localFaces = surf.localFaces();
        const edgeList& edges = surf.edges();

        label nSplit = 0;

        // Storage for new triangles.
        DynamicList<labelledTri> newTris(surf.size());

        // Whether face has been dealt with (either copied/split or deleted)
        boolList faceHandled(surf.size(), false);


        forAll(edges, edgeI)
        {
            const edge& e = edges[edgeI];

            // Detect if edge is in between collapseRegion and non-collapse face
            label regionI = edgeType(surf, collapseRegion, edgeI);

            if (regionI == -2)
            {
                // in between collapsed faces. nothing needs to be done.
            }
            else if (regionI == -1)
            {
                // edge in between uncollapsed faces. Handle these later on.
            }
            else
            {
                // some faces around edge are collapsed.

                // Find additional set of points on edge to be used to split
                // the remaining faces.

                labelList splitVerts;
                scalarField splitWeights;
                getSplitVerts
                (
                    surf,
                    regionI,
                    spanPoints[regionI],
                    orderedVertices[regionI],
                    orderedWeights[regionI],
                    edgeI,

                    splitVerts,
                    splitWeights
                );

                if (splitVerts.size())
                {
                    // Split edge using splitVerts. All non-collapsed triangles
                    // using edge will get split.

                    //{
                    //    const pointField& localPoints = surf.localPoints();
                    //    Pout<< "edge " << edgeI << ' ' << e
                    //        << "  points "
                    //        << localPoints[e[0]] << ' ' << localPoints[e[1]]
                    //        << " split into edges with extra points:"
                    //        << endl;
                    //    forAll(splitVerts, i)
                    //    {
                    //        Pout<< "    " << splitVerts[i] << " weight "
                    //            << splitWeights[i] << nl;
                    //    }
                    //}

                    const labelList& eFaces = surf.edgeFaces()[edgeI];

                    forAll(eFaces, i)
                    {
                        label facei = eFaces[i];

                        if (!faceHandled[facei] && faceToEdge[facei] == -1)
                        {
                            // Split face to use vertices.
                            splitTri
                            (
                                localFaces[facei],
                                e,
                                splitVerts,
                                newTris
                            );

                            faceHandled[facei] = true;

                            nSplit++;
                        }
                    }
                }
            }
        }

        // Copy all unsplit faces
        forAll(faceHandled, facei)
        {
            if (!faceHandled[facei] && faceToEdge[facei] == -1)
            {
                newTris.append(localFaces[facei]);
            }
        }

        Info<< "collapseBase : collapsing " << nSplit
            << " triangles by splitting their base edge."
            << endl;

        nTotalSplit += nSplit;

        if (nSplit == 0)
        {
            break;
        }

        // Pack the triangles
        newTris.shrink();

        //Pout<< "Resetting surface from " << surf.size() << " to "
        //    << newTris.size() << " triangles" << endl;
        surf = triSurface(newTris, surf.patches(), surf.localPoints());

        //{
        //    fileName fName("bla" + name(iter) + ".obj");
        //    Pout<< "Writing surf to " << fName << endl;
        //    surf.write(fName);
        //}

        iter++;
    }

    // Remove any unused vertices
    surf = triSurface(surf.localFaces(), surf.patches(), surf.localPoints());

    return nTotalSplit;
}


// ************************************************************************* //
