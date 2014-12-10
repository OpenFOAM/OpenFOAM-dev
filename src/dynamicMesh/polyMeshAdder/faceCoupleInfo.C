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

#include "faceCoupleInfo.H"
#include "polyMesh.H"
#include "matchPoints.H"
#include "indirectPrimitivePatch.H"
#include "meshTools.H"
#include "treeDataFace.H"
#include "indexedOctree.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(faceCoupleInfo, 0);

const scalar faceCoupleInfo::angleTol_ = 1e-3;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Write edges
void Foam::faceCoupleInfo::writeOBJ
(
    const fileName& fName,
    const edgeList& edges,
    const pointField& points,
    const bool compact
)
{
    OFstream str(fName);

    labelList pointMap(points.size(), -1);

    if (compact)
    {
        label newPointI = 0;

        forAll(edges, edgeI)
        {
            const edge& e = edges[edgeI];

            forAll(e, eI)
            {
                label pointI = e[eI];

                if (pointMap[pointI] == -1)
                {
                    pointMap[pointI] = newPointI++;

                    meshTools::writeOBJ(str, points[pointI]);
                }
            }
        }
    }
    else
    {
        forAll(points, pointI)
        {
            meshTools::writeOBJ(str, points[pointI]);
        }

        pointMap = identity(points.size());
    }

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        str<< "l " << pointMap[e[0]]+1 << ' ' << pointMap[e[1]]+1 << nl;
    }
}


//- Writes edges.
void Foam::faceCoupleInfo::writeOBJ
(
    const fileName& fName,
    const pointField& points0,
    const pointField& points1
)
{
    Pout<< "Writing connections as edges to " << fName << endl;

    OFstream str(fName);

    label vertI = 0;

    forAll(points0, i)
    {
        meshTools::writeOBJ(str, points0[i]);
        vertI++;
        meshTools::writeOBJ(str, points1[i]);
        vertI++;
        str << "l " << vertI-1 << ' ' << vertI << nl;
    }
}


//- Writes face and point connectivity as .obj files.
void Foam::faceCoupleInfo::writePointsFaces() const
{
    const indirectPrimitivePatch& m = masterPatch();
    const indirectPrimitivePatch& s = slavePatch();
    const primitiveFacePatch& c = cutFaces();

    // Patches
    {
        OFstream str("masterPatch.obj");
        Pout<< "Writing masterPatch to " << str.name() << endl;
        meshTools::writeOBJ(str, m.localFaces(), m.localPoints());
    }
    {
        OFstream str("slavePatch.obj");
        Pout<< "Writing slavePatch to " << str.name() << endl;
        meshTools::writeOBJ(str, s.localFaces(), s.localPoints());
    }
    {
        OFstream str("cutFaces.obj");
        Pout<< "Writing cutFaces to " << str.name() << endl;
        meshTools::writeOBJ(str, c.localFaces(), c.localPoints());
    }

    // Point connectivity
    {
        Pout<< "Writing cutToMasterPoints to cutToMasterPoints.obj" << endl;

        writeOBJ
        (
            "cutToMasterPoints.obj",
            m.localPoints(),
            pointField(c.localPoints(), masterToCutPoints_));
    }
    {
        Pout<< "Writing cutToSlavePoints to cutToSlavePoints.obj" << endl;

        writeOBJ
        (
            "cutToSlavePoints.obj",
            s.localPoints(),
            pointField(c.localPoints(), slaveToCutPoints_)
        );
    }

    // Face connectivity
    {
        Pout<< "Writing cutToMasterFaces to cutToMasterFaces.obj" << endl;

        pointField equivMasterFaces(c.size());

        forAll(cutToMasterFaces(), cutFaceI)
        {
            label masterFaceI = cutToMasterFaces()[cutFaceI];

            if (masterFaceI != -1)
            {
                equivMasterFaces[cutFaceI] = m[masterFaceI].centre(m.points());
            }
            else
            {
                WarningIn("writePointsFaces()")
                    << "No master face for cut face " << cutFaceI
                    << " at position " << c[cutFaceI].centre(c.points())
                    << endl;

                equivMasterFaces[cutFaceI] = vector::zero;
            }
        }

        writeOBJ
        (
            "cutToMasterFaces.obj",
            calcFaceCentres<List>(c, cutPoints(), 0, c.size()),
            equivMasterFaces
        );
    }

    {
        Pout<< "Writing cutToSlaveFaces to cutToSlaveFaces.obj" << endl;

        pointField equivSlaveFaces(c.size());

        forAll(cutToSlaveFaces(), cutFaceI)
        {
            label slaveFaceI = cutToSlaveFaces()[cutFaceI];

            equivSlaveFaces[cutFaceI] = s[slaveFaceI].centre(s.points());
        }

        writeOBJ
        (
            "cutToSlaveFaces.obj",
            calcFaceCentres<List>(c, cutPoints(), 0, c.size()),
            equivSlaveFaces
        );
    }

    Pout<< endl;
}


void Foam::faceCoupleInfo::writeEdges
(
    const labelList& cutToMasterEdges,
    const labelList& cutToSlaveEdges
) const
{
    const indirectPrimitivePatch& m = masterPatch();
    const indirectPrimitivePatch& s = slavePatch();
    const primitiveFacePatch& c = cutFaces();

    // Edge connectivity
    {
        OFstream str("cutToMasterEdges.obj");
        Pout<< "Writing cutToMasterEdges to " << str.name() << endl;

        label vertI = 0;

        forAll(cutToMasterEdges, cutEdgeI)
        {
            if (cutToMasterEdges[cutEdgeI] != -1)
            {
                const edge& masterEdge =
                    m.edges()[cutToMasterEdges[cutEdgeI]];
                const edge& cutEdge = c.edges()[cutEdgeI];

                meshTools::writeOBJ(str, m.localPoints()[masterEdge[0]]);
                vertI++;
                meshTools::writeOBJ(str, m.localPoints()[masterEdge[1]]);
                vertI++;
                meshTools::writeOBJ(str, c.localPoints()[cutEdge[0]]);
                vertI++;
                meshTools::writeOBJ(str, c.localPoints()[cutEdge[1]]);
                vertI++;
                str << "l " << vertI-3 << ' ' << vertI-2 << nl;
                str << "l " << vertI-3 << ' ' << vertI-1 << nl;
                str << "l " << vertI-3 << ' ' << vertI << nl;
                str << "l " << vertI-2 << ' ' << vertI-1 << nl;
                str << "l " << vertI-2 << ' ' << vertI << nl;
                str << "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }
    {
        OFstream str("cutToSlaveEdges.obj");
        Pout<< "Writing cutToSlaveEdges to " << str.name() << endl;

        label vertI = 0;

        labelList slaveToCut(invert(s.nEdges(), cutToSlaveEdges));

        forAll(slaveToCut, edgeI)
        {
            if (slaveToCut[edgeI] != -1)
            {
                const edge& slaveEdge = s.edges()[edgeI];
                const edge& cutEdge = c.edges()[slaveToCut[edgeI]];

                meshTools::writeOBJ(str, s.localPoints()[slaveEdge[0]]);
                vertI++;
                meshTools::writeOBJ(str, s.localPoints()[slaveEdge[1]]);
                vertI++;
                meshTools::writeOBJ(str, c.localPoints()[cutEdge[0]]);
                vertI++;
                meshTools::writeOBJ(str, c.localPoints()[cutEdge[1]]);
                vertI++;
                str << "l " << vertI-3 << ' ' << vertI-2 << nl;
                str << "l " << vertI-3 << ' ' << vertI-1 << nl;
                str << "l " << vertI-3 << ' ' << vertI << nl;
                str << "l " << vertI-2 << ' ' << vertI-1 << nl;
                str << "l " << vertI-2 << ' ' << vertI << nl;
                str << "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }

    Pout<< endl;
}


// Given an edgelist and a map for the points on the edges it tries to find
// the corresponding patch edges.
Foam::labelList Foam::faceCoupleInfo::findMappedEdges
(
    const edgeList& edges,
    const labelList& pointMap,
    const indirectPrimitivePatch& patch
)
{
    labelList toPatchEdges(edges.size());

    forAll(toPatchEdges, edgeI)
    {
        const edge& e = edges[edgeI];

        label v0 = pointMap[e[0]];
        label v1 = pointMap[e[1]];

        toPatchEdges[edgeI] =
            meshTools::findEdge
            (
                patch.edges(),
                patch.pointEdges()[v0],
                v0,
                v1
            );
    }
    return toPatchEdges;
}


// Detect a cut edge which originates from two boundary faces having different
// polyPatches.
bool Foam::faceCoupleInfo::regionEdge
(
    const polyMesh& slaveMesh,
    const label slaveEdgeI
) const
{
    const labelList& eFaces = slavePatch().edgeFaces()[slaveEdgeI];

    if (eFaces.size() == 1)
    {
        return true;
    }
    else
    {
        // Count how many different patches connected to this edge.

        label patch0 = -1;

        forAll(eFaces, i)
        {
            label faceI = eFaces[i];

            label meshFaceI = slavePatch().addressing()[faceI];

            label patchI = slaveMesh.boundaryMesh().whichPatch(meshFaceI);

            if (patch0 == -1)
            {
                patch0 = patchI;
            }
            else if (patchI != patch0)
            {
                // Found two different patches connected to this edge.
                return true;
            }
        }
        return false;
    }
}


// Find edge using pointI that is most aligned with vector between
// master points. Patchdivision tells us whether or not to use
// patch information to match edges.
Foam::label Foam::faceCoupleInfo::mostAlignedCutEdge
(
    const bool report,
    const polyMesh& slaveMesh,
    const bool patchDivision,
    const labelList& cutToMasterEdges,
    const labelList& cutToSlaveEdges,
    const label pointI,
    const label edgeStart,
    const label edgeEnd
) const
{
    const pointField& localPoints = cutFaces().localPoints();

    const labelList& pEdges = cutFaces().pointEdges()[pointI];

    if (report)
    {
        Pout<< "mostAlignedEdge : finding nearest edge among "
            << UIndirectList<edge>(cutFaces().edges(), pEdges)()
            << " connected to point " << pointI
            << " coord:" << localPoints[pointI]
            << " running between " << edgeStart << " coord:"
            << localPoints[edgeStart]
            << " and " << edgeEnd << " coord:"
            << localPoints[edgeEnd]
            << endl;
    }

    // Find the edge that gets us nearest end.

    label maxEdgeI = -1;
    scalar maxCos = -GREAT;

    forAll(pEdges, i)
    {
        label edgeI = pEdges[i];

        if
        (
           !(
                patchDivision
             && cutToMasterEdges[edgeI] == -1
            )
         || (
                patchDivision
             && regionEdge(slaveMesh, cutToSlaveEdges[edgeI])
            )
        )
        {
            const edge& e = cutFaces().edges()[edgeI];

            label otherPointI = e.otherVertex(pointI);

            if (otherPointI == edgeEnd)
            {
                // Shortcut: found edge end point.
                if (report)
                {
                    Pout<< "    mostAlignedEdge : found end point " << edgeEnd
                        << endl;
                }
                return edgeI;
            }

            // Get angle between edge and edge to masterEnd

            vector eVec(localPoints[otherPointI] - localPoints[pointI]);

            scalar magEVec = mag(eVec);

            if (magEVec < VSMALL)
            {
                WarningIn("faceCoupleInfo::mostAlignedEdge")
                    << "Crossing zero sized edge " << edgeI
                    << " coords:" << localPoints[otherPointI]
                    << localPoints[pointI]
                    << " when walking from " << localPoints[edgeStart]
                    << " to " << localPoints[edgeEnd]
                    << endl;
                return edgeI;
            }

            eVec /= magEVec;

            vector eToEndPoint(localPoints[edgeEnd] - localPoints[otherPointI]);
            eToEndPoint /= mag(eToEndPoint);

            scalar cosAngle = eVec & eToEndPoint;

            if (report)
            {
                Pout<< "    edge:" << e << " points:" << localPoints[pointI]
                    << localPoints[otherPointI]
                    << "  vec:" << eVec
                    << "  vecToEnd:" << eToEndPoint
                    << " cosAngle:" << cosAngle
                    << endl;
            }

            if (cosAngle > maxCos)
            {
                maxCos = cosAngle;
                maxEdgeI = edgeI;
            }
        }
    }

    if (maxCos > 1 - angleTol_)
    {
        return maxEdgeI;
    }
    else
    {
        return -1;
    }
}


// Construct points to split points map (in cut addressing)
void Foam::faceCoupleInfo::setCutEdgeToPoints(const labelList& cutToMasterEdges)
{
    labelListList masterToCutEdges
    (
        invertOneToMany
        (
            masterPatch().nEdges(),
            cutToMasterEdges
        )
    );

    const edgeList& cutEdges = cutFaces().edges();

    // Size extra big so searching is faster
    cutEdgeToPoints_.resize
    (
        masterPatch().nEdges()
      + slavePatch().nEdges()
      + cutEdges.size()
    );

    forAll(masterToCutEdges, masterEdgeI)
    {
        const edge& masterE = masterPatch().edges()[masterEdgeI];

        //Pout<< "Master:" << masterPatch().localPoints()[masterE[0]] << ' '
        //    << masterPatch().localPoints()[masterE[1]] << endl;

        const labelList& stringedEdges = masterToCutEdges[masterEdgeI];

        if (stringedEdges.empty())
        {
            FatalErrorIn
            (
                "faceCoupleInfo::setCutEdgeToPoints"
                "(const labelList&)"
            )   << "Did not match all of master edges to cutFace edges"
                << nl
                << "First unmatched edge:" << masterEdgeI << " endPoints:"
                << masterPatch().localPoints()[masterE[0]]
                << masterPatch().localPoints()[masterE[1]]
                << endl
                << "This usually means that the slave patch is not a"
                << " subdivision of the master patch"
                << abort(FatalError);
        }
        else if (stringedEdges.size() > 1)
        {
            // String up the edges between e[0] and e[1]. Store the points
            // inbetween e[0] and e[1] (all in cutFaces() labels)

            DynamicList<label> splitPoints(stringedEdges.size()-1);

            // Unsplit edge endpoints
            const edge unsplitEdge
            (
                masterToCutPoints_[masterE[0]],
                masterToCutPoints_[masterE[1]]
            );

            label startVertI = unsplitEdge[0];
            label startEdgeI = -1;

            while (startVertI != unsplitEdge[1])
            {
                // Loop over all string of edges. Update
                // - startVertI : previous vertex
                // - startEdgeI : previous edge
                // and insert any points into splitPoints

                // For checking
                label oldStart = startVertI;

                forAll(stringedEdges, i)
                {
                    label edgeI = stringedEdges[i];

                    if (edgeI != startEdgeI)
                    {
                        const edge& e = cutEdges[edgeI];

                        //Pout<< "    cut:" << e << " at:"
                        //    << cutFaces().localPoints()[e[0]]
                        //    << ' ' << cutFaces().localPoints()[e[1]] << endl;

                        if (e[0] == startVertI)
                        {
                            startEdgeI = edgeI;
                            startVertI = e[1];
                            if (e[1] != unsplitEdge[1])
                            {
                                splitPoints.append(e[1]);
                            }
                            break;
                        }
                        else if (e[1] == startVertI)
                        {
                            startEdgeI = edgeI;
                            startVertI = e[0];
                            if (e[0] != unsplitEdge[1])
                            {
                                splitPoints.append(e[0]);
                            }
                            break;
                        }
                    }
                }

                // Check
                if (oldStart == startVertI)
                {
                    FatalErrorIn
                    (
                        "faceCoupleInfo::setCutEdgeToPoints"
                        "(const labelList&)"
                    )   << " unsplitEdge:" << unsplitEdge
                        << " does not correspond to split edges "
                        << UIndirectList<edge>(cutEdges, stringedEdges)()
                        << abort(FatalError);
                }
            }

            //Pout<< "For master edge:"
            //    << unsplitEdge
            //    << " Found stringed points "
            //    <<  UIndirectList<point>
            //        (
            //            cutFaces().localPoints(),
            //            splitPoints.shrink()
            //        )()
            //    << endl;

            cutEdgeToPoints_.insert(unsplitEdge, splitPoints.shrink());
        }
    }
}


// Determines rotation for f1 to match up with f0, i.e. the index in f0 of
// the first point of f1.
Foam::label Foam::faceCoupleInfo::matchFaces
(
    const scalar absTol,
    const pointField& points0,
    const face& f0,
    const pointField& points1,
    const face& f1,
    const bool sameOrientation
)
{
    if (f0.size() != f1.size())
    {
        FatalErrorIn
        (
            "faceCoupleInfo::matchFaces"
            "(const scalar, const face&, const pointField&"
            ", const face&, const pointField&)"
        )   << "Different sizes for supposedly matching faces." << nl
            << "f0:" << f0 << " coords:" << UIndirectList<point>(points0, f0)()
            << nl
            << "f1:" << f1 << " coords:" << UIndirectList<point>(points1, f1)()
            << abort(FatalError);
    }

    const scalar absTolSqr = sqr(absTol);


    label matchFp = -1;

    forAll(f0, startFp)
    {
        // See -if starting from startFp on f0- the two faces match.
        bool fullMatch = true;

        label fp0 = startFp;
        label fp1 = 0;

        forAll(f1, i)
        {
            scalar distSqr = Foam::magSqr(points0[f0[fp0]] - points1[f1[fp1]]);

            if (distSqr > absTolSqr)
            {
                fullMatch = false;
                break;
            }

            fp0 = f0.fcIndex(fp0);  // walk forward

            if (sameOrientation)
            {
                fp1 = f1.fcIndex(fp1);
            }
            else
            {
                fp1 = f1.rcIndex(fp1);
            }
        }

        if (fullMatch)
        {
            matchFp = startFp;
            break;
        }
    }

    if (matchFp == -1)
    {
        FatalErrorIn
        (
            "faceCoupleInfo::matchFaces"
            "(const scalar, const face&, const pointField&"
            ", const face&, const pointField&)"
        )   << "No unique match between two faces" << nl
            << "Face " << f0 << " coords "
            << UIndirectList<point>(points0, f0)() << nl
            << "Face " << f1 << " coords "
            << UIndirectList<point>(points1, f1)()
            << "when using tolerance " << absTol
            << " and forwardMatching:" << sameOrientation
            << abort(FatalError);
    }

    return matchFp;
}


// Find correspondence from patch points to cut points. This might
// detect shared points so the output is a patch-to-cut point list
// and a compaction list for the cut points (which will always be equal or more
// connected than the patch).
// Returns true if there are any duplicates.
bool Foam::faceCoupleInfo::matchPointsThroughFaces
(
    const scalar absTol,
    const pointField& cutPoints,
    const faceList& cutFaces,
    const pointField& patchPoints,
    const faceList& patchFaces,
    const bool sameOrientation,

    labelList& patchToCutPoints,    // patch to (uncompacted) cut points
    labelList& cutToCompact,        // compaction list for cut points
    labelList& compactToCut         // inverse ,,
)
{

    // From slave to cut point
    patchToCutPoints.setSize(patchPoints.size());
    patchToCutPoints = -1;

    // Compaction list for cut points: either -1 or index into master which
    // gives the point to compact to.
    labelList cutPointRegion(cutPoints.size(), -1);
    DynamicList<label> cutPointRegionMaster;

    forAll(patchFaces, patchFaceI)
    {
        const face& patchF = patchFaces[patchFaceI];

        //const face& cutF = cutFaces[patchToCutFaces[patchFaceI]];
        const face& cutF = cutFaces[patchFaceI];

        // Do geometric matching to get position of cutF[0] in patchF
        label patchFp = matchFaces
        (
            absTol,
            patchPoints,
            patchF,
            cutPoints,
            cutF,
            sameOrientation        // orientation
        );

        forAll(cutF, cutFp)
        {
            label cutPointI = cutF[cutFp];
            label patchPointI = patchF[patchFp];

            //const point& cutPt = cutPoints[cutPointI];
            //const point& patchPt = patchPoints[patchPointI];
            //if (mag(cutPt - patchPt) > SMALL)
            //{
            //    FatalErrorIn("matchPointsThroughFaces")
            //    << "cutP:" << cutPt
            //    << " patchP:" << patchPt
            //    << abort(FatalError);
            //}

            if (patchToCutPoints[patchPointI] == -1)
            {
                patchToCutPoints[patchPointI] = cutPointI;
            }
            else if (patchToCutPoints[patchPointI] != cutPointI)
            {
                // Multiple cut points connecting to same patch.
                // Check if already have region & region master for this set
                label otherCutPointI = patchToCutPoints[patchPointI];

                //Pout<< "PatchPoint:" << patchPt
                //    << " matches to:" << cutPointI
                //    << " coord:" << cutPoints[cutPointI]
                //    << " and to:" << otherCutPointI
                //    << " coord:" << cutPoints[otherCutPointI]
                //    << endl;

                if (cutPointRegion[otherCutPointI] != -1)
                {
                    // Have region for this set. Copy.
                    label region = cutPointRegion[otherCutPointI];
                    cutPointRegion[cutPointI] = region;

                    // Update region master with min point label
                    cutPointRegionMaster[region] = min
                    (
                        cutPointRegionMaster[region],
                        cutPointI
                    );
                }
                else
                {
                    // Create new region.
                    label region = cutPointRegionMaster.size();
                    cutPointRegionMaster.append
                    (
                        min(cutPointI, otherCutPointI)
                    );
                    cutPointRegion[cutPointI] = region;
                    cutPointRegion[otherCutPointI] = region;
                }
            }

            if (sameOrientation)
            {
                patchFp = patchF.fcIndex(patchFp);
            }
            else
            {
                patchFp = patchF.rcIndex(patchFp);
            }
        }
    }

    // Rework region&master into compaction array
    compactToCut.setSize(cutPointRegion.size());
    cutToCompact.setSize(cutPointRegion.size());
    cutToCompact = -1;
    label compactPointI = 0;

    forAll(cutPointRegion, i)
    {
        if (cutPointRegion[i] == -1)
        {
            // Unduplicated point. Allocate new compacted point.
            cutToCompact[i] = compactPointI;
            compactToCut[compactPointI] = i;
            compactPointI++;
        }
        else
        {
            // Duplicate point. Get master.

            label masterPointI = cutPointRegionMaster[cutPointRegion[i]];

            if (cutToCompact[masterPointI] == -1)
            {
                cutToCompact[masterPointI] = compactPointI;
                compactToCut[compactPointI] = masterPointI;
                compactPointI++;
            }
            cutToCompact[i] = cutToCompact[masterPointI];
        }
    }
    compactToCut.setSize(compactPointI);

    return compactToCut.size() != cutToCompact.size();
}


// Return max distance from any point on cutF to masterF
Foam::scalar Foam::faceCoupleInfo::maxDistance
(
    const face& cutF,
    const pointField& cutPoints,
    const face& masterF,
    const pointField& masterPoints
)
{
    scalar maxDist = -GREAT;

    forAll(cutF, fp)
    {
        const point& cutPt = cutPoints[cutF[fp]];

        pointHit pHit = masterF.nearestPoint(cutPt, masterPoints);

        maxDist = max(maxDist, pHit.distance());
    }
    return maxDist;
}


void Foam::faceCoupleInfo::findPerfectMatchingFaces
(
    const primitiveMesh& mesh0,
    const primitiveMesh& mesh1,
    const scalar absTol,

    labelList& mesh0Faces,
    labelList& mesh1Faces
)
{
    // Face centres of external faces (without invoking
    // mesh.faceCentres since mesh might have been clearedOut)

    pointField fc0
    (
        calcFaceCentres<List>
        (
            mesh0.faces(),
            mesh0.points(),
            mesh0.nInternalFaces(),
            mesh0.nFaces() - mesh0.nInternalFaces()
        )
    );

    pointField fc1
    (
        calcFaceCentres<List>
        (
            mesh1.faces(),
            mesh1.points(),
            mesh1.nInternalFaces(),
            mesh1.nFaces() - mesh1.nInternalFaces()
        )
    );


    if (debug)
    {
        Pout<< "Face matching tolerance : " << absTol << endl;
    }


    // Match geometrically
    labelList from1To0;
    bool matchedAllFaces = matchPoints
    (
        fc1,
        fc0,
        scalarField(fc1.size(), absTol),
        false,
        from1To0
    );

    if (matchedAllFaces)
    {
        Warning
            << "faceCoupleInfo::faceCoupleInfo : "
            << "Matched ALL " << fc1.size()
            << " boundary faces of mesh0 to boundary faces of mesh1." << endl
            << "This is only valid if the mesh to add is fully"
            << " enclosed by the mesh it is added to." << endl;
    }


    // Collect matches.
    label nMatched = 0;

    mesh0Faces.setSize(fc0.size());
    mesh1Faces.setSize(fc1.size());

    forAll(from1To0, i)
    {
        if (from1To0[i] != -1)
        {
            mesh1Faces[nMatched] = i + mesh1.nInternalFaces();
            mesh0Faces[nMatched] = from1To0[i] + mesh0.nInternalFaces();

            nMatched++;
        }
    }

    mesh0Faces.setSize(nMatched);
    mesh1Faces.setSize(nMatched);
}


void Foam::faceCoupleInfo::findSlavesCoveringMaster
(
    const primitiveMesh& mesh0,
    const primitiveMesh& mesh1,
    const scalar absTol,

    labelList& mesh0Faces,
    labelList& mesh1Faces
)
{
    // Construct octree from all mesh0 boundary faces
    labelList bndFaces(mesh0.nFaces()-mesh0.nInternalFaces());
    forAll(bndFaces, i)
    {
        bndFaces[i] = mesh0.nInternalFaces() + i;
    }

    treeBoundBox overallBb(mesh0.points());

    Random rndGen(123456);

    indexedOctree<treeDataFace> tree
    (
        treeDataFace    // all information needed to search faces
        (
            false,                      // do not cache bb
            mesh0,
            bndFaces                    // boundary faces only
        ),
        overallBb.extend(rndGen, 1e-4), // overall search domain
        8,                              // maxLevel
        10,                             // leafsize
        3.0                             // duplicity
    );

    if (debug)
    {
        Pout<< "findSlavesCoveringMaster :"
            << " constructed octree for mesh0 boundary faces" << endl;
    }

    // Search nearest face for every mesh1 boundary face.

    labelHashSet mesh0Set(mesh0.nFaces() - mesh0.nInternalFaces());
    labelHashSet mesh1Set(mesh1.nFaces() - mesh1.nInternalFaces());

    for
    (
        label mesh1FaceI = mesh1.nInternalFaces();
        mesh1FaceI < mesh1.nFaces();
        mesh1FaceI++
    )
    {
        const face& f1 = mesh1.faces()[mesh1FaceI];

        // Generate face centre (prevent cellCentres() reconstruction)
        point fc(f1.centre(mesh1.points()));

        pointIndexHit nearInfo = tree.findNearest(fc, Foam::sqr(absTol));

        if (nearInfo.hit())
        {
            label mesh0FaceI = tree.shapes().faceLabels()[nearInfo.index()];

            // Check if points of f1 actually lie on top of mesh0 face
            // This is the bit that might fail since if f0 is severely warped
            // and f1's points are not present on f0 (since subdivision)
            // then the distance of the points to f0 might be quite large
            // and the test will fail. This all should in fact be some kind
            // of walk from known corresponding points/edges/faces.
            scalar dist =
                maxDistance
                (
                    f1,
                    mesh1.points(),
                    mesh0.faces()[mesh0FaceI],
                    mesh0.points()
                );

            if (dist < absTol)
            {
                mesh0Set.insert(mesh0FaceI);
                mesh1Set.insert(mesh1FaceI);
            }
        }
    }

    if (debug)
    {
        Pout<< "findSlavesCoveringMaster :"
            << " matched " << mesh1Set.size() << " mesh1 faces to "
            << mesh0Set.size() << " mesh0 faces" << endl;
    }

    mesh0Faces = mesh0Set.toc();
    mesh1Faces = mesh1Set.toc();
}


// Grow cutToMasterFace across 'internal' edges.
Foam::label Foam::faceCoupleInfo::growCutFaces
(
    const labelList& cutToMasterEdges,
    Map<labelList>& candidates
)
{
    if (debug)
    {
        Pout<< "growCutFaces :"
            << " growing cut faces to masterPatch" << endl;
    }

    label nTotChanged = 0;

    while (true)
    {
        const labelListList& cutFaceEdges = cutFaces().faceEdges();
        const labelListList& cutEdgeFaces = cutFaces().edgeFaces();

        label nChanged = 0;

        forAll(cutToMasterFaces_, cutFaceI)
        {
            const label masterFaceI = cutToMasterFaces_[cutFaceI];

            if (masterFaceI != -1)
            {
                // We now have a cutFace for which we have already found a
                // master face. Grow this masterface across any internal edge
                // (internal: no corresponding master edge)

                const labelList& fEdges = cutFaceEdges[cutFaceI];

                forAll(fEdges, i)
                {
                    const label cutEdgeI = fEdges[i];

                    if (cutToMasterEdges[cutEdgeI] == -1)
                    {
                        // So this edge:
                        // - internal to the cutPatch (since all region edges
                        //   marked before)
                        // - on cutFaceI which corresponds to masterFace.
                        // Mark all connected faces with this masterFace.

                        const labelList& eFaces = cutEdgeFaces[cutEdgeI];

                        forAll(eFaces, j)
                        {
                            const label faceI = eFaces[j];

                            if (cutToMasterFaces_[faceI] == -1)
                            {
                                cutToMasterFaces_[faceI] = masterFaceI;
                                candidates.erase(faceI);
                                nChanged++;
                            }
                            else if (cutToMasterFaces_[faceI] != masterFaceI)
                            {
                                const pointField& cutPoints =
                                    cutFaces().points();
                                const pointField& masterPoints =
                                    masterPatch().points();

                                const edge& e = cutFaces().edges()[cutEdgeI];

                                label myMaster = cutToMasterFaces_[faceI];
                                const face& myF = masterPatch()[myMaster];

                                const face& nbrF = masterPatch()[masterFaceI];

                                FatalErrorIn
                                (
                                    "faceCoupleInfo::growCutFaces"
                                    "(const labelList&, Map<labelList>&)"
                                )   << "Cut face "
                                    << cutFaces()[faceI].points(cutPoints)
                                    << " has master "
                                    << myMaster
                                    << " but also connects to nbr face "
                                    << cutFaces()[cutFaceI].points(cutPoints)
                                    << " with master " << masterFaceI
                                    << nl
                                    << "myMasterFace:"
                                    << myF.points(masterPoints)
                                    << "  nbrMasterFace:"
                                    << nbrF.points(masterPoints) << nl
                                    << "Across cut edge " << cutEdgeI
                                    << " coord:"
                                    << cutFaces().localPoints()[e[0]]
                                    << cutFaces().localPoints()[e[1]]
                                    << abort(FatalError);
                            }
                        }
                    }
                }
            }
        }

        if (debug)
        {
            Pout<< "growCutFaces : Grown an additional " << nChanged
                << " cut to master face" << " correspondence" << endl;
        }

        nTotChanged += nChanged;

        if (nChanged == 0)
        {
            break;
        }
    }

    return nTotChanged;
}


void Foam::faceCoupleInfo::checkMatch(const labelList& cutToMasterEdges) const
{
    const pointField& cutLocalPoints = cutFaces().localPoints();

    const pointField& masterLocalPoints = masterPatch().localPoints();
    const faceList& masterLocalFaces = masterPatch().localFaces();

    forAll(cutToMasterEdges, cutEdgeI)
    {
        const edge& e = cutFaces().edges()[cutEdgeI];

        if (cutToMasterEdges[cutEdgeI] == -1)
        {
            // Internal edge. Check that master face is same on both sides.
            const labelList& cutEFaces = cutFaces().edgeFaces()[cutEdgeI];

            label masterFaceI = -1;

            forAll(cutEFaces, i)
            {
                label cutFaceI = cutEFaces[i];

                if (cutToMasterFaces_[cutFaceI] != -1)
                {
                    if (masterFaceI == -1)
                    {
                        masterFaceI = cutToMasterFaces_[cutFaceI];
                    }
                    else if (masterFaceI != cutToMasterFaces_[cutFaceI])
                    {
                        label myMaster = cutToMasterFaces_[cutFaceI];
                        const face& myF = masterLocalFaces[myMaster];

                        const face& nbrF = masterLocalFaces[masterFaceI];

                        FatalErrorIn
                        (
                            "faceCoupleInfo::checkMatch(const labelList&) const"
                        )
                            << "Internal CutEdge " << e
                            << " coord:"
                            << cutLocalPoints[e[0]]
                            << cutLocalPoints[e[1]]
                            << " connects to master " << myMaster
                            << " and to master " << masterFaceI << nl
                            << "myMasterFace:"
                            << myF.points(masterLocalPoints)
                            << "  nbrMasterFace:"
                            << nbrF.points(masterLocalPoints)
                            << abort(FatalError);
                    }
                }
            }
        }
    }
}


// Extends matching information by elimination across cutFaces using more
// than one region edge. Updates cutToMasterFaces_ and sets candidates
// which is for every cutface on a region edge the possible master faces.
Foam::label Foam::faceCoupleInfo::matchEdgeFaces
(
    const labelList& cutToMasterEdges,
    Map<labelList>& candidates
)
{
    // For every unassigned cutFaceI the possible list of master faces.
    candidates.clear();
    candidates.resize(cutFaces().size());

    label nChanged = 0;

    forAll(cutToMasterEdges, cutEdgeI)
    {
        label masterEdgeI = cutToMasterEdges[cutEdgeI];

        if (masterEdgeI != -1)
        {
            // So cutEdgeI is matched to masterEdgeI. For all cut faces using
            // the cut edge store the master face as a candidate match.

            const labelList& cutEFaces = cutFaces().edgeFaces()[cutEdgeI];
            const labelList& masterEFaces =
                masterPatch().edgeFaces()[masterEdgeI];

            forAll(cutEFaces, i)
            {
                label cutFaceI = cutEFaces[i];

                if (cutToMasterFaces_[cutFaceI] == -1)
                {
                    // So this face (cutFaceI) is on an edge (cutEdgeI) that
                    // is used by some master faces (masterEFaces). Check if
                    // the cutFaceI and some of these masterEFaces share more
                    // than one edge (which uniquely defines face)

                    // Combine master faces with current set of candidate
                    // master faces.
                    Map<labelList>::iterator fnd = candidates.find(cutFaceI);

                    if (fnd == candidates.end())
                    {
                        // No info yet for cutFaceI. Add all master faces as
                        // candidates
                        candidates.insert(cutFaceI, masterEFaces);
                    }
                    else
                    {
                        // From some other cutEdgeI there are already some
                        // candidate master faces. Check the overlap with
                        // the current set of master faces.
                        const labelList& masterFaces = fnd();

                        DynamicList<label> newCandidates(masterFaces.size());

                        forAll(masterEFaces, j)
                        {
                            if (findIndex(masterFaces, masterEFaces[j]) != -1)
                            {
                                newCandidates.append(masterEFaces[j]);
                            }
                        }

                        if (newCandidates.size() == 1)
                        {
                            // We found a perfect match. Delete entry from
                            // candidates map.
                            cutToMasterFaces_[cutFaceI] = newCandidates[0];
                            candidates.erase(cutFaceI);
                            nChanged++;
                        }
                        else
                        {
                            // Should not happen?
                            //Pout<< "On edge:" << cutEdgeI
                            //    << " have connected masterFaces:"
                            //    << masterEFaces
                            //    << " and from previous edge we have"
                            //    << " connected masterFaces" << masterFaces
                            //    << " . Overlap:" << newCandidates << endl;

                            fnd() = newCandidates.shrink();
                        }
                    }
                }

            }
        }
    }

    if (debug)
    {
        Pout<< "matchEdgeFaces : Found " << nChanged
            << " faces where there was"
            << " only one remaining choice for cut-master correspondence"
            << endl;
    }

    return nChanged;
}


// Gets a list of cutFaces (that use a master edge) and the candidate
// master faces.
// Finds most aligned master face.
Foam::label Foam::faceCoupleInfo::geometricMatchEdgeFaces
(
    Map<labelList>& candidates
)
{
    const pointField& cutPoints = cutFaces().points();

    label nChanged = 0;

    // Mark all master faces that have been matched so far.

    labelListList masterToCutFaces
    (
        invertOneToMany
        (
            masterPatch().size(),
            cutToMasterFaces_
        )
    );

    forAllConstIter(Map<labelList>, candidates, iter)
    {
        label cutFaceI = iter.key();

        const face& cutF = cutFaces()[cutFaceI];

        if (cutToMasterFaces_[cutFaceI] == -1)
        {
            const labelList& masterFaces = iter();

            // Find the best matching master face.
            scalar minDist = GREAT;
            label minMasterFaceI = -1;

            forAll(masterFaces, i)
            {
                label masterFaceI = masterFaces[i];

                if (masterToCutFaces[masterFaces[i]].empty())
                {
                    scalar dist = maxDistance
                    (
                        cutF,
                        cutPoints,
                        masterPatch()[masterFaceI],
                        masterPatch().points()
                    );

                    if (dist < minDist)
                    {
                        minDist = dist;
                        minMasterFaceI = masterFaceI;
                    }
                }
            }

            if (minMasterFaceI != -1)
            {
                cutToMasterFaces_[cutFaceI] = minMasterFaceI;
                masterToCutFaces[minMasterFaceI] = cutFaceI;
                nChanged++;
            }
        }
    }

    // (inefficiently) force candidates to be uptodate.
    forAll(cutToMasterFaces_, cutFaceI)
    {
        if (cutToMasterFaces_[cutFaceI] != -1)
        {
            candidates.erase(cutFaceI);
        }
    }


    if (debug)
    {
        Pout<< "geometricMatchEdgeFaces : Found " << nChanged
            << " faces where there was"
            << " only one remaining choice for cut-master correspondence"
            << endl;
    }

    return nChanged;
}


// Calculate the set of cut faces inbetween master and slave patch
// assuming perfect match (and optional face ordering on slave)
void Foam::faceCoupleInfo::perfectPointMatch
(
    const scalar absTol,
    const bool slaveFacesOrdered
)
{
    if (debug)
    {
        Pout<< "perfectPointMatch :"
            << " Matching master and slave to cut."
            << " Master and slave faces are identical" << nl;

        if (slaveFacesOrdered)
        {
            Pout<< "and master and slave faces are ordered"
                << " (on coupled patches)" << endl;
        }
        else
        {
            Pout<< "and master and slave faces are not ordered"
                << " (on coupled patches)" << endl;
        }
    }

    cutToMasterFaces_ = identity(masterPatch().size());
    cutPoints_ = masterPatch().localPoints();
    cutFacesPtr_.reset
    (
        new primitiveFacePatch
        (
            masterPatch().localFaces(),
            cutPoints_
        )
    );
    masterToCutPoints_ = identity(cutPoints_.size());


    // Cut faces to slave patch.
    bool matchedAllFaces = false;

    if (slaveFacesOrdered)
    {
        cutToSlaveFaces_ = identity(cutFaces().size());
        matchedAllFaces = (cutFaces().size() == slavePatch().size());
    }
    else
    {
        // Faces do not have to be ordered (but all have
        // to match). Note: Faces will be already ordered if we enter here from
        // construct from meshes.
        matchedAllFaces = matchPoints
        (
            calcFaceCentres<List>
            (
                cutFaces(),
                cutPoints_,
                0,
                cutFaces().size()
            ),
            calcFaceCentres<IndirectList>
            (
                slavePatch(),
                slavePatch().points(),
                0,
                slavePatch().size()
            ),
            scalarField(slavePatch().size(), absTol),
            true,
            cutToSlaveFaces_
        );
    }


    if (!matchedAllFaces)
    {
        FatalErrorIn
        (
            "faceCoupleInfo::perfectPointMatch"
            "(const scalar, const bool)"
        )   << "Did not match all of the master faces to the slave faces"
            << endl
            << "This usually means that the slave patch and master patch"
            << " do not align to within " << absTol << " metre."
            << abort(FatalError);
    }


    // Find correspondence from slave points to cut points. This might
    // detect shared points so the output is a slave-to-cut point list
    // and a compaction list.

    labelList cutToCompact, compactToCut;
    matchPointsThroughFaces
    (
        absTol,
        cutFaces().localPoints(),
        reorder(cutToSlaveFaces_, cutFaces().localFaces()),
        slavePatch().localPoints(),
        slavePatch().localFaces(),
        false,                      // slave and cut have opposite orientation

        slaveToCutPoints_,          // slave to (uncompacted) cut points
        cutToCompact,               // compaction map: from cut to compacted
        compactToCut                // compaction map: from compacted to cut
    );


    // Use compaction lists to renumber cutPoints.
    cutPoints_ = UIndirectList<point>(cutPoints_, compactToCut)();
    {
        const faceList& cutLocalFaces = cutFaces().localFaces();

        faceList compactFaces(cutLocalFaces.size());
        forAll(cutLocalFaces, i)
        {
            compactFaces[i] = renumber(cutToCompact, cutLocalFaces[i]);
        }
        cutFacesPtr_.reset
        (
            new primitiveFacePatch
            (
                compactFaces,
                cutPoints_
            )
        );
    }
    inplaceRenumber(cutToCompact, slaveToCutPoints_);
    inplaceRenumber(cutToCompact, masterToCutPoints_);
}


// Calculate the set of cut faces inbetween master and slave patch
// assuming that slave patch is subdivision of masterPatch.
void Foam::faceCoupleInfo::subDivisionMatch
(
    const polyMesh& slaveMesh,
    const bool patchDivision,
    const scalar absTol
)
{
    if (debug)
    {
        Pout<< "subDivisionMatch :"
            << " Matching master and slave to cut."
            << " Slave can be subdivision of master but all master edges"
            << " have to be covered by slave edges." << endl;
    }

    // Assume slave patch is subdivision of the master patch so cutFaces is a
    // copy of the slave (but reversed since orientation has to be like master).
    cutPoints_ = slavePatch().localPoints();
    {
        faceList cutFaces(slavePatch().size());

        forAll(cutFaces, i)
        {
            cutFaces[i] = slavePatch().localFaces()[i].reverseFace();
        }
        cutFacesPtr_.reset(new primitiveFacePatch(cutFaces, cutPoints_));
    }

    // Cut is copy of slave so addressing to slave is trivial.
    cutToSlaveFaces_ = identity(cutFaces().size());
    slaveToCutPoints_ = identity(slavePatch().nPoints());

    // Cut edges to slave patch
    labelList cutToSlaveEdges
    (
        findMappedEdges
        (
            cutFaces().edges(),
            slaveToCutPoints_,  //note:should be cutToSlavePoints but since iden
            slavePatch()
        )
    );


    if (debug)
    {
        OFstream str("regionEdges.obj");

        Pout<< "subDivisionMatch :"
            << " addressing for slave patch fully done."
            << " Dumping region edges to " << str.name() << endl;

        label vertI = 0;

        forAll(slavePatch().edges(), slaveEdgeI)
        {
            if (regionEdge(slaveMesh, slaveEdgeI))
            {
                const edge& e = slavePatch().edges()[slaveEdgeI];

                meshTools::writeOBJ(str, slavePatch().localPoints()[e[0]]);
                vertI++;
                meshTools::writeOBJ(str, slavePatch().localPoints()[e[1]]);
                vertI++;
                str<< "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }


    // Addressing from cut to master.

    // Cut points to master patch
    // - determine master points to cut points. Has to be full!
    // - invert to get cut to master
    if (debug)
    {
        Pout<< "subDivisionMatch :"
            << " matching master points to cut points" << endl;
    }

    bool matchedAllPoints = matchPoints
    (
        masterPatch().localPoints(),
        cutPoints_,
        scalarField(masterPatch().nPoints(), absTol),
        true,
        masterToCutPoints_
    );

    if (!matchedAllPoints)
    {
        FatalErrorIn
        (
            "faceCoupleInfo::subDivisionMatch"
            "(const polyMesh&, const bool, const scalar)"
        )   << "Did not match all of the master points to the slave points"
            << endl
            << "This usually means that the slave patch is not a"
            << " subdivision of the master patch"
            << abort(FatalError);
    }


    // Do masterEdges to cutEdges :
    // - mark all edges between two masterEdge endpoints. (geometric test since
    //   this is the only distinction)
    // - this gives the borders inbetween which all cutfaces come from
    //   a single master face.
    if (debug)
    {
        Pout<< "subDivisionMatch :"
            << " matching cut edges to master edges" << endl;
    }

    const edgeList& masterEdges = masterPatch().edges();
    const pointField& masterPoints = masterPatch().localPoints();

    const edgeList& cutEdges = cutFaces().edges();

    labelList cutToMasterEdges(cutFaces().nEdges(), -1);

    forAll(masterEdges, masterEdgeI)
    {
        const edge& masterEdge = masterEdges[masterEdgeI];

        label cutPoint0 = masterToCutPoints_[masterEdge[0]];
        label cutPoint1 = masterToCutPoints_[masterEdge[1]];

        // Find edges between cutPoint0 and cutPoint1.

        label cutPointI = cutPoint0;
        do
        {
            // Find edge (starting at pointI on cut), aligned with master
            // edge.
            label cutEdgeI =
                mostAlignedCutEdge
                (
                    false,
                    slaveMesh,
                    patchDivision,
                    cutToMasterEdges,
                    cutToSlaveEdges,
                    cutPointI,
                    cutPoint0,
                    cutPoint1
                );

            if (cutEdgeI == -1)
            {
                // Problem. Report while matching to produce nice error message
                mostAlignedCutEdge
                (
                    true,
                    slaveMesh,
                    patchDivision,
                    cutToMasterEdges,
                    cutToSlaveEdges,
                    cutPointI,
                    cutPoint0,
                    cutPoint1
                );

                Pout<< "Dumping unmatched pointEdges to errorEdges.obj"
                    << endl;

                writeOBJ
                (
                    "errorEdges.obj",
                    edgeList
                    (
                        UIndirectList<edge>
                        (
                            cutFaces().edges(),
                            cutFaces().pointEdges()[cutPointI]
                        )
                    ),
                    cutFaces().localPoints(),
                    false
                );

                FatalErrorIn
                (
                    "faceCoupleInfo::subDivisionMatch"
                    "(const polyMesh&, const bool, const scalar)"
                )   << "Problem in finding cut edges corresponding to"
                    << " master edge " << masterEdge
                    << " points:" << masterPoints[masterEdge[0]]
                    << ' ' << masterPoints[masterEdge[1]]
                    << " corresponding cut edge: (" << cutPoint0
                    << ") " << cutPoint1
                    << abort(FatalError);
            }

            cutToMasterEdges[cutEdgeI] = masterEdgeI;

            cutPointI = cutEdges[cutEdgeI].otherVertex(cutPointI);

        } while (cutPointI != cutPoint1);
    }

    if (debug)
    {
        // Write all matched edges.
        writeEdges(cutToMasterEdges, cutToSlaveEdges);
    }

    // Rework cutToMasterEdges into list of points inbetween two endpoints
    // (cutEdgeToPoints_)
    setCutEdgeToPoints(cutToMasterEdges);


    // Now that we have marked the cut edges that lie on top of master edges
    // we can find cut faces that have two (or more) of these edges.
    // Note that there can be multiple faces having two or more matched edges
    // since the cut faces can form a non-manifold surface(?)
    // So the matching is done as an elimination process where for every
    // cutFace on a matched edge we store the possible master faces and
    // eliminate more and more until we only have one possible master face
    // left.

    if (debug)
    {
        Pout<< "subDivisionMatch :"
            << " matching (topological) cut faces to masterPatch" << endl;
    }

    // For every unassigned cutFaceI the possible list of master faces.
    Map<labelList> candidates(cutFaces().size());

    cutToMasterFaces_.setSize(cutFaces().size());
    cutToMasterFaces_ = -1;

    while (true)
    {
        // See if there are any edges left where there is only one remaining
        // candidate.
        label nChanged = matchEdgeFaces(cutToMasterEdges, candidates);

        checkMatch(cutToMasterEdges);


        // Now we should have found a face correspondence for every cutFace
        // that uses two or more cutEdges. Floodfill from these across
        // 'internal' cutedges (i.e. edges that do not have a master
        // equivalent) to determine some more cutToMasterFaces
        nChanged += growCutFaces(cutToMasterEdges, candidates);

        checkMatch(cutToMasterEdges);

        if (nChanged == 0)
        {
            break;
        }
    }

    if (debug)
    {
        Pout<< "subDivisionMatch :"
            << " matching (geometric) cut faces to masterPatch" << endl;
    }

    // Above should have matched all cases fully. Cannot prove this yet so
    // do any remaining unmatched edges by a geometric test.
    while (true)
    {
        // See if there are any edges left where there is only one remaining
        // candidate.
        label nChanged = geometricMatchEdgeFaces(candidates);

        checkMatch(cutToMasterEdges);

        nChanged += growCutFaces(cutToMasterEdges, candidates);

        checkMatch(cutToMasterEdges);

        if (nChanged == 0)
        {
            break;
        }
    }


    // All cut faces matched?
    forAll(cutToMasterFaces_, cutFaceI)
    {
        if (cutToMasterFaces_[cutFaceI] == -1)
        {
            const face& cutF = cutFaces()[cutFaceI];

            FatalErrorIn
            (
                "faceCoupleInfo::subDivisionMatch"
                "(const polyMesh&, const bool, const scalar)"
            )   << "Did not match all of cutFaces to a master face" << nl
                << "First unmatched cut face:" << cutFaceI << " with points:"
                << UIndirectList<point>(cutFaces().points(), cutF)()
                << nl
                << "This usually means that the slave patch is not a"
                << " subdivision of the master patch"
                << abort(FatalError);
        }
    }

    if (debug)
    {
        Pout<< "subDivisionMatch :"
            << " finished matching master and slave to cut" << endl;
    }

    if (debug)
    {
        writePointsFaces();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh data
Foam::faceCoupleInfo::faceCoupleInfo
(
    const polyMesh& masterMesh,
    const polyMesh& slaveMesh,
    const scalar absTol,
    const bool perfectMatch
)
:
    masterPatchPtr_(NULL),
    slavePatchPtr_(NULL),
    cutPoints_(0),
    cutFacesPtr_(NULL),
    cutToMasterFaces_(0),
    masterToCutPoints_(0),
    cutToSlaveFaces_(0),
    slaveToCutPoints_(0),
    cutEdgeToPoints_(0)
{
    // Get faces on both meshes that are aligned.
    // (not ordered i.e. masterToMesh[0] does
    // not couple to slaveToMesh[0]. This ordering is done later on)
    labelList masterToMesh;
    labelList slaveToMesh;

    if (perfectMatch)
    {
        // Identical faces so use tight face-centre comparison.
        findPerfectMatchingFaces
        (
            masterMesh,
            slaveMesh,
            absTol,
            masterToMesh,
            slaveToMesh
        );
    }
    else
    {
        // Slave subdivision of master so use 'nearest'. Bit dodgy, especially
        // with using absTol (which is quite small)
        findSlavesCoveringMaster
        (
            masterMesh,
            slaveMesh,
            absTol,
            masterToMesh,
            slaveToMesh
        );
    }

    // Construct addressing engines for both sides
    masterPatchPtr_.reset
    (
        new indirectPrimitivePatch
        (
            IndirectList<face>(masterMesh.faces(), masterToMesh),
            masterMesh.points()
        )
    );

    slavePatchPtr_.reset
    (
        new indirectPrimitivePatch
        (
            IndirectList<face>(slaveMesh.faces(), slaveToMesh),
            slaveMesh.points()
        )
    );


    if (perfectMatch)
    {
        // Faces are perfectly aligned but probably not ordered.
        perfectPointMatch(absTol, false);
    }
    else
    {
        // Slave faces are subdivision of master face. Faces not ordered.
        subDivisionMatch(slaveMesh, false, absTol);
    }

    if (debug)
    {
        writePointsFaces();
    }
}


// Slave is subdivision of master patch.
// (so -both cover the same area -all of master points are present in slave)
Foam::faceCoupleInfo::faceCoupleInfo
(
    const polyMesh& masterMesh,
    const labelList& masterAddressing,
    const polyMesh& slaveMesh,
    const labelList& slaveAddressing,
    const scalar absTol,
    const bool perfectMatch,
    const bool orderedFaces,
    const bool patchDivision
)
:
    masterPatchPtr_
    (
        new indirectPrimitivePatch
        (
            IndirectList<face>(masterMesh.faces(), masterAddressing),
            masterMesh.points()
        )
    ),
    slavePatchPtr_
    (
        new indirectPrimitivePatch
        (
            IndirectList<face>(slaveMesh.faces(), slaveAddressing),
            slaveMesh.points()
        )
    ),
    cutPoints_(0),
    cutFacesPtr_(NULL),
    cutToMasterFaces_(0),
    masterToCutPoints_(0),
    cutToSlaveFaces_(0),
    slaveToCutPoints_(0),
    cutEdgeToPoints_(0)
{
    if (perfectMatch && (masterAddressing.size() != slaveAddressing.size()))
    {
        FatalErrorIn
        (
            "faceCoupleInfo::faceCoupleInfo(const primitiveMesh&"
            ", const primitiveMesh&, const scalar, const bool"
        )   << "Perfect match specified but number of master and slave faces"
            << " differ." << endl
            << "master:" << masterAddressing.size()
            << "  slave:" << slaveAddressing.size()
            << abort(FatalError);
    }

    if
    (
        masterAddressing.size()
     && min(masterAddressing) < masterMesh.nInternalFaces()
    )
    {
        FatalErrorIn
        (
            "faceCoupleInfo::faceCoupleInfo(const primitiveMesh&"
            ", const primitiveMesh&, const scalar, const bool"
        )   << "Supplied internal face on master mesh to couple." << nl
            << "Faces to be coupled have to be boundary faces."
            << abort(FatalError);
    }
    if
    (
        slaveAddressing.size()
     && min(slaveAddressing) < slaveMesh.nInternalFaces()
    )
    {
        FatalErrorIn
        (
            "faceCoupleInfo::faceCoupleInfo(const primitiveMesh&"
            ", const primitiveMesh&, const scalar, const bool"
        )   << "Supplied internal face on slave mesh to couple." << nl
            << "Faces to be coupled have to be boundary faces."
            << abort(FatalError);
    }


    if (perfectMatch)
    {
        perfectPointMatch(absTol, orderedFaces);
    }
    else
    {
        // Slave faces are subdivision of master face. Faces not ordered.
        subDivisionMatch(slaveMesh, patchDivision, absTol);
    }

    if (debug)
    {
        writePointsFaces();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faceCoupleInfo::~faceCoupleInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::faceCoupleInfo::faceLabels(const polyPatch& pp)
{
    labelList faces(pp.size());

    label faceI = pp.start();

    forAll(pp, i)
    {
        faces[i] = faceI++;
    }
    return faces;
}


Foam::Map<Foam::label> Foam::faceCoupleInfo::makeMap(const labelList& lst)
{
    Map<label> map(lst.size());

    forAll(lst, i)
    {
        if (lst[i] != -1)
        {
            map.insert(i, lst[i]);
        }
    }
    return map;
}


Foam::Map<Foam::labelList> Foam::faceCoupleInfo::makeMap
(
    const labelListList& lst
)
{
    Map<labelList> map(lst.size());

    forAll(lst, i)
    {
        if (lst[i].size())
        {
            map.insert(i, lst[i]);
        }
    }
    return map;
}


// ************************************************************************* //
