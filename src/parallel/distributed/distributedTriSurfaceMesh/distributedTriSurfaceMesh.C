/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "distributedTriSurfaceMesh.H"
#include "mapDistribute.H"
#include "Random.H"
#include "addToRunTimeSelectionTable.H"
#include "triangleFuncs.H"
#include "matchPoints.H"
#include "globalIndex.H"
#include "Time.H"

#include "IFstream.H"
#include "decompositionMethod.H"
#include "geomDecomp.H"
#include "vectorList.H"
#include "PackedBoolList.H"
#include "PatchTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(distributedTriSurfaceMesh, 0);
    addToRunTimeSelectionTable
    (
        searchableSurface,
        distributedTriSurfaceMesh,
        dict
    );

    template<>
    const char* Foam::NamedEnum
    <
        Foam::distributedTriSurfaceMesh::distributionType,
        3
    >::names[] =
    {
        "follow",
        "independent",
        "frozen"
    };
}


const Foam::NamedEnum<Foam::distributedTriSurfaceMesh::distributionType, 3>
    Foam::distributedTriSurfaceMesh::distributionTypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Read my additional data from the dictionary
bool Foam::distributedTriSurfaceMesh::read()
{
    // Get bb of all domains.
    procBb_.setSize(Pstream::nProcs());

    procBb_[Pstream::myProcNo()] = List<treeBoundBox>(dict_.lookup("bounds"));
    Pstream::gatherList(procBb_);
    Pstream::scatterList(procBb_);

    // Distribution type
    distType_ = distributionTypeNames_.read(dict_.lookup("distributionType"));

    // Merge distance
    mergeDist_ = dict_.lookup<scalar>("mergeDistance");

    return true;
}


// Is segment fully local?
bool Foam::distributedTriSurfaceMesh::isLocal
(
    const List<treeBoundBox>& myBbs,
    const point& start,
    const point& end
)
{
    forAll(myBbs, bbI)
    {
        if (myBbs[bbI].contains(start) && myBbs[bbI].contains(end))
        {
            return true;
        }
    }
    return false;
}


//void Foam::distributedTriSurfaceMesh::splitSegment
//(
//    const label segmentI,
//    const point& start,
//    const point& end,
//    const treeBoundBox& bb,
//
//    DynamicList<segment>& allSegments,
//    DynamicList<label>& allSegmentMap,
//    DynamicList<label> sendMap
//) const
//{
//    // Work points
//    point clipPt0, clipPt1;
//
//    if (bb.contains(start))
//    {
//        // start within, trim end to bb
//        bool clipped = bb.intersects(end, start, clipPt0);
//
//        if (clipped)
//        {
//            // segment from start to clippedStart passes
//            // through proc.
//            sendMap[proci].append(allSegments.size());
//            allSegmentMap.append(segmentI);
//            allSegments.append(segment(start, clipPt0));
//        }
//    }
//    else if (bb.contains(end))
//    {
//        // end within, trim start to bb
//        bool clipped = bb.intersects(start, end, clipPt0);
//
//        if (clipped)
//        {
//            sendMap[proci].append(allSegments.size());
//            allSegmentMap.append(segmentI);
//            allSegments.append(segment(clipPt0, end));
//        }
//    }
//    else
//    {
//        // trim both
//        bool clippedStart = bb.intersects(start, end, clipPt0);
//
//        if (clippedStart)
//        {
//            bool clippedEnd = bb.intersects(end, clipPt0, clipPt1);
//
//            if (clippedEnd)
//            {
//                // middle part of segment passes through proc.
//                sendMap[proci].append(allSegments.size());
//                allSegmentMap.append(segmentI);
//                allSegments.append(segment(clipPt0, clipPt1));
//            }
//        }
//    }
//}


void Foam::distributedTriSurfaceMesh::distributeSegment
(
    const label segmentI,
    const point& start,
    const point& end,

    DynamicList<segment>& allSegments,
    DynamicList<label>& allSegmentMap,
    List<DynamicList<label>>& sendMap
) const
{
    // 1. Fully local already handled outside. Note: retest is cheap.
    if (isLocal(procBb_[Pstream::myProcNo()], start, end))
    {
        return;
    }


    // 2. If fully inside one other processor, then only need to send
    // to that one processor even if it intersects another. Rare occurrence
    // but cheap to test.
    forAll(procBb_, proci)
    {
        if (proci != Pstream::myProcNo())
        {
            const List<treeBoundBox>& bbs = procBb_[proci];

            if (isLocal(bbs, start, end))
            {
                sendMap[proci].append(allSegments.size());
                allSegmentMap.append(segmentI);
                allSegments.append(segment(start, end));
                return;
            }
        }
    }


    // 3. If not contained in single processor send to all intersecting
    // processors.
    forAll(procBb_, proci)
    {
        const List<treeBoundBox>& bbs = procBb_[proci];

        forAll(bbs, bbI)
        {
            const treeBoundBox& bb = bbs[bbI];

            // Scheme a: any processor that intersects the segment gets
            // the segment.

            // Intersection point
            point clipPt;

            if (bb.intersects(start, end, clipPt))
            {
                sendMap[proci].append(allSegments.size());
                allSegmentMap.append(segmentI);
                allSegments.append(segment(start, end));
            }

            // Alternative: any processor only gets clipped bit of
            // segment. This gives small problems with additional
            // truncation errors.
            // splitSegment
            //(
            //    segmentI,
            //    start,
            //    end,
            //    bb,
            //
            //    allSegments,
            //    allSegmentMap,
            //   sendMap[proci]
            //);
        }
    }
}


Foam::autoPtr<Foam::mapDistribute>
Foam::distributedTriSurfaceMesh::distributeSegments
(
    const pointField& start,
    const pointField& end,

    List<segment>& allSegments,
    labelList& allSegmentMap
) const
{
    // Determine send map
    // ~~~~~~~~~~~~~~~~~~

    labelListList sendMap(Pstream::nProcs());

    {
        // Since intersection test is quite expensive compared to memory
        // allocation we use DynamicList to immediately store the segment
        // in the correct bin.

        // Segments to test
        DynamicList<segment> dynAllSegments(start.size());
        // Original index of segment
        DynamicList<label> dynAllSegmentMap(start.size());
        // Per processor indices into allSegments to send
        List<DynamicList<label>> dynSendMap(Pstream::nProcs());

        forAll(start, segmentI)
        {
            distributeSegment
            (
                segmentI,
                start[segmentI],
                end[segmentI],

                dynAllSegments,
                dynAllSegmentMap,
                dynSendMap
            );
        }

        // Convert dynamicList to labelList
        sendMap.setSize(Pstream::nProcs());
        forAll(sendMap, proci)
        {
            dynSendMap[proci].shrink();
            sendMap[proci].transfer(dynSendMap[proci]);
        }

        allSegments.transfer(dynAllSegments.shrink());
        allSegmentMap.transfer(dynAllSegmentMap.shrink());
    }


    // Send over how many I need to receive.
    labelListList sendSizes(Pstream::nProcs());
    sendSizes[Pstream::myProcNo()].setSize(Pstream::nProcs());
    forAll(sendMap, proci)
    {
        sendSizes[Pstream::myProcNo()][proci] = sendMap[proci].size();
    }
    Pstream::gatherList(sendSizes);
    Pstream::scatterList(sendSizes);


    // Determine order of receiving
    labelListList constructMap(Pstream::nProcs());

    // My local segments first
    constructMap[Pstream::myProcNo()] = identity
    (
        sendMap[Pstream::myProcNo()].size()
    );

    label segmentI = constructMap[Pstream::myProcNo()].size();
    forAll(constructMap, proci)
    {
        if (proci != Pstream::myProcNo())
        {
            // What I need to receive is what other processor is sending to me.
            label nRecv = sendSizes[proci][Pstream::myProcNo()];
            constructMap[proci].setSize(nRecv);

            for (label i = 0; i < nRecv; i++)
            {
                constructMap[proci][i] = segmentI++;
            }
        }
    }

    return autoPtr<mapDistribute>
    (
        new mapDistribute
        (
            segmentI,       // size after construction
            move(sendMap),
            move(constructMap)
        )
    );
}


void Foam::distributedTriSurfaceMesh::findLine
(
    const bool nearestIntersection,
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    const indexedOctree<treeDataTriSurface>& octree = tree();

    // Initialise
    info.setSize(start.size());
    forAll(info, i)
    {
        info[i].setMiss();
    }

    if (!Pstream::parRun())
    {
        forAll(start, i)
        {
            if (nearestIntersection)
            {
                info[i] = octree.findLine(start[i], end[i]);
            }
            else
            {
                info[i] = octree.findLineAny(start[i], end[i]);
            }
        }
    }
    else
    {
        // Important:force synchronised construction of indexing
        const globalIndex& triIndexer = globalTris();


        // Do any local queries
        // ~~~~~~~~~~~~~~~~~~~~

        label nLocal = 0;

        forAll(start, i)
        {
            if (isLocal(procBb_[Pstream::myProcNo()], start[i], end[i]))
            {
                if (nearestIntersection)
                {
                    info[i] = octree.findLine(start[i], end[i]);
                }
                else
                {
                    info[i] = octree.findLineAny(start[i], end[i]);
                }

                if (info[i].hit())
                {
                    info[i].setIndex(triIndexer.toGlobal(info[i].index()));
                }
                nLocal++;
            }
        }


        if
        (
            returnReduce(nLocal, sumOp<label>())
          < returnReduce(start.size(), sumOp<label>())
        )
        {
            // Not all can be resolved locally. Build segments and map,
            // send over segments, do intersections, send back and merge.


            // Construct queries (segments)
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            // Segments to test
            List<segment> allSegments(start.size());
            // Original index of segment
            labelList allSegmentMap(start.size());

            const autoPtr<mapDistribute> mapPtr
            (
                distributeSegments
                (
                    start,
                    end,
                    allSegments,
                    allSegmentMap
                )
            );
            const mapDistribute& map = mapPtr();

            label nOldAllSegments = allSegments.size();


            // Exchange the segments
            // ~~~~~~~~~~~~~~~~~~~~~

            map.distribute(allSegments);


            // Do tests I need to do
            // ~~~~~~~~~~~~~~~~~~~~~

            // Intersections
            List<pointIndexHit> intersections(allSegments.size());

            forAll(allSegments, i)
            {
                if (nearestIntersection)
                {
                    intersections[i] = octree.findLine
                    (
                        allSegments[i].first(),
                        allSegments[i].second()
                    );
                }
                else
                {
                    intersections[i] = octree.findLineAny
                    (
                        allSegments[i].first(),
                        allSegments[i].second()
                    );
                }

                // Convert triangle index to global numbering
                if (intersections[i].hit())
                {
                    intersections[i].setIndex
                    (
                        triIndexer.toGlobal(intersections[i].index())
                    );
                }
            }


            // Exchange the intersections (opposite to segments)
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            map.reverseDistribute(nOldAllSegments, intersections);


            // Extract the hits
            // ~~~~~~~~~~~~~~~~

            forAll(intersections, i)
            {
                const pointIndexHit& allInfo = intersections[i];
                label segmentI = allSegmentMap[i];
                pointIndexHit& hitInfo = info[segmentI];

                if (allInfo.hit())
                {
                    if (!hitInfo.hit())
                    {
                        // No intersection yet so take this one
                        hitInfo = allInfo;
                    }
                    else if (nearestIntersection)
                    {
                        // Nearest intersection
                        if
                        (
                            magSqr(allInfo.hitPoint()-start[segmentI])
                          < magSqr(hitInfo.hitPoint()-start[segmentI])
                        )
                        {
                            hitInfo = allInfo;
                        }
                    }
                }
            }
        }
    }
}


// Exchanges indices to the processor they come from.
// - calculates exchange map
// - uses map to calculate local triangle index
Foam::autoPtr<Foam::mapDistribute>
Foam::distributedTriSurfaceMesh::calcLocalQueries
(
    const List<pointIndexHit>& info,
    labelList& triangleIndex
) const
{
    triangleIndex.setSize(info.size());

    const globalIndex& triIndexer = globalTris();


    // Determine send map
    // ~~~~~~~~~~~~~~~~~~

    // Since determining which processor the query should go to is
    // cheap we do a multi-pass algorithm to save some memory temporarily.

    // 1. Count
    labelList nSend(Pstream::nProcs(), 0);

    forAll(info, i)
    {
        if (info[i].hit())
        {
            label proci = triIndexer.whichProcID(info[i].index());
            nSend[proci]++;
        }
    }

    // 2. Size sendMap
    labelListList sendMap(Pstream::nProcs());
    forAll(nSend, proci)
    {
        sendMap[proci].setSize(nSend[proci]);
        nSend[proci] = 0;
    }

    // 3. Fill sendMap
    forAll(info, i)
    {
        if (info[i].hit())
        {
            label proci = triIndexer.whichProcID(info[i].index());
            triangleIndex[i] = triIndexer.toLocal(proci, info[i].index());
            sendMap[proci][nSend[proci]++] = i;
        }
        else
        {
            triangleIndex[i] = -1;
        }
    }


    // Send over how many I need to receive
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelListList sendSizes(Pstream::nProcs());
    sendSizes[Pstream::myProcNo()].setSize(Pstream::nProcs());
    forAll(sendMap, proci)
    {
        sendSizes[Pstream::myProcNo()][proci] = sendMap[proci].size();
    }
    Pstream::gatherList(sendSizes);
    Pstream::scatterList(sendSizes);


    // Determine receive map
    // ~~~~~~~~~~~~~~~~~~~~~

    labelListList constructMap(Pstream::nProcs());

    // My local segments first
    constructMap[Pstream::myProcNo()] = identity
    (
        sendMap[Pstream::myProcNo()].size()
    );

    label segmentI = constructMap[Pstream::myProcNo()].size();
    forAll(constructMap, proci)
    {
        if (proci != Pstream::myProcNo())
        {
            // What I need to receive is what other processor is sending to me.
            label nRecv = sendSizes[proci][Pstream::myProcNo()];
            constructMap[proci].setSize(nRecv);

            for (label i = 0; i < nRecv; i++)
            {
                constructMap[proci][i] = segmentI++;
            }
        }
    }


    // Pack into distribution map
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    autoPtr<mapDistribute> mapPtr
    (
        new mapDistribute
        (
            segmentI,       // size after construction
            move(sendMap),
            move(constructMap)
        )
    );
    const mapDistribute& map = mapPtr();


    // Send over queries
    // ~~~~~~~~~~~~~~~~~

    map.distribute(triangleIndex);


    return mapPtr;
}


Foam::label Foam::distributedTriSurfaceMesh::calcOverlappingProcs
(
    const point& centre,
    const scalar radiusSqr,
    boolList& overlaps
) const
{
    overlaps = false;
    label nOverlaps = 0;

    forAll(procBb_, proci)
    {
        const List<treeBoundBox>& bbs = procBb_[proci];

        forAll(bbs, bbI)
        {
            if (bbs[bbI].overlaps(centre, radiusSqr))
            {
                overlaps[proci] = true;
                nOverlaps++;
                break;
            }
        }
    }
    return nOverlaps;
}


// Generate queries for parallel distance calculation
// - calculates exchange map
// - uses map to exchange points and radius
Foam::autoPtr<Foam::mapDistribute>
Foam::distributedTriSurfaceMesh::calcLocalQueries
(
    const pointField& centres,
    const scalarField& radiusSqr,

    pointField& allCentres,
    scalarField& allRadiusSqr,
    labelList& allSegmentMap
) const
{
    // Determine queries
    // ~~~~~~~~~~~~~~~~~

    labelListList sendMap(Pstream::nProcs());

    {
        // Queries
        DynamicList<point> dynAllCentres(centres.size());
        DynamicList<scalar> dynAllRadiusSqr(centres.size());
        // Original index of segment
        DynamicList<label> dynAllSegmentMap(centres.size());
        // Per processor indices into allSegments to send
        List<DynamicList<label>> dynSendMap(Pstream::nProcs());

        // Work array - whether processor bb overlaps the bounding sphere.
        boolList procBbOverlaps(Pstream::nProcs());

        forAll(centres, centreI)
        {
            // Find the processor this sample+radius overlaps.
            calcOverlappingProcs
            (
                centres[centreI],
                radiusSqr[centreI],
                procBbOverlaps
            );

            forAll(procBbOverlaps, proci)
            {
                if (proci != Pstream::myProcNo() && procBbOverlaps[proci])
                {
                    dynSendMap[proci].append(dynAllCentres.size());
                    dynAllSegmentMap.append(centreI);
                    dynAllCentres.append(centres[centreI]);
                    dynAllRadiusSqr.append(radiusSqr[centreI]);
                }
            }
        }

        // Convert dynamicList to labelList
        sendMap.setSize(Pstream::nProcs());
        forAll(sendMap, proci)
        {
            dynSendMap[proci].shrink();
            sendMap[proci].transfer(dynSendMap[proci]);
        }

        allCentres.transfer(dynAllCentres.shrink());
        allRadiusSqr.transfer(dynAllRadiusSqr.shrink());
        allSegmentMap.transfer(dynAllSegmentMap.shrink());
    }


    // Send over how many I need to receive.
    labelListList sendSizes(Pstream::nProcs());
    sendSizes[Pstream::myProcNo()].setSize(Pstream::nProcs());
    forAll(sendMap, proci)
    {
        sendSizes[Pstream::myProcNo()][proci] = sendMap[proci].size();
    }
    Pstream::gatherList(sendSizes);
    Pstream::scatterList(sendSizes);


    // Determine order of receiving
    labelListList constructMap(Pstream::nProcs());

    // My local segments first
    constructMap[Pstream::myProcNo()] = identity
    (
        sendMap[Pstream::myProcNo()].size()
    );

    label segmentI = constructMap[Pstream::myProcNo()].size();
    forAll(constructMap, proci)
    {
        if (proci != Pstream::myProcNo())
        {
            // What I need to receive is what other processor is sending to me.
            label nRecv = sendSizes[proci][Pstream::myProcNo()];
            constructMap[proci].setSize(nRecv);

            for (label i = 0; i < nRecv; i++)
            {
                constructMap[proci][i] = segmentI++;
            }
        }
    }

    autoPtr<mapDistribute> mapPtr
    (
        new mapDistribute
        (
            segmentI,       // size after construction
            move(sendMap),
            move(constructMap)
        )
    );
    return mapPtr;
}


// Find bounding boxes that guarantee a more or less uniform distribution
// of triangles. Decomposition in here is only used to get the bounding
// boxes, actual decomposition is done later on.
// Returns a per processor a list of bounding boxes that most accurately
// describe the shape. For now just a single bounding box per processor but
// optimisation might be to determine a better fitting shape.
Foam::List<Foam::List<Foam::treeBoundBox>>
Foam::distributedTriSurfaceMesh::independentlyDistributedBbs
(
    const triSurface& s
)
{
    if (!decomposer_.valid())
    {
        // Use current decomposer.
        // Note: or always use hierarchical?
        IOdictionary decomposeDict
        (
            IOobject
            (
                "decomposeParDict",
                searchableSurface::time().system(),
                searchableSurface::time(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        );
        decomposer_ = decompositionMethod::New(decomposeDict);

        if (!decomposer_().parallelAware())
        {
            FatalErrorInFunction
                << "The decomposition method " << decomposer_().typeName
                << " does not decompose in parallel."
                << " Please choose one that does." << exit(FatalError);
        }

        if (!isA<geomDecomp>(decomposer_()))
        {
            FatalErrorInFunction
                << "The decomposition method " << decomposer_().typeName
                << " is not a geometric decomposition method." << endl
                << "Only geometric decomposition methods are currently"
                << " supported."
                << exit(FatalError);
        }
    }

    // Do decomposition according to triangle centre
    pointField triCentres(s.size());
    forAll(s, triI)
    {
        triCentres[triI] = s[triI].centre(s.points());
    }


    geomDecomp& decomposer = refCast<geomDecomp>(decomposer_());

    // Do the actual decomposition
    labelList distribution(decomposer.decompose(triCentres));

    // Find bounding box for all triangles on new distribution.

    // Initialise to inverted box (vGreat, -vGreat)
    List<List<treeBoundBox>> bbs(Pstream::nProcs());
    forAll(bbs, proci)
    {
        bbs[proci].setSize(1);
        // bbs[proci][0] = boundBox::invertedBox;
        bbs[proci][0].min() = point( vGreat,  vGreat,  vGreat);
        bbs[proci][0].max() = point(-vGreat, -vGreat, -vGreat);
    }

    forAll(s, triI)
    {
        point& bbMin = bbs[distribution[triI]][0].min();
        point& bbMax = bbs[distribution[triI]][0].max();

        const triSurface::FaceType& f = s[triI];
        forAll(f, fp)
        {
            const point& pt = s.points()[f[fp]];
            bbMin = ::Foam::min(bbMin, pt);
            bbMax = ::Foam::max(bbMax, pt);
        }
    }

    // Now combine for all processors and convert to correct format.
    forAll(bbs, proci)
    {
        forAll(bbs[proci], i)
        {
            reduce(bbs[proci][i].min(), minOp<point>());
            reduce(bbs[proci][i].max(), maxOp<point>());
        }
    }
    return bbs;
}


// Does any part of triangle overlap bb.
bool Foam::distributedTriSurfaceMesh::overlaps
(
    const List<treeBoundBox>& bbs,
    const point& p0,
    const point& p1,
    const point& p2
)
{
    forAll(bbs, bbI)
    {
        const treeBoundBox& bb = bbs[bbI];

        treeBoundBox triBb(p0, p0);
        triBb.min() = min(triBb.min(), p1);
        triBb.min() = min(triBb.min(), p2);

        triBb.max() = max(triBb.max(), p1);
        triBb.max() = max(triBb.max(), p2);

        // Exact test of triangle intersecting bb

        // Quick rejection. If whole bounding box of tri is outside cubeBb then
        // there will be no intersection.
        if (bb.overlaps(triBb))
        {
            // Check if one or more triangle point inside
            if (bb.contains(p0) || bb.contains(p1) || bb.contains(p2))
            {
                // One or more points inside
                return true;
            }

            // Now we have the difficult case: all points are outside but
            // connecting edges might go through cube. Use fast intersection
            // of bounding box.
            bool intersect = triangleFuncs::intersectBb(p0, p1, p2, bb);

            if (intersect)
            {
                return true;
            }
        }
    }
    return false;
}


void Foam::distributedTriSurfaceMesh::subsetMeshMap
(
    const triSurface& s,
    const boolList& include,
    const label nIncluded,
    labelList& newToOldPoints,
    labelList& oldToNewPoints,
    labelList& newToOldFaces
)
{
    newToOldFaces.setSize(nIncluded);
    newToOldPoints.setSize(s.points().size());
    oldToNewPoints.setSize(s.points().size());
    oldToNewPoints = -1;
    {
        label facei = 0;
        label pointi = 0;

        forAll(include, oldFacei)
        {
            if (include[oldFacei])
            {
                // Store new faces compact
                newToOldFaces[facei++] = oldFacei;

                // Renumber labels for face
                const triSurface::FaceType& f = s[oldFacei];

                forAll(f, fp)
                {
                    label oldPointi = f[fp];

                    if (oldToNewPoints[oldPointi] == -1)
                    {
                        oldToNewPoints[oldPointi] = pointi;
                        newToOldPoints[pointi++] = oldPointi;
                    }
                }
            }
        }
        newToOldPoints.setSize(pointi);
    }
}


Foam::triSurface Foam::distributedTriSurfaceMesh::subsetMesh
(
    const triSurface& s,
    const labelList& newToOldPoints,
    const labelList& oldToNewPoints,
    const labelList& newToOldFaces
)
{
    // Extract points
    pointField newPoints(newToOldPoints.size());
    forAll(newToOldPoints, i)
    {
        newPoints[i] = s.points()[newToOldPoints[i]];
    }
    // Extract faces
    List<labelledTri> newTriangles(newToOldFaces.size());

    forAll(newToOldFaces, i)
    {
        // Get old vertex labels
        const labelledTri& tri = s[newToOldFaces[i]];

        newTriangles[i][0] = oldToNewPoints[tri[0]];
        newTriangles[i][1] = oldToNewPoints[tri[1]];
        newTriangles[i][2] = oldToNewPoints[tri[2]];
        newTriangles[i].region() = tri.region();
    }

    // Reuse storage.
    return triSurface(newTriangles, s.patches(), newPoints, true);
}


Foam::triSurface Foam::distributedTriSurfaceMesh::subsetMesh
(
    const triSurface& s,
    const boolList& include,
    labelList& newToOldPoints,
    labelList& newToOldFaces
)
{
    label n = 0;

    forAll(include, i)
    {
        if (include[i])
        {
            n++;
        }
    }

    labelList oldToNewPoints;
    subsetMeshMap
    (
        s,
        include,
        n,
        newToOldPoints,
        oldToNewPoints,
        newToOldFaces
    );

    return subsetMesh
    (
        s,
        newToOldPoints,
        oldToNewPoints,
        newToOldFaces
    );
}


Foam::triSurface Foam::distributedTriSurfaceMesh::subsetMesh
(
    const triSurface& s,
    const labelList& newToOldFaces,
    labelList& newToOldPoints
)
{
    const boolList include
    (
        createWithValues<boolList>
        (
            s.size(),
            false,
            newToOldFaces,
            true
        )
    );

    newToOldPoints.setSize(s.points().size());
    labelList oldToNewPoints(s.points().size(), -1);
    {
        label pointi = 0;

        forAll(include, oldFacei)
        {
            if (include[oldFacei])
            {
                // Renumber labels for face
                const triSurface::FaceType& f = s[oldFacei];

                forAll(f, fp)
                {
                    label oldPointi = f[fp];

                    if (oldToNewPoints[oldPointi] == -1)
                    {
                        oldToNewPoints[oldPointi] = pointi;
                        newToOldPoints[pointi++] = oldPointi;
                    }
                }
            }
        }
        newToOldPoints.setSize(pointi);
    }

    return subsetMesh
    (
        s,
        newToOldPoints,
        oldToNewPoints,
        newToOldFaces
    );
}


Foam::label Foam::distributedTriSurfaceMesh::findTriangle
(
    const List<labelledTri>& allFaces,
    const labelListList& allPointFaces,
    const labelledTri& otherF
)
{
    // allFaces connected to otherF[0]
    const labelList& pFaces = allPointFaces[otherF[0]];

    forAll(pFaces, i)
    {
        const labelledTri& f = allFaces[pFaces[i]];

        if (f.region() == otherF.region())
        {
            // Find index of otherF[0]
            label fp0 = findIndex(f, otherF[0]);
            // Check rest of triangle in same order
            label fp1 = f.fcIndex(fp0);
            label fp2 = f.fcIndex(fp1);

            if (f[fp1] == otherF[1] && f[fp2] == otherF[2])
            {
                return pFaces[i];
            }
        }
    }
    return -1;
}


// Merge into allSurf.
void Foam::distributedTriSurfaceMesh::merge
(
    const scalar mergeDist,
    const List<labelledTri>& subTris,
    const pointField& subPoints,

    List<labelledTri>& allTris,
    pointField& allPoints,

    labelList& faceConstructMap,
    labelList& pointConstructMap
)
{
    labelList subToAll;
    matchPoints
    (
        subPoints,
        allPoints,
        scalarField(subPoints.size(), mergeDist),   // match distance
        false,                                      // verbose
        pointConstructMap
    );

    label nOldAllPoints = allPoints.size();


    // Add all unmatched points
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    label allPointi = nOldAllPoints;
    forAll(pointConstructMap, pointi)
    {
        if (pointConstructMap[pointi] == -1)
        {
            pointConstructMap[pointi] = allPointi++;
        }
    }

    if (allPointi > nOldAllPoints)
    {
        allPoints.setSize(allPointi);

        forAll(pointConstructMap, pointi)
        {
            if (pointConstructMap[pointi] >= nOldAllPoints)
            {
                allPoints[pointConstructMap[pointi]] = subPoints[pointi];
            }
        }
    }


    // To check whether triangles are same we use pointFaces.
    labelListList allPointFaces;
    invertManyToMany(nOldAllPoints, allTris, allPointFaces);


    // Add all unmatched triangles
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    label allTriI = allTris.size();
    allTris.setSize(allTriI + subTris.size());

    faceConstructMap.setSize(subTris.size());

    forAll(subTris, triI)
    {
        const labelledTri& subTri = subTris[triI];

        // Get triangle in new numbering
        labelledTri mappedTri
        (
            pointConstructMap[subTri[0]],
            pointConstructMap[subTri[1]],
            pointConstructMap[subTri[2]],
            subTri.region()
        );


        // Check if all points were already in surface
        bool fullMatch = true;

        forAll(mappedTri, fp)
        {
            if (mappedTri[fp] >= nOldAllPoints)
            {
                fullMatch = false;
                break;
            }
        }

        if (fullMatch)
        {
            // All three points are mapped to old points. See if same
            // triangle.
            label i = findTriangle
            (
                allTris,
                allPointFaces,
                mappedTri
            );

            if (i == -1)
            {
                // Add
                faceConstructMap[triI] = allTriI;
                allTris[allTriI] = mappedTri;
                allTriI++;
            }
            else
            {
                faceConstructMap[triI] = i;
            }
        }
        else
        {
            // Add
            faceConstructMap[triI] = allTriI;
            allTris[allTriI] = mappedTri;
            allTriI++;
        }
    }
    allTris.setSize(allTriI);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributedTriSurfaceMesh::distributedTriSurfaceMesh
(
    const IOobject& io,
    const triSurface& s,
    const dictionary& dict
)
:
    triSurfaceMesh(io, s),
    dict_
    (
        IOobject
        (
            searchableSurface::name() + "Dict",
            searchableSurface::instance(),
            searchableSurface::local(),
            searchableSurface::db(),
            searchableSurface::NO_READ,
            searchableSurface::writeOpt(),
            searchableSurface::registerObject()
        ),
        dict
    )
{
    read();

    reduce(bounds().min(), minOp<point>());
    reduce(bounds().max(), maxOp<point>());

    if (debug)
    {
        InfoInFunction << "Constructed from triSurface:" << endl;
        writeStats(Info);

        labelList nTris(Pstream::nProcs());
        nTris[Pstream::myProcNo()] = triSurface::size();
        Pstream::gatherList(nTris);
        Pstream::scatterList(nTris);

        Info<< endl<< "\tproc\ttris\tbb" << endl;
        forAll(nTris, proci)
        {
            Info<< '\t' << proci << '\t' << nTris[proci]
                 << '\t' << procBb_[proci] << endl;
        }
        Info<< endl;
    }
}


Foam::distributedTriSurfaceMesh::distributedTriSurfaceMesh(const IOobject& io)
:
    triSurfaceMesh
    (
        IOobject
        (
            io.name(),
            io.time().findInstance(io.local(), word::null),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            io.registerObject()
        ),
        false
    ),
    dict_
    (
        IOobject
        (
            searchableSurface::name() + "Dict",
            searchableSurface::instance(),
            searchableSurface::local(),
            searchableSurface::db(),
            searchableSurface::readOpt(),
            searchableSurface::writeOpt(),
            searchableSurface::registerObject()
        )
    )
{
    if
    (
        Pstream::parRun()
     && (
            dict_.readOpt() == IOobject::MUST_READ
         || dict_.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     && (
            regIOobject::fileModificationChecking == timeStampMaster
         || regIOobject::fileModificationChecking == inotifyMaster
        )
    )
    {
        FatalErrorInFunction
            << "    using 'timeStampMaster' or 'inotifyMaster.'\n"
            << "    Modify the entry fileModificationChecking\n"
            << "    in the etc/controlDict.\n"
            << "    Use 'timeStamp' instead."
            << exit(FatalError);
    }

    read();

    reduce(bounds().min(), minOp<point>());
    reduce(bounds().max(), maxOp<point>());

    if (debug)
    {
        InfoInFunction << "Read distributedTriSurface from " << io.objectPath()
            << ':' << endl;
        writeStats(Info);

        labelList nTris(Pstream::nProcs());
        nTris[Pstream::myProcNo()] = triSurface::size();
        Pstream::gatherList(nTris);
        Pstream::scatterList(nTris);

        Info<< endl<< "\tproc\ttris\tbb" << endl;
        forAll(nTris, proci)
        {
            Info<< '\t' << proci << '\t' << nTris[proci]
                 << '\t' << procBb_[proci] << endl;
        }
        Info<< endl;
    }
}


Foam::distributedTriSurfaceMesh::distributedTriSurfaceMesh
(
    const IOobject& io,
    const dictionary& dict
)
:
    // triSurfaceMesh(io, dict),
    triSurfaceMesh
    (
        IOobject
        (
            io.name(),
            io.time().findInstance(io.local(), word::null),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            io.registerObject()
        ),
        dict,
        false
    ),
    dict_
    (
        IOobject
        (
            searchableSurface::name() + "Dict",
            searchableSurface::instance(),
            searchableSurface::local(),
            searchableSurface::db(),
            searchableSurface::readOpt(),
            searchableSurface::writeOpt(),
            searchableSurface::registerObject()
        )
    )
{
    if
    (
        Pstream::parRun()
     && (
            dict_.readOpt() == IOobject::MUST_READ
         || dict_.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     && (
            regIOobject::fileModificationChecking == timeStampMaster
         || regIOobject::fileModificationChecking == inotifyMaster
        )
    )
    {
        FatalErrorInFunction
            << "    using 'timeStampMaster' or 'inotifyMaster.'\n"
            << "    Modify the entry fileModificationChecking\n"
            << "    in the etc/controlDict.\n"
            << "    Use 'timeStamp' instead."
            << exit(FatalError);
    }

    read();

    reduce(bounds().min(), minOp<point>());
    reduce(bounds().max(), maxOp<point>());

    if (debug)
    {
        InfoInFunction << "Read distributedTriSurface from " << io.objectPath()
            << " and dictionary:" << endl;
        writeStats(Info);

        labelList nTris(Pstream::nProcs());
        nTris[Pstream::myProcNo()] = triSurface::size();
        Pstream::gatherList(nTris);
        Pstream::scatterList(nTris);

        Info<< endl<< "\tproc\ttris\tbb" << endl;
        forAll(nTris, proci)
        {
            Info<< '\t' << proci << '\t' << nTris[proci]
                 << '\t' << procBb_[proci] << endl;
        }
        Info<< endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributedTriSurfaceMesh::~distributedTriSurfaceMesh()
{
    clearOut();
}


void Foam::distributedTriSurfaceMesh::clearOut()
{
    globalTris_.clear();
    triSurfaceMesh::clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::globalIndex& Foam::distributedTriSurfaceMesh::globalTris() const
{
    if (!globalTris_.valid())
    {
        globalTris_.reset(new globalIndex(triSurface::size()));
    }
    return globalTris_;
}


void Foam::distributedTriSurfaceMesh::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    const indexedOctree<treeDataTriSurface>& octree = tree();

    // Important:force synchronised construction of indexing
    const globalIndex& triIndexer = globalTris();


    // Initialise
    // ~~~~~~~~~~

    info.setSize(samples.size());
    forAll(info, i)
    {
        info[i].setMiss();
    }



    // Do any local queries
    // ~~~~~~~~~~~~~~~~~~~~

    label nLocal = 0;

    {
        // Work array - whether processor bb overlaps the bounding sphere.
        boolList procBbOverlaps(Pstream::nProcs());

        forAll(samples, i)
        {
            // Find the processor this sample+radius overlaps.
            label nProcs = calcOverlappingProcs
            (
                samples[i],
                nearestDistSqr[i],
                procBbOverlaps
            );

            // Overlaps local processor?
            if (procBbOverlaps[Pstream::myProcNo()])
            {
                info[i] = octree.findNearest(samples[i], nearestDistSqr[i]);
                if (info[i].hit())
                {
                    info[i].setIndex(triIndexer.toGlobal(info[i].index()));
                }
                if (nProcs == 1)
                {
                    // Fully local
                    nLocal++;
                }
            }
        }
    }


    if
    (
        Pstream::parRun()
     && (
            returnReduce(nLocal, sumOp<label>())
          < returnReduce(samples.size(), sumOp<label>())
        )
    )
    {
        // Not all can be resolved locally. Build queries and map, send over
        // queries, do intersections, send back and merge.

        // Calculate queries and exchange map
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        pointField allCentres;
        scalarField allRadiusSqr;
        labelList allSegmentMap;
        autoPtr<mapDistribute> mapPtr
        (
            calcLocalQueries
            (
                samples,
                nearestDistSqr,

                allCentres,
                allRadiusSqr,
                allSegmentMap
            )
        );
        const mapDistribute& map = mapPtr();


        // swap samples to local processor
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        map.distribute(allCentres);
        map.distribute(allRadiusSqr);


        // Do my tests
        // ~~~~~~~~~~~

        List<pointIndexHit> allInfo(allCentres.size());
        forAll(allInfo, i)
        {
            allInfo[i] = octree.findNearest
            (
                allCentres[i],
                allRadiusSqr[i]
            );
            if (allInfo[i].hit())
            {
                allInfo[i].setIndex(triIndexer.toGlobal(allInfo[i].index()));
            }
        }


        // Send back results
        // ~~~~~~~~~~~~~~~~~

        map.reverseDistribute(allSegmentMap.size(), allInfo);


        // Extract information
        // ~~~~~~~~~~~~~~~~~~~

        forAll(allInfo, i)
        {
            if (allInfo[i].hit())
            {
                label pointi = allSegmentMap[i];

                if (!info[pointi].hit())
                {
                    // No intersection yet so take this one
                    info[pointi] = allInfo[i];
                }
                else
                {
                    // Nearest intersection
                    if
                    (
                        magSqr(allInfo[i].hitPoint()-samples[pointi])
                      < magSqr(info[pointi].hitPoint()-samples[pointi])
                    )
                    {
                        info[pointi] = allInfo[i];
                    }
                }
            }
        }
    }
}


void Foam::distributedTriSurfaceMesh::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    findLine
    (
        true,   // nearestIntersection
        start,
        end,
        info
    );
}


void Foam::distributedTriSurfaceMesh::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    findLine
    (
        true,   // nearestIntersection
        start,
        end,
        info
    );
}


void Foam::distributedTriSurfaceMesh::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit>>& info
) const
{
    // Reuse fineLine. We could modify all of findLine to do multiple
    // intersections but this would mean a lot of data transferred so
    // for now we just find nearest intersection and retest from that
    // intersection onwards.

    // Work array.
    List<pointIndexHit> hitInfo(start.size());

    findLine
    (
        true,   // nearestIntersection
        start,
        end,
        hitInfo
    );

    // Tolerances:
    // To find all intersections we add a small vector to the last intersection
    // This is chosen such that
    // - it is significant (small is smallest representative relative tolerance;
    //   we need something bigger since we're doing calculations)
    // - if the start-end vector is zero we still progress
    const vectorField dirVec(end-start);
    const scalarField magSqrDirVec(magSqr(dirVec));
    const vectorField smallVec
    (
        rootSmall*dirVec
      + vector(rootVSmall,rootVSmall,rootVSmall)
    );

    // Copy to input and compact any hits
    labelList pointMap(start.size());
    pointField e0(start.size());
    pointField e1(start.size());
    label compactI = 0;

    info.setSize(hitInfo.size());
    forAll(hitInfo, pointi)
    {
        if (hitInfo[pointi].hit())
        {
            info[pointi].setSize(1);
            info[pointi][0] = hitInfo[pointi];

            point pt = hitInfo[pointi].hitPoint() + smallVec[pointi];

            if (((pt-start[pointi])&dirVec[pointi]) <= magSqrDirVec[pointi])
            {
                e0[compactI] = pt;
                e1[compactI] = end[pointi];
                pointMap[compactI] = pointi;
                compactI++;
            }
        }
        else
        {
            info[pointi].clear();
        }
    }

    e0.setSize(compactI);
    e1.setSize(compactI);
    pointMap.setSize(compactI);

    while (returnReduce(e0.size(), sumOp<label>()) > 0)
    {
        findLine
        (
            true,   // nearestIntersection
            e0,
            e1,
            hitInfo
        );


        // Extract
        label compactI = 0;
        forAll(hitInfo, i)
        {
            if (hitInfo[i].hit())
            {
                label pointi = pointMap[i];

                label sz = info[pointi].size();
                info[pointi].setSize(sz+1);
                info[pointi][sz] = hitInfo[i];

                point pt = hitInfo[i].hitPoint() + smallVec[pointi];

                if (((pt-start[pointi])&dirVec[pointi]) <= magSqrDirVec[pointi])
                {
                    e0[compactI] = pt;
                    e1[compactI] = end[pointi];
                    pointMap[compactI] = pointi;
                    compactI++;
                }
            }
        }

        // Trim
        e0.setSize(compactI);
        e1.setSize(compactI);
        pointMap.setSize(compactI);
    }
}


void Foam::distributedTriSurfaceMesh::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    if (!Pstream::parRun())
    {
        region.setSize(info.size());
        forAll(info, i)
        {
            if (info[i].hit())
            {
                region[i] = triSurface::operator[](info[i].index()).region();
            }
            else
            {
                region[i] = -1;
            }
        }
        return;
    }


    // Get query data (= local index of triangle)
    // ~~~~~~~~~~~~~~

    labelList triangleIndex(info.size());
    autoPtr<mapDistribute> mapPtr
    (
        calcLocalQueries
        (
            info,
            triangleIndex
        )
    );
    const mapDistribute& map = mapPtr();


    // Do my tests
    // ~~~~~~~~~~~

    const triSurface& s = static_cast<const triSurface&>(*this);

    region.setSize(triangleIndex.size());

    forAll(triangleIndex, i)
    {
        label triI = triangleIndex[i];
        region[i] = s[triI].region();
    }


    // Send back results
    // ~~~~~~~~~~~~~~~~~

    map.reverseDistribute(info.size(), region);
}


void Foam::distributedTriSurfaceMesh::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    if (!Pstream::parRun())
    {
        triSurfaceMesh::getNormal(info, normal);
        return;
    }


    // Get query data (= local index of triangle)
    // ~~~~~~~~~~~~~~

    labelList triangleIndex(info.size());
    autoPtr<mapDistribute> mapPtr
    (
        calcLocalQueries
        (
            info,
            triangleIndex
        )
    );
    const mapDistribute& map = mapPtr();


    // Do my tests
    // ~~~~~~~~~~~

    const triSurface& s = static_cast<const triSurface&>(*this);

    normal.setSize(triangleIndex.size());

    forAll(triangleIndex, i)
    {
        label triI = triangleIndex[i];
        normal[i] = s[triI].normal(s.points());
    }


    // Send back results
    // ~~~~~~~~~~~~~~~~~

    map.reverseDistribute(info.size(), normal);
}


void Foam::distributedTriSurfaceMesh::getField
(
    const List<pointIndexHit>& info,
    labelList& values
) const
{
    if (!Pstream::parRun())
    {
        triSurfaceMesh::getField(info, values);
        return;
    }

    if (foundObject<triSurfaceLabelField>("values"))
    {
        const triSurfaceLabelField& fld = lookupObject<triSurfaceLabelField>
        (
            "values"
        );


        // Get query data (= local index of triangle)
        // ~~~~~~~~~~~~~~

        labelList triangleIndex(info.size());
        autoPtr<mapDistribute> mapPtr
        (
            calcLocalQueries
            (
                info,
                triangleIndex
            )
        );
        const mapDistribute& map = mapPtr();


        // Do my tests
        // ~~~~~~~~~~~

        values.setSize(triangleIndex.size());

        forAll(triangleIndex, i)
        {
            label triI = triangleIndex[i];
            values[i] = fld[triI];
        }


        // Send back results
        // ~~~~~~~~~~~~~~~~~

        map.reverseDistribute(info.size(), values);
    }
}


void Foam::distributedTriSurfaceMesh::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    FatalErrorInFunction
        << "Volume type not supported for distributed surfaces."
        << exit(FatalError);
}


// Subset the part of surface that is overlapping bb.
Foam::triSurface Foam::distributedTriSurfaceMesh::overlappingSurface
(
    const triSurface& s,
    const List<treeBoundBox>& bbs,

    labelList& subPointMap,
    labelList& subFaceMap
)
{
    // Determine what triangles to keep.
    boolList includedFace(s.size(), false);

    // Create a slightly larger bounding box.
    List<treeBoundBox> bbsX(bbs.size());
    const scalar eps = 1.0e-4;
    forAll(bbs, i)
    {
        const point mid = 0.5*(bbs[i].min() + bbs[i].max());
        const vector halfSpan = (1.0+eps)*(bbs[i].max() - mid);

        bbsX[i].min() = mid - halfSpan;
        bbsX[i].max() = mid + halfSpan;
    }

    forAll(s, triI)
    {
        const labelledTri& f = s[triI];
        const point& p0 = s.points()[f[0]];
        const point& p1 = s.points()[f[1]];
        const point& p2 = s.points()[f[2]];

        if (overlaps(bbsX, p0, p1, p2))
        {
            includedFace[triI] = true;
        }
    }

    return subsetMesh(s, includedFace, subPointMap, subFaceMap);
}


void Foam::distributedTriSurfaceMesh::distribute
(
    const List<treeBoundBox>& bbs,
    const bool keepNonLocal,
    autoPtr<mapDistribute>& faceMap,
    autoPtr<mapDistribute>& pointMap
)
{
    // Get bbs of all domains
    // ~~~~~~~~~~~~~~~~~~~~~~

    {
        List<List<treeBoundBox>> newProcBb(Pstream::nProcs());

        switch(distType_)
        {
            case FOLLOW:
                newProcBb[Pstream::myProcNo()].setSize(bbs.size());
                forAll(bbs, i)
                {
                    newProcBb[Pstream::myProcNo()][i] = bbs[i];
                }
                Pstream::gatherList(newProcBb);
                Pstream::scatterList(newProcBb);
            break;

            case INDEPENDENT:
                newProcBb = independentlyDistributedBbs(*this);
            break;

            case FROZEN:
                return;
            break;

            default:
                FatalErrorInFunction
                    << "Unsupported distribution type." << exit(FatalError);
            break;
        }

        if (newProcBb == procBb_)
        {
            return;
        }
        else
        {
            procBb_.transfer(newProcBb);
            dict_.set("bounds", procBb_[Pstream::myProcNo()]);
        }
    }


    // Debug information
    if (debug)
    {
        labelList nTris(Pstream::nProcs());
        nTris[Pstream::myProcNo()] = triSurface::size();
        Pstream::gatherList(nTris);
        Pstream::scatterList(nTris);

        InfoInFunction
            << "before distribution:" << endl << "\tproc\ttris" << endl;

        forAll(nTris, proci)
        {
            Info<< '\t' << proci << '\t' << nTris[proci] << endl;
        }
        Info<< endl;
    }


    // Use procBbs to determine which faces go where
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelListList faceSendMap(Pstream::nProcs());
    labelListList pointSendMap(Pstream::nProcs());

    forAll(procBb_, proci)
    {
        overlappingSurface
        (
            *this,
            procBb_[proci],
            pointSendMap[proci],
            faceSendMap[proci]
        );

        if (debug)
        {
            // Pout<< "Overlapping with proc " << proci
            //    << " faces:" << faceSendMap[proci].size()
            //    << " points:" << pointSendMap[proci].size() << endl << endl;
        }
    }

    if (keepNonLocal)
    {
        // Include in faceSendMap/pointSendMap the triangles that are
        // not mapped to any processor so they stay local.

        const triSurface& s = static_cast<const triSurface&>(*this);

        boolList includedFace(s.size(), true);

        forAll(faceSendMap, proci)
        {
            if (proci != Pstream::myProcNo())
            {
                forAll(faceSendMap[proci], i)
                {
                    includedFace[faceSendMap[proci][i]] = false;
                }
            }
        }

        // Combine includedFace (all triangles that are not on any neighbour)
        // with overlap.

        forAll(faceSendMap[Pstream::myProcNo()], i)
        {
            includedFace[faceSendMap[Pstream::myProcNo()][i]] = true;
        }

        subsetMesh
        (
            s,
            includedFace,
            pointSendMap[Pstream::myProcNo()],
            faceSendMap[Pstream::myProcNo()]
        );
    }


    // Send over how many faces/points I need to receive
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelListList faceSendSizes(Pstream::nProcs());
    faceSendSizes[Pstream::myProcNo()].setSize(Pstream::nProcs());
    forAll(faceSendMap, proci)
    {
        faceSendSizes[Pstream::myProcNo()][proci] = faceSendMap[proci].size();
    }
    Pstream::gatherList(faceSendSizes);
    Pstream::scatterList(faceSendSizes);



    // Exchange surfaces
    // ~~~~~~~~~~~~~~~~~

    // Storage for resulting surface
    List<labelledTri> allTris;
    pointField allPoints;

    labelListList faceConstructMap(Pstream::nProcs());
    labelListList pointConstructMap(Pstream::nProcs());


    // My own surface first
    // ~~~~~~~~~~~~~~~~~~~~

    {
        labelList pointMap;
        triSurface subSurface
        (
            subsetMesh
            (
                *this,
                faceSendMap[Pstream::myProcNo()],
                pointMap
            )
        );

        allTris = subSurface;
        allPoints = subSurface.points();

        faceConstructMap[Pstream::myProcNo()] = identity
        (
            faceSendMap[Pstream::myProcNo()].size()
        );
        pointConstructMap[Pstream::myProcNo()] = identity
        (
            pointSendMap[Pstream::myProcNo()].size()
        );
    }



    // Send all
    // ~~~~~~~~

    forAll(faceSendSizes, proci)
    {
        if (proci != Pstream::myProcNo())
        {
            if (faceSendSizes[Pstream::myProcNo()][proci] > 0)
            {
                OPstream str(Pstream::commsTypes::blocking, proci);

                labelList pointMap;
                triSurface subSurface
                (
                    subsetMesh
                    (
                        *this,
                        faceSendMap[proci],
                        pointMap
                    )
                );

                // if (debug)
                //{
                //    Pout<< "Sending to " << proci
                //        << " faces:" << faceSendMap[proci].size()
                //        << " points:" << subSurface.points().size() << endl
                //        << endl;
                //}

                str << subSurface;
           }
        }
    }


    // Receive and merge all
    // ~~~~~~~~~~~~~~~~~~~~~

    forAll(faceSendSizes, proci)
    {
        if (proci != Pstream::myProcNo())
        {
            if (faceSendSizes[proci][Pstream::myProcNo()] > 0)
            {
                IPstream str(Pstream::commsTypes::blocking, proci);

                // Receive
                triSurface subSurface(str);

                // if (debug)
                //{
                //    Pout<< "Received from " << proci
                //        << " faces:" << subSurface.size()
                //        << " points:" << subSurface.points().size() << endl
                //        << endl;
                //}

                // Merge into allSurf
                merge
                (
                    mergeDist_,
                    subSurface,
                    subSurface.points(),

                    allTris,
                    allPoints,
                    faceConstructMap[proci],
                    pointConstructMap[proci]
                );

                // if (debug)
                //{
                //    Pout<< "Current merged surface : faces:" << allTris.size()
                //        << " points:" << allPoints.size() << endl << endl;
                //}
           }
        }
    }


    faceMap.reset
    (
        new mapDistribute
        (
            allTris.size(),
            move(faceSendMap),
            move(faceConstructMap)
        )
    );
    pointMap.reset
    (
        new mapDistribute
        (
            allPoints.size(),
            move(pointSendMap),
            move(pointConstructMap)
        )
    );

    // Construct triSurface. Reuse storage.
    triSurface::operator=(triSurface(allTris, patches(), allPoints, true));

    clearOut();

    // Set the bounds() value to the boundBox of the undecomposed surface
    triSurfaceMesh::bounds() = boundBox(points());

    reduce(bounds().min(), minOp<point>());
    reduce(bounds().max(), maxOp<point>());

    // Regions stays same
    // Volume type stays same.

    distributeFields<label>(faceMap());
    distributeFields<scalar>(faceMap());
    distributeFields<vector>(faceMap());
    distributeFields<sphericalTensor>(faceMap());
    distributeFields<symmTensor>(faceMap());
    distributeFields<tensor>(faceMap());

    if (debug)
    {
        labelList nTris(Pstream::nProcs());
        nTris[Pstream::myProcNo()] = triSurface::size();
        Pstream::gatherList(nTris);
        Pstream::scatterList(nTris);

        InfoInFunction
            << "after distribution:" << endl << "\tproc\ttris" << endl;

        forAll(nTris, proci)
        {
            Info<< '\t' << proci << '\t' << nTris[proci] << endl;
        }
        Info<< endl;
    }
}


bool Foam::distributedTriSurfaceMesh::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    // Make sure dictionary goes to same directory as surface
    const_cast<fileName&>(dict_.instance()) = searchableSurface::instance();

    // Copy of triSurfaceMesh::writeObject except for the sorting of
    // triangles by region. This is done so we preserve region names,
    // even if locally we have zero triangles.
    {
        fileName fullPath(searchableSurface::objectPath());

        if (!mkDir(fullPath.path()))
        {
            return false;
        }

        // Important: preserve any zero-sized patches
        triSurface::write(fullPath, true);

        if (!isFile(fullPath))
        {
            return false;
        }
    }

    // Dictionary needs to be written in ascii - binary output not supported.
    return dict_.writeObject(IOstream::ASCII, ver, cmp, true);
}


void Foam::distributedTriSurfaceMesh::writeStats(Ostream& os) const
{
    boundBox bb;
    label nPoints;
    PatchTools::calcBounds(static_cast<const triSurface&>(*this), bb, nPoints);
    reduce(bb.min(), minOp<point>());
    reduce(bb.max(), maxOp<point>());

    os  << "Triangles    : " << returnReduce(triSurface::size(), sumOp<label>())
        << endl
        << "Vertices     : " << returnReduce(nPoints, sumOp<label>()) << endl
        << "Bounding Box : " << bb << endl;
}


// ************************************************************************* //
