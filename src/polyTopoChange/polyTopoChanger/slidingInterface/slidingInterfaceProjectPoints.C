/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "slidingInterface.H"
#include "polyMesh.H"
#include "line.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::slidingInterface::pointMergeTolDefault_ = 0.05;
const Foam::scalar Foam::slidingInterface::edgeMergeTolDefault_ = 0.01;
const Foam::label Foam::slidingInterface::nFacesPerSlaveEdgeDefault_ = 5;
const Foam::label Foam::slidingInterface::edgeFaceEscapeLimitDefault_ = 10;

const Foam::scalar Foam::slidingInterface::integralAdjTolDefault_ = 0.05;
const Foam::scalar
    Foam::slidingInterface::edgeMasterCatchFractionDefault_ = 0.4;
const Foam::scalar Foam::slidingInterface::edgeEndCutoffTolDefault_ = 0.0001;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Index of debug signs:
// a - integral match adjustment: point adjusted
// n - integral match adjustment: point not adjusted
// . - loop of the edge-to-face interaction detection
// x - reversal of direction in edge-to-face interaction detection
// + - complete edge-to-face interaction detection
// z - incomplete edge-to-face interaction detection.  This may be OK for edges
//     crossing from one to the other side of multiply connected master patch
// * - colinear triangle: adjusting projection with slave face normal
// m - master point inserted into the edge

bool Foam::slidingInterface::projectPoints() const
{
    if (debug)
    {
        Pout<< "bool slidingInterface::projectPoints() : "
            << " for object " << name() << " : "
            << "Projecting slave points onto master surface." << endl;
    }

    // Algorithm:
    // 1) Go through all the points of the master and slave patch and calculate
    //    minimum edge length coming from the point.  Calculate the point
    //    merge tolerance as the fraction of minimum edge length.
    // 2) Project all the slave points onto the master patch
    //    in the normal direction.
    // 3) If some points have missed and the match is integral, the
    //    points in question will be adjusted.  Find the nearest face for
    //    those points and continue with the procedure.
    // 4) For every point, check all the points of the face it has hit.
    //    For every pair of points find if their distance is smaller than
    //    both the master and slave merge tolerance.  If so, the slave point
    //    is moved to the location of the master point.  Remember the master
    //    point index.
    // 5) For every unmerged slave point, check its distance to all the
    //    edges of the face it has hit.  If the distance is smaller than the
    //    edge merge tolerance, the point will be moved onto the edge.
    //    Remember the master edge index.
    // 6) The remaining slave points will be projected into faces.  Remember the
    //    master face index.
    // 7) For the points that miss the master patch, grab the nearest face
    //    on the master and leave the slave point where it started
    //    from and the miss is recorded.

    const polyMesh& mesh = this->mesh();

    const primitiveFacePatch& masterPatch =
        mesh.faceZones()[masterFaceZoneIndex_.index()]();

    const primitiveFacePatch& slavePatch =
        mesh.faceZones()[slaveFaceZoneIndex_.index()]();

    // Get references to local points, local edges and local faces
    // for master and slave patch
//     const labelList& masterMeshPoints = masterPatch.meshPoints();
    const pointField& masterLocalPoints = masterPatch.localPoints();
    const faceList& masterLocalFaces = masterPatch.localFaces();
    const edgeList& masterEdges = masterPatch.edges();
    const labelListList& masterEdgeFaces = masterPatch.edgeFaces();
    const labelListList& masterFaceEdges = masterPatch.faceEdges();
    const labelListList& masterFaceFaces = masterPatch.faceFaces();
//     Pout<< "Master patch.  Local points: " << masterLocalPoints << nl
//         << "Master patch.  Mesh points: " << masterPatch.meshPoints() << nl
//         << "Local faces: " << masterLocalFaces << nl
//         << "local edges: " << masterEdges << endl;

//     const labelList& slaveMeshPoints = slavePatch.meshPoints();
    const pointField& slaveLocalPoints = slavePatch.localPoints();
    const edgeList& slaveEdges = slavePatch.edges();
    const labelListList& slaveEdgeFaces = slavePatch.edgeFaces();
    const vectorField& slavePointNormals = slavePatch.pointNormals();
//     const vectorField& slaveFaceNormals = slavePatch.faceNormals();
//     Pout<< "Slave patch.  Local points: " << slaveLocalPoints << nl
//         << "Slave patch.  Mesh points: " << slavePatch.meshPoints() << nl
//         << "Local faces: " << slavePatch.localFaces() << nl
//         << "local edges: " << slaveEdges << endl;

    // Calculate min edge distance for points and faces

    // Calculate min edge length for the points and faces of master patch
    scalarField minMasterPointLength(masterLocalPoints.size(), great);
    scalarField minMasterFaceLength(masterPatch.size(), great);

    forAll(masterEdges, edgeI)
    {
        const edge& curEdge = masterEdges[edgeI];

        const scalar curLength =
            masterEdges[edgeI].mag(masterLocalPoints);

        // Do points
        minMasterPointLength[curEdge.start()] =
            min
            (
                minMasterPointLength[curEdge.start()],
                curLength
            );

        minMasterPointLength[curEdge.end()] =
            min
            (
                minMasterPointLength[curEdge.end()],
                curLength
            );

        // Do faces
        const labelList& curFaces = masterEdgeFaces[edgeI];

        forAll(curFaces, facei)
        {
            minMasterFaceLength[curFaces[facei]] =
                min
                (
                    minMasterFaceLength[curFaces[facei]],
                    curLength
                );
        }
    }

//     Pout<< "min length for master points: " << minMasterPointLength << endl
//         << "min length for master faces: " << minMasterFaceLength << endl;

    // Calculate min edge length for the points and faces of slave patch
    scalarField minSlavePointLength(slaveLocalPoints.size(), great);
    scalarField minSlaveFaceLength(slavePatch.size(), great);

    forAll(slaveEdges, edgeI)
    {
        const edge& curEdge = slaveEdges[edgeI];

        const scalar curLength =
            slaveEdges[edgeI].mag(slaveLocalPoints);

        // Do points
        minSlavePointLength[curEdge.start()] =
            min
            (
                minSlavePointLength[curEdge.start()],
                curLength
            );

        minSlavePointLength[curEdge.end()] =
            min
            (
                minSlavePointLength[curEdge.end()],
                curLength
            );

        // Do faces
        const labelList& curFaces = slaveEdgeFaces[edgeI];

        forAll(curFaces, facei)
        {
            minSlaveFaceLength[curFaces[facei]] =
                min
                (
                    minSlaveFaceLength[curFaces[facei]],
                    curLength
                );
        }
    }

//     Pout<< "min length for slave points: " << minSlavePointLength << endl
//         << "min length for slave faces: " << minSlaveFaceLength << endl;

    // Project slave points onto the master patch

    // Face hit by the slave point
    List<objectHit> slavePointFaceHits =
        slavePatch.projectPoints
        (
            masterPatch,
            slavePointNormals,
            projectionAlgo_
        );

//     Pout<< "USING N-SQAURED!!!" << endl;
//     List<objectHit> slavePointFaceHits =
//         projectPointsNSquared<face, List, const pointField&>
//         (
//             slavePatch,
//             masterPatch,
//             slavePointNormals,
//             projectionAlgo_
//         );

    if (debug)
    {
        label nHits = 0;

        forAll(slavePointFaceHits, pointi)
        {
            if (slavePointFaceHits[pointi].hit())
            {
                nHits++;
            }
        }

        Pout<< "Number of hits in point projection: " << nHits
            << " out of " << slavePointFaceHits.size() << " points."
            << endl;
    }

    // Projected slave points are stored for mesh motion correction
    if (projectedSlavePointsPtr_) delete projectedSlavePointsPtr_;

    projectedSlavePointsPtr_ =
        new pointField(slavePointFaceHits.size(), Zero);
    pointField& projectedSlavePoints = *projectedSlavePointsPtr_;

    // Adjust projection to type of match

    label nAdjustedPoints = 0;

    // If the match is integral and some points have missed,
    // find nearest master face
    if (matchType_ == INTEGRAL)
    {
        if (debug)
        {
            Pout<< "bool slidingInterface::projectPoints() for object "
                << name() << " : "
                << "Adjusting point projection for integral match: ";
        }

        forAll(slavePointFaceHits, pointi)
        {
            if (slavePointFaceHits[pointi].hit())
            {
                // Grab the hit point
                projectedSlavePoints[pointi] =
                    masterLocalFaces
                        [slavePointFaceHits[pointi].hitObject()].ray
                        (
                            slaveLocalPoints[pointi],
                            slavePointNormals[pointi],
                            masterLocalPoints,
                            projectionAlgo_
                        ).hitPoint();
            }
            else
            {
                // Grab the nearest point on the face (edge)
                pointHit missAdjust =
                    masterLocalFaces[slavePointFaceHits[pointi].hitObject()].ray
                    (
                        slaveLocalPoints[pointi],
                        slavePointNormals[pointi],
                        masterLocalPoints,
                        projectionAlgo_
                    );

                const point nearPoint = missAdjust.missPoint();
                const point missPoint =
                    slaveLocalPoints[pointi]
                  + slavePointNormals[pointi]*missAdjust.distance();

                // Calculate the tolerance
                const scalar mergeTol =
                    integralAdjTol_*minSlavePointLength[pointi];

                // Adjust the hit
                if (mag(nearPoint - missPoint) < mergeTol)
                {
                    if (debug)
                    {
                        Pout<< "a";
                    }

//                     Pout<< "Moving slave point in integral adjustment "
//                         << pointi << " miss point: " << missPoint
//                         << " near point: " << nearPoint
//                         << " mergeTol: " << mergeTol
//                         << " dist: " << mag(nearPoint - missPoint) << endl;

                    projectedSlavePoints[pointi] = nearPoint;

                    slavePointFaceHits[pointi] =
                        objectHit(true, slavePointFaceHits[pointi].hitObject());

                    nAdjustedPoints++;
                }
                else
                {
                    projectedSlavePoints[pointi] = slaveLocalPoints[pointi];

                    if (debug)
                    {
                        Pout<< "n";
                    }
                }
            }
        }

        if (debug)
        {
            Pout<< " done." << endl;
        }
    }
    else if (matchType_ == PARTIAL)
    {
        forAll(slavePointFaceHits, pointi)
        {
            if (slavePointFaceHits[pointi].hit())
            {
                // Grab the hit point
                projectedSlavePoints[pointi] =
                    masterLocalFaces
                        [slavePointFaceHits[pointi].hitObject()].ray
                        (
                            slaveLocalPoints[pointi],
                            slavePointNormals[pointi],
                            masterLocalPoints,
                            projectionAlgo_
                        ).hitPoint();
            }
            else
            {
                // The point remains where it started from
                projectedSlavePoints[pointi] = slaveLocalPoints[pointi];
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << " for object " << name()
            << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "Number of adjusted points in projection: "
            << nAdjustedPoints << endl;

        // Check for zero-length edges in slave projection
        scalar minEdgeLength = great;
        scalar el = 0;
        label nShortEdges = 0;

        forAll(slaveEdges, edgeI)
        {
            el = slaveEdges[edgeI].mag(projectedSlavePoints);

            if (el < small)
            {
                Pout<< "Point projection problems for edge: "
                    << slaveEdges[edgeI] << ". Length = " << el
                    << endl;

                nShortEdges++;
            }

            minEdgeLength = min(minEdgeLength, el);
        }

        if (nShortEdges > 0)
        {
            FatalErrorInFunction
                << " short projected edges "
                << "after adjustment for object " << name()
                << abort(FatalError);
        }
        else
        {
            Pout<< " ... projection OK." << endl;
        }
    }
//     scalarField magDiffs(mag(slaveLocalPoints - projectedSlavePoints));
//     Pout<< "Slave point face hits: " << slavePointFaceHits << nl
//         << "slave points: " << slaveLocalPoints << nl
//         << "Projected slave points: " << projectedSlavePoints
//         << "diffs: " << magDiffs << endl;

    // Merge projected points against master points

    labelList slavePointPointHits(slaveLocalPoints.size(), -1);
    labelList masterPointPointHits(masterLocalPoints.size(), -1);

    // Go through all the slave points and compare them against all the points
    // in the master face they hit.  If the distance between the projected point
    // and any of the master face points is smaller than both tolerances,
    // merge the projected point by:
    // 1) adjusting the projected point to coincide with the
    //    master point it merges with
    // 2) remembering the hit master point index in slavePointPointHits
    // 3) resetting the hit face to -1
    // 4) If a master point has been hit directly, it cannot be merged
    // into the edge. Mark it as used in the reverse list

    label nMergedPoints = 0;

    forAll(projectedSlavePoints, pointi)
    {
        if (slavePointFaceHits[pointi].hit())
        {
            // Taking a non-const reference so the point can be adjusted
            point& curPoint = projectedSlavePoints[pointi];

            // Get the hit face
            const face& hitFace =
                masterLocalFaces[slavePointFaceHits[pointi].hitObject()];

            label mergePoint = -1;
            scalar mergeDist = great;

            // Try all point before deciding on best fit.
            forAll(hitFace, hitPointi)
            {
                scalar dist =
                    mag(masterLocalPoints[hitFace[hitPointi]] - curPoint);

                // Calculate the tolerance
                const scalar mergeTol =
                    pointMergeTol_*
                    min
                    (
                        minSlavePointLength[pointi],
                        minMasterPointLength[hitFace[hitPointi]]
                    );

                if (dist < mergeTol && dist < mergeDist)
                {
                    mergePoint = hitFace[hitPointi];
                    mergeDist = dist;

//                     Pout<< "Merging slave point "
//                         << slavePatch.meshPoints()[pointi] << " at "
//                         << slaveLocalPoints[pointi] << " with master "
//                         << masterPatch.meshPoints()[mergePoint] << " at "
//                         << masterLocalPoints[mergePoint]
//                         << ". dist: " << mergeDist
//                         << " mergeTol: " << mergeTol << endl;
                }
            }

            if (mergePoint > -1)
            {
                // Point is to be merged with master point
                nMergedPoints++;

                slavePointPointHits[pointi] = mergePoint;
                curPoint = masterLocalPoints[mergePoint];
                masterPointPointHits[mergePoint] = pointi;
            }
        }
    }

//     Pout<< "slavePointPointHits: " << slavePointPointHits << nl
//         << "masterPointPointHits: " << masterPointPointHits << endl;

    if (debug)
    {
        // Check for zero-length edges in slave projection
        scalar minEdgeLength = great;
        scalar el = 0;

        forAll(slaveEdges, edgeI)
        {
            el = slaveEdges[edgeI].mag(projectedSlavePoints);

            if (el < small)
            {
                Pout<< "Point projection problems for edge: "
                    << slaveEdges[edgeI] << ". Length = " << el
                    << endl;
            }

            minEdgeLength = min(minEdgeLength, el);
        }

        if (minEdgeLength < small)
        {
            FatalErrorInFunction
                << " after point merge for object " << name()
                << abort(FatalError);
        }
        else
        {
            Pout<< " ... point merge OK." << endl;
        }
    }

    // Merge projected points against master edges

    labelList slavePointEdgeHits(slaveLocalPoints.size(), -1);

    label nMovedPoints = 0;

    forAll(projectedSlavePoints, pointi)
    {
        // Eliminate the points merged into points
        if (slavePointPointHits[pointi] < 0)
        {
            // Get current point position
            point& curPoint = projectedSlavePoints[pointi];

            // Get the hit face
            const labelList& hitFaceEdges =
                masterFaceEdges[slavePointFaceHits[pointi].hitObject()];

            // Calculate the tolerance
            const scalar mergeLength =
                min
                (
                    minSlavePointLength[pointi],
                    minMasterFaceLength[slavePointFaceHits[pointi].hitObject()]
                );

            const scalar mergeTol = pointMergeTol_*mergeLength;

            scalar minDistance = great;

            forAll(hitFaceEdges, edgeI)
            {
                const edge& curEdge = masterEdges[hitFaceEdges[edgeI]];

                pointHit edgeHit =
                    curEdge.line(masterLocalPoints).nearestDist(curPoint);

                if (edgeHit.hit())
                {
                    scalar dist =
                        mag(edgeHit.hitPoint() - projectedSlavePoints[pointi]);

                    if (dist < mergeTol && dist < minDistance)
                    {
                        // Point is to be moved onto master edge
                        nMovedPoints++;

                        slavePointEdgeHits[pointi] = hitFaceEdges[edgeI];
                        projectedSlavePoints[pointi] = edgeHit.hitPoint();

                        minDistance = dist;

//                         Pout<< "Moving slave point "
//                             << slavePatch.meshPoints()[pointi]
//                             << " (" << pointi
//                             << ") at " << slaveLocalPoints[pointi]
//                             << " onto master edge " << hitFaceEdges[edgeI]
//                             << " or ("
//                             << masterLocalPoints[curEdge.start()]
//                             << masterLocalPoints[curEdge.end()]
//                             << ") hit: " << edgeHit.hitPoint()
//                             << ". dist: " << dist
//                             << " mergeTol: " << mergeTol << endl;
                    }
                }
            } // end of hit face edges

            if (slavePointEdgeHits[pointi] > -1)
            {
                // Projected slave point has moved.  Re-attempt merge with
                // master points of the edge
                point& curPoint = projectedSlavePoints[pointi];

                const edge& hitMasterEdge =
                    masterEdges[slavePointEdgeHits[pointi]];

                label mergePoint = -1;
                scalar mergeDist = great;

                forAll(hitMasterEdge, hmeI)
                {
                    scalar hmeDist =
                        mag(masterLocalPoints[hitMasterEdge[hmeI]] - curPoint);

                    // Calculate the tolerance
                    const scalar mergeTol =
                        pointMergeTol_*
                        min
                        (
                            minSlavePointLength[pointi],
                            minMasterPointLength[hitMasterEdge[hmeI]]
                    );

                    if (hmeDist < mergeTol && hmeDist < mergeDist)
                    {
                        mergePoint = hitMasterEdge[hmeI];
                        mergeDist = hmeDist;

//                     Pout<< "Merging slave point; SECOND TRY "
//                         << slavePatch.meshPoints()[pointi] << " local "
//                         << pointi << " at "
//                         << slaveLocalPoints[pointi] << " with master "
//                         << masterPatch.meshPoints()[mergePoint] << " at "
//                         << masterLocalPoints[mergePoint]
//                         << ". hmeDist: " << mergeDist
//                         << " mergeTol: " << mergeTol << endl;
                    }
                }

                if (mergePoint > -1)
                {
                    // Point is to be merged with master point
                    slavePointPointHits[pointi] = mergePoint;
                    curPoint = masterLocalPoints[mergePoint];
                    masterPointPointHits[mergePoint] = pointi;

                    slavePointFaceHits[pointi] =
                        objectHit(true, slavePointFaceHits[pointi].hitObject());


                    // Disable edge merge
                    slavePointEdgeHits[pointi] = -1;
                }
            }
        }
    }

    if (debug)
    {
        Pout<< "Number of merged master points: " << nMergedPoints << nl
            << "Number of adjusted slave points: " << nMovedPoints << endl;

        // Check for zero-length edges in slave projection
        scalar minEdgeLength = great;
        scalar el = 0;

        forAll(slaveEdges, edgeI)
        {
            el = slaveEdges[edgeI].mag(projectedSlavePoints);

            if (el < small)
            {
                Pout<< "Point projection problems for edge: "
                    << slaveEdges[edgeI] << ". Length = " << el
                    << endl;
            }

            minEdgeLength = min(minEdgeLength, el);
        }

        if (minEdgeLength < small)
        {
            FatalErrorInFunction
            << " after correction for object " << name()
            << abort(FatalError);
        }
    }

//     Pout<< "slavePointEdgeHits: " << slavePointEdgeHits << endl;

    // Insert the master points into closest slave edge if appropriate

    // Algorithm:
    //    The face potentially interacts with all the points of the
    //    faces covering the path between its two ends.  Since we are
    //    examining an arbitrary shadow of the edge on a non-Euclidian
    //    surface, it is typically quite hard to do a geometric point
    //    search (under a shadow of a straight line).  Therefore, the
    //    search will be done topologically.
    //
    // I) Point collection
    // -------------------
    // 1) Grab the master faces to which the end points of the edge
    //    have been projected.
    // 2) Starting from the face containing the edge start, grab all
    //    its points and put them into the point lookup map.  Put the
    //    face onto the face lookup map.
    // 3) If the face of the end point is on the face lookup, complete
    //    the point collection step (this is checked during insertion.
    // 4) Start a new round of insertion.  Visit all faces in the face
    //    lookup and add their neighbours which are not already on the
    //    map.  Every time the new neighbour is found, add its points to
    //    the point lookup.  If the face of the end point is inserted,
    //    continue with the current roundof insertion and stop the
    //    algorithm.
    //
    // II) Point check
    // ---------------
    //    Grab all the points from the point collection and check them
    //    against the current edge.  If the edge-to-point distance is
    //    smaller than the smallest tolerance in the game (min of
    //    master point tolerance and left and right slave face around
    //    the edge tolerance) and the nearest point hits the edge, the
    //    master point will break the slave edge.  Check the actual
    //    distance of the original master position from the edge.  If
    //    the distance is smaller than a fraction of slave edge
    //    length, the hit is considered valid.  Store the slave edge
    //    index for every master point.

    labelList masterPointEdgeHits(masterLocalPoints.size(), -1);
    scalarField masterPointEdgeDist(masterLocalPoints.size(), great);

    // Note.  "Processing slave edges" code is repeated twice in the
    // sliding interface coupling in order to allow the point
    // projection to be done separately from the actual cutting.
    // Please change consistently with coupleSlidingInterface.C
    //

    if (debug)
    {
        Pout<< "Processing slave edges " << endl;
    }

    // Create a map of faces the edge can interact with
    labelHashSet curFaceMap
    (
        nFacesPerSlaveEdge_*primitiveMesh::edgesPerFace_
    );

    labelHashSet addedFaces(2*primitiveMesh::edgesPerFace_);

    forAll(slaveEdges, edgeI)
    {
        const edge& curEdge = slaveEdges[edgeI];

        {
            // Clear the maps
            curFaceMap.clear();
            addedFaces.clear();

            // Grab the faces for start and end points
            const label startFace =
                slavePointFaceHits[curEdge.start()].hitObject();
            const label endFace = slavePointFaceHits[curEdge.end()].hitObject();

            // Insert the start face into the list
            curFaceMap.insert(startFace);
            addedFaces.insert(startFace);

            // Pout<< "Doing edge " << edgeI
            //     << " or " << curEdge
            //     << " start: "
            //     << slavePointFaceHits[curEdge.start()].hitObject()
            //     << " end "
            //     << slavePointFaceHits[curEdge.end()].hitObject()
            //     << endl;

            // If the end face is on the list, the face collection is finished
            label nSweeps = 0;
            bool completed = false;

            while (nSweeps < edgeFaceEscapeLimit_)
            {
                nSweeps++;

                if (addedFaces.found(endFace))
                {
                    completed = true;
                }

                // Add all face neighbours of face in the map
                const labelList cf = addedFaces.toc();
                addedFaces.clear();

                forAll(cf, cfI)
                {
                    const labelList& curNbrs = masterFaceFaces[cf[cfI]];

                    forAll(curNbrs, nbrI)
                    {
                        if (!curFaceMap.found(curNbrs[nbrI]))
                        {
                            curFaceMap.insert(curNbrs[nbrI]);
                            addedFaces.insert(curNbrs[nbrI]);
                        }
                    }
                }

                if (completed) break;

                if (debug)
                {
                    Pout<< ".";
                }
            }

            if (!completed)
            {
                if (debug)
                {
                    Pout<< "x";
                }

                // It is impossible to reach the end from the start, probably
                // due to disconnected domain.  Do search in opposite direction

                label nReverseSweeps = 0;

                addedFaces.clear();
                curFaceMap.insert(endFace);
                addedFaces.insert(endFace);

                while (nReverseSweeps < edgeFaceEscapeLimit_)
                {
                    nReverseSweeps++;

                    if (addedFaces.found(startFace))
                    {
                        completed = true;
                    }

                    // Add all face neighbours of face in the map
                    const labelList cf = addedFaces.toc();
                    addedFaces.clear();

                    forAll(cf, cfI)
                    {
                        const labelList& curNbrs = masterFaceFaces[cf[cfI]];

                        forAll(curNbrs, nbrI)
                        {
                            if (!curFaceMap.found(curNbrs[nbrI]))
                            {
                                curFaceMap.insert(curNbrs[nbrI]);
                                addedFaces.insert(curNbrs[nbrI]);
                            }
                        }
                    }

                    if (completed) break;

                    if (debug)
                    {
                        Pout<< ".";
                    }
                }
            }

            if (completed)
            {
                if (debug)
                {
                    Pout<< "+ ";
                }
            }
            else
            {
                if (debug)
                {
                    Pout<< "z ";
                }
            }

            // Collect the points

            // Create a map of points the edge can interact with
            labelHashSet curPointMap
            (
                nFacesPerSlaveEdge_*primitiveMesh::pointsPerFace_
            );

            const labelList curFaces = curFaceMap.toc();
//             Pout<< "curFaces: " << curFaces << endl;
            forAll(curFaces, facei)
            {
                const face& f = masterLocalFaces[curFaces[facei]];

                forAll(f, pointi)
                {
                    curPointMap.insert(f[pointi]);
                }
            }

            const labelList curMasterPoints = curPointMap.toc();

            // Check all the points against the edge.

            linePointRef edgeLine = curEdge.line(projectedSlavePoints);

            const vector edgeVec = edgeLine.vec();
            const scalar edgeMag = edgeLine.mag();

            // Calculate actual distance involved in projection.  This
            // is used to reject master points out of reach.
            // Calculated as a combination of travel distance in projection and
            // edge length
            scalar slaveCatchDist =
                edgeMasterCatchFraction_*edgeMag
              + 0.5*
                (
                    mag
                    (
                        projectedSlavePoints[curEdge.start()]
                      - slaveLocalPoints[curEdge.start()]
                    )
                  + mag
                    (
                        projectedSlavePoints[curEdge.end()]
                      - slaveLocalPoints[curEdge.end()]
                    )
                );

            // The point merge distance needs to be measured in the
            // plane of the slave edge.  The unit vector is calculated
            // as a cross product of the edge vector and the edge
            // projection direction.  When checking for the distance
            // in plane, a minimum of the master-to-edge and
            // projected-master-to-edge distance is used, to avoid
            // problems with badly defined master planes.  HJ,
            // 17/Oct/2004
            vector edgeNormalInPlane =
                edgeVec
              ^ (
                    slavePointNormals[curEdge.start()]
                  + slavePointNormals[curEdge.end()]
                );

            edgeNormalInPlane /= mag(edgeNormalInPlane);

            forAll(curMasterPoints, pointi)
            {
                const label cmp = curMasterPoints[pointi];

                // Skip the current point if the edge start or end has
                // been adjusted onto in
                if
                (
                    slavePointPointHits[curEdge.start()] == cmp
                 || slavePointPointHits[curEdge.end()] == cmp
                 || masterPointPointHits[cmp] > -1
                )
                {
// Pout<< "Edge already snapped to point.  Skipping." << endl;
                    continue;
                }

                // Check if the point actually hits the edge within bounds
                pointHit edgeLineHit =
                    edgeLine.nearestDist(masterLocalPoints[cmp]);

                if (edgeLineHit.hit())
                {
                    // If the distance to the line is smaller than
                    // the tolerance the master point needs to be
                    // inserted into the edge

                    // Strict checking of slave cut to avoid capturing
                    // end points.
                    scalar cutOnSlave =
                        ((edgeLineHit.hitPoint() - edgeLine.start()) & edgeVec)
                        /sqr(edgeMag);

                    scalar distInEdgePlane =
                        min
                        (
                            edgeLineHit.distance(),
                            mag
                            (
                                (
                                    masterLocalPoints[cmp]
                                  - edgeLineHit.hitPoint()
                                )
                              & edgeNormalInPlane
                            )
                        );
//                     Pout<< "master point: " << cmp
//                         << " cutOnSlave " << cutOnSlave
//                         << " distInEdgePlane: " << distInEdgePlane
//                         << " tol1: " << pointMergeTol_*edgeMag
//                         << " hitDist: " << edgeLineHit.distance()
//                         << " tol2: " <<
//                         min
//                         (
//                             slaveCatchDist,
//                             masterPointEdgeDist[cmp]
//                         ) << endl;

                    // Not a point hit, check for edge
                    if
                    (
                        cutOnSlave > edgeEndCutoffTol_
                     && cutOnSlave < 1.0 - edgeEndCutoffTol_ // check edge cut
                     && distInEdgePlane < edgeMergeTol_*edgeMag // merge plane
                     && edgeLineHit.distance()
                      < min
                        (
                            slaveCatchDist,
                            masterPointEdgeDist[cmp]
                        )
                    )
                    {
                        if (debug)
                        {
                            if (masterPointEdgeHits[cmp] == -1)
                            {
                                // First hit
                                Pout<< "m";
                            }
                            else
                            {
                                // Repeat hit
                                Pout<< "M";
                            }
                        }

                        // Snap to point onto edge
                        masterPointEdgeHits[cmp] = edgeI;
                        masterPointEdgeDist[cmp] = edgeLineHit.distance();

//                         Pout<< "Inserting master point "
//                             << masterPatch.meshPoints()[cmp]
//                             << " (" << cmp
//                             << ") at " << masterLocalPoints[cmp]
//                             << " into slave edge " << edgeI
//                             << " " << curEdge
//                             << " cutOnSlave: " << cutOnSlave
//                             << " distInEdgePlane: " << distInEdgePlane
//                             << ". dist: " << edgeLineHit.distance()
//                             << " mergeTol: "
//                             << edgeMergeTol_*edgeMag
//                             << " other tol: " <<
//                             min
//                             (
//                                 slaveCatchDist,
//                                 masterPointEdgeDist[cmp]
//                             )
//                             << endl;
                    }
                }
            }
        } // End if both ends missing
    } // End all slave edges

    if (debug)
    {
        Pout<< endl;
    }
//     Pout<< "masterPointEdgeHits: " << masterPointEdgeHits << endl;

    if (debug)
    {
        Pout<< "bool slidingInterface::projectPoints() for object "
            << name() << " : "
            << "Finished projecting  points. Topology = ";
    }

//     Pout<< "slavePointPointHits: " << slavePointPointHits << nl
//         << "slavePointEdgeHits: " << slavePointEdgeHits << nl
//         << "slavePointFaceHits: " << slavePointFaceHits << nl
//         << "masterPointEdgeHits: " << masterPointEdgeHits << endl;

    // The four lists:
    // - slavePointPointHits
    // - slavePointEdgeHits
    // - slavePointFaceHits
    // - masterPointEdgeHits
    //   define how the two patches will be merged topologically.
    //   If the lists have not changed since the last merge, the
    //   sliding interface changes only geometrically and simple mesh
    //   motion will suffice.  Otherwise, a topological change is
    //   required.

    // Compare the result with the current state
    if (!attached_)
    {
        // The mesh needs to change topologically
        trigger_ = true;

        // Store the addressing arrays and projected points
        deleteDemandDrivenData(slavePointPointHitsPtr_);
        slavePointPointHitsPtr_ = new labelList(slavePointPointHits);

        deleteDemandDrivenData(slavePointEdgeHitsPtr_);
        slavePointEdgeHitsPtr_ = new labelList(slavePointEdgeHits);

        deleteDemandDrivenData(slavePointFaceHitsPtr_);
        slavePointFaceHitsPtr_ = new List<objectHit>(slavePointFaceHits);

        deleteDemandDrivenData(masterPointEdgeHitsPtr_);
        masterPointEdgeHitsPtr_ = new labelList(masterPointEdgeHits);

        if (debug)
        {
            Pout<< "(Detached interface) changing." << endl;
        }
    }
    else
    {
        // Compare the lists against the stored lists.  If any of them
        // has changed, topological change will be executed.
        trigger_ = false;

        if
        (
            !slavePointPointHitsPtr_
         || !slavePointEdgeHitsPtr_
         || !slavePointFaceHitsPtr_
         || !masterPointEdgeHitsPtr_
        )
        {
            // Must be restart.  Force topological change.
            deleteDemandDrivenData(slavePointPointHitsPtr_);
            slavePointPointHitsPtr_ = new labelList(slavePointPointHits);

            deleteDemandDrivenData(slavePointEdgeHitsPtr_);
            slavePointEdgeHitsPtr_ = new labelList(slavePointEdgeHits);

            deleteDemandDrivenData(slavePointFaceHitsPtr_);
            slavePointFaceHitsPtr_ = new List<objectHit>(slavePointFaceHits);

            deleteDemandDrivenData(masterPointEdgeHitsPtr_);
            masterPointEdgeHitsPtr_ = new labelList(masterPointEdgeHits);

            if (debug)
            {
                Pout<< "(Attached interface restart) changing." << endl;
            }

            trigger_ = true;
            return trigger_;
        }

        if (slavePointPointHits != (*slavePointPointHitsPtr_))
        {
            if (debug)
            {
                Pout<< "(Point projection) ";
            }

            trigger_ = true;
        }

        if (slavePointEdgeHits != (*slavePointEdgeHitsPtr_))
        {
            if (debug)
            {
                Pout<< "(Edge projection) ";
            }

            trigger_ = true;
        }

        // Comparison for face hits needs to exclude points that have hit
        // another point or edge
        bool faceHitsDifferent = false;

        const List<objectHit>& oldPointFaceHits = *slavePointFaceHitsPtr_;

        forAll(slavePointFaceHits, pointi)
        {
            if
            (
                slavePointPointHits[pointi] < 0
             && slavePointEdgeHits[pointi] < 0
            )
            {
                // This is a straight face hit
                if (slavePointFaceHits[pointi] != oldPointFaceHits[pointi])
                {
                    // Two lists are different
                    faceHitsDifferent = true;
                    break;
                }
            }
        }

        if (faceHitsDifferent)
        {
            if (debug)
            {
                Pout<< "(Face projection) ";
            }

            trigger_ = true;

        }

        if (masterPointEdgeHits != (*masterPointEdgeHitsPtr_))
        {
            if (debug)
            {
                Pout<< "(Master point projection) ";
            }

            trigger_ = true;
        }


        if (trigger_)
        {
            // Grab new data
            deleteDemandDrivenData(slavePointPointHitsPtr_);
            slavePointPointHitsPtr_ = new labelList(slavePointPointHits);

            deleteDemandDrivenData(slavePointEdgeHitsPtr_);
            slavePointEdgeHitsPtr_ = new labelList(slavePointEdgeHits);

            deleteDemandDrivenData(slavePointFaceHitsPtr_);
            slavePointFaceHitsPtr_ = new List<objectHit>(slavePointFaceHits);

            deleteDemandDrivenData(masterPointEdgeHitsPtr_);
            masterPointEdgeHitsPtr_ = new labelList(masterPointEdgeHits);

            if (debug)
            {
                Pout<< "changing." << endl;
            }
        }
        else
        {
            if (debug)
            {
                Pout<< "preserved." << endl;
            }
        }
    }

    return trigger_;
}


void Foam::slidingInterface::clearPointProjection() const
{
    deleteDemandDrivenData(slavePointPointHitsPtr_);
    deleteDemandDrivenData(slavePointEdgeHitsPtr_);
    deleteDemandDrivenData(slavePointFaceHitsPtr_);
    deleteDemandDrivenData(masterPointEdgeHitsPtr_);

    deleteDemandDrivenData(projectedSlavePointsPtr_);
}


// ************************************************************************* //
