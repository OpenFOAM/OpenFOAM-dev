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

#include "slidingInterface.H"
#include "polyTopoChange.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "enrichedPatch.H"
#include "DynamicList.H"
#include "pointHit.H"
#include "triPointRef.H"
#include "plane.H"
#include "polyTopoChanger.H"
#include "polyAddPoint.H"
#include "polyRemovePoint.H"
#include "polyAddFace.H"
#include "polyModifyPoint.H"
#include "polyModifyFace.H"
#include "polyRemoveFace.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::slidingInterface::edgeCoPlanarTolDefault_ = 0.8;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Index of debug signs:
// Index of debug signs:
// . - loop of the edge-to-face interaction detection
// x - reversal of direction in edge-to-face interaction detection
// + - complete edge-to-face interaction detection
// z - incomplete edge-to-face interaction detection.  This may be OK for edges
//     crossing from one to the other side of multiply connected master patch

// e - adding a point on edge
// f - adding a point on face
// . - collecting edges off another face for edge-to-edge cut
// + - finished collection of edges
// * - cut both master and slave
// n - cutting new edge
// t - intersection exists but it is outside of tolerance
// x - missed slave edge in cut
// - - missed master edge in cut
// u - edge already used in cutting

void Foam::slidingInterface::coupleInterface(polyTopoChange& ref) const
{
    if (debug)
    {
        Pout<< "void slidingInterface::coupleInterface"
            << "(polyTopoChange& ref) : "
            << "Coupling sliding interface " << name() << endl;
    }

    const polyMesh& mesh = topoChanger().mesh();

    const pointField& points = mesh.points();
    const faceList& faces = mesh.faces();

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const faceZoneMesh& faceZones = mesh.faceZones();

    const primitiveFacePatch& masterPatch =
        faceZones[masterFaceZoneID_.index()]();

    const labelList& masterPatchAddr = faceZones[masterFaceZoneID_.index()];

    const boolList& masterPatchFlip =
        faceZones[masterFaceZoneID_.index()].flipMap();

    const primitiveFacePatch& slavePatch =
        faceZones[slaveFaceZoneID_.index()]();

    const labelList& slavePatchAddr = faceZones[slaveFaceZoneID_.index()];

    const boolList& slavePatchFlip =
        faceZones[slaveFaceZoneID_.index()].flipMap();

    const edgeList& masterEdges = masterPatch.edges();
    const labelListList& masterPointEdges = masterPatch.pointEdges();
    const labelList& masterMeshPoints = masterPatch.meshPoints();
    const pointField& masterLocalPoints = masterPatch.localPoints();
    const labelListList& masterFaceFaces = masterPatch.faceFaces();
    const labelListList& masterFaceEdges = masterPatch.faceEdges();
    const Map<label>& masterMeshPointMap = masterPatch.meshPointMap();

    const edgeList& slaveEdges = slavePatch.edges();
    const labelListList& slavePointEdges = slavePatch.pointEdges();
    const labelList& slaveMeshPoints = slavePatch.meshPoints();
    const pointField& slaveLocalPoints = slavePatch.localPoints();
    const Map<label>& slaveMeshPointMap = slavePatch.meshPointMap();
    const vectorField& slavePointNormals = slavePatch.pointNormals();

    // Collect projection addressing
    if
    (
        !(
            slavePointPointHitsPtr_
         && slavePointEdgeHitsPtr_
         && slavePointFaceHitsPtr_
         && masterPointEdgeHitsPtr_
        )
    )
    {
        FatalErrorInFunction
            << "Point projection addressing not available."
            << abort(FatalError);
    }

    const labelList& slavePointPointHits = *slavePointPointHitsPtr_;
    const labelList& slavePointEdgeHits = *slavePointEdgeHitsPtr_;
    const List<objectHit>& slavePointFaceHits = *slavePointFaceHitsPtr_;
    const labelList& masterPointEdgeHits = *masterPointEdgeHitsPtr_;
    const pointField& projectedSlavePoints = *projectedSlavePointsPtr_;

    // Create enriched patch
    enrichedPatch cutPatch
    (
        masterPatch,
        slavePatch,
        slavePointPointHits,
        slavePointEdgeHits,
        slavePointFaceHits
    );

    // Get reference to list of added point for the enriched patch.
    // This will be used for point addition
    Map<point>& pointMap = cutPatch.pointMap();

    // Get reference to the list of merged points
    Map<label>& pointMergeMap = cutPatch.pointMergeMap();

    // Create mapping for every merged point of the slave patch
    forAll(slavePointPointHits, pointi)
    {
        if (slavePointPointHits[pointi] >= 0)
        {
            // Pout<< "Inserting point merge pair: " << slaveMeshPoints[pointi]
            //     << " : " << masterMeshPoints[slavePointPointHits[pointi]]
            //     << endl;

            pointMergeMap.insert
            (
                slaveMeshPoints[pointi],
                masterMeshPoints[slavePointPointHits[pointi]]
            );
        }
    }

    // Collect the list of used edges for every slave edge

    List<labelHashSet> usedMasterEdges(slaveEdges.size());

    // Collect of slave point hits
    forAll(slavePointPointHits, pointi)
    {
        // For point hits, add all point-edges from master side to all point
        const labelList& curSlaveEdges = slavePointEdges[pointi];

        if (slavePointPointHits[pointi] > -1)
        {
            const labelList& curMasterEdges =
                masterPointEdges[slavePointPointHits[pointi]];

            // Mark all current master edges as used for all the current slave
            // edges
            forAll(curSlaveEdges, slaveEdgeI)
            {
                labelHashSet& sm = usedMasterEdges[curSlaveEdges[slaveEdgeI]];

                forAll(curMasterEdges, masterEdgeI)
                {
                    sm.insert(curMasterEdges[masterEdgeI]);
                }
            }
        }
        else if (slavePointEdgeHits[pointi] > -1)
        {
            // For edge hits, add the master edge
            forAll(curSlaveEdges, slaveEdgeI)
            {
                usedMasterEdges[curSlaveEdges[slaveEdgeI]].insert
                (
                    slavePointEdgeHits[pointi]
                );
            }
        }
    }

    // Collect off master point hits
    // For every master point that has hit an edge, add all edges coming from
    // that point to the slave edge that has been hit
    forAll(masterPointEdgeHits, masterPointi)
    {
        if (masterPointEdgeHits[masterPointi] > -1)
        {
            const labelList& curMasterEdges = masterPointEdges[masterPointi];

            labelHashSet& sm =
                usedMasterEdges[masterPointEdgeHits[masterPointi]];

            forAll(curMasterEdges, masterEdgeI)
            {
                sm.insert(curMasterEdges[masterEdgeI]);
            }
        }
    }

    // Pout<< "used edges: " << endl;
    // forAll(usedMasterEdges, edgeI)
    // {
    //     Pout<< "edge: " << edgeI
    //         << " used: " << usedMasterEdges[edgeI].toc()
    //         << endl;
    // }

    // For every master and slave edge make a list of points to be added into
    // that edge.
    List<DynamicList<label>> pointsIntoMasterEdges(masterEdges.size());
    List<DynamicList<label>> pointsIntoSlaveEdges(slaveEdges.size());

    // Add all points from the slave patch that have hit the edge
    forAll(slavePointEdgeHits, pointi)
    {
        if (slavePointEdgeHits[pointi] > -1)
        {
            // Create a new point on the master edge

            point edgeCutPoint =
                masterEdges[slavePointEdgeHits[pointi]].line
                (
                    masterLocalPoints
                ).nearestDist(projectedSlavePoints[pointi]).hitPoint();

            label newPoint =
                ref.setAction
                (
                    polyAddPoint
                    (
                        edgeCutPoint,             // point
                        slaveMeshPoints[pointi],  // master point
                        cutPointZoneID_.index(),  // zone for point
                        true                      // supports a cell
                    )
                );

            // Pout<< "Inserting merge pair off edge: "
            //     << slaveMeshPoints[pointi] << " " << newPoint
            //     << " cut point: " << edgeCutPoint
            //     << " orig: " << slaveLocalPoints[pointi]
            //     << " proj: " << projectedSlavePoints[pointi]
            //     << endl;

            // Add the new edge point into the merge map
            pointMergeMap.insert(slaveMeshPoints[pointi], newPoint);

            pointsIntoMasterEdges[slavePointEdgeHits[pointi]].append
            (
                newPoint
            );

            // Add the point into the enriched patch map
            pointMap.insert
            (
                newPoint,
                edgeCutPoint
            );

            if (debug)
            {
                Pout<< "e";
                // Pout<< newPoint << " = " << edgeCutPoint << endl;
            }
        }
    }

    if (debug)
    {
        Pout<< endl;
    }

    // Add all points from the slave patch that have hit a face
    forAll(slavePointFaceHits, pointi)
    {
        if
        (
            slavePointPointHits[pointi] < 0
         && slavePointEdgeHits[pointi] < 0
         && slavePointFaceHits[pointi].hit()
        )
        {
            label newPoint =
                ref.setAction
                (
                    polyAddPoint
                    (
                        projectedSlavePoints[pointi],   // point
                        slaveMeshPoints[pointi],        // master point
                        cutPointZoneID_.index(),        // zone for point
                        true                            // supports a cell
                    )
                );

            // Pout<< "Inserting merge pair off face: "
            //     << slaveMeshPoints[pointi]
            //     << " " << newPoint
            //     << endl;

            // Add the new edge point into the merge map
            pointMergeMap.insert(slaveMeshPoints[pointi], newPoint);

            // Add the point into the enriched patch map
            pointMap.insert
            (
                newPoint,
                projectedSlavePoints[pointi]
            );

            if (debug)
            {
                Pout<< "f: " << newPoint << " = "
                    << projectedSlavePoints[pointi] << endl;
            }
        }
    }

    forAll(masterPointEdgeHits, pointi)
    {
        if (masterPointEdgeHits[pointi] > -1)
        {
            pointsIntoSlaveEdges[masterPointEdgeHits[pointi]].append
            (
                masterMeshPoints[pointi]
            );
        }
    }

    // Cut all slave edges.
    // Collect all master edges the slave edge interacts with.  Skip
    // all the edges already marked as used.  For every unused edge,
    // calculate the cut and insert the new point into the master and
    // slave edge.
    // For the edge selection algorithm, see, comment in
    // slidingInterfaceProjectPoints.C.
    // Edge cutting algorithm:
    // As the master patch defines the cutting surface, the newly
    // inserted point needs to be on the master edge.  Also, in 3-D
    // the pair of edges generally misses each other rather than
    // intersect.  Therefore, the intersection is calculated using the
    // plane the slave edge defines during projection.  The plane is
    // defined by the centrepoint of the edge in the original
    // configuration and the projected end points.  In case of flat
    // geometries (when the three points are colinear), the plane is
    // defined by the two projected end-points and the slave edge
    // normal used as the in-plane vector.  When the intersection
    // point is created, it is added as a new point for both the
    // master and the slave edge.
    // In order to be able to re-create the points on edges in mesh
    // motion without the topology change, the edge pair used to
    // create the cut point needs to be recorded in
    // cutPointEdgePairMap

    // Note.  "Processing slave edges" code is repeated twice in the
    // sliding intergace coupling in order to allow the point
    // projection to be done separately from the actual cutting.
    // Please change consistently with slidingInterfaceProjectPoints.C
    //
    if (debug)
    {
        Pout<< "Processing slave edges " << endl;
    }

    if (!cutPointEdgePairMapPtr_)
    {
        FatalErrorInFunction
            << "Cut point edge pair map pointer not set."
            << abort(FatalError);
    }

    Map<Pair<edge>>& addToCpepm = *cutPointEdgePairMapPtr_;

    // Clear the old map
    addToCpepm.clear();

    // Create a map of faces the edge can interact with
    labelHashSet curFaceMap
    (
        nFacesPerSlaveEdge_*primitiveMesh::edgesPerFace_
    );

    labelHashSet addedFaces(2*primitiveMesh::edgesPerFace_);

    forAll(slaveEdges, edgeI)
    {
        const edge& curEdge = slaveEdges[edgeI];

        if
        (
            slavePointFaceHits[curEdge.start()].hit()
         || slavePointFaceHits[curEdge.end()].hit()
        )
        {
            labelHashSet& curUme = usedMasterEdges[edgeI];

            // Pout<< "Doing edge " << edgeI
            //     << " curEdge: " << curEdge
            //     << " curUme: " << curUme
            //     << endl;

            // Clear the maps
            curFaceMap.clear();
            addedFaces.clear();

            // Grab the faces for start and end points.
            const label startFace =
                slavePointFaceHits[curEdge.start()].hitObject();
            const label endFace = slavePointFaceHits[curEdge.end()].hitObject();

            // Pout<< "startFace: " << slavePointFaceHits[curEdge.start()]
            //     << " endFace: " << slavePointFaceHits[curEdge.end()]
            //     << endl;

            // Insert the start face into the list
            curFaceMap.insert(startFace);
            addedFaces.insert(startFace);

            // Pout<< "curFaceMap: " << curFaceMap.toc() << endl;

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

            // Collect the edges

            // Create a map of edges the edge can interact with
            labelHashSet curMasterEdgesMap
            (
                nFacesPerSlaveEdge_*primitiveMesh::edgesPerFace_
            );

            const labelList curFaces = curFaceMap.toc();

            // Pout<< "curFaces: " << curFaces << endl;

            forAll(curFaces, facei)
            {
                // Pout<< "face: " << curFaces[facei] << " "
                //     << masterPatch[curFaces[facei]]
                //     << " local: "
                //     << masterPatch.localFaces()[curFaces[facei]]
                //     << endl;

                const labelList& me = masterFaceEdges[curFaces[facei]];

                forAll(me, meI)
                {
                    curMasterEdgesMap.insert(me[meI]);
                }
            }

            const labelList curMasterEdges = curMasterEdgesMap.toc();

            // For all master edges to intersect, skip the ones
            // already used and cut the rest with a cutting plane.  If
            // the intersection point, falls inside of both edges, it
            // is valid.

            // Note.
            // The edge cutting code is repeated in
            // slidingInterface::modifyMotionPoints.  This is done for
            // efficiency reasons and avoids multiple creation of cutting
            // planes.  Please update both simultaneously.

            const point& a = projectedSlavePoints[curEdge.start()];
            const point& b = projectedSlavePoints[curEdge.end()];

            point c =
                0.5*
                (
                    slaveLocalPoints[curEdge.start()]
                  + slavePointNormals[curEdge.start()] // Add start normal
                  + slaveLocalPoints[curEdge.end()]
                  + slavePointNormals[curEdge.end()] // Add end normal
                );

            // Create the plane
            plane cutPlane(a, b, c);

            // Pout<< "a: " << a
            //     << " b: " << b
            //     << " c: " << c
            //     << " plane: " << cutPlane
            //     << endl;

            linePointRef curSlaveLine = curEdge.line(projectedSlavePoints);
            const scalar curSlaveLineMag = curSlaveLine.mag();

            // Pout<< "curSlaveLine: " << curSlaveLine << endl;

            forAll(curMasterEdges, masterEdgeI)
            {
                if (!curUme.found(curMasterEdges[masterEdgeI]))
                {
                    // New edge
                    if (debug)
                    {
                        Pout<< "n";
                    }

                    const label cmeIndex = curMasterEdges[masterEdgeI];
                    const edge& cme = masterEdges[cmeIndex];

                    // Pout<< "Edge " << cmeIndex
                    //     << " cme: " << cme
                    //     << " line: " << cme.line(masterLocalPoints)
                    //     << endl;

                    scalar cutOnMaster =
                        cutPlane.lineIntersect
                        (
                            cme.line(masterLocalPoints)
                        );

                    if
                    (
                        cutOnMaster > edgeEndCutoffTol_
                     && cutOnMaster < 1.0 - edgeEndCutoffTol_
                    )
                    {
                        // Master is cut, check the slave
                        point masterCutPoint =
                            masterLocalPoints[cme.start()]
                          + cutOnMaster*cme.vec(masterLocalPoints);

                        pointHit slaveCut =
                            curSlaveLine.nearestDist(masterCutPoint);

                        if (slaveCut.hit())
                        {
                            // Strict checking of slave cut to avoid capturing
                            // end points.
                            scalar cutOnSlave =
                                (
                                    (
                                        slaveCut.hitPoint()
                                      - curSlaveLine.start()
                                    ) & curSlaveLine.vec()
                                )/sqr(curSlaveLineMag);

                            // Calculate merge tolerance from the
                            // target edge length
                            scalar mergeTol = edgeCoPlanarTol_*mag(b - a);

                            // Pout<< "cutOnMaster: " << cutOnMaster
                            //     << " masterCutPoint: " << masterCutPoint
                            //     << " slaveCutPoint: " << slaveCut.hitPoint()
                            //     << " slaveCut.distance(): "
                            //     << slaveCut.distance()
                            //     << " slave length: " << mag(b - a)
                            //     << " mergeTol: " << mergeTol
                            //     << " 1: " << mag(b - a)
                            //     << " 2: "
                            //     << cme.line(masterLocalPoints).mag()
                            //     << endl;

                            if
                            (
                                cutOnSlave > edgeEndCutoffTol_
                             && cutOnSlave < 1.0 - edgeEndCutoffTol_
                             && slaveCut.distance() < mergeTol
                            )
                            {
                                // Cut both master and slave.  Add point
                                // to edge points The point is nominally
                                // added from the start of the master edge
                                // and added to the cut point zone
                                label newPoint =
                                    ref.setAction
                                    (
                                        polyAddPoint
                                        (
                                            masterCutPoint,           // point
                                            masterMeshPoints[cme.start()],// m p
                                            cutPointZoneID_.index(),  // zone
                                            true                      // active
                                        )
                                    );

                                // Pout<< "Inserting point: " << newPoint
                                //     << " as edge to edge intersection.  "
                                //     << "Slave edge: "
                                //     << edgeI << " " << curEdge
                                //     << " master edge: "
                                //     << cmeIndex << " " << cme
                                //     << endl;

                                pointsIntoSlaveEdges[edgeI].append(newPoint);
                                pointsIntoMasterEdges[cmeIndex].append
                                (
                                    newPoint
                                );

                                // Add the point into the enriched patch map
                                pointMap.insert
                                (
                                    newPoint,
                                    masterCutPoint
                                );

                                // Record which two edges intersect to
                                // create cut point
                                addToCpepm.insert
                                (
                                    newPoint,    // Cut point index
                                    Pair<edge>
                                    (
                                        edge
                                        (
                                            masterMeshPoints[cme.start()],
                                            masterMeshPoints[cme.end()]
                                        ),    // Master edge
                                        edge
                                        (
                                            slaveMeshPoints[curEdge.start()],
                                            slaveMeshPoints[curEdge.end()]
                                        )// Slave edge
                                    )
                                );

                                if (debug)
                                {
                                    Pout<< " " << newPoint << " = "
                                        << masterCutPoint << " ";
                                }
                            }
                            else
                            {
                                if (debug)
                                {
                                    // Intersection exists but it is too far
                                    Pout<< "t";
                                }
                            }
                        }
                        else
                        {
                            if (debug)
                            {
                                // Missed slave edge
                                Pout<< "x";
                            }
                        }
                    }
                    else
                    {
                        if (debug)
                        {
                            // Missed master edge
                            Pout<< "-";
                        }
                    }
                }
                else
                {
                    if (debug)
                    {
                        Pout<< "u";
                    }
                }
            }

            if (debug)
            {
                Pout<< endl;
            }
        } // End if both ends missing
    } // End for all slave edges

//     Pout<< "pointsIntoMasterEdges: " << pointsIntoMasterEdges << endl;
//     Pout<< "pointsIntoSlaveEdges: " << pointsIntoSlaveEdges << endl;

    // Re-pack the points into edges lists
    labelListList pime(pointsIntoMasterEdges.size());

    forAll(pointsIntoMasterEdges, i)
    {
        pime[i].transfer(pointsIntoMasterEdges[i]);
    }

    labelListList pise(pointsIntoSlaveEdges.size());

    forAll(pointsIntoSlaveEdges, i)
    {
        pise[i].transfer(pointsIntoSlaveEdges[i]);
    }

    // Prepare the enriched faces
    cutPatch.calcEnrichedFaces
    (
        pime,
        pise,
        projectedSlavePoints
    );

    // Demand driven calculate the cut faces. Apart from the
    // cutFaces/cutFaceMaster/cutFaceSlave no information from the cutPatch
    // is used anymore!
    const faceList& cutFaces = cutPatch.cutFaces();
    const labelList& cutFaceMaster = cutPatch.cutFaceMaster();
    const labelList& cutFaceSlave = cutPatch.cutFaceSlave();

    const labelList& masterFc = masterFaceCells();
    const labelList& slaveFc = slaveFaceCells();

    // Couple the interface

    // Algorithm:
    // 1) Go through the cut faces and check if the cut face is the same as the
    //    defining master or slave.  If so, modify the appropriate
    //    face and mark the other for relegation into the face zone.
    // 2) If different, mark both sides for relegation and insert the new face


    boolList orphanedMaster(masterPatch.size(), false);
    boolList orphanedSlave(slavePatch.size(), false);

    forAll(cutFaces, facei)
    {
        const face& curCutFace = cutFaces[facei];
        const label curMaster = cutFaceMaster[facei];
        const label curSlave = cutFaceSlave[facei];

//         Pout<< "Doing insertion of face " << facei << ": ";

        // Check if the face has changed topologically
        bool insertedFace = false;

        if (curMaster >= 0)
        {
            // Face has got a master
            if (curCutFace == masterPatch[curMaster])
            {
                // Face is equal to master.  Modify master face.
//                 Pout<< "Face is equal to master and is ";

                // If the face has got both master and slave, it is an
                // internal face; otherwise it is a patch face in the
                // master patch.  Keep it in the master face zone.

                if (curSlave >= 0)
                {
//                     Pout<< "internal" << endl;
                    if (masterFc[curMaster] < slaveFc[curSlave])
                    {
                        // Cut face should point into slave.
                        // Be careful about flips in zone!
                        ref.setAction
                        (
                            polyModifyFace
                            (
                                curCutFace,                  // new face
                                masterPatchAddr[curMaster],  // master face id
                                masterFc[curMaster],         // owner
                                slaveFc[curSlave],           // neighbour
                                false,                       // flux flip
                                -1,                          // patch ID
                                false,                       // remove from zone
                                masterFaceZoneID_.index(),   // zone ID
                                masterPatchFlip[curMaster]   // zone flip
                            )
                        );

                        // Pout<< "modifying master face. Old master: "
                        //     << masterPatch[curMaster]
                        //     << " new face: " << curCutFace.reverseFace()
                        //     << " own: " << masterFc[curMaster]
                        //     << " nei: " << slaveFc[curSlave] << endl;
                    }
                    else
                    {
                        // Cut face should point into master.  Flip required.
                        // Be careful about flips in zone!
                        ref.setAction
                        (
                            polyModifyFace
                            (
                                curCutFace.reverseFace(),    // new face
                                masterPatchAddr[curMaster],  // master face id
                                slaveFc[curSlave],           // owner
                                masterFc[curMaster],         // neighbour
                                true,                        // flux flip
                                -1,                          // patch ID
                                false,                       // remove from zone
                                masterFaceZoneID_.index(),   // zone ID
                                !masterPatchFlip[curMaster]  // zone flip
                            )
                        );
                    }

                    // Orphan the slave
                    orphanedSlave[curSlave] = true;
                }
                else
                {
//                     Pout<< "master boundary" << endl;
                    ref.setAction
                    (
                        polyModifyFace
                        (
                            curCutFace,                  // new face
                            masterPatchAddr[curMaster],  // master face index
                            masterFc[curMaster],         // owner
                            -1,                          // neighbour
                            false,                       // flux flip
                            masterPatchID_.index(),      // patch ID
                            false,                       // remove from zone
                            masterFaceZoneID_.index(),   // zone ID
                            masterPatchFlip[curMaster]   // zone flip
                        )
                    );
                }

                insertedFace = true;
            }
        }
        else if (curSlave >= 0)
        {
            // Face has got a slave

            // Renumber the slave face into merged labels
            face rsf(slavePatch[curSlave]);

            forAll(rsf, i)
            {
                Map<label>::const_iterator mpIter = pointMergeMap.find(rsf[i]);

                if (mpIter != pointMergeMap.end())
                {
                    rsf[i] = mpIter();
                }
            }

            if (curCutFace == rsf)
            {
                // Face is equal to slave.  Modify slave face.
                // Pout<< "Face is equal to slave and is ";

                if (curMaster >= 0)
                {
                    // Pout<< "regular internal" << endl;
                    if (masterFc[curMaster] < slaveFc[curSlave])
                    {
                        ref.setAction
                        (
                            polyModifyFace
                            (
                                curCutFace,                  // new face
                                slavePatchAddr[curSlave],    // master face id
                                masterFc[curMaster],         // owner
                                slaveFc[curSlave],           // neighbour
                                true,                        // flux flip
                                -1,                          // patch ID
                                false,                       // remove from zone
                                slaveFaceZoneID_.index(),    // zone ID
                                !slavePatchFlip[curMaster]   // zone flip
                            )
                        );
                    }
                    else
                    {
                        // Cut face should point into master.
                        // Be careful about flips in zone!
                        // Pout<< "flipped internal" << endl;
                        ref.setAction
                        (
                            polyModifyFace
                            (
                                curCutFace.reverseFace(),    // new face
                                slavePatchAddr[curSlave],    // master face id
                                slaveFc[curSlave],           // owner
                                masterFc[curMaster],         // neighbour
                                false,                       // flux flip
                                -1,                          // patch ID
                                false,                       // remove from zone
                                slaveFaceZoneID_.index(),    // zone ID
                                slavePatchFlip[curSlave]     // zone flip
                            )
                        );
                    }

                    // Orphan the master
                    orphanedMaster[curMaster] = true;
                }
                else
                {
                    // Pout<< "slave boundary" << endl;
                    ref.setAction
                    (
                        polyModifyFace
                        (
                            curCutFace,                  // new face
                            slavePatchAddr[curSlave],    // master face index
                            slaveFc[curSlave],           // owner
                            -1,                          // neighbour
                            false,                       // flux flip
                            slavePatchID_.index(),       // patch ID
                            false,                       // remove from zone
                            slaveFaceZoneID_.index(),    // zone ID
                            slavePatchFlip[curSlave]     // zone flip
                        )
                    );
                }

                insertedFace = true;
            }
        }
        else
        {
            FatalErrorInFunction
                << "Face " << facei << " in cut faces has neither a master "
                << "nor a slave.  Error in the cutting algorithm on modify."
                << abort(FatalError);
        }

        if (!insertedFace)
        {
            // Face is different from both master and slave
            // Pout<< "Face different from both master and slave" << endl;

            if (curMaster >= 0)
            {
                if (curSlave >= 0)
                {
                    // Add internal face
                    if (masterFc[curMaster] < slaveFc[curSlave])
                    {
                        // Pout<< "Adding internal face " << curCutFace
                        //     << " owner: " << masterFc[curMaster]
                        //     << " slave: " << slaveFc[curSlave]
                        //     << " master face: " << masterPatchAddr[curMaster]
                        //     << endl;

                        // Cut face should point into slave.
                        ref.setAction
                        (
                            polyAddFace
                            (
                                curCutFace,                  // new face
                                masterFc[curMaster],         // owner
                                slaveFc[curSlave],           // neighbour
                                -1,                          // master point
                                -1,                          // master edge
                                masterPatchAddr[curMaster],  // master face id
                                false,                       // flux flip
                                -1,                          // patch ID
                                cutFaceZoneID_.index(),      // zone ID
                                false                        // zone flip
                            )
                        );
                    }
                    else
                    {
                        // Cut face should point into master.  Flip required.
                        ref.setAction
                        (
                            polyAddFace
                            (
                                curCutFace.reverseFace(),    // new face
                                slaveFc[curSlave],           // owner
                                masterFc[curMaster],         // neighbour
                                -1,                          // master point
                                -1,                          // master edge
                                masterPatchAddr[curMaster],  // master face id
                                true,                        // flux flip
                                -1,                          // patch ID
                                cutFaceZoneID_.index(),      // zone ID
                                true                         // zone flip
                            )
                        );
                    }

                    // Orphan slave
                    orphanedSlave[curSlave] = true;
                }
                else
                {
                    // Pout<< "Adding solo master face " << curCutFace
                    //     << " owner: " << masterFc[curMaster]
                    //     << " master face: " << masterPatchAddr[curMaster]
                    //     << endl;

                    // Add master patch face
                    ref.setAction
                    (
                        polyAddFace
                        (
                            curCutFace,                  // new face
                            masterFc[curMaster],         // owner
                            -1,                          // neighbour
                            -1,                          // master point
                            -1,                          // master edge
                            masterPatchAddr[curMaster],  // master face index
                            false,                       // flux flip
                            masterPatchID_.index(),      // patch ID
                            cutFaceZoneID_.index(),      // zone ID
                            false                        // zone flip
                        )
                    );
                }

                // Orphan master
                orphanedMaster[curMaster] = true;
            }
            else if (curSlave >= 0)
            {
                // Pout<< "Adding solo slave face " << curCutFace
                //     << " owner: " << slaveFc[curSlave]
                //     << " master face: " << slavePatchAddr[curSlave]
                //     << endl;

                // Add slave patch face
                ref.setAction
                (
                    polyAddFace
                    (
                        curCutFace,                  // new face
                        slaveFc[curSlave],           // owner
                        -1,                          // neighbour
                        -1,                          // master point
                        -1,                          // master edge
                        slavePatchAddr[curSlave],    // master face index
                        false,                       // flux flip
                        slavePatchID_.index(),       // patch ID
                        cutFaceZoneID_.index(),      // zone ID
                        false                        // zone flip
                    )
                );

                // Orphan slave
                orphanedSlave[curSlave] = true;
            }
            else
            {
                FatalErrorInFunction
                    << "Face " << facei << " in cut faces has neither a master "
                    << "nor a slave.  Error in the cutting algorithm on add."
                    << abort(FatalError);
            }
        }
    }

    // Move the orphaned faces into the face zone
    // Pout<< "Orphaned master faces: " << orphanedMaster << endl;
    // Pout<< "Orphaned slave faces: " << orphanedSlave << endl;

    label nOrphanedMasters = 0;

    forAll(orphanedMaster, facei)
    {
        if (orphanedMaster[facei])
        {
            nOrphanedMasters++;

            //// Recover original orientation
            //ref.setAction
            //(
            //    polyModifyFace
            //    (
            //        masterPatch[facei],                 // new face
            //        masterPatchAddr[facei],             // master face index
            //        -1,                                 // owner
            //        -1,                                 // neighbour
            //        false,                              // flux flip
            //        -1,                                 // patch ID
            //        false,                              // remove from zone
            //        masterFaceZoneID_.index(),          // zone ID
            //        false                               // zone flip
            //    )
            //);

            //Pout<< "**MJ:deleting master face " << masterPatchAddr[facei]
            //    << " old verts:" << masterPatch[facei] << endl;
            ref.setAction(polyRemoveFace(masterPatchAddr[facei]));
        }
    }

    label nOrphanedSlaves = 0;

    forAll(orphanedSlave, facei)
    {
        if (orphanedSlave[facei])
        {
            nOrphanedSlaves++;

            //// Recover original orientation
            //ref.setAction
            //(
            //    polyModifyFace
            //    (
            //        slavePatch[facei],                // new face
            //        slavePatchAddr[facei],            // slave face index
            //        -1,                               // owner
            //        -1,                               // neighbour
            //        false,                            // flux flip
            //        -1,                               // patch ID
            //        false,                            // remove from zone
            //        slaveFaceZoneID_.index(),         // zone ID
            //        false                             // zone flip
            //    )
            //);

            //Pout<< "**MJ:deleting slave face " << slavePatchAddr[facei]
            //    << " old verts:" << slavePatch[facei] << endl;
            ref.setAction(polyRemoveFace(slavePatchAddr[facei]));
        }
    }

    if (debug)
    {
        Pout<< "Number of orphaned faces: "
            << "master = " << nOrphanedMasters << " out of "
            << orphanedMaster.size()
            << " slave = " << nOrphanedSlaves << " out of "
            << orphanedSlave.size() << endl;
    }

    // Finished coupling the plane of sliding interface.

    // Modify faces influenced by the sliding interface
    // These are the faces belonging to the master and slave cell
    // layer which have not already been modified.
    // Modification comes in three types:
    // 1) eliminate the points that have been removed when the old sliding
    //    interface has been removed
    // 2) Merge the slave points that have hit points on the master patch
    // 3) Introduce new points resulting from the new sliding interface cut

    // Collect all affected faces

    // Master side

    // Grab the list of faces in the layer
    const labelList& masterStickOuts = masterStickOutFaces();

    // Pout<< "masterStickOuts: " << masterStickOuts << endl;

    // Re-create the master stick-out faces
    forAll(masterStickOuts, facei)
    {
        // Renumber the face and remove additional points

        const label curFaceID = masterStickOuts[facei];

        const face& oldRichFace = faces[curFaceID];

        bool changed = false;

        // Remove removed points from the face
        face oldFace(oldRichFace.size());
        label nOldFace = 0;

        forAll(oldRichFace, pointi)
        {
            if (ref.pointRemoved(oldRichFace[pointi]))
            {
                changed = true;
            }
            else
            {
                // Point off patch
                oldFace[nOldFace] = oldRichFace[pointi];
                nOldFace++;
            }
        }

        oldFace.setSize(nOldFace);

        // Pout<< "old rich master face: " << oldRichFace
        //     << " old face: " << oldFace
        //     << endl;

        DynamicList<label> newFaceLabels(2*oldFace.size());

        forAll(oldFace, pointi)
        {
            if (masterMeshPointMap.found(oldFace[pointi]))
            {
                // Point is in master patch. Add it

                // If the point is a direct hit, grab its label; otherwise
                // keep the original
                if (pointMergeMap.found(oldFace[pointi]))
                {
                    changed = true;
                    newFaceLabels.append
                    (
                        pointMergeMap.find(oldFace[pointi])()
                    );
                }
                else
                {
                    newFaceLabels.append(oldFace[pointi]);
                }

                // Find if there are additional points inserted onto the edge
                // between current point and the next point
                // Algorithm:
                // 1) Find all the edges in the master patch coming
                //    out of the current point.
                // 2) If the next point in the face to pick the right edge
                const label localFirstLabel =
                    masterMeshPointMap.find(oldFace[pointi])();

                const labelList& curEdges = masterPointEdges[localFirstLabel];

                const label  nextLabel = oldFace.nextLabel(pointi);

                Map<label>::const_iterator mmpmIter =
                    masterMeshPointMap.find(nextLabel);

                if (mmpmIter != masterMeshPointMap.end())
                {
                    // Pout<< "found label pair " << oldFace[pointi]
                    //     << " and " << nextLabel;
                    // Find the points on the edge between them
                    const label localNextLabel = mmpmIter();

                    forAll(curEdges, curEdgeI)
                    {
                        if
                        (
                            masterEdges[curEdges[curEdgeI]].otherVertex
                            (
                                localFirstLabel
                            )
                         == localNextLabel
                        )
                        {
                            // Pout<< " found edge: " << curEdges[curEdgeI]
                            //     << endl;

                            // Get points on current edge
                            const labelList& curPime = pime[curEdges[curEdgeI]];

                            if (curPime.size())
                            {
                                changed = true;
                                // Pout<< "curPime: " << curPime << endl;
                                // Insert the edge points into the face
                                // in the correct order
                                const point& startPoint =
                                    masterLocalPoints[localFirstLabel];

                                vector e =
                                    masterLocalPoints[localNextLabel]
                                  - startPoint;

                                e /= magSqr(e);

                                scalarField edgePointWeights(curPime.size());

                                forAll(curPime, curPimeI)
                                {
                                    edgePointWeights[curPimeI] =
                                        (
                                            e
                                          & (
                                              pointMap.find(curPime[curPimeI])()
                                            - startPoint
                                            )
                                        );
                                }

                                if (debug)
                                {
                                    if
                                    (
                                        min(edgePointWeights) < 0
                                     || max(edgePointWeights) > 1
                                    )
                                    {
                                        FatalErrorInFunction
                                            << "Error in master stick-out edge "
                                            << "point collection."
                                            << abort(FatalError);
                                    }
                                }

                                // Go through the points and collect
                                // them based on weights from lower to
                                // higher.  This gives the correct
                                // order of points along the edge.
                                for
                                (
                                    label passI = 0;
                                    passI < edgePointWeights.size();
                                    passI++
                                )
                                {
                                    // Max weight can only be one, so
                                    // the sorting is done by
                                    // elimination.
                                    label nextPoint = -1;
                                    scalar dist = 2;

                                    forAll(edgePointWeights, wI)
                                    {
                                        if (edgePointWeights[wI] < dist)
                                        {
                                            dist = edgePointWeights[wI];
                                            nextPoint = wI;
                                        }
                                    }

                                    // Insert the next point and reset
                                    // its weight to exclude it from
                                    // future picks
                                    newFaceLabels.append(curPime[nextPoint]);
                                    edgePointWeights[nextPoint] = great;
                                }
                            }

                            break;
                        } // End of found edge
                    } // End of looking through current edges
                }
            }
            else
            {
                newFaceLabels.append(oldFace[pointi]);
            }
        }

        if (changed)
        {
            if (newFaceLabels.size() < 3)
            {
                FatalErrorInFunction
                    << "Face " << curFaceID << " reduced to less than "
                    << "3 points.  Topological/cutting error A." << nl
                    << "Old face: " << oldFace << " new face: " << newFaceLabels
                    << abort(FatalError);
            }

            // Get face zone and its flip
            label modifiedFaceZone = faceZones.whichZone(curFaceID);
            bool modifiedFaceZoneFlip = false;

            if (modifiedFaceZone >= 0)
            {
                modifiedFaceZoneFlip =
                    faceZones[modifiedFaceZone].flipMap()
                    [
                        faceZones[modifiedFaceZone].whichFace(curFaceID)
                    ];
            }

            face newFace;
            newFace.transfer(newFaceLabels);

            // Pout<< "Modifying master stick-out face " << curFaceID
            //     << " old face: " << oldFace
            //     << " new face: " << newFace
            //     << endl;

            // Modify the face
            if (mesh.isInternalFace(curFaceID))
            {
                ref.setAction
                (
                    polyModifyFace
                    (
                        newFace,                // modified face
                        curFaceID,              // label of face being modified
                        own[curFaceID],         // owner
                        nei[curFaceID],         // neighbour
                        false,                  // face flip
                        -1,                     // patch for face
                        false,                  // remove from zone
                        modifiedFaceZone,       // zone for face
                        modifiedFaceZoneFlip    // face flip in zone
                    )
                );
            }
            else
            {
                ref.setAction
                (
                    polyModifyFace
                    (
                        newFace,                // modified face
                        curFaceID,              // label of face being modified
                        own[curFaceID],         // owner
                        -1,                     // neighbour
                        false,                  // face flip
                        mesh.boundaryMesh().whichPatch(curFaceID),
                                                // patch for face
                        false,                  // remove from zone
                        modifiedFaceZone,       // zone for face
                        modifiedFaceZoneFlip    // face flip in zone
                    )
                );
            }
        }
    }

    // Pout<< "Finished master side" << endl;

    // Slave side

    // Grab the list of faces in the layer
    const labelList& slaveStickOuts = slaveStickOutFaces();

    // Pout<< "slaveStickOuts: " << slaveStickOuts << endl;

    const Map<label>& rpm = retiredPointMap();

    // Re-create the slave stick-out faces

    forAll(slaveStickOuts, facei)
    {
        // Renumber the face and remove additional points
        const label curFaceID = slaveStickOuts[facei];

        const face& oldRichFace = faces[curFaceID];

        bool changed = false;

        // Remove removed points from the face
        face oldFace(oldRichFace.size());
        label nOldFace = 0;

        forAll(oldRichFace, pointi)
        {
            if
            (
                rpm.found(oldRichFace[pointi])
             || slaveMeshPointMap.found(oldRichFace[pointi])
            )
            {
                // Point definitely live. Add it
                oldFace[nOldFace] = oldRichFace[pointi];
                nOldFace++;
            }
            else if
            (
                ref.pointRemoved(oldRichFace[pointi])
             || masterMeshPointMap.found(oldRichFace[pointi])
            )
            {
                // Point removed and not on slave patch
                // (first if takes care of that!) or
                // point belonging to master patch
                changed = true;
            }
            else
            {
                // Point off patch
                oldFace[nOldFace] = oldRichFace[pointi];
                nOldFace++;
            }
        }

        oldFace.setSize(nOldFace);

        DynamicList<label> newFaceLabels(2*oldFace.size());

        // Pout<< "old rich slave face: " << oldRichFace
        //     << " old face: " << oldFace
        //     << endl;

        forAll(oldFace, pointi)
        {
            // Try to find the point in retired points
            label curP = oldFace[pointi];

            Map<label>::const_iterator rpmIter = rpm.find(oldFace[pointi]);

            if (rpmIter != rpm.end())
            {
                changed = true;
                curP = rpmIter();
            }

            if (slaveMeshPointMap.found(curP))
            {
                // Point is in slave patch. Add it

                // If the point is a direct hit, grab its label; otherwise
                // keep the original
                if (pointMergeMap.found(curP))
                {
                    changed = true;
                    newFaceLabels.append
                    (
                        pointMergeMap.find(curP)()
                    );
                }
                else
                {
                    newFaceLabels.append(curP);
                }

                // Find if there are additional points inserted onto the edge
                // between current point and the next point
                // Algorithm:
                // 1) Find all the edges in the slave patch coming
                //    out of the current point.
                // 2) Use the next point in the face to pick the right edge

                const label localFirstLabel =
                    slaveMeshPointMap.find(curP)();

                const labelList& curEdges = slavePointEdges[localFirstLabel];

                label nextLabel = oldFace.nextLabel(pointi);

                Map<label>::const_iterator rpmNextIter =
                    rpm.find(nextLabel);

                if (rpmNextIter != rpm.end())
                {
                    nextLabel = rpmNextIter();
                }

                Map<label>::const_iterator mmpmIter =
                    slaveMeshPointMap.find(nextLabel);

                if (mmpmIter != slaveMeshPointMap.end())
                {
                    // Both points on the slave patch.
                    // Find the points on the edge between them
                    const label localNextLabel = mmpmIter();

                    forAll(curEdges, curEdgeI)
                    {
                        if
                        (
                            slaveEdges[curEdges[curEdgeI]].otherVertex
                            (
                                localFirstLabel
                            )
                         == localNextLabel
                        )
                        {
                            // Pout<< " found edge: " << curEdges[curEdgeI]
                            //     << endl;

                            // Get points on current edge
                            const labelList& curPise = pise[curEdges[curEdgeI]];

                            if (curPise.size())
                            {
                                changed = true;
                                // Pout<< "curPise: " << curPise << endl;
                                // Insert the edge points into the face
                                // in the correct order
                                const point& startPoint =
                                    projectedSlavePoints[localFirstLabel];

                                vector e =
                                    projectedSlavePoints[localNextLabel]
                                  - startPoint;

                                e /= magSqr(e);

                                scalarField edgePointWeights(curPise.size());

                                forAll(curPise, curPiseI)
                                {
                                    edgePointWeights[curPiseI] =
                                    (
                                        e
                                      & (
                                            pointMap.find(curPise[curPiseI])()
                                          - startPoint
                                        )
                                    );
                                }

                                if (debug)
                                {
                                    if
                                    (
                                        min(edgePointWeights) < 0
                                     || max(edgePointWeights) > 1
                                    )
                                    {
                                        FatalErrorInFunction
                                            << "Error in slave stick-out edge "
                                            << "point collection."
                                            << abort(FatalError);
                                        }
                                    }

                                // Go through the points and collect
                                // them based on weights from lower to
                                // higher.  This gives the correct
                                // order of points along the edge.
                                for
                                (
                                    label passI = 0;
                                    passI < edgePointWeights.size();
                                    passI++
                                )
                                {
                                    // Max weight can only be one, so
                                    // the sorting is done by
                                    // elimination.
                                    label nextPoint = -1;
                                    scalar dist = 2;

                                    forAll(edgePointWeights, wI)
                                    {
                                        if (edgePointWeights[wI] < dist)
                                        {
                                            dist = edgePointWeights[wI];
                                            nextPoint = wI;
                                        }
                                    }

                                    // Insert the next point and reset
                                    // its weight to exclude it from
                                    // future picks
                                    newFaceLabels.append(curPise[nextPoint]);
                                    edgePointWeights[nextPoint] = great;
                                }
                            }

                            break;
                        }
                    } // End of found edge
                } // End of looking through current edges
            }
            else
            {
                newFaceLabels.append(oldFace[pointi]);
            }
        }

        if (changed)
        {
            if (newFaceLabels.size() < 3)
            {
                FatalErrorInFunction
                    << "Face " << curFaceID << " reduced to less than "
                    << "3 points.  Topological/cutting error B." << nl
                    << "Old face: " << oldFace << " new face: " << newFaceLabels
                    << abort(FatalError);
            }

            // Get face zone and its flip
            label modifiedFaceZone = faceZones.whichZone(curFaceID);
            bool modifiedFaceZoneFlip = false;

            if (modifiedFaceZone >= 0)
            {
                modifiedFaceZoneFlip =
                    faceZones[modifiedFaceZone].flipMap()
                    [
                        faceZones[modifiedFaceZone].whichFace(curFaceID)
                    ];
            }

            face newFace;
            newFace.transfer(newFaceLabels);

            // Pout<< "Modifying slave stick-out face " << curFaceID
            //     << " old face: " << oldFace
            //     << " new face: " << newFace
            //     << endl;

            // Modify the face
            if (mesh.isInternalFace(curFaceID))
            {
                ref.setAction
                (
                    polyModifyFace
                    (
                        newFace,                // modified face
                        curFaceID,              // label of face being modified
                        own[curFaceID],         // owner
                        nei[curFaceID],         // neighbour
                        false,                  // face flip
                        -1,                     // patch for face
                        false,                  // remove from zone
                        modifiedFaceZone,       // zone for face
                        modifiedFaceZoneFlip    // face flip in zone
                    )
                );
            }
            else
            {
                ref.setAction
                (
                    polyModifyFace
                    (
                        newFace,                // modified face
                        curFaceID,              // label of face being modified
                        own[curFaceID],         // owner
                        -1,                     // neighbour
                        false,                  // face flip
                        mesh.boundaryMesh().whichPatch(curFaceID),
                                                // patch for face
                        false,                  // remove from zone
                        modifiedFaceZone,       // zone for face
                        modifiedFaceZoneFlip    // face flip in zone
                    )
                );
            }
        }
    }

    // Activate and retire slave patch points
    // This needs to be done last, so that the map of removed points
    // does not get damaged by point modifications

    if (!retiredPointMapPtr_)
    {
        FatalErrorInFunction
            << "Retired point map pointer not set."
            << abort(FatalError);
    }

    Map<label>& addToRpm = *retiredPointMapPtr_;

    // Clear the old map
    addToRpm.clear();

    label nRetiredPoints = 0;

    forAll(slaveMeshPoints, pointi)
    {
        if (pointMergeMap.found(slaveMeshPoints[pointi]))
        {
            // Retire the point - only used for supporting the detached
            // slave patch
            nRetiredPoints++;

            // ref.setAction
            // (
            //    polyModifyPoint
            //    (
            //        slaveMeshPoints[pointi],             // point ID
            //        points[slaveMeshPoints[pointi]],     // point
            //        false,                               // remove from zone
            //        mesh.pointZones().whichZone(slaveMeshPoints[pointi]),
            //                                             // zone
            //        false                                // in a cell
            //    )
            // );
            //Pout<< "MJ retire slave point " << slaveMeshPoints[pointi]
            //    << " coord " << points[slaveMeshPoints[pointi]]
            //    << endl;
            ref.setAction
            (
                polyRemovePoint
                (
                    slaveMeshPoints[pointi]
                )
            );

            addToRpm.insert
            (
                pointMergeMap.find(slaveMeshPoints[pointi])(),
                slaveMeshPoints[pointi]
            );
        }
        else
        {
            ref.setAction
            (
                polyModifyPoint
                (
                    slaveMeshPoints[pointi],             // point ID
                    points[slaveMeshPoints[pointi]],     // point
                    false,                               // remove from zone
                    mesh.pointZones().whichZone(slaveMeshPoints[pointi]),// zone
                    true                                 // in a cell
                )
            );
        }
    }

    if (debug)
    {
        Pout<< "Retired " << nRetiredPoints << " out of "
            << slaveMeshPoints.size() << " points." << endl;
    }

    // Grab cut face master and slave addressing
    if (cutFaceMasterPtr_) deleteDemandDrivenData(cutFaceMasterPtr_);
    cutFaceMasterPtr_ = new labelList(cutPatch.cutFaceMaster());

    if (cutFaceSlavePtr_) deleteDemandDrivenData(cutFaceSlavePtr_);
    cutFaceSlavePtr_ = new labelList(cutPatch.cutFaceSlave());

    // Finished coupling
    attached_ = true;

    if (debug)
    {
        Pout<< "void slidingInterface::coupleInterface("
            << "polyTopoChange& ref) : "
            << "Finished coupling sliding interface " << name() << endl;
    }
}


// ************************************************************************* //
