/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Description
    Create coupled match faces and add them to the cells

\*---------------------------------------------------------------------------*/

#include "starMesh.H"
#include "boolList.H"
#include "pointHit.H"
#include "IOmanip.H"
#include "boundBox.H"
#include "Map.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void starMesh::createCoupleMatches()
{
    // Loop through all couples and create intersection faces. Add all points
    // of intersection faces to the couple points lists. The numbering of
    // the list is set such that the list can be appended to the
    // existing points list

    // Estimate the number of cells affected by couple matches
    const label cellMapSize = min
    (
        cellShapes_.size()/10,
        couples_.size()*2
    );

    // Store newly created faces for each cell
    Map<SLList<face> > cellAddedFaces(cellMapSize);

    Map<SLList<label> > cellRemovedFaces(cellMapSize);

    // In order to remove often allocation, remember the number of live points.
    // If you run out of space in point creation, increase it by the number of
    // couples (good scale) and resize at the end;
    label nLivePoints = points_.size();

    const label infoJump = max(1000, couples_.size()/20);

    forAll(couples_, coupleI)
    {
        if (coupleI % infoJump == 0)
        {
            Info<< "Doing couple " << coupleI << ". STAR couple ID: "
                << couples_[coupleI].coupleID() << endl;
        }

        // Initialise cell edges for master and slave cells
        const coupledFacePair& fp = couples_[coupleI];
        const face& masterFace = cellFaces_[fp.masterCell()][fp.masterFace()];
        const face& slaveFace = cellFaces_[fp.slaveCell()][fp.slaveFace()];

#       ifdef DEBUG_COUPLE
        Info<< "coupleI: " << coupleI << endl
            << "masterFace: " << masterFace << endl
            << "master points: " << masterFace.points(points_) << endl
            << "slaveFace: " << slaveFace << endl
            << "slave points: " << slaveFace.points(points_)
            << endl << endl;
#       endif

        // check the angle of face area vectors
         scalar faceAreaAngle =
             mag
             (
                 -(masterFace.normal(points_) & slaveFace.normal(points_))/
                 (masterFace.mag(points_)*slaveFace.mag(points_) + VSMALL)
             );

         if (faceAreaAngle < 0.94)
         {
             Info<< "Couple direction mismatch in the couple match "
                 << coupleI << ". STAR couple ID: "
                 << couples_[coupleI].coupleID() << endl
                 << "The angle between face normals is "
                 << radToDeg(Foam::acos(faceAreaAngle))
                 << " deg." << endl
                 << "master cell: " << fp.masterCell()
                 << " STAR number: " << starCellID_[fp.masterCell()]
                 << " type: " << cellShapes_[fp.masterCell()].model().name()
                 << " face: " << fp.masterFace() << endl
                 << "slave cell : " << fp.slaveCell()
                 << " STAR number: " << starCellID_[fp.slaveCell()]
                 << " type: " << cellShapes_[fp.slaveCell()].model().name()
                 << " face: " << fp.slaveFace() << endl;
         }

        // Deal with integral patches
        if (fp.integralMatch())
        {
            // Master face is replaced by a set of slave faces

            Map<SLList<label> >::iterator crfIter =
                cellRemovedFaces.find(fp.masterCell());

            if (crfIter == cellRemovedFaces.end())
            {
                cellRemovedFaces.insert
                (
                    fp.masterCell(),
                    SLList<label>(fp.masterFace())
                );
            }
            else
            {
                crfIter().append(fp.masterFace());
            }

            Map<SLList<face> >::iterator cafIter =
                cellAddedFaces.find(fp.masterCell());
            if (cafIter == cellAddedFaces.end())
            {
                cellAddedFaces.insert
                (
                    fp.masterCell(),
                    SLList<face>(slaveFace.reverseFace())
                );
            }
            else
            {
                cafIter().append(slaveFace.reverseFace());
            }
        }
        else
        {
            // Create cut faces, which replace both master and slave faces

            // Store newly created points
            SLList<point> coupleFacePoints;

            // Master data
            edgeList masterEdges = masterFace.edges();
            List<SLList<label> > masterEdgePoints(masterEdges.size());

            // Slave data
            edgeList slaveEdges = slaveFace.edges();
            List<SLList<label> > slaveEdgePoints(slaveEdges.size());

            // Find common plane
            vector n = masterFace.normal(points_);
            n /= mag(n) + VSMALL;

            // Loop through all edges of the master face. For every edge,
            // intersect it with all edges of the cutting face.
            forAll(masterEdges, masterEdgeI)
            {
                const edge& curMasterEdge = masterEdges[masterEdgeI];

                point P = points_[curMasterEdge.start()];

                // get d and return it into plane
                vector d = curMasterEdge.vec(points_);
                d -= n*(n & d);

#               ifdef DEBUG_COUPLE_INTERSECTION
                Info<< "curMasterEdge: " << curMasterEdge << endl
                    << "P: " << P << endl << "d: " << d << endl;
#               endif

                // go through all slave edges and try to get an intersection.
                // The point is created along the original master edge rather
                // than its corrected direction.
                forAll(slaveEdges, slaveEdgeI)
                {
                    const edge& curSlaveEdge = slaveEdges[slaveEdgeI];

                    point S = points_[curSlaveEdge.start()];

                    // get e and return it into plane
                    vector e = curSlaveEdge.vec(points_);
                    e -= n*(n & e);
                    scalar det = -(e & (n ^ d));

#                   ifdef DEBUG_COUPLE_INTERSECTION
                    Info<< "curSlaveEdge: " << curSlaveEdge << endl
                        << "S: " << S << endl
                        << "e: " << e << endl;
#                   endif

                    if (mag(det) > SMALL)
                    {
                        // non-singular matrix. Look for intersection
                        scalar beta = ((S - P) & (n ^ d))/det;

#                       ifdef DEBUG_COUPLE_INTERSECTION
                        Info<< " beta: " << beta << endl;
#                       endif

                        if (beta > -smallMergeTol_ && beta < 1 + smallMergeTol_)
                        {
                            // slave intersection OK. Try master intersection
                            scalar alpha =
                                (((S - P) & d) + beta*(d & e))/magSqr(d);

#                           ifdef DEBUG_COUPLE_INTERSECTION
                            Info<< " alpha: " << alpha << endl;
#                           endif

                            if
                            (
                                alpha > -smallMergeTol_
                             && alpha < 1 + smallMergeTol_
                            )
                            {
                                // intersection of non-parallel edges
#                               ifdef DEBUG_COUPLE_INTERSECTION
                                Info<< "intersection of non-parallel edges"
                                    << endl;
#                               endif


                                // check for insertion of start-end
                                // points in the middle of the other
                                // edge
                                if (alpha < smallMergeTol_)
                                {
                                    // inserting the start of master edge
                                    if
                                    (
                                        beta > smallMergeTol_
                                     && beta < 1 - smallMergeTol_
                                    )
                                    {
                                        slaveEdgePoints[slaveEdgeI].append
                                        (
                                            curMasterEdge.start()
                                        );
                                    }
                                }
                                else if (alpha > 1 - smallMergeTol_)
                                {
                                    // inserting the end of master edge
                                    if
                                    (
                                        beta > smallMergeTol_
                                     && beta < 1 - smallMergeTol_
                                    )
                                    {
                                        slaveEdgePoints[slaveEdgeI].append
                                        (
                                            curMasterEdge.end()
                                        );
                                    }
                                }
                                else if (beta < smallMergeTol_)
                                {
                                    // inserting the start of the slave edge
                                    if
                                    (
                                        alpha > smallMergeTol_
                                     && alpha < 1 - smallMergeTol_
                                    )
                                    {
                                        masterEdgePoints[masterEdgeI].append
                                        (
                                            curSlaveEdge.start()
                                        );
                                    }
                                }
                                else if (beta > 1 - smallMergeTol_)
                                {
                                    // inserting the start of the slave edge
                                    if
                                    (
                                        alpha > smallMergeTol_
                                     && alpha < 1 - smallMergeTol_
                                    )
                                    {
                                        masterEdgePoints[masterEdgeI].append
                                        (
                                            curSlaveEdge.end()
                                        );
                                    }
                                }
                                else
                                {
                                    masterEdgePoints[masterEdgeI].append
                                    (
                                        nLivePoints + coupleFacePoints.size()
                                    );

                                    slaveEdgePoints[slaveEdgeI].append
                                    (
                                        nLivePoints + coupleFacePoints.size()
                                    );

#                                   ifdef DEBUG_COUPLE_INTERSECTION
                                    Info<< "regular intersection. "
                                        << "Adding point: "
                                        << coupleFacePoints.size()
                                        << " which is "
                                        << P + alpha*curMasterEdge.vec(points_)
                                        << endl;
#                                   endif

                                    // A new point is created. Warning:
                                    // using original edge for accuracy.
                                    //
                                    coupleFacePoints.append
                                        (P + alpha*curMasterEdge.vec(points_));
                                }
                            }
                        }
                    }
                    else
                    {
                        // Add special cases, for intersection of two
                        // parallel line Warning. Here, typically, no new
                        // points will be created. Either one or two of
                        // the slave edge points need to be added to the
                        // master edge and vice versa. The problem is that
                        // no symmetry exists, i.e. both operations needs
                        // to be done separately for both master and slave
                        // side.

                        // Master side
                        // check if the first or second point of slave edge is
                        // on the master edge
                        vector ps = S - P;

                        bool colinear = false;

                        if (mag(ps) < SMALL)
                        {
                            // colinear because P and S are the same point
                            colinear = true;
                        }
                        else if
                        (
                            (ps & d)/(mag(ps)*mag(d)) > 1.0 - smallMergeTol_
                        )
                        {
                            // colinear because ps and d are parallel
                            colinear = true;
                        }

                        if (colinear)
                        {
                            scalar alpha1 = (ps & d)/magSqr(d);

                            if
                            (
                                alpha1 > -smallMergeTol_
                             && alpha1 < 1 + smallMergeTol_
                            )
                            {
#                               ifdef DEBUG_COUPLE_INTERSECTION
                                Info<< "adding irregular master "
                                    << "intersection1: "
                                    << points_[slaveEdges[slaveEdgeI].start()]
                                    << endl;
#                               endif

                                masterEdgePoints[masterEdgeI].append
                                (
                                    slaveEdges[slaveEdgeI].start()
                                );
                            }

                             scalar alpha2 = ((ps + e) & d)/magSqr(d);

                            if
                            (
                                alpha2 > -smallMergeTol_
                             && alpha2 < 1 + smallMergeTol_
                            )
                            {
#                               ifdef DEBUG_COUPLE_INTERSECTION
                                Info<< "adding irregular master "
                                    << "intersection2: "
                                    << points_[slaveEdges[slaveEdgeI].end()]
                                    << endl;
#                               endif

                                masterEdgePoints[masterEdgeI].append
                                (
                                    slaveEdges[slaveEdgeI].end()
                                );
                            }

                            // Slave side
                            // check if the first or second point of
                            // master edge is on the slave edge

                            vector sp = P - S;

                            scalar beta1 = (sp & e)/magSqr(e);

#                           ifdef DEBUG_COUPLE_INTERSECTION
                            Info<< "P: " << P << " S: " << S << " d: " << d
                                << " e: " << e << " sp: " << sp
                                << " beta1: " << beta1 << endl;
#                           endif

                            if
                            (
                                beta1 > -smallMergeTol_
                             && beta1 < 1 + smallMergeTol_
                            )
                            {
#                               ifdef DEBUG_COUPLE_INTERSECTION
                                Info<< "adding irregular slave "
                                    << "intersection1: "
                                    << points_[masterEdges[masterEdgeI].start()]
                                    << endl;
#                               endif

                                slaveEdgePoints[slaveEdgeI].append
                                (
                                    masterEdges[masterEdgeI].start()
                                );
                            }

                            scalar beta2 = ((sp + d) & e)/magSqr(e);

                            if
                            (
                                beta2 > -smallMergeTol_
                             && beta2 < 1 + smallMergeTol_
                            )
                            {
#                               ifdef DEBUG_COUPLE_INTERSECTION
                                Info<< "adding irregular slave "
                                    << "intersection2: "
                                    << points_[masterEdges[masterEdgeI].end()]
                                    << endl;
#                               endif

                                slaveEdgePoints[slaveEdgeI].append
                                (
                                    masterEdges[masterEdgeI].end()
                                );
                            }
                        } // end of colinear
                    } // end of singular intersection
                } // end of slave edges
            } // end of master edges

#           ifdef DEBUG_COUPLE_INTERSECTION
            Info<< "additional slave edge points: " << endl;
            forAll(slaveEdgePoints, edgeI)
            {
                Info<< "edge: " << edgeI << ": " << slaveEdgePoints[edgeI]
                    << endl;
            }
#           endif

            // Add new points
            if (nLivePoints + coupleFacePoints.size() >= points_.size())
            {
                // increase the size of the points list
                Info<< "Resizing points list" << endl;
                points_.setSize(points_.size() + couples_.size());
            }

            for
            (
                SLList<point>::iterator coupleFacePointsIter =
                    coupleFacePoints.begin();
                coupleFacePointsIter != coupleFacePoints.end();
                ++coupleFacePointsIter
            )
            {
                points_[nLivePoints] = coupleFacePointsIter();
                nLivePoints++;
            }

            // edge intersection finished

            // Creating new master side

            // count the number of additional points for face
            label nAdditionalMasterPoints = 0;

            forAll(masterEdgePoints, edgeI)
            {
                nAdditionalMasterPoints += masterEdgePoints[edgeI].size();
            }

            face tmpMasterFace
            (
                masterFace.size()
              + nAdditionalMasterPoints
            );
            label nTmpMasterLabels = 0;

#           ifdef DEBUG_COUPLE_INTERSECTION
            Info<< "masterFace: " << masterFace << endl
                << "nAdditionalMasterPoints: " << nAdditionalMasterPoints
                << endl;
#           endif

            forAll(masterEdges, masterEdgeI)
            {
                // Insert the starting point of the edge
                tmpMasterFace[nTmpMasterLabels] =
                    masterEdges[masterEdgeI].start();
                nTmpMasterLabels++;

                // get reference to added points of current edge
                const SLList<label>& curMEdgePoints =
                    masterEdgePoints[masterEdgeI];

                // create a markup list of points that have been used
                boolList usedMasterPoint(curMEdgePoints.size(), false);

                vector edgeVector = masterEdges[masterEdgeI].vec(points_);

#               ifdef DEBUG_FACE_ORDERING
                Info<< "edgeVector: " << edgeVector << endl
                    << "curMEdgePoints.size(): " << curMEdgePoints.size()
                    << endl;
#               endif

                // renormalise
                edgeVector /= magSqr(edgeVector);

                point edgeStartPoint =
                    points_[masterEdges[masterEdgeI].start()];

                // loop until the next label to add is -1
                for (;;)
                {
                    label nextPointLabel = -1;
                    label usedI = -1;
                    scalar minAlpha = GREAT;

                    label i = 0;

                    for
                    (
                        SLList<label>::const_iterator curMEdgePointsIter =
                            curMEdgePoints.begin();
                        curMEdgePointsIter != curMEdgePoints.end();
                        ++curMEdgePointsIter
                    )
                    {
                        if (!usedMasterPoint[i])
                        {
                            scalar alpha =
                                edgeVector
                              & (
                                    points_[curMEdgePointsIter()]
                                  - edgeStartPoint
                                );

#                           ifdef DEBUG_FACE_ORDERING
                            Info<< " edgeStartPoint: " << edgeStartPoint
                                << " edgeEndPoint: "
                                << points_[masterEdges[masterEdgeI].end()]
                                << " other point: "
                                << points_[curMEdgePointsIter()]
                                << " alpha: " << alpha << endl;
#                           endif

                            if (alpha < minAlpha)
                            {
                                minAlpha = alpha;
                                usedI = i;
                                nextPointLabel = curMEdgePointsIter();
                            }
                        }

#                       ifdef DEBUG_FACE_ORDERING
                        Info<< "nextPointLabel: " << nextPointLabel << endl;
#                       endif

                        i++;
                    }

                    if (nextPointLabel > -1)
                    {
#                       ifdef DEBUG_FACE_ORDERING
                        Info<< "added nextPointLabel: " << nextPointLabel
                            << " nTmpMasterLabels: " << nTmpMasterLabels
                            << " to place " << nTmpMasterLabels << endl;
#                       endif

                        usedMasterPoint[usedI] = true;
                        // add the next point
                        tmpMasterFace[nTmpMasterLabels] =
                            nextPointLabel;
                        nTmpMasterLabels++;
                    }
                    else
                    {
                        break;
                    }
                }
            }

            // reset the size of master
            tmpMasterFace.setSize(nTmpMasterLabels);

#           ifdef DEBUG_FACE_ORDERING
            Info<< "tmpMasterFace: " << tmpMasterFace << endl;
#           endif

            // Eliminate all zero-length edges
            face newMasterFace(labelList(tmpMasterFace.size(), labelMax));

            // insert first point by hand. Careful: the first one is
            // used for comparison to allow the edge collapse across
            // point zero.
            //
            newMasterFace[0] = tmpMasterFace[0];
            label nMaster = 0;

            edgeList mstEdgesToCollapse = tmpMasterFace.edges();

            scalar masterTol =
                cpMergePointTol_*boundBox(tmpMasterFace.points(points_)).mag();

            forAll(mstEdgesToCollapse, edgeI)
            {
#               ifdef DEBUG_FACE_ORDERING
                Info<< "edgeI: " << edgeI << " curEdge: "
                    << mstEdgesToCollapse[edgeI] << endl
                    << "master edge " << edgeI << ", "
                    << mstEdgesToCollapse[edgeI].mag(points_) << endl;
#               endif

                // Edge merge tolerance = masterTol
                if (mstEdgesToCollapse[edgeI].mag(points_) < masterTol)
                {
                    newMasterFace[nMaster] =
                        min
                        (
                            newMasterFace[nMaster],
                            mstEdgesToCollapse[edgeI].end()
                        );

#                   ifdef DEBUG_FACE_ORDERING
                    Info<< "Collapsed: nMaster: " << nMaster
                        << " label: " << newMasterFace[nMaster] << endl;
#                   endif

                }
                else
                {
                    nMaster++;

                    if (edgeI < mstEdgesToCollapse.size() - 1)
                    {
                        // last edge does not add the point
#                   ifdef DEBUG_FACE_ORDERING
                        Info<< "Added: nMaster: " << nMaster
                            << " label: " << mstEdgesToCollapse[edgeI].end()
                            << endl;
#                   endif

                        newMasterFace[nMaster] =
                            mstEdgesToCollapse[edgeI].end();
                    }
                }
            }

            newMasterFace.setSize(nMaster);

#           ifdef DEBUG_COUPLE
            Info<< "newMasterFace: " << newMasterFace << endl
                << "points: " << newMasterFace.points(points_) << endl;
#           endif

            // Creating new slave side

            // count the number of additional points for face
            label nAdditionalSlavePoints = 0;

            forAll(slaveEdgePoints, edgeI)
            {
                nAdditionalSlavePoints += slaveEdgePoints[edgeI].size();
            }

            face tmpSlaveFace
            (
                slaveFace.size()
              + nAdditionalSlavePoints
            );
            label nTmpSlaveLabels = 0;

#           ifdef DEBUG_COUPLE_INTERSECTION
            Info<< "slaveFace: " << slaveFace << endl
                << "nAdditionalSlavePoints: " << nAdditionalSlavePoints << endl;
#           endif

            forAll(slaveEdges, slaveEdgeI)
            {
                // Insert the starting point of the edge
                tmpSlaveFace[nTmpSlaveLabels] =
                    slaveEdges[slaveEdgeI].start();
                nTmpSlaveLabels++;

                // get reference to added points of current edge
                const SLList<label>& curSEdgePoints =
                    slaveEdgePoints[slaveEdgeI];

                // create a markup list of points that have been used
                boolList usedSlavePoint(curSEdgePoints.size(), false);

                vector edgeVector = slaveEdges[slaveEdgeI].vec(points_);

#               ifdef DEBUG_FACE_ORDERING
                Info<< "curSEdgePoints.size(): "
                    << curSEdgePoints.size() << endl
                    << "edgeVector: " << edgeVector << endl;
#               endif

                // renormalise
                edgeVector /= magSqr(edgeVector);

                point edgeStartPoint =
                    points_[slaveEdges[slaveEdgeI].start()];

                // loop until the next label to add is -1
                for (;;)
                {
                    label nextPointLabel = -1;
                    label usedI = -1;
                    scalar minAlpha = GREAT;

                    label i = 0;

                    for
                    (
                        SLList<label>::const_iterator curSEdgePointsIter =
                            curSEdgePoints.begin();
                        curSEdgePointsIter != curSEdgePoints.end();
                        ++curSEdgePointsIter
                    )
                    {
                        if (!usedSlavePoint[i])
                        {
                            scalar alpha =
                                edgeVector
                              & (
                                    points_[curSEdgePointsIter()]
                                  - edgeStartPoint
                                );

#                           ifdef DEBUG_FACE_ORDERING
                            Info<< " edgeStartPoint: " << edgeStartPoint
                                << " edgeEndPoint: "
                                << points_[slaveEdges[slaveEdgeI].end()]
                                << " other point: "
                                << points_[curSEdgePointsIter()]
                                << " alpha: " << alpha << endl;
#                           endif

                            if (alpha < minAlpha)
                            {
                                minAlpha = alpha;
                                usedI = i;
                                nextPointLabel = curSEdgePointsIter();
                            }
                        }

#                       ifdef DEBUG_FACE_ORDERING
                        Info<< "nextPointLabel: " << nextPointLabel << endl;
#                       endif

                        i++;
                    }

                    if (nextPointLabel > -1)
                    {
#                       ifdef DEBUG_FACE_ORDERING
                        Info<< "added nextPointLabel: " << nextPointLabel
                            << " nTmpSlaveLabels: " << nTmpSlaveLabels
                            << " to place " << nTmpSlaveLabels << endl;
#                       endif

                        usedSlavePoint[usedI] = true;
                        // add the next point
                        tmpSlaveFace[nTmpSlaveLabels] =
                            nextPointLabel;
                        nTmpSlaveLabels++;
                    }
                    else
                    {
                        break;
                    }
                }
            }

            // reset the size of slave
            tmpSlaveFace.setSize(nTmpSlaveLabels);

#           ifdef DEBUG_FACE_ORDERING
            Info<< "tmpSlaveFace: " << tmpSlaveFace << endl;
#           endif

            // Eliminate all zero-length edges
            face newSlaveFace(labelList(tmpSlaveFace.size(), labelMax));

            // insert first point by hand. Careful: the first one is
            // used for comparison to allow the edge collapse across
            // point zero.
            //
            newSlaveFace[0] = tmpSlaveFace[0];
            label nSlave = 0;

            edgeList slvEdgesToCollapse = tmpSlaveFace.edges();

            scalar slaveTol =
                cpMergePointTol_*boundBox(tmpSlaveFace.points(points_)).mag();

            forAll(slvEdgesToCollapse, edgeI)
            {
#               ifdef DEBUG_FACE_ORDERING
                Info<< "slave edge length: " << edgeI << ", "
                    << slvEdgesToCollapse[edgeI].mag(points_)<< endl;
#               endif

                 // edge merge tolerance = slaveTol
                if (slvEdgesToCollapse[edgeI].mag(points_) < slaveTol)
                {
                    newSlaveFace[nSlave] =
                        min
                        (
                            newSlaveFace[nSlave],
                            slvEdgesToCollapse[edgeI].end()
                        );
                }
                else
                {
                    nSlave++;
                    if (edgeI < slvEdgesToCollapse.size() - 1)
                    {
                        // last edge does not add the point
                        newSlaveFace[nSlave] = slvEdgesToCollapse[edgeI].end();
                    }
                }
            }

            newSlaveFace.setSize(nSlave);

#           ifdef DEBUG_COUPLE
            Info<< "newSlaveFace: " << newSlaveFace << endl
                << "points: " << newSlaveFace.points(points_) << endl << endl;
#           endif

            // Create the intersection face

            // Algorithm:
            // Loop through
            // points of the master and try to find one which falls
            // within the slave.  If not found, look through all
            // edges of the slave and find one which falls within the
            // master.  This point will be the starting location for
            // the cut face.

            edgeList newMasterEdges = newMasterFace.edges();
            edgeList newSlaveEdges = newSlaveFace.edges();

#           ifdef DEBUG_RIGHT_HAND_WALK
            Info<< "newMasterEdges: " << newMasterEdges << endl
                << "newSlaveEdges: " << newSlaveEdges << endl;
#           endif

            edge startEdge(-1, -1);

            // Remember where the start edge was found:
            // 0 for not found
            // 1 for master
            // 2 for slave
            label startEdgeFound = 0;

            vector masterProjDir = -newMasterFace.normal(points_);

            forAll(newSlaveEdges, edgeI)
            {
                // Take the slave edge points and project into the master.
                // In order to create a good intersection, move the
                // point away from the master in the direction of its
                // normal.
                point pointStart = points_[newSlaveEdges[edgeI].start()];

                point pointEnd = points_[newSlaveEdges[edgeI].end()];

                if
                (
                    newMasterFace.ray
                    (
                        pointStart,
                        masterProjDir,
                        points_,
                        intersection::FULL_RAY
                    ).hit()
                 && newMasterFace.ray
                    (
                        pointEnd,
                        masterProjDir,
                        points_,
                        intersection::FULL_RAY
                    ).hit()
                )
                {
                    startEdge = newSlaveEdges[edgeI];
                    startEdgeFound = 2;

#                   ifdef DEBUG_RIGHT_HAND_WALK
                    Info<< "slave edge found" << endl;
#                   endif

                    break;
                }
            }

            if (startEdgeFound == 0)
            {
                vector slaveProjDir = -newSlaveFace.normal(points_);

                forAll(newMasterEdges, edgeI)
                {
                    // Take the edge master points and project into the slave.
                    // In order to create a good intersection, move the
                    // point away from the slave in the direction of its
                    // normal.
                    point pointStart = points_[newMasterEdges[edgeI].start()];

                    point pointEnd = points_[newMasterEdges[edgeI].end()];

                    if
                    (
                        newSlaveFace.ray
                        (
                            pointStart,
                            slaveProjDir,
                            points_,
                        intersection::FULL_RAY
                        ).hit()
                     && newSlaveFace.ray
                        (
                            pointEnd,
                            slaveProjDir,
                            points_,
                            intersection::FULL_RAY
                        ).hit()
                    )
                    {
                        startEdge = newMasterEdges[edgeI];
                        startEdgeFound = 1;

#                       ifdef DEBUG_RIGHT_HAND_WALK
                        Info<< "master edge found" << endl;
#                       endif

                        break;
                    }
                }
            }

            // create the intersected face using right-hand walk rule
            face intersectedFace
            (
                labelList(newMasterFace.size() + newSlaveFace.size(), -1)
            );

            if (startEdgeFound > 0)
            {
#               ifdef DEBUG_RIGHT_HAND_WALK
                Info<< "start edge: " << startEdge << endl;
#               endif

                // Loop through both faces and add all edges
                // containing the current point and add them to the
                // list of edges to consider.  Make sure all edges are
                // added such that the current point is their start.
                // Loop through all edges to consider and find the one
                // which produces the buggest right-hand-turn.  This
                // is the next edge to be added to the face.  If its
                // end is the same as the starting point, the face is
                // complete; resize it to the number of active points
                // and exit.

                vector planeNormal = newMasterFace.normal(points_);
                planeNormal /= mag(planeNormal) + VSMALL;

#               ifdef DEBUG_RIGHT_HAND_WALK
                Info<< "planeNormal: " << planeNormal << endl;
#               endif

                // Do a check to control the right-hand turn.  This is
                // based on the triple product of the edge start
                // vector to face centre, the edge vector and the
                // plane normal.  If the triple product is negative,
                // the edge needs to be reversed to allow the
                // right-hand-turn rule to work.

                vector faceCentre;

                if (startEdgeFound == 1)
                {
                    faceCentre = newMasterFace.centre(points_);
                }
                else
                {
                    faceCentre = newSlaveFace.centre(points_);
                }

                scalar tripleProduct =
                    (
                        (faceCentre - points_[startEdge.start()])
                      ^ startEdge.vec(points_)
                    ) & planeNormal;

                if (tripleProduct < 0)
                {
#                   ifdef DEBUG_RIGHT_HAND_WALK
                    Info<< "Turning edge for right-hand turn rule" << endl;
#                   endif
                    startEdge.flip();
                }

                // prepare the loop for the right-hand walk
                intersectedFace[0] = startEdge.start();
                intersectedFace[1] = startEdge.end();
                label nIntFacePoints = 2;

                edge curEdge = startEdge;

                bool completedFace = false;

                do
                {
                    SLList<edge> edgesToConsider;

                    // collect master edges
                    forAll(newMasterEdges, edgeI)
                    {
                        const edge& cme = newMasterEdges[edgeI];

                        if (cme != curEdge)
                        {
                            if (cme.start() == curEdge.end())
                            {
                                edgesToConsider.append(cme);
                            }
                            else if (cme.end() == curEdge.end())
                            {
                                edgesToConsider.append(cme.reverseEdge());
                            }
                            // otherwise, it does not have the current point
                        }
                    }

                    // collect slave edges
                    forAll(newSlaveEdges, edgeI)
                    {
                        const edge& cse = newSlaveEdges[edgeI];

                        if (cse != curEdge)
                        {
                            if (cse.start() == curEdge.end())
                            {
                                edgesToConsider.append(cse);
                            }
                            else if (cse.end() == curEdge.end())
                            {
                                edgesToConsider.append(cse.reverseEdge());
                            }
                            // otherwise, it does not have the current point
                        }
                    }

#                   ifdef DEBUG_RIGHT_HAND_WALK
                    Info<< "number of edges to consider: "
                        << edgesToConsider.size() << endl
                        << "edges to consider: " << edgesToConsider << endl;
#                   endif

                    if (edgesToConsider.empty())
                    {
                        FatalErrorIn("void starMesh::createCoupleMatches()")
                            << setprecision(12)
                            << "void starMesh::createCoupleMatches() : "
                            << endl << "error in face intersection: "
                            << "no edges to consider for closing the loop"
                            << coupleI << ". STAR couple ID: "
                            << couples_[coupleI].coupleID() << endl
                            << "Cut Master Face: " << newMasterFace << endl
                            << "points: " << newMasterFace.points(points_)
                            << endl
                            << "Cut Slave Face: " << newSlaveFace << endl
                            << "points: " << newSlaveFace.points(points_)
                            << endl << "intersected face: "
                            << intersectedFace
                            << abort(FatalError);
                    }

                    // vector along the edge
                    vector ahead = curEdge.vec(points_);
                    ahead -= planeNormal*(planeNormal & ahead);
                    ahead /= mag(ahead) + VSMALL;

                    // vector pointing right
                    vector right = ahead ^ planeNormal;
                    right /= mag(right) + VSMALL;

                    // first edge taken for reference
                    edge nextEdge = edgesToConsider.first();
                    vector nextEdgeVec = nextEdge.vec(points_);
                    nextEdgeVec -= planeNormal*(planeNormal & nextEdgeVec);
                    nextEdgeVec /= mag(nextEdgeVec) + VSMALL;

                    scalar rightTurn = nextEdgeVec & right;
                    scalar goStraight = nextEdgeVec & ahead;

#                   ifdef DEBUG_RIGHT_HAND_WALK
                    Info<< "rightTurn: " << rightTurn
                        << " goStraight: " << goStraight << endl;
#                   endif

                    for
                    (
                        SLList<edge>::iterator etcIter =
                            edgesToConsider.begin();
                        etcIter != edgesToConsider.end();
                        ++etcIter
                    )
                    {
                        // right-hand walk rule
                        vector newDir = etcIter().vec(points_);
                        newDir -= planeNormal*(planeNormal & newDir);
                        newDir /= mag(newDir) + VSMALL;

                        scalar curRightTurn = newDir & right;
                        scalar curGoStraight = newDir & ahead;

#                       ifdef DEBUG_RIGHT_HAND_WALK
                        Info<< "curRightTurn: " << curRightTurn
                            << " curGoStraight: " << curGoStraight << endl;
#                       endif

                        if (rightTurn < 0) // old edge turning left
                        {
                            if (curRightTurn < 0) // new edge turning left
                            {
                                // both go left. Grab one with greater ahead
                                if (curGoStraight > goStraight)
                                {
#                                   ifdef DEBUG_RIGHT_HAND_WALK
                                    Info<< "a" << endl;
#                                   endif

                                    // Good edge, turning left less than before
                                    nextEdge = etcIter();
                                    rightTurn = curRightTurn;
                                    goStraight = curGoStraight;
                                }
                            }
                            else // new edge turning right
                            {
#                               ifdef DEBUG_RIGHT_HAND_WALK
                                Info<< "b" << endl;
#                               endif

                                // good edge, turning right
                                nextEdge = etcIter();
                                rightTurn = curRightTurn;
                                goStraight = curGoStraight;
                            }
                        }
                        else // old edge turning right
                        {
                            // new edge turning left rejected
                            if (curRightTurn >= 0) // new edge turning right
                            {
                                // grab one with smaller ahead
                                if (curGoStraight < goStraight)
                                {
#                                   ifdef DEBUG_RIGHT_HAND_WALK
                                    Info<< "c" << endl;
#                                   endif

                                    // good edge, turning right more than before
                                    nextEdge = etcIter();
                                    rightTurn = curRightTurn;
                                    goStraight = curGoStraight;
                                }
                            }
                        }
                    }

                    // check if the loop is completed
                    if (nextEdge.end() == intersectedFace[0])
                    {
                        // loop is completed. No point to add
                        completedFace = true;
                    }
                    else
                    {
                        // Check if there is room for more points
                        if (nIntFacePoints >= intersectedFace.size())
                        {
                            FatalErrorIn("void starMesh::createCoupleMatches()")
                                << setprecision(12)
                                << "void starMesh::createCoupleMatches() : "
                                << endl << "error in intersected face: "
                                << "lost thread for intersection of couple "
                                << coupleI << ". STAR couple ID: "
                                << couples_[coupleI].coupleID() << endl
                                << "Cut Master Face: " << newMasterFace << endl
                                << "points: " << newMasterFace.points(points_)
                                << endl
                                << "Cut Slave Face: " << newSlaveFace << endl
                                << "points: " << newSlaveFace.points(points_)
                                << endl << "intersected face: "
                                << intersectedFace
                                << abort(FatalError);
                        }

                        // insert the point
                        intersectedFace[nIntFacePoints] = nextEdge.end();
                        nIntFacePoints++;

                        // grab the current point and the current edge
                        curEdge = nextEdge;

#                       ifdef DEBUG_RIGHT_HAND_WALK
                        Info<< "inserted point " << nextEdge.end() << endl
                            << "curEdge: " << curEdge << endl;
#                       endif
                    }
                }
                while (!completedFace);

                // resize the face
                intersectedFace.setSize(nIntFacePoints);

#               ifdef DEBUG_COUPLE
                Info<< "intersectedFace: " << intersectedFace << endl;
#               endif

                // check the intersection face for duplicate points
                forAll(intersectedFace, checkI)
                {
                    for
                    (
                        label checkJ = checkI + 1;
                        checkJ < intersectedFace.size();
                        checkJ++
                    )
                    {
                        if (intersectedFace[checkI] == intersectedFace[checkJ])
                        {
                            FatalErrorIn("void starMesh::createCoupleMatches()")
                                << setprecision(12)
                                << "void starMesh::createCoupleMatches() : "
                                << endl << "error in intersected face: "
                                << "duplicate point in intersected face "
                                << "for couple no " << coupleI
                                << ". STAR couple ID: "
                                << couples_[coupleI].coupleID() << endl
                                << "Duplicate point label: "
                                << intersectedFace[checkI] << endl
                                << "Cut Master Face: " << newMasterFace << endl
                                << "points: " << newMasterFace.points(points_)
                                << endl
                                << "Cut Slave Face: " << newSlaveFace << endl
                                << "points: " << newSlaveFace.points(points_)
                                << endl << "intersected face: "
                                << intersectedFace
                                << abort(FatalError);
                        }
                    }
                }
            }
            else
            {
                FatalErrorIn("void starMesh::createCoupleMatches()")
                    << setprecision(12)
                    << "void starMesh::createCoupleMatches() : " << endl
                    << "could not find start edge for intersection of couple "
                    << coupleI << ". STAR couple ID: "
                    << couples_[coupleI].coupleID() << endl
                    << "Cut Master Face: " << newMasterFace << endl
                    << "points: " << newMasterFace.points(points_) << endl
                    << "Cut Slave Face: " << newSlaveFace << endl
                    << "points: " << newSlaveFace.points(points_)
                    << abort(FatalError);
            }

            // Project all points of the intersected face
            // onto the master face to ensure closedness
            vector pointProjectionNormal = -masterFace.normal(points_);

            forAll(intersectedFace, intPointI)
            {
#               ifdef DEBUG_COUPLE_PROJECTION
                Info<< "Proj: old point: "
                    << points_[intersectedFace[intPointI]] << endl;
#               endif

                pointHit projHit =
                    masterFace.ray
                    (
                        points_[intersectedFace[intPointI]],
                        pointProjectionNormal,
                        points_,
                        intersection::FULL_RAY
                    );

                if (projHit.hit())
                {
                    points_[intersectedFace[intPointI]] =
                        projHit.hitPoint();

#                   ifdef DEBUG_COUPLE_PROJECTION
                    Info<< "      new point: "
                        << points_[intersectedFace[intPointI]] << endl;
#                   endif
                }
            }

            // Check the direction of the intersection face
            if
            (
                (
                    masterFace.normal(points_)
                  & intersectedFace.normal(points_)
                ) < VSMALL
            )
            {
                intersectedFace.flip();
            }

            // Add the new face to both master and slave

            // Master face is replaced by a set of slave faces
            Map<SLList<label> >::iterator crfMasterIter =
                cellRemovedFaces.find(fp.masterCell());

            if (crfMasterIter == cellRemovedFaces.end())
            {
                cellRemovedFaces.insert
                (
                    fp.masterCell(),
                    SLList<label>(fp.masterFace())
                );
            }
            else
            {
                crfMasterIter().append(fp.masterFace());
            }

            Map<SLList<label> >::iterator crfSlaveIter =
                cellRemovedFaces.find(fp.slaveCell());

            if (crfSlaveIter == cellRemovedFaces.end())
            {
                cellRemovedFaces.insert
                (
                    fp.slaveCell(),
                    SLList<label>(fp.slaveFace())
                );
            }
            else
            {
                crfSlaveIter().append(fp.slaveFace());
            }

            Map<SLList<face> >::iterator cafMasterIter =
                cellAddedFaces.find(fp.masterCell());
            if (cafMasterIter == cellAddedFaces.end())
            {
                cellAddedFaces.insert
                (
                    fp.masterCell(),
                    SLList<face>(intersectedFace)
                );
            }
            else
            {
                cafMasterIter().append(intersectedFace);
            }

            Map<SLList<face> >::iterator cafSlaveIter =
                cellAddedFaces.find(fp.slaveCell());
            if (cafSlaveIter == cellAddedFaces.end())
            {
                cellAddedFaces.insert
                (
                    fp.slaveCell(),
                    SLList<face>(intersectedFace.reverseFace())
                );
            }
            else
            {
                cafSlaveIter().append(intersectedFace.reverseFace());
            }
        } // end of arbitrary match
    }

    if (couples_.size())
    {
        // Loop through all cells and reset faces for removal to zero size
        const labelList crfToc = cellRemovedFaces.toc();

        forAll(crfToc, cellI)
        {
            const label curCell = crfToc[cellI];

            const SLList<label>& curRemovedFaces = cellRemovedFaces[curCell];

            for
            (
                SLList<label>::const_iterator curRemovedFacesIter =
                    curRemovedFaces.begin();
                curRemovedFacesIter != curRemovedFaces.end();
                ++curRemovedFacesIter
            )
            {
                cellFaces_[curCell][curRemovedFacesIter()].setSize(0);
            }

            if (curRemovedFaces.size())
            {
                // reset the shape pointer to unknown
                cellShapes_[curCell] = cellShape(*unknownPtr_, labelList(0));
            }
        }

        const labelList cafToc = cellAddedFaces.toc();

        // Insert the new faces into the list
        forAll(cafToc, cellI)
        {
            const label curCell = cafToc[cellI];

            const SLList<face>& curAddedFaces = cellAddedFaces[curCell];

            faceList oldFaces = cellFaces_[curCell];

            faceList& newFaces = cellFaces_[curCell];

            newFaces.setSize(oldFaces.size() + curAddedFaces.size());
            label nNewFaces = 0;

            // copy original faces that have not been removed
            forAll(oldFaces, faceI)
            {
                if (oldFaces[faceI].size())
                {
                    newFaces[nNewFaces] = oldFaces[faceI];
                    nNewFaces++;
                }
            }

            // add new faces
            for
            (
                SLList<face>::const_iterator curAddedFacesIter =
                    curAddedFaces.begin();
                curAddedFacesIter != curAddedFaces.end();
                ++curAddedFacesIter
            )
            {
                newFaces[nNewFaces] = curAddedFacesIter();
                nNewFaces++;
            }

            // reset the size of the face list
            newFaces.setSize(nNewFaces);

            if (curAddedFaces.size())
            {
                // reset the shape pointer to unknown
                cellShapes_[curCell] = cellShape(*unknownPtr_, labelList(0));
            }
        }

        // Resize the point list to the number of created points
        points_.setSize(nLivePoints);

        // Finished
        Info<< "Finished doing couples" << endl;
    }
}


// ************************************************************************* //
