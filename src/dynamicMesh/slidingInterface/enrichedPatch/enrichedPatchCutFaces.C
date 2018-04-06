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

Description
    Calculating cut faces of the enriched patch, together with the addressing
    into master and slave patch.

\*---------------------------------------------------------------------------*/

#include "enrichedPatch.H"
#include "boolList.H"
#include "DynamicList.H"
#include "labelPair.H"
#include "primitiveMesh.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// If the cut face gets more than this number of points, it will be checked
const Foam::label Foam::enrichedPatch::maxFaceSizeDebug_ = 100;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Index of debug signs:
// x - skip a point
// l - left turn
// r - right turn

void Foam::enrichedPatch::calcCutFaces() const
{
    if (cutFacesPtr_ || cutFaceMasterPtr_ || cutFaceSlavePtr_)
    {
        FatalErrorInFunction
            << "Cut faces addressing already calculated."
            << abort(FatalError);
    }

    const faceList& enFaces = enrichedFaces();
    const labelList& mp = meshPoints();
    const faceList& lf = localFaces();
    const pointField& lp = localPoints();
    const labelListList& pp = pointPoints();
    // Pout<< "enFaces: " << enFaces << endl;
    // Pout<< "lf: " << lf << endl;
    // Pout<< "lp: " << lp << endl;
    // Pout<< "pp: " << pp << endl;
    const Map<labelList>& masterPointFaceAddr = masterPointFaces();

    // Prepare the storage
    DynamicList<face> cf(2*lf.size());
    DynamicList<label> cfMaster(2*lf.size());
    DynamicList<label> cfSlave(2*lf.size());

    // Algorithm
    // Go through all the faces
    // 1) For each face, start at point zero and grab the edge zero.
    //    Both points are added into the cut face.  Mark the first edge
    //    as used and position the current point as the end of the first
    //    edge and previous point as point zero.
    // 2) Grab all possible points from the current point.  Excluding
    //    the previous point find the point giving the biggest right
    //    turn. Add the point to the face and mark the edges as used.
    //    Continue doing this until the face is complete, i.e. the next point
    //    to add is the first point of the face.
    // 3) Once the facet is completed, register its owner from the face
    //    it has been created from (remember that the first lot of faces
    //    inserted into the enriched faces list are the slaves, to allow
    //    matching of the other side).
    // 4) If the facet is created from an original slave face, go
    //    through the master patch and try to identify the master face
    //    this facet belongs to.  This is recognised by the fact that the
    //    faces consists exclusively of the points of the master face and
    //    the points projected onto the face.

    // Create a set of edge usage parameters
    HashSet<edge, Hash<edge>> edgesUsedOnce(pp.size());
    HashSet<edge, Hash<edge>> edgesUsedTwice
        (pp.size()*primitiveMesh::edgesPerPoint_);


    forAll(lf, facei)
    {
        const face& curLocalFace = lf[facei];
        const face& curGlobalFace = enFaces[facei];

        // Pout<< "Doing face " << facei
        //     << " local: " << curLocalFace
        //     << " or " << curGlobalFace
        //     << endl;

        // if (facei < slavePatch_.size())
        // {
        //     Pout<< "original slave: " << slavePatch_[facei]
        //         << " local: " << slavePatch_.localFaces()[facei] << endl;
        // }
        // else
        // {
        //     Pout<< "original master: "
        //         << masterPatch_[facei - slavePatch_.size()] << " "
        //         << masterPatch_.localFaces()[facei - slavePatch_.size()]
        //         << endl;
        // }
        // {
        //     pointField facePoints = curLocalFace.points(lp);
        //     forAll(curLocalFace, pointi)
        //     {
        //         Pout<< "v " << facePoints[pointi].x() << " "
        //             << facePoints[pointi].y() << " "
        //             << facePoints[pointi].z() << endl;
        //     }
        // }

        // Track the usage of face edges.  When all edges are used, the
        // face decomposition is complete.
        // The seed edges include all the edges of the original face + all edges
        // of other faces that have been used in the construction of the
        // facet.  Edges from other faces can be considered as
        // internal to the current face if used only once.

        // Track the edge usage to avoid duplicate faces and reset it to unused
        boolList usedFaceEdges(curLocalFace.size(), false);

        SLList<edge> edgeSeeds;

        // Insert the edges of current face into the seed list.
        edgeList cfe = curLocalFace.edges();
        forAll(curLocalFace, edgeI)
        {
            edgeSeeds.append(cfe[edgeI]);
        }

        // Grab face normal
        const vector normal = curLocalFace.normal(lp);

        while (edgeSeeds.size())
        {
            // Pout<< "edgeSeeds.size(): "
            //     << edgeSeeds.size()
            //     << endl;

            const edge curEdge = edgeSeeds.removeHead();

            // Locate the edge in current face
            const label curEdgeWhich = curLocalFace.which(curEdge.start());

            // Check if the edge is in current face and if it has been
            // used already.  If so, skip it
            if
            (
                curEdgeWhich > -1
             && curLocalFace.nextLabel(curEdgeWhich) == curEdge.end()
            )
            {
                // Edge is in the starting face.
                // If unused, mark it as used; if used, skip it
                if (usedFaceEdges[curEdgeWhich]) continue;

                usedFaceEdges[curEdgeWhich] = true;
            }

            // If the edge has already been used twice, skip it
            if (edgesUsedTwice.found(curEdge)) continue;

            // Pout<< "Trying new edge (" << mp[curEdge.start()]
            //     << ", " << mp[curEdge.end()]
            //     << ") seed: " << curEdge
            //     << " used: " << edgesUsedTwice.found(curEdge)
            //     << endl;

            // Estimate the size of cut face as twice the size of original face
            DynamicList<label> cutFaceGlobalPoints(2*curLocalFace.size());
            DynamicList<label> cutFaceLocalPoints(2*curLocalFace.size());

            // Found unused edge.
            label prevPointLabel = curEdge.start();
            cutFaceGlobalPoints.append(mp[prevPointLabel]);
            cutFaceLocalPoints.append(prevPointLabel);
            // Pout<< "prevPointLabel: " << mp[prevPointLabel] << endl;
            // Grab current point and append it to the list
            label curPointLabel = curEdge.end();
            point curPoint = lp[curPointLabel];
            cutFaceGlobalPoints.append(mp[curPointLabel]);
            cutFaceLocalPoints.append(curPointLabel);

            bool completedCutFace = false;

            label faceSizeDebug = maxFaceSizeDebug_;

            do
            {
                // Grab the next point options

                // Pout<< "curPointLabel: " << mp[curPointLabel] << endl;

                const labelList& nextPoints = pp[curPointLabel];

                // Pout<< "nextPoints: "
                //     << UIndirectList<label>(mp, nextPoints)
                //     << endl;

                // Get the vector along the edge and the right vector
                vector ahead = curPoint - lp[prevPointLabel];
                ahead -= normal*(normal & ahead);
                ahead /= mag(ahead);

                vector right = normal ^ ahead;
                right /= mag(right);

                // Pout<< "normal: " << normal
                //     << " ahead: " << ahead
                //     << " right: " << right
                //     << endl;

                scalar atanTurn = -great;
                label bestAtanPoint = -1;

                forAll(nextPoints, nextI)
                {
                    // Exclude the point we are coming from; there will always
                    // be more than one edge, so this is safe
                    if (nextPoints[nextI] != prevPointLabel)
                    {
                        // Pout<< "cur point: " << curPoint
                        //     << " trying for point: "
                        //     << mp[nextPoints[nextI]]
                        //     << " " << lp[nextPoints[nextI]];
                        vector newDir = lp[nextPoints[nextI]] - curPoint;
                        // Pout<< " newDir: " << newDir
                        //     << " mag: " << mag(newDir) << flush;
                        newDir -= normal*(normal & newDir);
                        scalar magNewDir = mag(newDir);
                        // Pout<< " corrected: " << newDir
                        //     << " mag: " << mag(newDir) << flush;

                        if (magNewDir < small)
                        {
                            FatalErrorInFunction
                                << "projection error: slave patch probably "
                                << "does not project onto master.  "
                                << "Please switch on "
                                << "enriched patch debug for more info"
                                << abort(FatalError);
                        }

                        newDir /= magNewDir;

                        scalar curAtanTurn =
                            atan2(newDir & right, newDir & ahead);

                        // Pout<< " atan: " << curAtanTurn << endl;

                        if (curAtanTurn > atanTurn)
                        {
                            bestAtanPoint = nextPoints[nextI];
                            atanTurn = curAtanTurn;
                        }
                    } // end of prev point skip
                } // end of next point selection

                // Pout<< "   bestAtanPoint: " << bestAtanPoint << " or "
                //     << mp[bestAtanPoint]
                //     << endl;

                // Selected next best point.

                // Pout<< "cutFaceGlobalPoints: "
                //     << cutFaceGlobalPoints
                //     << endl;

                // Check if the edge about to be added has been used
                // in the current face or twice in other faces.  If
                // so, the face is bad.
                const label whichNextPoint = curLocalFace.which(curPointLabel);

                if
                (
                    whichNextPoint > -1
                 && curLocalFace.nextLabel(whichNextPoint) == bestAtanPoint
                 && usedFaceEdges[whichNextPoint]
                )
                {
                    // This edge is already used in current face
                    // face cannot be good; start on a new one

                    // Pout<< "Double usage in current face, cannot be good"
                    //     << endl;

                    completedCutFace = true;
                }


                if (edgesUsedTwice.found(edge(curPointLabel, bestAtanPoint)))
                {
                    // This edge is already used -
                    // face cannot be good; start on a new one
                    completedCutFace = true;

                    // Pout<< "Double usage elsewhere, cannot be good" << endl;
                }

                if (completedCutFace) continue;

                // Insert the next best point into the list
                if (mp[bestAtanPoint] == cutFaceGlobalPoints[0])
                {
                    // Face is completed.  Save it and move on to the next
                    // available edge
                    completedCutFace = true;

                    if (debug)
                    {
                        Pout<< " local: " << cutFaceLocalPoints
                            << " one side: " << facei;
                    }

                    // Append the face
                    face cutFaceGlobal;
                    cutFaceGlobal.transfer(cutFaceGlobalPoints);

                    cf.append(cutFaceGlobal);

                    face cutFaceLocal;
                    cutFaceLocal.transfer(cutFaceLocalPoints);

                    // Pout<< "\ncutFaceLocal: "
                    //     << cutFaceLocal.points(lp)
                    //     << endl;

                    // Go through all edges of the cut faces.
                    // If the edge corresponds to a starting face edge,
                    // mark the starting face edge as true

                    forAll(cutFaceLocal, cutI)
                    {
                        const edge curCutFaceEdge
                        (
                            cutFaceLocal[cutI],
                            cutFaceLocal.nextLabel(cutI)
                        );

                        // Increment the usage count using two hash sets
                        HashSet<edge, Hash<edge>>::iterator euoIter =
                            edgesUsedOnce.find(curCutFaceEdge);

                        if (euoIter == edgesUsedOnce.end())
                        {
                            // Pout<< "Found edge not used before: "
                            //     << curCutFaceEdge
                            //     << endl;
                            edgesUsedOnce.insert(curCutFaceEdge);
                        }
                        else
                        {
                            // Pout<< "Found edge used once: "
                            //     << curCutFaceEdge
                            //     << endl;
                            edgesUsedOnce.erase(euoIter);
                            edgesUsedTwice.insert(curCutFaceEdge);
                        }

                        const label curCutFaceEdgeWhich = curLocalFace.which
                        (
                            curCutFaceEdge.start()
                        );

                        if
                        (
                            curCutFaceEdgeWhich > -1
                         && curLocalFace.nextLabel(curCutFaceEdgeWhich)
                         == curCutFaceEdge.end()
                        )
                        {
                            // Found edge in original face

                            // Pout<< "Found edge in orig face: "
                            //     << curCutFaceEdge << ": "
                            //     << curCutFaceEdgeWhich
                            //     << endl;

                            usedFaceEdges[curCutFaceEdgeWhich] = true;
                        }
                        else
                        {
                            // Edge not in original face.  Add it to seeds

                            // Pout<< "Found new edge seed: "
                            //     << curCutFaceEdge
                            //     << endl;

                            edgeSeeds.append(curCutFaceEdge.reverseEdge());
                        }
                    }

                    // Find out what the other side is

                    // Algorithm

                    // If the face is in the slave half of the
                    // enrichedFaces lists, it may be matched against
                    // the master face.  It can be recognised by the
                    // fact that all its points belong to one master
                    // face.  For this purpose:
                    // 1) Grab the master faces around the point zero
                    // of the cut face and collect all master faces
                    // around it.
                    // 2) Loop through the rest of cut face points and
                    // if the point refers to one of the faces marked
                    // by point zero, increment its count.
                    // 3) When completed, the face whose count is
                    // equal to the number of points in the cut face
                    // is the other side.  If this is not the case, there is no
                    // face on the other side.

                    if (facei < slavePatch_.size())
                    {
                        Map<labelList>::const_iterator mpfAddrIter =
                            masterPointFaceAddr.find(cutFaceGlobal[0]);

                        bool otherSideFound = false;

                        if (mpfAddrIter != masterPointFaceAddr.end())
                        {
                            bool miss = false;

                            // Create the label count using point zero
                            const labelList& masterFacesOfPZero = mpfAddrIter();

                            labelList hits(masterFacesOfPZero.size(), 1);

                            for
                            (
                                label pointi = 1;
                                pointi < cutFaceGlobal.size();
                                pointi++
                            )
                            {
                                Map<labelList>::const_iterator
                                    mpfAddrPointIter =
                                        masterPointFaceAddr.find
                                        (
                                            cutFaceGlobal[pointi]
                                        );

                                if
                                (
                                    mpfAddrPointIter
                                 == masterPointFaceAddr.end()
                                )
                                {
                                    // Point is off the master patch. Skip
                                    miss = true;
                                    break;
                                }

                                const labelList& curMasterFaces =
                                    mpfAddrPointIter();

                                // For every current face, try to find it in the
                                // zero-list
                                forAll(curMasterFaces, i)
                                {
                                    forAll(masterFacesOfPZero, j)
                                    {
                                        if
                                        (
                                            curMasterFaces[i]
                                         == masterFacesOfPZero[j]
                                        )
                                        {
                                            hits[j]++;
                                            break;
                                        }
                                    }
                                }
                            }

                            // If all point are found attempt matching
                            if (!miss)
                            {
                                forAll(hits, pointi)
                                {
                                    if (hits[pointi] == cutFaceGlobal.size())
                                    {
                                        // Found other side.
                                        otherSideFound = true;

                                        cfMaster.append
                                        (
                                            masterFacesOfPZero[pointi]
                                        );

                                        cfSlave.append(facei);

                                        // Reverse the face such that it
                                        // points out of the master patch
                                        cf.last().flip();

                                        if (debug)
                                        {
                                            Pout<< " other side: "
                                                << masterFacesOfPZero[pointi]
                                                << endl;
                                        }
                                    } // end of hits
                                } // end of for all hits

                            } // end of not miss

                            if (!otherSideFound || miss)
                            {
                                if (debug)
                                {
                                    Pout<< " solo slave A" << endl;
                                }

                                cfMaster.append(-1);
                                cfSlave.append(facei);
                            }
                        }
                        else
                        {
                            // First point not in master patch
                            if (debug)
                            {
                                Pout<< " solo slave B" << endl;
                            }

                            cfMaster.append(-1);
                            cfSlave.append(facei);
                        }
                    }
                    else
                    {
                        if (debug)
                        {
                            Pout<< " master side" << endl;
                        }

                        cfMaster.append(facei - slavePatch_.size());
                        cfSlave.append(-1);
                    }
                }
                else
                {
                    // Face not completed.  Prepare for the next point search
                    prevPointLabel = curPointLabel;
                    curPointLabel = bestAtanPoint;
                    curPoint = lp[curPointLabel];
                    cutFaceGlobalPoints.append(mp[curPointLabel]);
                    cutFaceLocalPoints.append(curPointLabel);

                    if (debug || cutFaceGlobalPoints.size() > faceSizeDebug)
                    {
                        faceSizeDebug *= 2;

                        // Check for duplicate points in the face
                        forAll(cutFaceGlobalPoints, checkI)
                        {
                            for
                            (
                                label checkJ = checkI + 1;
                                checkJ < cutFaceGlobalPoints.size();
                                checkJ++
                            )
                            {
                                if
                                (
                                    cutFaceGlobalPoints[checkI]
                                 == cutFaceGlobalPoints[checkJ]
                                )
                                {
                                    // Shrink local points for debugging
                                    cutFaceLocalPoints.shrink();

                                    face origFace;
                                    face origFaceLocal;
                                    if (facei < slavePatch_.size())
                                    {
                                        origFace = slavePatch_[facei];
                                        origFaceLocal =
                                            slavePatch_.localFaces()[facei];
                                    }
                                    else
                                    {
                                        origFace =
                                            masterPatch_
                                            [facei - slavePatch_.size()];

                                        origFaceLocal =
                                            masterPatch_.localFaces()
                                            [facei - slavePatch_.size()];
                                    }

                                    FatalErrorInFunction
                                        << "Duplicate point found in cut face. "
                                        << "Error in the face cutting "
                                        << "algorithm for global face "
                                        << origFace << " local face "
                                        << origFaceLocal << nl
                                        << "Slave size: " << slavePatch_.size()
                                        << " Master size: "
                                        << masterPatch_.size()
                                        << " index: " << facei << ".\n"
                                        << "Face: " << curGlobalFace << nl
                                        << "Cut face: " << cutFaceGlobalPoints
                                        << " local: " << cutFaceLocalPoints
                                        << nl << "Points: "
                                        << face(cutFaceLocalPoints).points(lp)
                                        << abort(FatalError);
                                }
                            }
                        }
                    } // end of debug
                }
            } while (!completedCutFace);
        } // end of used edges

        if (debug)
        {
            Pout<< " Finished face " << facei << endl;
        }

    } // end of local faces

    // Re-pack the list into compact storage
    cutFacesPtr_ = new faceList();
    cutFacesPtr_->transfer(cf);

    cutFaceMasterPtr_ = new labelList();
    cutFaceMasterPtr_->transfer(cfMaster);

    cutFaceSlavePtr_ = new labelList();
    cutFaceSlavePtr_->transfer(cfSlave);
}


void Foam::enrichedPatch::clearCutFaces()
{
    deleteDemandDrivenData(cutFacesPtr_);
    deleteDemandDrivenData(cutFaceMasterPtr_);
    deleteDemandDrivenData(cutFaceSlavePtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::faceList& Foam::enrichedPatch::cutFaces() const
{
    if (!cutFacesPtr_)
    {
        calcCutFaces();
    }

    return *cutFacesPtr_;
}


const Foam::labelList& Foam::enrichedPatch::cutFaceMaster() const
{
    if (!cutFaceMasterPtr_)
    {
        calcCutFaces();
    }

    return *cutFaceMasterPtr_;
}


const Foam::labelList& Foam::enrichedPatch::cutFaceSlave() const
{
    if (!cutFaceSlavePtr_)
    {
        calcCutFaces();
    }

    return *cutFaceSlavePtr_;
}


// ************************************************************************* //
