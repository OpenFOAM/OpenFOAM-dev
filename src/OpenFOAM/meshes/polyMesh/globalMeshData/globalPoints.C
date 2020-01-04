/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "globalPoints.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "polyMesh.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(globalPoints, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::globalPoints::countPatchPoints
(
    const polyBoundaryMesh& patches
)
{
    label nTotPoints = 0;

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];
        if (pp.coupled())
        {
            nTotPoints += pp.nPoints();
        }
    }
    return nTotPoints;
}


Foam::label Foam::globalPoints::findSamePoint
(
    const labelPairList& allInfo,
    const labelPair& info
) const
{
    const label proci = globalTransforms_.processor(info);
    const label index = globalTransforms_.index(info);

    forAll(allInfo, i)
    {
        if
        (
            globalTransforms_.processor(allInfo[i]) == proci
         && globalTransforms_.index(allInfo[i]) == index
        )
        {
            return i;
        }
    }
    return -1;
}


Foam::labelPairList Foam::globalPoints::addSendTransform
(
    const label patchi,
    const labelPairList& info
) const
{
    scalar tol = refCast<const coupledPolyPatch>
    (
        mesh_.boundaryMesh()[patchi]
    ).matchTolerance();

    labelPairList sendInfo(info.size());

    forAll(info, i)
    {
        // Pout<< "    adding send transform to" << nl
        //    << "    proc:" << globalTransforms_.processor(info[i])
        //    << nl
        //    << "    index:" << globalTransforms_.index(info[i]) << nl
        //    << "    trafo:"
        //    <<  globalTransforms_.decodeTransformIndex
        //        (globalTransforms_.transformIndex(info[i]))
        //    << endl;

        sendInfo[i] = globalTransforms_.encode
        (
            globalTransforms_.processor(info[i]),
            globalTransforms_.index(info[i]),
            globalTransforms_.addToTransformIndex
            (
                globalTransforms_.transformIndex(info[i]),
                patchi,
                true,           // patchi is sending side
                tol             // tolerance for comparison
            )
        );
    }
    return sendInfo;
}


void Foam::globalPoints::addToSend
(
    const polyPatch& pp,
    const label patchPointi,
    const labelPairList& knownInfo,

    DynamicList<label>& patchFaces,
    DynamicList<label>& indexInFace,
    DynamicList<labelPairList>& allInfo
) const
{
    // Collect all topological information about a point on a patch.  (this
    // information is the patch faces using the point and the relative position
    // of the point in the face)

    label meshPointi = pp.meshPoints()[patchPointi];

    // Add all faces using the point so we are sure we find it on the
    // other side.
    const labelList& pFaces = pp.pointFaces()[patchPointi];

    forAll(pFaces, i)
    {
        label patchFacei = pFaces[i];

        const face& f = pp[patchFacei];

        patchFaces.append(patchFacei);
        indexInFace.append(findIndex(f, meshPointi));

        // Add patch transformation
        allInfo.append(addSendTransform(pp.index(), knownInfo));
    }
}


bool Foam::globalPoints::mergeInfo
(
    const labelPairList& nbrInfo,
    const label localPointi,
    labelPairList& myInfo
) const
{
    // Add nbrInfo to myInfo. Return true if anything changed.  nbrInfo is for a
    // point a list of all the global points using it

    bool anyChanged = false;

    // Extend to make space for the nbrInfo (trimmed later)
    labelPairList newInfo(myInfo);
    label newI = newInfo.size();
    newInfo.setSize(newI + nbrInfo.size());

    forAll(nbrInfo, i)
    {
        // Check if already have information about nbr point. There are two
        // possibilities:
        // - information found about same point but different transform.
        //   Combine transforms
        // - information not found.

        label index = findSamePoint(myInfo, nbrInfo[i]);

        if (index == -1)
        {
            // New point
            newInfo[newI++] = nbrInfo[i];
            anyChanged = true;
        }
        else
        {
            // Same point. So we already have a connection between localPointi
            // and the nbrIndex. Two situations:
            // - same transform
            // - one transform takes two steps, the other just a single.
            if (myInfo[index] == nbrInfo[i])
            {
                // Everything same (so also transform). Nothing changed.
            }
            else
            {
                label myTransform = globalTransforms_.transformIndex
                (
                    myInfo[index]
                );
                label nbrTransform = globalTransforms_.transformIndex
                (
                    nbrInfo[i]
                );

                // Different transform. See which is 'simplest'.
                label minTransform = globalTransforms_.minimumTransformIndex
                (
                    myTransform,
                    nbrTransform
                );

                if (minTransform != myTransform)
                {
                    // Use nbr info.
                    newInfo[index] = nbrInfo[i];
                    anyChanged = true;
                }
            }
        }
    }

    newInfo.setSize(newI);
    myInfo.transfer(newInfo);

    return anyChanged;
}


Foam::label Foam::globalPoints::meshToLocalPoint
(
    const Map<label>& meshToPatchPoint, // from mesh point to local numbering
    const label meshPointi
)
{
    return
    (
        meshToPatchPoint.size() == 0
      ? meshPointi
      : meshToPatchPoint[meshPointi]
    );
}


Foam::label Foam::globalPoints::localToMeshPoint
(
    const labelList& patchToMeshPoint,
    const label localPointi
)
{
    return
    (
        patchToMeshPoint.size() == 0
      ? localPointi
      : patchToMeshPoint[localPointi]
    );
}


bool Foam::globalPoints::mergeInfo
(
    const labelPairList& nbrInfo,
    const label localPointi
)
{
    // Updates database of current information on meshpoints with nbrInfo.  Uses
    // mergeInfo above. Returns true if data kept for meshPointi changed.

    label infoChanged = false;

    // Get the index into the procPoints list.
    Map<label>::iterator iter = meshToProcPoint_.find(localPointi);

    if (iter != meshToProcPoint_.end())
    {
        if (mergeInfo(nbrInfo, localPointi, procPoints_[iter()]))
        {
            infoChanged = true;
        }
    }
    else
    {
        // Construct local index for point
        labelPairList knownInfo
        (
            1,
            globalTransforms_.encode
            (
                Pstream::myProcNo(),
                localPointi,
                globalTransforms_.nullTransformIndex()
            )
        );

        if (mergeInfo(nbrInfo, localPointi, knownInfo))
        {
            // Update addressing from into procPoints
            meshToProcPoint_.insert(localPointi, procPoints_.size());
            // Insert into list of equivalences.
            procPoints_.append(knownInfo);

            infoChanged = true;
        }
    }
    return infoChanged;
}


bool Foam::globalPoints::storeInitialInfo
(
    const labelPairList& nbrInfo,
    const label localPointi
)
{
    // Updates database of current information on meshpoints with nbrInfo.  Uses
    // mergeInfo above. Returns true if data kept for meshPointi changed.

    label infoChanged = false;

    // Get the index into the procPoints list.
    Map<label>::iterator iter = meshToProcPoint_.find(localPointi);

    if (iter != meshToProcPoint_.end())
    {
        if (mergeInfo(nbrInfo, localPointi, procPoints_[iter()]))
        {
            infoChanged = true;
        }
    }
    else
    {
        // Update addressing into procPoints
        meshToProcPoint_.insert(localPointi, procPoints_.size());
        // Insert into list of equivalences.
        procPoints_.append(nbrInfo);

        infoChanged = true;
    }
    return infoChanged;
}


void Foam::globalPoints::printProcPoint
(
    const labelList& patchToMeshPoint,
    const labelPair& pointInfo
) const
{
    label proci = globalTransforms_.processor(pointInfo);
    label index = globalTransforms_.index(pointInfo);
    label trafoI = globalTransforms_.transformIndex(pointInfo);

    Pout<< "    proc:" << proci;
    Pout<< " localpoint:";
    Pout<< index;
    Pout<< " through transform:"
        << trafoI << " bits:"
        << globalTransforms_.decodeTransformIndex(trafoI);

    if (proci == Pstream::myProcNo())
    {
        label meshPointi = localToMeshPoint(patchToMeshPoint, index);
        Pout<< " at:" <<  mesh_.points()[meshPointi];
    }
}


void Foam::globalPoints::printProcPoints
(
    const labelList& patchToMeshPoint,
    const labelPairList& pointInfo
) const
{
    forAll(pointInfo, i)
    {
        printProcPoint(patchToMeshPoint, pointInfo[i]);
        Pout<< endl;
    }
}


void Foam::globalPoints::initOwnPoints
(
    const Map<label>& meshToPatchPoint,
    const bool allPoints,
    labelHashSet& changedPoints
)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled())
        {
            const labelList& meshPoints = pp.meshPoints();

            if (allPoints)
            {
                // All points on patch
                forAll(meshPoints, patchPointi)
                {
                    label meshPointi = meshPoints[patchPointi];
                    label localPointi = meshToLocalPoint
                    (
                        meshToPatchPoint,
                        meshPointi
                    );

                    labelPairList knownInfo
                    (
                        1,
                        globalTransforms_.encode
                        (
                            Pstream::myProcNo(),
                            localPointi,
                            globalTransforms_.nullTransformIndex()
                        )
                    );

                    // Pout<< "For point "<< pp.points()[meshPointi]
                    //    << " inserting info " << knownInfo
                    //    << endl;

                    // Update changedpoints info.
                    if (storeInitialInfo(knownInfo, localPointi))
                    {
                        changedPoints.insert(localPointi);
                    }
                }
            }
            else
            {
                // Boundary points only
                const labelList& boundaryPoints = pp.boundaryPoints();

                forAll(boundaryPoints, i)
                {
                    label meshPointi = meshPoints[boundaryPoints[i]];
                    label localPointi = meshToLocalPoint
                    (
                        meshToPatchPoint,
                        meshPointi
                    );

                    labelPairList knownInfo
                    (
                        1,
                        globalTransforms_.encode
                        (
                            Pstream::myProcNo(),
                            localPointi,
                            globalTransforms_.nullTransformIndex()
                        )
                    );

                    if (storeInitialInfo(knownInfo, localPointi))
                    {
                        changedPoints.insert(localPointi);
                    }
                }
            }
        }
    }
}


void Foam::globalPoints::sendPatchPoints
(
    const bool mergeSeparated,
    const Map<label>& meshToPatchPoint,
    PstreamBuffers& pBufs,
    const labelHashSet& changedPoints
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const labelPairList& patchInfo = globalTransforms_.patchTransformSign();

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        // mergeSeparated=true : send from all processor patches
        //               =false: send from ones without transform

        if
        (
            (Pstream::parRun() && isA<processorPolyPatch>(pp))
         && (mergeSeparated || patchInfo[patchi].first() == -1)
        )
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pp);

            // Information to send:
            // patch face
            DynamicList<label> patchFaces(pp.nPoints());
            // index in patch face
            DynamicList<label> indexInFace(pp.nPoints());
            // all information I currently hold about this patchPoint
            DynamicList<labelPairList> allInfo(pp.nPoints());


            // Now collect information on all points mentioned in
            // changedPoints. Note that these points only should occur on
            // processorPatches (or rather this is a limitation!).

            const labelList& meshPoints = pp.meshPoints();

            forAll(meshPoints, patchPointi)
            {
                label meshPointi = meshPoints[patchPointi];
                label localPointi = meshToLocalPoint
                (
                    meshToPatchPoint,
                    meshPointi
                );

                if (changedPoints.found(localPointi))
                {
                    label index = meshToProcPoint_[localPointi];

                    const labelPairList& knownInfo = procPoints_[index];

                    // Add my information about localPointi to the
                    // send buffers. Encode the transformation
                    addToSend
                    (
                        pp,
                        patchPointi,
                        knownInfo,

                        patchFaces,
                        indexInFace,
                        allInfo
                    );
                }
            }

            // Send to neighbour
            if (debug)
            {
                Pout<< " Sending from " << pp.name() << " to "
                    << procPatch.neighbProcNo() << "   point information:"
                    << patchFaces.size() << endl;
            }

            UOPstream toNeighbour(procPatch.neighbProcNo(), pBufs);
            toNeighbour << patchFaces << indexInFace << allInfo;
        }
    }
}


void Foam::globalPoints::receivePatchPoints
(
    const bool mergeSeparated,
    const Map<label>& meshToPatchPoint,
    const labelList& patchToMeshPoint,
    PstreamBuffers& pBufs,
    labelHashSet& changedPoints
)
{
    // Receive all my neighbours' information and merge with mine.
    // After finishing will have updated
    // - procPoints_ : all neighbour information merged in.
    // - meshToProcPoint_
    // - changedPoints: all points for which something changed.

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const labelPairList& patchInfo = globalTransforms_.patchTransformSign();

    // Reset changed points
    changedPoints.clear();

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if
        (
            (Pstream::parRun() && isA<processorPolyPatch>(pp))
         && (mergeSeparated || patchInfo[patchi].first() == -1)
        )
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pp);

            labelList patchFaces;
            labelList indexInFace;
            List<labelPairList> nbrInfo;

            {
                UIPstream fromNeighbour(procPatch.neighbProcNo(), pBufs);
                fromNeighbour >> patchFaces >> indexInFace >> nbrInfo;
            }

            if (debug)
            {
                Pout<< " On " << pp.name()
                    << " Received from "
                    << procPatch.neighbProcNo() << "   point information:"
                    << patchFaces.size() << endl;
            }

            forAll(patchFaces, i)
            {
                const face& f = pp[patchFaces[i]];

                // Get index in this face from index on face on other side.
                label index = (f.size() - indexInFace[i]) % f.size();

                // Get the meshpoint on my side
                label meshPointi = f[index];

                label localPointi = meshToLocalPoint
                (
                    meshToPatchPoint,
                    meshPointi
                );

                if (mergeInfo(nbrInfo[i], localPointi))
                {
                    changedPoints.insert(localPointi);
                }
            }
        }
        else if
        (
            (
                isA<cyclicPolyPatch>(pp)
             && refCast<const cyclicPolyPatch>(pp).owner()
            )
         && (mergeSeparated || patchInfo[patchi].first() == -1)
        )
        {
            // Handle cyclics: send lower half to upper half and vice versa.
            // Or since they both are in memory just do it point by point.

            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(pp);

            // Pout<< "Patch:" << patchi << " name:" << pp.name() << endl;

            const labelList& meshPoints = pp.meshPoints();
            const labelList coupledMeshPoints(reverseMeshPoints(cycPatch));

            forAll(meshPoints, i)
            {
                label meshPointA = meshPoints[i];
                label meshPointB = coupledMeshPoints[i];

                if (meshPointA != meshPointB)
                {
                    // Pout<< "Connection between point " << meshPointA
                    //    << " at " << mesh_.points()[meshPointA]
                    //    << " and " << meshPointB
                    //    << " at " << mesh_.points()[meshPointB] << endl;

                    label localA = meshToLocalPoint
                    (
                        meshToPatchPoint,
                        meshPointA
                    );
                    label localB = meshToLocalPoint
                    (
                        meshToPatchPoint,
                        meshPointB
                    );


                    // Do we have information on pointA?
                    Map<label>::iterator procPointA =
                        meshToProcPoint_.find(localA);

                    if (procPointA != meshToProcPoint_.end())
                    {
                        const labelPairList infoA = addSendTransform
                        (
                            cycPatch.index(),
                            procPoints_[procPointA()]
                        );

                        if (mergeInfo(infoA, localB))
                        {
                            changedPoints.insert(localB);
                        }
                    }

                    // Same for info on pointB
                    Map<label>::iterator procPointB =
                        meshToProcPoint_.find(localB);

                    if (procPointB != meshToProcPoint_.end())
                    {
                        const labelPairList infoB = addSendTransform
                        (
                            cycPatch.nbrPatchID(),
                            procPoints_[procPointB()]
                        );

                        if (mergeInfo(infoB, localA))
                        {
                            changedPoints.insert(localA);
                        }
                    }
                }
            }
        }
    }
}


void Foam::globalPoints::remove
(
    const labelList& patchToMeshPoint,
    const Map<label>& directNeighbours
)
{
    // Remove entries which are handled by normal face-face communication. I.e.
    // those points where the equivalence list is only me and my (face)neighbour

    // Save old ones.
    Map<label> oldMeshToProcPoint(move(meshToProcPoint_));
    meshToProcPoint_.resize(oldMeshToProcPoint.size());
    DynamicList<labelPairList> oldProcPoints(move(procPoints_));
    procPoints_.setCapacity(oldProcPoints.size());

    // Go through all equivalences
    forAllConstIter(Map<label>, oldMeshToProcPoint, iter)
    {
        label localPointi = iter.key();
        const labelPairList& pointInfo = oldProcPoints[iter()];

        if (pointInfo.size() == 2)
        {
            // I will be in this equivalence list.
            // Check whether my direct (=face) neighbour
            // is in it. This would be an ordinary connection and can be
            // handled by normal face-face connectivity.

            label proc0 = globalTransforms_.processor(pointInfo[0]);
            label proc1 = globalTransforms_.processor(pointInfo[1]);

            if
            (
                (
                    proc0 == Pstream::myProcNo()
                 && directNeighbours.found
                    (
                        globalTransforms_.index(pointInfo[0])
                    )
                )
             || (
                    proc1 == Pstream::myProcNo()
                 && directNeighbours.found
                    (
                        globalTransforms_.index(pointInfo[1])
                    )
                )
            )
            {
                // Normal faceNeighbours
                if (proc0 == Pstream::myProcNo())
                {
                    // Pout<< "Removing direct neighbour:"
                    //    << mesh_.points()
                    //       [globalTransforms_.index(pointInfo[0])]
                    //    << endl;
                }
                else if (proc1 == Pstream::myProcNo())
                {
                    // Pout<< "Removing direct neighbour:"
                    //    << mesh_.points()
                    //       [globalTransforms_.index(pointInfo[1])]
                    //    << endl;
                }
            }
            else
            {
                // This condition will be very rare: points are used by
                // two processors which are not face-face connected.
                // e.g.
                // +------+------+
                // | wall |  B   |
                // +------+------+
                // |   A  | wall |
                // +------+------+
                // Processor A and B share a point. Note that this only will
                // be found if the two domains are face connected at all
                // (not shown in the picture)

                meshToProcPoint_.insert(localPointi, procPoints_.size());
                procPoints_.append(pointInfo);
            }
        }
        else if (pointInfo.size() == 1)
        {
            // This happens for 'wedge' like cyclics where the two halves
            // come together in the same point so share the same meshPoint.
            // So this meshPoint will have info of size one only.
            if
            (
                globalTransforms_.processor(pointInfo[0])
             != Pstream::myProcNo()
             || !directNeighbours.found
                (
                    globalTransforms_.index(pointInfo[0])
                )
            )
            {
                meshToProcPoint_.insert(localPointi, procPoints_.size());
                procPoints_.append(pointInfo);
            }
        }
        else
        {
            meshToProcPoint_.insert(localPointi, procPoints_.size());
            procPoints_.append(pointInfo);
        }
    }

    procPoints_.shrink();
    meshToProcPoint_.resize(2*procPoints_.size());
}


Foam::labelList Foam::globalPoints::reverseMeshPoints
(
    const cyclicPolyPatch& pp
)
{
    const cyclicPolyPatch& nbrPatch = pp.nbrPatch();

    faceList masterFaces(nbrPatch.size());

    forAll(nbrPatch, facei)
    {
        masterFaces[facei] = nbrPatch[facei].reverseFace();
    }

    return primitiveFacePatch
    (
        masterFaces,
        nbrPatch.points()
    ).meshPoints();
}


void Foam::globalPoints::calculateSharedPoints
(
    const Map<label>& meshToPatchPoint, // from mesh point to local numbering
    const labelList& patchToMeshPoint,  // from local numbering to mesh point
    const bool keepAllPoints,
    const bool mergeSeparated
)
{
    if (debug)
    {
        Pout<< "globalPoints::calculateSharedPoints(..) : "
            << "doing processor to processor communication to get sharedPoints"
            << endl
            << "    keepAllPoints :" << keepAllPoints << endl
            << "    mergeSeparated:" << mergeSeparated << endl
            << endl;
    }


    labelHashSet changedPoints(2*nPatchPoints_);

    // Initialize procPoints with my patch points. Keep track of points
    // inserted (in changedPoints)
    // There are two possible forms of this:
    // - initialize with all patch points (allPoints = true). This causes all
    //   patch points to be exchanged so a lot of information gets stored and
    //   transferred. This all gets filtered out later when removing the
    //   equivalence lists of size 2.
    // - initialize with boundary points of patches only (allPoints = false).
    //   This should work for all decompositions except extreme ones where a
    //   shared point is not on the boundary of any processor patches using it.
    //   This would happen if a domain was pinched such that two patches share
    //   a point or edge.
    initOwnPoints(meshToPatchPoint, true, changedPoints);

    // Do one exchange iteration to get neighbour points.
    {
        // Note: to use 'scheduled' would have to intersperse send and receive.
        // So for now just use nonBlocking. Also globalPoints itself gets
        // constructed by mesh.globalData().patchSchedule() so creates a loop.
        PstreamBuffers pBufs
        (
            (
                Pstream::defaultCommsType == Pstream::commsTypes::scheduled
              ? Pstream::commsTypes::nonBlocking
              : Pstream::defaultCommsType
            )
        );
        sendPatchPoints
        (
            mergeSeparated,
            meshToPatchPoint,
            pBufs,
            changedPoints
        );
        pBufs.finishedSends();
        receivePatchPoints
        (
            mergeSeparated,
            meshToPatchPoint,
            patchToMeshPoint,
            pBufs,
            changedPoints
        );
    }

    // Save neighbours reachable through face-face communication.
    Map<label> neighbourList;
    if (!keepAllPoints)
    {
        neighbourList = meshToProcPoint_;
    }

    // Exchange until nothing changes on all processors.
    bool changed = false;

    do
    {
        PstreamBuffers pBufs
        (
            (
                Pstream::defaultCommsType == Pstream::commsTypes::scheduled
              ? Pstream::commsTypes::nonBlocking
              : Pstream::defaultCommsType
            )
        );
        sendPatchPoints
        (
            mergeSeparated,
            meshToPatchPoint,
            pBufs,
            changedPoints
        );
        pBufs.finishedSends();
        receivePatchPoints
        (
            mergeSeparated,
            meshToPatchPoint,
            patchToMeshPoint,
            pBufs,
            changedPoints
        );

        changed = changedPoints.size() > 0;
        reduce(changed, orOp<bool>());

    } while (changed);


    // Pout<< "**ALL** connected points:" << endl;
    // forAllConstIter(Map<label>, meshToProcPoint_, iter)
    //{
    //    label localI = iter.key();
    //    const labelPairList& pointInfo = procPoints_[iter()];
    //    Pout<< "pointi:" << localI << " index:" << iter()
    //        << " coord:"
    //        << mesh_.points()[localToMeshPoint(patchToMeshPoint, localI)]
    //        << endl;
    //    printProcPoints(patchToMeshPoint, pointInfo);
    //    Pout<< endl;
    //}


    // Remove direct neighbours from point equivalences.
    if (!keepAllPoints)
    {
        remove(patchToMeshPoint, neighbourList);
    }


    // Sort procPoints in incremental order. This will make
    // the master the first element on all processors.
    // Note: why not sort in decreasing order? Give more work to higher
    //       processors.
    forAllConstIter(Map<label>, meshToProcPoint_, iter)
    {
        labelPairList& pointInfo = procPoints_[iter()];
        sort(pointInfo, globalIndexAndTransform::less(globalTransforms_));
    }


    // We now have - in procPoints_ - a list of points which are shared between
    // multiple processors. Filter into non-transformed and transformed
    // connections.

    pointPoints_.setSize(globalIndices_.localSize());
    List<labelPairList> transformedPoints(globalIndices_.localSize());
    forAllConstIter(Map<label>, meshToProcPoint_, iter)
    {
        const labelPairList& pointInfo = procPoints_[iter()];

        if (pointInfo.size() >= 2)
        {
            // Since sorted master point is the first element
            const labelPair& masterInfo = pointInfo[0];

            if
            (
                (
                    globalTransforms_.processor(masterInfo)
                 == Pstream::myProcNo()
                )
             && (globalTransforms_.index(masterInfo) == iter.key())
            )
            {
                labelList& pPoints = pointPoints_[iter.key()];
                pPoints.setSize(pointInfo.size()-1);

                labelPairList& trafoPPoints = transformedPoints[iter.key()];
                trafoPPoints.setSize(pointInfo.size()-1);

                label nonTransformI = 0;
                label transformI = 0;

                for (label i = 1; i < pointInfo.size(); i++)
                {
                    const labelPair& info = pointInfo[i];
                    label proci = globalTransforms_.processor(info);
                    label index = globalTransforms_.index(info);
                    label transform = globalTransforms_.transformIndex
                    (
                        info
                    );

                    if (transform == globalTransforms_.nullTransformIndex())
                    {
                        pPoints[nonTransformI++] = globalIndices_.toGlobal
                        (
                            proci,
                            index
                        );
                    }
                    else
                    {
                        trafoPPoints[transformI++] = info;
                    }
                }

                pPoints.setSize(nonTransformI);
                trafoPPoints.setSize(transformI);
            }
        }
    }


    List<Map<label>> compactMap;
    map_.reset
    (
        new mapDistribute
        (
            globalIndices_,
            pointPoints_,

            globalTransforms_,
            transformedPoints,
            transformedPointPoints_,

            compactMap
        )
    );

    if (debug)
    {
        Pout<< "globalPoints::calculateSharedPoints(..) : "
            << "Finished global points" << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::globalPoints::globalPoints
(
    const polyMesh& mesh,
    const bool keepAllPoints,
    const bool mergeSeparated
)
:
    mesh_(mesh),
    globalIndices_(mesh_.nPoints()),
    globalTransforms_(mesh),
    nPatchPoints_(countPatchPoints(mesh.boundaryMesh())),
    procPoints_(nPatchPoints_),
    meshToProcPoint_(nPatchPoints_)
{
    // Empty patch maps to signal storing mesh point labels
    Map<label> meshToPatchPoint(0);
    labelList patchToMeshPoint(0);

    calculateSharedPoints
    (
        meshToPatchPoint,
        patchToMeshPoint,
        keepAllPoints,
        mergeSeparated
    );
}


Foam::globalPoints::globalPoints
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& coupledPatch,
    const bool keepAllPoints,
    const bool mergeSeparated
)
:
    mesh_(mesh),
    globalIndices_(coupledPatch.nPoints()),
    globalTransforms_(mesh),
    nPatchPoints_(coupledPatch.nPoints()),
    procPoints_(nPatchPoints_),
    meshToProcPoint_(nPatchPoints_)
{
    calculateSharedPoints
    (
        coupledPatch.meshPointMap(),
        coupledPatch.meshPoints(),
        keepAllPoints,
        mergeSeparated
    );
}


// ************************************************************************* //
