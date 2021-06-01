/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "syncTools.H"
#include "polyMesh.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "globalMeshData.H"
#include "contiguous.H"
#include "transform.H"
#include "SubField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T, class CombineOp>
void Foam::syncTools::combine
(
    Map<T>& pointValues,
    const CombineOp& cop,
    const label index,
    const T& val
)
{
    typename Map<T>::iterator iter = pointValues.find(index);

    if (iter != pointValues.end())
    {
        cop(iter(), val);
    }
    else
    {
        pointValues.insert(index, val);
    }
}


template<class T, class CombineOp>
void Foam::syncTools::combine
(
    EdgeMap<T>& edgeValues,
    const CombineOp& cop,
    const edge& index,
    const T& val
)
{
    typename EdgeMap<T>::iterator iter = edgeValues.find(index);

    if (iter != edgeValues.end())
    {
        cop(iter(), val);
    }
    else
    {
        edgeValues.insert(index, val);
    }
}


template<class T, class CombineOp, class TransformOp>
void Foam::syncTools::syncPointMap
(
    const polyMesh& mesh,
    Map<T>& pointValues,        // from mesh point label to value
    const CombineOp& cop,
    const TransformOp& top
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Synchronise multiple shared points.
    const globalMeshData& pd = mesh.globalData();

    // Values on shared points. Keyed on global shared index.
    Map<T> sharedPointValues(0);

    if (pd.nGlobalPoints() > 0)
    {
        // meshPoint per local index
        const labelList& sharedPtLabels = pd.sharedPointLabels();
        // global shared index per local index
        const labelList& sharedPtAddr = pd.sharedPointAddr();

        sharedPointValues.resize(sharedPtAddr.size());

        // Fill my entries in the shared points
        forAll(sharedPtLabels, i)
        {
            label meshPointi = sharedPtLabels[i];

            typename Map<T>::const_iterator fnd =
                pointValues.find(meshPointi);

            if (fnd != pointValues.end())
            {
                combine
                (
                    sharedPointValues,
                    cop,
                    sharedPtAddr[i],    // index
                    fnd()               // value
                );
            }
        }
    }


    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send

        forAll(patches, patchi)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchi])
             && patches[patchi].nPoints() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchi]);

                // Get data per patchPoint in neighbouring point numbers.

                const labelList& meshPts = procPatch.meshPoints();
                const labelList& nbrPts = procPatch.nbrPoints();

                // Extract local values. Create map from nbrPoint to value.
                // Note: how small initial size?
                Map<T> patchInfo(meshPts.size() / 20);

                forAll(meshPts, i)
                {
                    typename Map<T>::const_iterator iter =
                        pointValues.find(meshPts[i]);

                    if (iter != pointValues.end())
                    {
                        patchInfo.insert(nbrPts[i], iter());
                    }
                }

                UOPstream toNeighb(procPatch.neighbProcNo(), pBufs);
                toNeighb << patchInfo;
            }
        }

        pBufs.finishedSends();

        // Receive and combine.

        forAll(patches, patchi)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchi])
             && patches[patchi].nPoints() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchi]);

                UIPstream fromNb(procPatch.neighbProcNo(), pBufs);
                Map<T> nbrPatchInfo(fromNb);

                // Transform
                top(procPatch, nbrPatchInfo);

                const labelList& meshPts = procPatch.meshPoints();

                // Only update those values which come from neighbour
                forAllConstIter(typename Map<T>, nbrPatchInfo, nbrIter)
                {
                    combine
                    (
                        pointValues,
                        cop,
                        meshPts[nbrIter.key()],
                        nbrIter()
                    );
                }
            }
        }
    }

    // Do the cyclics.
    forAll(patches, patchi)
    {
        if (isA<cyclicPolyPatch>(patches[patchi]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchi]);

            if (cycPatch.owner())
            {
                // Owner does all.

                const cyclicPolyPatch& nbrPatch = cycPatch.nbrPatch();
                const edgeList& coupledPoints = cycPatch.coupledPoints();
                const labelList& meshPtsA = cycPatch.meshPoints();
                const labelList& meshPtsB = nbrPatch.meshPoints();

                // Extract local values. Create map from coupled-edge to value.
                Map<T> half0Values(meshPtsA.size() / 20);
                Map<T> half1Values(half0Values.size());

                forAll(coupledPoints, i)
                {
                    const edge& e = coupledPoints[i];

                    typename Map<T>::const_iterator point0Fnd =
                        pointValues.find(meshPtsA[e[0]]);

                    if (point0Fnd != pointValues.end())
                    {
                        half0Values.insert(i, point0Fnd());
                    }

                    typename Map<T>::const_iterator point1Fnd =
                        pointValues.find(meshPtsB[e[1]]);

                    if (point1Fnd != pointValues.end())
                    {
                        half1Values.insert(i, point1Fnd());
                    }
                }

                // Transform to receiving side
                top(cycPatch, half1Values);
                top(nbrPatch, half0Values);

                forAll(coupledPoints, i)
                {
                    const edge& e = coupledPoints[i];

                    typename Map<T>::const_iterator half0Fnd =
                        half0Values.find(i);

                    if (half0Fnd != half0Values.end())
                    {
                        combine
                        (
                            pointValues,
                            cop,
                            meshPtsB[e[1]],
                            half0Fnd()
                        );
                    }

                    typename Map<T>::const_iterator half1Fnd =
                        half1Values.find(i);

                    if (half1Fnd != half1Values.end())
                    {
                        combine
                        (
                            pointValues,
                            cop,
                            meshPtsA[e[0]],
                            half1Fnd()
                        );
                    }
                }
            }
        }
    }

    // Synchronise multiple shared points.
    if (pd.nGlobalPoints() > 0)
    {
        // meshPoint per local index
        const labelList& sharedPtLabels = pd.sharedPointLabels();
        // global shared index per local index
        const labelList& sharedPtAddr = pd.sharedPointAddr();

        // Reduce on master.

        if (Pstream::parRun())
        {
            if (Pstream::master())
            {
                // Receive the edges using shared points from the slave.
                for
                (
                    int slave=Pstream::firstSlave();
                    slave<=Pstream::lastSlave();
                    slave++
                )
                {
                    IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
                    Map<T> nbrValues(fromSlave);

                    // Merge neighbouring values with my values
                    forAllConstIter(typename Map<T>, nbrValues, iter)
                    {
                        combine
                        (
                            sharedPointValues,
                            cop,
                            iter.key(), // edge
                            iter()      // value
                        );
                    }
                }

                // Send back
                for
                (
                    int slave=Pstream::firstSlave();
                    slave<=Pstream::lastSlave();
                    slave++
                )
                {
                    OPstream toSlave(Pstream::commsTypes::scheduled, slave);
                    toSlave << sharedPointValues;
                }
            }
            else
            {
                // Slave: send to master
                {
                    OPstream toMaster
                    (
                        Pstream::commsTypes::scheduled,
                        Pstream::masterNo()
                    );
                    toMaster << sharedPointValues;
                }
                // Receive merged values
                {
                    IPstream fromMaster
                    (
                        Pstream::commsTypes::scheduled,
                        Pstream::masterNo()
                    );
                    fromMaster >> sharedPointValues;
                }
            }
        }


        // Merge sharedPointValues (keyed on sharedPointAddr) into
        // pointValues (keyed on mesh points).

        // Map from global shared index to meshpoint
        Map<label> sharedToMeshPoint(2*sharedPtAddr.size());
        forAll(sharedPtAddr, i)
        {
            sharedToMeshPoint.insert(sharedPtAddr[i], sharedPtLabels[i]);
        }

        forAllConstIter(Map<label>, sharedToMeshPoint, iter)
        {
            // Do I have a value for my shared point
            typename Map<T>::const_iterator sharedFnd =
                sharedPointValues.find(iter.key());

            if (sharedFnd != sharedPointValues.end())
            {
                pointValues.set(iter(), sharedFnd());
            }
        }
    }
}


template<class T, class CombineOp, class TransformOp>
void Foam::syncTools::syncEdgeMap
(
    const polyMesh& mesh,
    EdgeMap<T>& edgeValues,
    const CombineOp& cop,
    const TransformOp& top
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();


    // Do synchronisation without constructing globalEdge addressing
    // (since this constructs mesh edge addressing)


    // Swap proc patch info
    // ~~~~~~~~~~~~~~~~~~~~

    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send

        forAll(patches, patchi)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchi])
             && patches[patchi].nEdges() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchi]);


                // Get data per patch edge in neighbouring edge.

                const edgeList& edges = procPatch.edges();
                const labelList& meshPts = procPatch.meshPoints();
                const labelList& nbrPts = procPatch.nbrPoints();

                EdgeMap<T> patchInfo(edges.size() / 20);

                forAll(edges, i)
                {
                    const edge& e = edges[i];
                    const edge meshEdge(meshPts[e[0]], meshPts[e[1]]);

                    typename EdgeMap<T>::const_iterator iter =
                        edgeValues.find(meshEdge);

                    if (iter != edgeValues.end())
                    {
                        const edge nbrEdge(nbrPts[e[0]], nbrPts[e[1]]);
                        patchInfo.insert(nbrEdge, iter());
                    }
                }

                UOPstream toNeighb(procPatch.neighbProcNo(), pBufs);
                toNeighb << patchInfo;
            }
        }

        pBufs.finishedSends();

        // Receive and combine.

        forAll(patches, patchi)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchi])
             && patches[patchi].nEdges() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchi]);

                EdgeMap<T> nbrPatchInfo;
                {
                    UIPstream fromNbr(procPatch.neighbProcNo(), pBufs);
                    fromNbr >> nbrPatchInfo;
                }

                // Apply transform to convert to this side properties.
                top(procPatch, nbrPatchInfo);


                // Only update those values which come from neighbour
                const labelList& meshPts = procPatch.meshPoints();

                forAllConstIter(typename EdgeMap<T>, nbrPatchInfo, nbrIter)
                {
                    const edge& e = nbrIter.key();
                    const edge meshEdge(meshPts[e[0]], meshPts[e[1]]);

                    combine
                    (
                        edgeValues,
                        cop,
                        meshEdge,   // edge
                        nbrIter()   // value
                    );
                }
            }
        }
    }


    // Swap cyclic info
    // ~~~~~~~~~~~~~~~~

    forAll(patches, patchi)
    {
        if (isA<cyclicPolyPatch>(patches[patchi]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchi]);

            if (cycPatch.owner())
            {
                // Owner does all.

                const edgeList& coupledEdges = cycPatch.coupledEdges();
                const labelList& meshPtsA = cycPatch.meshPoints();
                const edgeList& edgesA = cycPatch.edges();
                const cyclicPolyPatch& nbrPatch = cycPatch.nbrPatch();
                const labelList& meshPtsB = nbrPatch.meshPoints();
                const edgeList& edgesB = nbrPatch.edges();

                // Extract local values. Create map from edge to value.
                Map<T> half0Values(edgesA.size() / 20);
                Map<T> half1Values(half0Values.size());

                forAll(coupledEdges, i)
                {
                    const edge& twoEdges = coupledEdges[i];

                    {
                        const edge& e0 = edgesA[twoEdges[0]];
                        const edge meshEdge0(meshPtsA[e0[0]], meshPtsA[e0[1]]);

                        typename EdgeMap<T>::const_iterator iter =
                            edgeValues.find(meshEdge0);

                        if (iter != edgeValues.end())
                        {
                            half0Values.insert(i, iter());
                        }
                    }
                    {
                        const edge& e1 = edgesB[twoEdges[1]];
                        const edge meshEdge1(meshPtsB[e1[0]], meshPtsB[e1[1]]);

                        typename EdgeMap<T>::const_iterator iter =
                            edgeValues.find(meshEdge1);

                        if (iter != edgeValues.end())
                        {
                            half1Values.insert(i, iter());
                        }
                    }
                }

                // Transform to this side
                top(cycPatch, half1Values);
                top(nbrPatch, half0Values);


                // Extract and combine information

                forAll(coupledEdges, i)
                {
                    const edge& twoEdges = coupledEdges[i];

                    typename Map<T>::const_iterator half1Fnd =
                        half1Values.find(i);

                    if (half1Fnd != half1Values.end())
                    {
                        const edge& e0 = edgesA[twoEdges[0]];
                        const edge meshEdge0(meshPtsA[e0[0]], meshPtsA[e0[1]]);

                        combine
                        (
                            edgeValues,
                            cop,
                            meshEdge0,  // edge
                            half1Fnd()  // value
                        );
                    }

                    typename Map<T>::const_iterator half0Fnd =
                        half0Values.find(i);
                    if (half0Fnd != half0Values.end())
                    {
                        const edge& e1 = edgesB[twoEdges[1]];
                        const edge meshEdge1(meshPtsB[e1[0]], meshPtsB[e1[1]]);

                        combine
                        (
                            edgeValues,
                            cop,
                            meshEdge1,  // edge
                            half0Fnd()  // value
                        );
                    }
                }
            }
        }
    }

    // Synchronise multiple shared points.
    // Problem is that we don't want to construct shared edges so basically
    // we do here like globalMeshData but then using sparse edge representation
    // (EdgeMap instead of mesh.edges())

    const globalMeshData& pd = mesh.globalData();
    const labelList& sharedPtAddr = pd.sharedPointAddr();
    const labelList& sharedPtLabels = pd.sharedPointLabels();

    // 1. Create map from meshPoint to globalShared index.
    Map<label> meshToShared(2*sharedPtLabels.size());
    forAll(sharedPtLabels, i)
    {
        meshToShared.insert(sharedPtLabels[i], sharedPtAddr[i]);
    }

    // Values on shared points. From two sharedPtAddr indices to a value.
    EdgeMap<T> sharedEdgeValues(meshToShared.size());

    // From shared edge to mesh edge. Used for merging later on.
    EdgeMap<edge> potentialSharedEdge(meshToShared.size());

    // 2. Find any edges using two global shared points. These will always be
    // on the outside of the mesh. (though might not be on coupled patch
    // if is single edge and not on coupled face)
    // Store value (if any) on sharedEdgeValues
    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        const face& f = mesh.faces()[facei];

        forAll(f, fp)
        {
            label v0 = f[fp];
            label v1 = f[f.fcIndex(fp)];

            Map<label>::const_iterator v0Fnd = meshToShared.find(v0);

            if (v0Fnd != meshToShared.end())
            {
                Map<label>::const_iterator v1Fnd = meshToShared.find(v1);

                if (v1Fnd != meshToShared.end())
                {
                    const edge meshEdge(v0, v1);

                    // edge in shared point labels
                    const edge sharedEdge(v0Fnd(), v1Fnd());

                    // Store mesh edge as a potential shared edge.
                    potentialSharedEdge.insert(sharedEdge, meshEdge);

                    typename EdgeMap<T>::const_iterator edgeFnd =
                        edgeValues.find(meshEdge);

                    if (edgeFnd != edgeValues.end())
                    {
                        // edge exists in edgeValues. See if already in map
                        // (since on same processor, e.g. cyclic)
                        combine
                        (
                            sharedEdgeValues,
                            cop,
                            sharedEdge, // edge
                            edgeFnd()   // value
                        );
                    }
                }
            }
        }
    }


    // Now sharedEdgeValues will contain per potential sharedEdge the value.
    // (potential since an edge having two shared points is not necessary a
    //  shared edge).
    // Reduce this on the master.

    if (Pstream::parRun())
    {
        if (Pstream::master())
        {
            // Receive the edges using shared points from the slave.
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
                EdgeMap<T> nbrValues(fromSlave);

                // Merge neighbouring values with my values
                forAllConstIter(typename EdgeMap<T>, nbrValues, iter)
                {
                    combine
                    (
                        sharedEdgeValues,
                        cop,
                        iter.key(), // edge
                        iter()      // value
                    );
                }
            }

            // Send back
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {

                OPstream toSlave(Pstream::commsTypes::scheduled, slave);
                toSlave << sharedEdgeValues;
            }
        }
        else
        {
            // Send to master
            {
                OPstream toMaster
                (
                    Pstream::commsTypes::scheduled,
                    Pstream::masterNo()
                );
                toMaster << sharedEdgeValues;
            }
            // Receive merged values
            {
                IPstream fromMaster
                (
                    Pstream::commsTypes::scheduled,
                    Pstream::masterNo()
                );
                fromMaster >> sharedEdgeValues;
            }
        }
    }


    // Merge sharedEdgeValues (keyed on sharedPointAddr) into edgeValues
    // (keyed on mesh points).

    // Loop over all my shared edges.
    forAllConstIter(typename EdgeMap<edge>, potentialSharedEdge, iter)
    {
        const edge& sharedEdge = iter.key();
        const edge& meshEdge = iter();

        // Do I have a value for the shared edge?
        typename EdgeMap<T>::const_iterator sharedFnd =
            sharedEdgeValues.find(sharedEdge);

        if (sharedFnd != sharedEdgeValues.end())
        {
            combine
            (
                edgeValues,
                cop,
                meshEdge,       // edge
                sharedFnd()     // value
            );
        }
    }
}


//template<class T, class CombineOp, class TransformOp>
//void Foam::syncTools::syncPointList
//(
//    const polyMesh& mesh,
//    List<T>& pointValues,
//    const CombineOp& cop,
//    const T& nullValue,
//    const TransformOp& top
//)
//{
//    if (pointValues.size() != mesh.nPoints())
//    {
//        FatalErrorInFunction
//            << "Number of values " << pointValues.size()
//            << " is not equal to the number of points in the mesh "
//            << mesh.nPoints() << abort(FatalError);
//    }
//
//    const polyBoundaryMesh& patches = mesh.boundaryMesh();
//
//    // Synchronise multiple shared points.
//    const globalMeshData& pd = mesh.globalData();
//
//    // Values on shared points.
//    Field<T> sharedPts(0);
//    if (pd.nGlobalPoints() > 0)
//    {
//        // Values on shared points.
//        sharedPts.setSize(pd.nGlobalPoints(), nullValue);
//
//        forAll(pd.sharedPointLabels(), i)
//        {
//            label meshPointi = pd.sharedPointLabels()[i];
//            // Fill my entries in the shared points
//            sharedPts[pd.sharedPointAddr()[i]] = pointValues[meshPointi];
//        }
//    }
//
//    if (Pstream::parRun())
//    {
//        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
//
//        // Send
//
//        forAll(patches, patchi)
//        {
//            if
//            (
//                isA<processorPolyPatch>(patches[patchi])
//             && patches[patchi].nPoints() > 0
//            )
//            {
//                const processorPolyPatch& procPatch =
//                    refCast<const processorPolyPatch>(patches[patchi]);
//
//                // Get data per patchPoint in neighbouring point numbers.
//                Field<T> patchInfo(procPatch.nPoints());
//
//                const labelList& meshPts = procPatch.meshPoints();
//                const labelList& nbrPts = procPatch.nbrPoints();
//
//                forAll(nbrPts, pointi)
//                {
//                    label nbrPointi = nbrPts[pointi];
//                    patchInfo[nbrPointi] = pointValues[meshPts[pointi]];
//                }
//
//                UOPstream toNbr(procPatch.neighbProcNo(), pBufs);
//                toNbr << patchInfo;
//            }
//        }
//
//        pBufs.finishedSends();
//
//        // Receive and combine.
//
//        forAll(patches, patchi)
//        {
//            if
//            (
//                isA<processorPolyPatch>(patches[patchi])
//             && patches[patchi].nPoints() > 0
//            )
//            {
//                const processorPolyPatch& procPatch =
//                    refCast<const processorPolyPatch>(patches[patchi]);
//
//                Field<T> nbrPatchInfo(procPatch.nPoints());
//                {
//                    UIPstream fromNbr(procPatch.neighbProcNo(), pBufs);
//                    fromNbr >> nbrPatchInfo;
//                }
//
//                // Transform to this side
//                top(procPatch, nbrPatchInfo);
//
//                const labelList& meshPts = procPatch.meshPoints();
//
//                forAll(meshPts, pointi)
//                {
//                    label meshPointi = meshPts[pointi];
//                    cop(pointValues[meshPointi], nbrPatchInfo[pointi]);
//                }
//            }
//        }
//    }
//
//    // Do the cyclics.
//    forAll(patches, patchi)
//    {
//        if (isA<cyclicPolyPatch>(patches[patchi]))
//        {
//            const cyclicPolyPatch& cycPatch =
//                refCast<const cyclicPolyPatch>(patches[patchi]);
//
//            if (cycPatch.owner())
//            {
//                // Owner does all.
//
//                const edgeList& coupledPoints = cycPatch.coupledPoints();
//                const labelList& meshPts = cycPatch.meshPoints();
//                const cyclicPolyPatch& nbrPatch = cycPatch.nbrPatch();
//                const labelList& nbrMeshPoints = nbrPatch.meshPoints();
//
//                Field<T> half0Values(coupledPoints.size());
//                Field<T> half1Values(coupledPoints.size());
//
//                forAll(coupledPoints, i)
//                {
//                    const edge& e = coupledPoints[i];
//                    half0Values[i] = pointValues[meshPts[e[0]]];
//                    half1Values[i] = pointValues[nbrMeshPoints[e[1]]];
//                }
//
//                // SubField<T> slice0(half0Values, half0Values.size());
//                // SubField<T> slice1(half1Values, half1Values.size());
//                // top(cycPatch, reinterpret_cast<Field<T>&>(slice1));
//                // top(nbrPatch, reinterpret_cast<Field<T>&>(slice0));
//
//                top(cycPatch, half1Values);
//                top(nbrPatch, half0Values);
//
//                forAll(coupledPoints, i)
//                {
//                    const edge& e = coupledPoints[i];
//                    cop(pointValues[meshPts[e[0]]], half1Values[i]);
//                    cop(pointValues[nbrMeshPoints[e[1]]], half0Values[i]);
//                }
//            }
//        }
//    }
//
//    // Synchronise multiple shared points.
//    const globalMeshData& pd = mesh.globalData();
//
//    if (pd.nGlobalPoints() > 0)
//    {
//        // Combine on master.
//        Pstream::listCombineGather(sharedPts, cop);
//        Pstream::listCombineScatter(sharedPts);
//
//        // Now we will all have the same information. Merge it back with
//        // my local information.
//        forAll(pd.sharedPointLabels(), i)
//        {
//            label meshPointi = pd.sharedPointLabels()[i];
//            pointValues[meshPointi] = sharedPts[pd.sharedPointAddr()[i]];
//        }
//    }
//}


//template<class T, class CombineOp, class TransformOp>
//void Foam::syncTools::syncPointList
//(
//    const polyMesh& mesh,
//    const labelList& meshPoints,
//    List<T>& pointValues,
//    const CombineOp& cop,
//    const T& nullValue,
//    const TransformOp& top
//)
//{
//    if (pointValues.size() != meshPoints.size())
//    {
//        FatalErrorInFunction
//            << "Number of values " << pointValues.size()
//            << " is not equal to the number of points "
//            << meshPoints.size() << abort(FatalError);
//    }
//
//    if (!hasCouples(mesh.boundaryMesh()))
//    {
//        return;
//    }
//
//    Field<T> meshValues(mesh.nPoints(), nullValue);
//
//    forAll(meshPoints, i)
//    {
//        meshValues[meshPoints[i]] = pointValues[i];
//    }
//
//    syncTools::syncPointList
//    (
//        mesh,
//        meshValues,
//        cop,            // combine op
//        nullValue,      // null value
//        top             // position or field
//    );
//
//    forAll(meshPoints, i)
//    {
//        pointValues[i] = meshValues[meshPoints[i]];
//    }
//}

template<class T, class CombineOp, class TransformOp>
void Foam::syncTools::syncPointList
(
    const polyMesh& mesh,
    List<T>& pointValues,
    const CombineOp& cop,
    const T& nullValue,
    const TransformOp& top
)
{
    if (pointValues.size() != mesh.nPoints())
    {
        FatalErrorInFunction
            << "Number of values " << pointValues.size()
            << " is not equal to the number of points in the mesh "
            << mesh.nPoints() << abort(FatalError);
    }

    mesh.globalData().syncPointData(pointValues, cop, top);
}


//template<class CombineOp>
//void Foam::syncTools::syncPointPositions
//(
//    const polyMesh& mesh,
//    List<point>& pointValues,
//    const CombineOp& cop,
//    const point& nullValue
//)
//{
//    if (pointValues.size() != mesh.nPoints())
//    {
//        FatalErrorInFunction
//            << "Number of values " << pointValues.size()
//            << " is not equal to the number of points in the mesh "
//            << mesh.nPoints() << abort(FatalError);
//    }
//
//    mesh.globalData().syncPointData(pointValues, cop, true);
//}


template<class T, class CombineOp, class TransformOp>
void Foam::syncTools::syncPointList
(
    const polyMesh& mesh,
    const labelList& meshPoints,
    List<T>& pointValues,
    const CombineOp& cop,
    const T& nullValue,
    const TransformOp& top
)
{
    if (pointValues.size() != meshPoints.size())
    {
        FatalErrorInFunction
            << "Number of values " << pointValues.size()
            << " is not equal to the number of meshPoints "
            << meshPoints.size() << abort(FatalError);
    }
    const globalMeshData& gd = mesh.globalData();
    const indirectPrimitivePatch& cpp = gd.coupledPatch();
    const Map<label>& mpm = cpp.meshPointMap();

    List<T> cppFld(cpp.nPoints(), nullValue);

    forAll(meshPoints, i)
    {
        label pointi = meshPoints[i];
        Map<label>::const_iterator iter = mpm.find(pointi);
        if (iter != mpm.end())
        {
            cppFld[iter()] = pointValues[i];
        }
    }

    globalMeshData::syncData
    (
        cppFld,
        gd.globalPointSlaves(),
        gd.globalPointTransformedSlaves(),
        gd.globalPointSlavesMap(),
        gd.globalTransforms(),
        cop,
        top
    );

    forAll(meshPoints, i)
    {
        label pointi = meshPoints[i];
        Map<label>::const_iterator iter = mpm.find(pointi);
        if (iter != mpm.end())
        {
            pointValues[i] = cppFld[iter()];
        }
    }
}


//template<class CombineOp>
//void Foam::syncTools::syncPointPositions
//(
//    const polyMesh& mesh,
//    const labelList& meshPoints,
//    List<point>& pointValues,
//    const CombineOp& cop,
//    const point& nullValue
//)
//{
//    if (pointValues.size() != meshPoints.size())
//    {
//        FatalErrorInFunction
//            << "Number of values " << pointValues.size()
//            << " is not equal to the number of meshPoints "
//            << meshPoints.size() << abort(FatalError);
//    }
//    const globalMeshData& gd = mesh.globalData();
//    const indirectPrimitivePatch& cpp = gd.coupledPatch();
//    const Map<label>& mpm = cpp.meshPointMap();
//
//    List<point> cppFld(cpp.nPoints(), nullValue);
//
//    forAll(meshPoints, i)
//    {
//        label pointi = meshPoints[i];
//        Map<label>::const_iterator iter = mpm.find(pointi);
//        if (iter != mpm.end())
//        {
//            cppFld[iter()] = pointValues[i];
//        }
//    }
//
//    globalMeshData::syncData
//    (
//        cppFld,
//        gd.globalPointSlaves(),
//        gd.globalPointTransformedSlaves(),
//        gd.globalPointSlavesMap(),
//        gd.globalTransforms(),
//        cop,
//        true,   // position?
//        mapDistribute::transform()  // not used
//    );
//
//    forAll(meshPoints, i)
//    {
//        label pointi = meshPoints[i];
//        Map<label>::const_iterator iter = mpm.find(pointi);
//        if (iter != mpm.end())
//        {
//            pointValues[i] = cppFld[iter()];
//        }
//    }
//}


template<class T, class CombineOp, class TransformOp>
void Foam::syncTools::syncEdgeList
(
    const polyMesh& mesh,
    List<T>& edgeValues,
    const CombineOp& cop,
    const T& nullValue,
    const TransformOp& top
)
{
    if (edgeValues.size() != mesh.nEdges())
    {
        FatalErrorInFunction
            << "Number of values " << edgeValues.size()
            << " is not equal to the number of edges in the mesh "
            << mesh.nEdges() << abort(FatalError);
    }

    const globalMeshData& gd = mesh.globalData();
    const labelList& meshEdges = gd.coupledPatchMeshEdges();
    const globalIndexAndTransform& git = gd.globalTransforms();
    const mapDistribute& edgeMap = gd.globalEdgeSlavesMap();

    List<T> cppFld(UIndirectList<T>(edgeValues, meshEdges));

    globalMeshData::syncData
    (
        cppFld,
        gd.globalEdgeSlaves(),
        gd.globalEdgeTransformedSlaves(),
        edgeMap,
        git,
        cop,
        top
    );

    // Extract back onto mesh
    forAll(meshEdges, i)
    {
        edgeValues[meshEdges[i]] = cppFld[i];
    }
}


//template<class CombineOp>
//void Foam::syncTools::syncEdgePositions
//(
//    const polyMesh& mesh,
//    List<point>& edgeValues,
//    const CombineOp& cop,
//    const point& nullValue
//)
//{
//    if (edgeValues.size() != mesh.nEdges())
//    {
//        FatalErrorInFunction
//            << "Number of values " << edgeValues.size()
//            << " is not equal to the number of edges in the mesh "
//            << mesh.nEdges() << abort(FatalError);
//    }
//
//    const globalMeshData& gd = mesh.globalData();
//    const labelList& meshEdges = gd.coupledPatchMeshEdges();
//    const globalIndexAndTransform& git = gd.globalTransforms();
//    const mapDistribute& map = gd.globalEdgeSlavesMap();
//
//    List<point> cppFld(UIndirectList<point>(edgeValues, meshEdges));
//
//    globalMeshData::syncData
//    (
//        cppFld,
//        gd.globalEdgeSlaves(),
//        gd.globalEdgeTransformedSlaves(),
//        map,
//        git,
//        cop,
//        true,       // position?
//        mapDistribute::transform()  // not used
//    );
//
//    // Extract back onto mesh
//    forAll(meshEdges, i)
//    {
//        edgeValues[meshEdges[i]] = cppFld[i];
//    }
//}


template<class T, class CombineOp, class TransformOp>
void Foam::syncTools::syncEdgeList
(
    const polyMesh& mesh,
    const labelList& meshEdges,
    List<T>& edgeValues,
    const CombineOp& cop,
    const T& nullValue,
    const TransformOp& top
)
{
    if (edgeValues.size() != meshEdges.size())
    {
        FatalErrorInFunction
            << "Number of values " << edgeValues.size()
            << " is not equal to the number of meshEdges "
            << meshEdges.size() << abort(FatalError);
    }
    const globalMeshData& gd = mesh.globalData();
    const indirectPrimitivePatch& cpp = gd.coupledPatch();
    const Map<label>& mpm = gd.coupledPatchMeshEdgeMap();

    List<T> cppFld(cpp.nEdges(), nullValue);

    forAll(meshEdges, i)
    {
        label edgeI = meshEdges[i];
        Map<label>::const_iterator iter = mpm.find(edgeI);
        if (iter != mpm.end())
        {
            cppFld[iter()] = edgeValues[i];
        }
    }

    globalMeshData::syncData
    (
        cppFld,
        gd.globalEdgeSlaves(),
        gd.globalEdgeTransformedSlaves(),
        gd.globalEdgeSlavesMap(),
        gd.globalTransforms(),
        cop,
        top
    );

    forAll(meshEdges, i)
    {
        label edgeI = meshEdges[i];
        Map<label>::const_iterator iter = mpm.find(edgeI);
        if (iter != mpm.end())
        {
            edgeValues[i] = cppFld[iter()];
        }
    }
}

template<class T, class CombineOp, class TransformOp>
void Foam::syncTools::syncBoundaryFaceList
(
    const polyMesh& mesh,
    UList<T>& faceValues,
    const CombineOp& cop,
    const TransformOp& top,
    const bool parRun
)
{
    const label nBFaces = mesh.nFaces() - mesh.nInternalFaces();

    if (faceValues.size() != nBFaces)
    {
        FatalErrorInFunction
            << "Number of values " << faceValues.size()
            << " is not equal to the number of boundary faces in the mesh "
            << nBFaces << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (parRun)
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send

        forAll(patches, patchi)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchi])
             && patches[patchi].size() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchi]);

                label patchStart = procPatch.start()-mesh.nInternalFaces();

                UOPstream toNbr(procPatch.neighbProcNo(), pBufs);
                toNbr << SubField<T>(faceValues, procPatch.size(), patchStart);
            }
        }


        pBufs.finishedSends();


        // Receive and combine.

        forAll(patches, patchi)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchi])
             && patches[patchi].size() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchi]);

                Field<T> nbrPatchInfo(procPatch.size());

                UIPstream fromNeighb(procPatch.neighbProcNo(), pBufs);
                fromNeighb >> nbrPatchInfo;

                top(procPatch, nbrPatchInfo);

                label bFacei = procPatch.start()-mesh.nInternalFaces();

                forAll(nbrPatchInfo, i)
                {
                    cop(faceValues[bFacei++], nbrPatchInfo[i]);
                }
            }
        }
    }

    // Do the cyclics.
    forAll(patches, patchi)
    {
        if (isA<cyclicPolyPatch>(patches[patchi]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchi]);

            if (cycPatch.owner())
            {
                // Owner does all.
                const cyclicPolyPatch& nbrPatch = cycPatch.nbrPatch();
                label ownStart = cycPatch.start()-mesh.nInternalFaces();
                label nbrStart = nbrPatch.start()-mesh.nInternalFaces();

                label sz = cycPatch.size();

                // Transform (copy of) data on both sides
                Field<T> ownVals(SubField<T>(faceValues, sz, ownStart));
                top(nbrPatch, ownVals);

                Field<T> nbrVals(SubField<T>(faceValues, sz, nbrStart));
                top(cycPatch, nbrVals);

                label i0 = ownStart;
                forAll(nbrVals, i)
                {
                    cop(faceValues[i0++], nbrVals[i]);
                }

                label i1 = nbrStart;
                forAll(ownVals, i)
                {
                    cop(faceValues[i1++], ownVals[i]);
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<unsigned nBits, class CombineOp>
void Foam::syncTools::syncFaceList
(
    const polyMesh& mesh,
    PackedList<nBits>& faceValues,
    const CombineOp& cop,
    const bool parRun
)
{
    if (faceValues.size() != mesh.nFaces())
    {
        FatalErrorInFunction
            << "Number of values " << faceValues.size()
            << " is not equal to the number of faces in the mesh "
            << mesh.nFaces() << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (parRun)
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send

        forAll(patches, patchi)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchi])
             && patches[patchi].size() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchi]);

                List<unsigned int> patchInfo(procPatch.size());
                forAll(procPatch, i)
                {
                    patchInfo[i] = faceValues[procPatch.start()+i];
                }

                UOPstream toNbr(procPatch.neighbProcNo(), pBufs);
                toNbr << patchInfo;
            }
        }


        pBufs.finishedSends();

        // Receive and combine.

        forAll(patches, patchi)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchi])
             && patches[patchi].size() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchi]);

                List<unsigned int> patchInfo(procPatch.size());
                {
                    UIPstream fromNbr(procPatch.neighbProcNo(), pBufs);
                    fromNbr >> patchInfo;
                }

                // Combine (bitwise)
                forAll(procPatch, i)
                {
                    unsigned int patchVal = patchInfo[i];
                    label meshFacei = procPatch.start()+i;
                    unsigned int faceVal = faceValues[meshFacei];
                    cop(faceVal, patchVal);
                    faceValues[meshFacei] = faceVal;
                }
            }
        }
    }

    // Do the cyclics.
    forAll(patches, patchi)
    {
        if (isA<cyclicPolyPatch>(patches[patchi]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchi]);

            if (cycPatch.owner())
            {
                // Owner does all.
                const cyclicPolyPatch& nbrPatch = cycPatch.nbrPatch();

                for (label i = 0; i < cycPatch.size(); i++)
                {
                    label meshFace0 = cycPatch.start()+i;
                    unsigned int val0 = faceValues[meshFace0];
                    label meshFace1 = nbrPatch.start()+i;
                    unsigned int val1 = faceValues[meshFace1];

                    unsigned int t = val0;
                    cop(t, val1);
                    faceValues[meshFace0] = t;

                    cop(val1, val0);
                    faceValues[meshFace1] = val1;
                }
            }
        }
    }
}


template<class T>
void Foam::syncTools::swapBoundaryCellList
(
    const polyMesh& mesh,
    const UList<T>& cellData,
    List<T>& neighbourCellData
)
{
    if (cellData.size() != mesh.nCells())
    {
        FatalErrorInFunction
            << "Number of cell values " << cellData.size()
            << " is not equal to the number of cells in the mesh "
            << mesh.nCells() << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    label nBnd = mesh.nFaces()-mesh.nInternalFaces();

    neighbourCellData.setSize(nBnd);

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];
        const labelUList& faceCells = pp.faceCells();
        forAll(faceCells, i)
        {
            label bFacei = pp.start()+i-mesh.nInternalFaces();
            neighbourCellData[bFacei] = cellData[faceCells[i]];
        }
    }
    syncTools::swapBoundaryFaceList(mesh, neighbourCellData);
}


template<unsigned nBits>
void Foam::syncTools::swapFaceList
(
    const polyMesh& mesh,
    PackedList<nBits>& faceValues
)
{
    syncFaceList(mesh, faceValues, eqOp<unsigned int>());
}


template<unsigned nBits, class CombineOp>
void Foam::syncTools::syncPointList
(
    const polyMesh& mesh,
    PackedList<nBits>& pointValues,
    const CombineOp& cop,
    const unsigned int nullValue
)
{
    if (pointValues.size() != mesh.nPoints())
    {
        FatalErrorInFunction
            << "Number of values " << pointValues.size()
            << " is not equal to the number of points in the mesh "
            << mesh.nPoints() << abort(FatalError);
    }

    const globalMeshData& gd = mesh.globalData();
    const labelList& meshPoints = gd.coupledPatch().meshPoints();

    List<unsigned int> cppFld(gd.globalPointSlavesMap().constructSize());
    forAll(meshPoints, i)
    {
        cppFld[i] = pointValues[meshPoints[i]];
    }

    globalMeshData::syncData
    (
        cppFld,
        gd.globalPointSlaves(),
        gd.globalPointTransformedSlaves(),
        gd.globalPointSlavesMap(),
        cop
    );

    // Extract back to mesh
    forAll(meshPoints, i)
    {
        pointValues[meshPoints[i]] = cppFld[i];
    }
}


template<unsigned nBits, class CombineOp>
void Foam::syncTools::syncEdgeList
(
    const polyMesh& mesh,
    PackedList<nBits>& edgeValues,
    const CombineOp& cop,
    const unsigned int nullValue
)
{
    if (edgeValues.size() != mesh.nEdges())
    {
        FatalErrorInFunction
            << "Number of values " << edgeValues.size()
            << " is not equal to the number of edges in the mesh "
            << mesh.nEdges() << abort(FatalError);
    }

    const globalMeshData& gd = mesh.globalData();
    const labelList& meshEdges = gd.coupledPatchMeshEdges();

    List<unsigned int> cppFld(gd.globalEdgeSlavesMap().constructSize());
    forAll(meshEdges, i)
    {
        cppFld[i] = edgeValues[meshEdges[i]];
    }

    globalMeshData::syncData
    (
        cppFld,
        gd.globalEdgeSlaves(),
        gd.globalEdgeTransformedSlaves(),
        gd.globalEdgeSlavesMap(),
        cop
    );

    // Extract back to mesh
    forAll(meshEdges, i)
    {
        edgeValues[meshEdges[i]] = cppFld[i];
    }
}


// ************************************************************************* //
