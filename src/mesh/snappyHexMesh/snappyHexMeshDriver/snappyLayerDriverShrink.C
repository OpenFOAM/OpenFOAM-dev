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
    Shrinking mesh (part of adding cell layers)

\*----------------------------------------------------------------------------*/

#include "snappyLayerDriver.H"
#include "fvMesh.H"
#include "Time.H"
#include "pointFields.H"
#include "motionSmoother.H"
#include "pointData.H"
#include "PointEdgeWave.H"
#include "OBJstream.H"
#include "meshTools.H"
#include "PatchTools.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Calculate inverse sum of edge weights (currently always 1.0)
void Foam::snappyLayerDriver::sumWeights
(
    const PackedBoolList& isMasterEdge,
    const labelList& meshEdges,
    const labelList& meshPoints,
    const edgeList& edges,
    scalarField& invSumWeight
) const
{
    const pointField& pts = meshRefiner_.mesh().points();

    invSumWeight = 0;

    forAll(edges, edgeI)
    {
        if (isMasterEdge.get(meshEdges[edgeI]) == 1)
        {
            const edge& e = edges[edgeI];
            // scalar eWeight = edgeWeights[edgeI];
            // scalar eWeight = 1.0;

            scalar eMag = max
            (
                vSmall,
                mag
                (
                    pts[meshPoints[e[1]]]
                  - pts[meshPoints[e[0]]]
                )
            );
            scalar eWeight = 1.0/eMag;

            invSumWeight[e[0]] += eWeight;
            invSumWeight[e[1]] += eWeight;
        }
    }

    syncTools::syncPointList
    (
        meshRefiner_.mesh(),
        meshPoints,
        invSumWeight,
        plusEqOp<scalar>(),
        scalar(0)         // null value
    );

    forAll(invSumWeight, pointi)
    {
        scalar w = invSumWeight[pointi];

        if (w > 0.0)
        {
            invSumWeight[pointi] = 1.0/w;
        }
    }
}


// Smooth field on moving patch
void Foam::snappyLayerDriver::smoothField
(
    const motionSmoother& meshMover,
    const PackedBoolList& isMasterPoint,
    const PackedBoolList& isMasterEdge,
    const labelList& meshEdges,
    const scalarField& fieldMin,
    const label nSmoothDisp,
    scalarField& field
) const
{
    const indirectPrimitivePatch& pp = meshMover.patch();
    const edgeList& edges = pp.edges();
    const labelList& meshPoints = pp.meshPoints();

    scalarField invSumWeight(pp.nPoints());
    sumWeights
    (
        isMasterEdge,
        meshEdges,
        meshPoints,
        edges,
        invSumWeight
    );

    // Get smoothly varying patch field.
    Info<< "shrinkMeshDistance : Smoothing field ..." << endl;

    for (label iter = 0; iter < nSmoothDisp; iter++)
    {
        scalarField average(pp.nPoints());
        averageNeighbours
        (
            meshMover.mesh(),
            isMasterEdge,
            meshEdges,
            meshPoints,
            edges,
            invSumWeight,
            field,
            average
        );

        // Transfer to field
        forAll(field, pointi)
        {
            // full smoothing neighbours + point value
            average[pointi] = 0.5*(field[pointi]+average[pointi]);

            // perform monotonic smoothing
            if
            (
                average[pointi] < field[pointi]
             && average[pointi] >= fieldMin[pointi]
            )
            {
                field[pointi] = average[pointi];
            }
        }

        // Do residual calculation every so often.
        if ((iter % 10) == 0)
        {
            scalar resid = meshRefinement::gAverage
            (
                meshMover.mesh(),
                isMasterPoint,
                meshPoints,
                mag(field-average)()
            );
            Info<< "    Iteration " << iter << "   residual " << resid << endl;
        }
    }
}
//XXXXXXXXX
//void Foam::snappyLayerDriver::smoothField
//(
//    const motionSmoother& meshMover,
//    const PackedBoolList& isMasterEdge,
//    const labelList& meshEdges,
//    const scalarField& fieldMin,
//    const label nSmoothDisp,
//    scalarField& field
//) const
//{
//    const indirectPrimitivePatch& pp = meshMover.patch();
//    const edgeList& edges = pp.edges();
//    const labelList& meshPoints = pp.meshPoints();
//
//    scalarField invSumWeight(pp.nPoints());
//    sumWeights
//    (
//        isMasterEdge,
//        meshEdges,
//        meshPoints,
//        edges,
//        invSumWeight
//    );
//
//    // Get smoothly varying patch field.
//    Info<< "shrinkMeshDistance : (lambda-mu) Smoothing field ..." << endl;
//
//
//    const scalar lambda = 0.33;
//    const scalar mu = -0.34;
//
//    for (label iter = 0; iter < 90; iter++)
//    {
//        scalarField average(pp.nPoints());
//
//        // Calculate average of field
//        averageNeighbours
//        (
//            meshMover.mesh(),
//            isMasterEdge,
//            meshEdges,
//            meshPoints,
//            pp.edges(),
//            invSumWeight,
//            field,
//            average
//        );
//
//        forAll(field, i)
//        {
//            if (field[i] >= fieldMin[i])
//            {
//                field[i] = (1-lambda)*field[i]+lambda*average[i];
//            }
//        }
//
//
//        // Calculate average of field
//        averageNeighbours
//        (
//            meshMover.mesh(),
//            isMasterEdge,
//            meshEdges,
//            meshPoints,
//            pp.edges(),
//            invSumWeight,
//            field,
//            average
//        );
//
//        forAll(field, i)
//        {
//            if (field[i] >= fieldMin[i])
//            {
//                field[i] = (1-mu)*field[i]+mu*average[i];
//            }
//        }
//
//
//        // Do residual calculation every so often.
//        if ((iter % 10) == 0)
//        {
//            Info<< "    Iteration " << iter << "   residual "
//                <<  gSum(mag(field-average))
//                   /returnReduce(average.size(), sumOp<label>())
//                << endl;
//        }
//    }
//}
//XXXXXXXXX

// Smooth normals on moving patch.
void Foam::snappyLayerDriver::smoothPatchNormals
(
    const motionSmoother& meshMover,
    const PackedBoolList& isMasterPoint,
    const PackedBoolList& isMasterEdge,
    const labelList& meshEdges,
    const label nSmoothDisp,
    pointField& normals
) const
{
    const indirectPrimitivePatch& pp = meshMover.patch();
    const edgeList& edges = pp.edges();
    const labelList& meshPoints = pp.meshPoints();

    // Calculate inverse sum of weights

    scalarField invSumWeight(pp.nPoints());
    sumWeights
    (
        isMasterEdge,
        meshEdges,
        meshPoints,
        edges,
        invSumWeight
    );

    // Get smoothly varying internal normals field.
    Info<< "shrinkMeshDistance : Smoothing normals ..." << endl;

    for (label iter = 0; iter < nSmoothDisp; iter++)
    {
        vectorField average(pp.nPoints());
        averageNeighbours
        (
            meshMover.mesh(),
            isMasterEdge,
            meshEdges,
            meshPoints,
            edges,
            invSumWeight,
            normals,
            average
        );

        // Do residual calculation every so often.
        if ((iter % 10) == 0)
        {
            scalar resid = meshRefinement::gAverage
            (
                meshMover.mesh(),
                isMasterPoint,
                meshPoints,
                mag(normals-average)()
            );
            Info<< "    Iteration " << iter << "   residual " << resid << endl;
        }

        // Transfer to normals vector field
        forAll(average, pointi)
        {
            // full smoothing neighbours + point value
            average[pointi] = 0.5*(normals[pointi]+average[pointi]);
            normals[pointi] = average[pointi];
            normals[pointi] /= mag(normals[pointi]) + vSmall;
        }
    }
}


// Smooth normals in interior.
void Foam::snappyLayerDriver::smoothNormals
(
    const label nSmoothDisp,
    const PackedBoolList& isMasterPoint,
    const PackedBoolList& isMasterEdge,
    const labelList& fixedPoints,
    pointVectorField& normals
) const
{
    // Get smoothly varying internal normals field.
    Info<< "shrinkMeshDistance : Smoothing normals in interior ..." << endl;

    const fvMesh& mesh = meshRefiner_.mesh();
    const edgeList& edges = mesh.edges();

    // Points that do not change.
    PackedBoolList isFixedPoint(mesh.nPoints());

    // Internal points that are fixed
    forAll(fixedPoints, i)
    {
        label meshPointi = fixedPoints[i];
        isFixedPoint.set(meshPointi, 1);
    }

    // Make sure that points that are coupled to meshPoints but not on a patch
    // are fixed as well
    syncTools::syncPointList(mesh, isFixedPoint, maxEqOp<unsigned int>(), 0);


    // Correspondence between local edges/points and mesh edges/points
    const labelList meshEdges(identity(mesh.nEdges()));
    const labelList meshPoints(identity(mesh.nPoints()));

    // Calculate inverse sum of weights

    scalarField invSumWeight(mesh.nPoints(), 0);
    sumWeights
    (
        isMasterEdge,
        meshEdges,
        meshPoints,
        edges,
        invSumWeight
    );

    for (label iter = 0; iter < nSmoothDisp; iter++)
    {
        vectorField average(mesh.nPoints());
        averageNeighbours
        (
            mesh,
            isMasterEdge,
            meshEdges,
            meshPoints,
            edges,
            invSumWeight,
            normals,
            average
        );

        // Do residual calculation every so often.
        if ((iter % 10) == 0)
        {
            scalar resid = meshRefinement::gAverage
            (
                mesh,
                isMasterPoint,
                mag(normals-average)()
            );
            Info<< "    Iteration " << iter << "   residual " << resid << endl;
        }


        // Transfer to normals vector field
        forAll(average, pointi)
        {
            if (isFixedPoint.get(pointi) == 0)
            {
                // full smoothing neighbours + point value
                average[pointi] = 0.5*(normals[pointi]+average[pointi]);
                normals[pointi] = average[pointi];
                normals[pointi] /= mag(normals[pointi]) + vSmall;
            }
        }
    }
}


// Tries and find a medial axis point. Done by comparing vectors to nearest
// wall point for both vertices of edge.
bool Foam::snappyLayerDriver::isMaxEdge
(
    const List<pointData>& pointWallDist,
    const label edgeI,
    const scalar minCos
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const pointField& points = mesh.points();

    // Do not mark edges with one side on moving wall.

    const edge& e = mesh.edges()[edgeI];

    vector v0(points[e[0]] - pointWallDist[e[0]].origin());
    scalar magV0(mag(v0));

    if (magV0 < small)
    {
        return false;
    }

    vector v1(points[e[1]] - pointWallDist[e[1]].origin());
    scalar magV1(mag(v1));

    if (magV1 < small)
    {
        return false;
    }


    //- Detect based on vector to nearest point differing for both endpoints
    // v0 /= magV0;
    // v1 /= magV1;
    //
    //// Test angle.
    // if ((v0 & v1) < minCos)
    //{
    //    return true;
    //}
    // else
    //{
    //    return false;
    //}

    //- Detect based on extrusion vector differing for both endpoints
    //  the idea is that e.g. a sawtooth wall can still be extruded
    //  successfully as long as it is done all to the same direction.
    if ((pointWallDist[e[0]].v() & pointWallDist[e[1]].v()) < minCos)
    {
        return true;
    }
    else
    {
        return false;
    }
}


// Stop layer growth where mesh wraps around edge with a
// large feature angle
void Foam::snappyLayerDriver::handleFeatureAngleLayerTerminations
(
    const scalar minCos,
    const PackedBoolList& isMasterPoint,
    const indirectPrimitivePatch& pp,
    const labelList& meshEdges,
    List<extrudeMode>& extrudeStatus,
    pointField& patchDisp,
    labelList& patchNLayers,
    label& nPointCounter
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    // Mark faces that have all points extruded
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    boolList extrudedFaces(pp.size(), true);

    forAll(pp.localFaces(), facei)
    {
        const face& f = pp.localFaces()[facei];

        forAll(f, fp)
        {
            if (extrudeStatus[f[fp]] == NOEXTRUDE)
            {
                extrudedFaces[facei] = false;
                break;
            }
        }
    }



    // label nOldPointCounter = nPointCounter;

    // Detect situation where two featureedge-neighbouring faces are partly or
    // not extruded and the edge itself is extruded. In this case unmark the
    // edge for extrusion.


    List<List<point>> edgeFaceNormals(pp.nEdges());
    List<List<bool>> edgeFaceExtrude(pp.nEdges());

    const labelListList& edgeFaces = pp.edgeFaces();
    const vectorField& faceNormals = pp.faceNormals();
    const labelList& meshPoints = pp.meshPoints();

    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        edgeFaceNormals[edgeI].setSize(eFaces.size());
        edgeFaceExtrude[edgeI].setSize(eFaces.size());
        forAll(eFaces, i)
        {
            label facei = eFaces[i];
            edgeFaceNormals[edgeI][i] = faceNormals[facei];
            edgeFaceExtrude[edgeI][i] = extrudedFaces[facei];
        }
    }

    syncTools::syncEdgeList
    (
        mesh,
        meshEdges,
        edgeFaceNormals,
        globalMeshData::ListPlusEqOp<List<point>>(),   // combine operator
        List<point>()               // null value
    );

    syncTools::syncEdgeList
    (
        mesh,
        meshEdges,
        edgeFaceExtrude,
        globalMeshData::ListPlusEqOp<List<bool>>(),    // combine operator
        List<bool>()                // null value
    );


    forAll(edgeFaceNormals, edgeI)
    {
        const List<point>& eFaceNormals = edgeFaceNormals[edgeI];
        const List<bool>& eFaceExtrude = edgeFaceExtrude[edgeI];

        if (eFaceNormals.size() == 2)
        {
            const edge& e = pp.edges()[edgeI];
            label v0 = e[0];
            label v1 = e[1];

            if
            (
                extrudeStatus[v0] != NOEXTRUDE
             || extrudeStatus[v1] != NOEXTRUDE
            )
            {
                if (!eFaceExtrude[0] || !eFaceExtrude[1])
                {
                    const vector& n0 = eFaceNormals[0];
                    const vector& n1 = eFaceNormals[1];

                    if ((n0 & n1) < minCos)
                    {
                        if
                        (
                            unmarkExtrusion
                            (
                                v0,
                                patchDisp,
                                patchNLayers,
                                extrudeStatus
                            )
                        )
                        {
                            if (isMasterPoint[meshPoints[v0]])
                            {
                                nPointCounter++;
                            }
                        }
                        if
                        (
                            unmarkExtrusion
                            (
                                v1,
                                patchDisp,
                                patchNLayers,
                                extrudeStatus
                            )
                        )
                        {
                            if (isMasterPoint[meshPoints[v1]])
                            {
                                nPointCounter++;
                            }
                        }
                    }
                }
            }
        }
    }

    // Info<< "Added "
    //    << returnReduce(nPointCounter-nOldPointCounter, sumOp<label>())
    //    << " point not to extrude." << endl;
}


// Find isolated islands (points, edges and faces and layer terminations)
// in the layer mesh and stop any layer growth at these points.
void Foam::snappyLayerDriver::findIsolatedRegions
(
    const scalar minCosLayerTermination,
    const PackedBoolList& isMasterPoint,
    const PackedBoolList& isMasterEdge,
    const indirectPrimitivePatch& pp,
    const labelList& meshEdges,
    const scalarField& minThickness,
    List<extrudeMode>& extrudeStatus,
    pointField& patchDisp,
    labelList& patchNLayers
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    Info<< "shrinkMeshDistance : Removing isolated regions ..." << endl;

    // Keep count of number of points unextruded
    label nPointCounter = 0;

    while (true)
    {
        // Stop layer growth where mesh wraps around edge with a
        // large feature angle
        handleFeatureAngleLayerTerminations
        (
            minCosLayerTermination,
            isMasterPoint,
            pp,
            meshEdges,

            extrudeStatus,
            patchDisp,
            patchNLayers,
            nPointCounter
        );

        syncPatchDisplacement
        (
            pp,
            minThickness,
            patchDisp,
            patchNLayers,
            extrudeStatus
        );


        // Do not extrude from point where all neighbouring
        // faces are not grown
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        boolList extrudedFaces(pp.size(), true);
        forAll(pp.localFaces(), facei)
        {
            const face& f = pp.localFaces()[facei];
            forAll(f, fp)
            {
                if (extrudeStatus[f[fp]] == NOEXTRUDE)
                {
                    extrudedFaces[facei] = false;
                    break;
                }
            }
        }

        const labelListList& pointFaces = pp.pointFaces();

        boolList keptPoints(pp.nPoints(), false);
        forAll(keptPoints, patchPointi)
        {
            const labelList& pFaces = pointFaces[patchPointi];

            forAll(pFaces, i)
            {
                label facei = pFaces[i];
                if (extrudedFaces[facei])
                {
                    keptPoints[patchPointi] = true;
                    break;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            keptPoints,
            orEqOp<bool>(),
            false               // null value
        );

        label nChanged = 0;

        forAll(keptPoints, patchPointi)
        {
            if (!keptPoints[patchPointi])
            {
                if
                (
                    unmarkExtrusion
                    (
                        patchPointi,
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    )
                )
                {
                   nPointCounter++;
                   nChanged++;
                }
            }
        }


        if (returnReduce(nChanged, sumOp<label>()) == 0)
        {
            break;
        }
    }

    const edgeList& edges = pp.edges();


    // Count number of mesh edges using a point
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList isolatedPoint(pp.nPoints(),0);

    forAll(edges, edgeI)
    {
        if (isMasterEdge.get(meshEdges[edgeI]) == 1)
        {
            const edge& e = edges[edgeI];

            label v0 = e[0];
            label v1 = e[1];

            if (extrudeStatus[v1] != NOEXTRUDE)
            {
                isolatedPoint[v0] += 1;
            }
            if (extrudeStatus[v0] != NOEXTRUDE)
            {
                isolatedPoint[v1] += 1;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        isolatedPoint,
        plusEqOp<label>(),
        label(0)        // null value
    );

    // stop layer growth on isolated faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    forAll(pp, facei)
    {
        const face& f = pp.localFaces()[facei];
        bool failed = false;
        forAll(f, fp)
        {
            if (isolatedPoint[f[fp]] > 2)
            {
                failed = true;
                break;
            }
        }
        bool allPointsExtruded = true;
        if (!failed)
        {
            forAll(f, fp)
            {
                if (extrudeStatus[f[fp]] == NOEXTRUDE)
                {
                    allPointsExtruded = false;
                    break;
                }
            }

            if (allPointsExtruded)
            {
                forAll(f, fp)
                {
                    if
                    (
                        unmarkExtrusion
                        (
                            f[fp],
                            patchDisp,
                            patchNLayers,
                            extrudeStatus
                        )
                    )
                    {
                        nPointCounter++;
                    }
                }
            }
        }
    }

    reduce(nPointCounter, sumOp<label>());
    Info<< "Number isolated points extrusion stopped : "<< nPointCounter
        << endl;
}


// Calculates medial axis fields:
// dispVec     : normal of nearest wall. Where this normal changes direction
//               defines the medial axis
// medialDist  : distance to medial axis
// medialRatio : ratio of medial distance to wall distance.
//               (1 at wall, 0 at medial axis)
void Foam::snappyLayerDriver::medialAxisSmoothingInfo
(
    const motionSmoother& meshMover,
    const label nSmoothNormals,
    const label nSmoothSurfaceNormals,
    const scalar minMedialAxisAngleCos,
    const scalar featureAngle,

    pointVectorField& dispVec,
    pointScalarField& medialRatio,
    pointScalarField& medialDist,
    pointVectorField& medialVec
) const
{

    Info<< "medialAxisSmoothingInfo :"
        << " Calculate distance to Medial Axis ..." << endl;

    const polyMesh& mesh = meshMover.mesh();
    const pointField& points = mesh.points();

    const indirectPrimitivePatch& pp = meshMover.patch();
    const labelList& meshPoints = pp.meshPoints();

    // Predetermine mesh edges
    // ~~~~~~~~~~~~~~~~~~~~~~~

    // Precalulate master point/edge (only relevant for shared points/edges)
    const PackedBoolList isMasterPoint(syncTools::getMasterPoints(mesh));
    const PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh));
    // Precalculate meshEdge per pp edge
    const labelList meshEdges
    (
        pp.meshEdges
        (
            mesh.edges(),
            mesh.pointEdges()
        )
    );


    // Determine pointNormal
    // ~~~~~~~~~~~~~~~~~~~~~

    pointField pointNormals(PatchTools::pointNormals(mesh, pp));

    // pointNormals
    if (debug&meshRefinement::MESH || debug&meshRefinement::LAYERINFO)
    {
        pointField meshPointNormals(mesh.nPoints(), point(1, 0, 0));
        UIndirectList<point>(meshPointNormals, meshPoints) = pointNormals;
        meshRefinement::testSyncPointList
        (
            "pointNormals",
            mesh,
            meshPointNormals
        );
    }

    // Smooth patch normal vectors
    smoothPatchNormals
    (
        meshMover,
        isMasterPoint,
        isMasterEdge,
        meshEdges,
        nSmoothSurfaceNormals,
        pointNormals
    );

    // smoothed pointNormals
    if (debug&meshRefinement::MESH || debug&meshRefinement::LAYERINFO)
    {
        pointField meshPointNormals(mesh.nPoints(), point(1, 0, 0));
        UIndirectList<point>(meshPointNormals, meshPoints) = pointNormals;
        meshRefinement::testSyncPointList
        (
            "smoothed pointNormals",
            mesh,
            meshPointNormals
        );
    }

    // Calculate distance to pp points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Distance to wall
    List<pointData> pointWallDist(mesh.nPoints());

    // Dummy additional info for PointEdgeWave
    int dummyTrackData = 0;


    // 1. Calculate distance to points where displacement is specified.
    {
        // Seed data.
        List<pointData> wallInfo(meshPoints.size());

        forAll(meshPoints, patchPointi)
        {
            label pointi = meshPoints[patchPointi];
            wallInfo[patchPointi] = pointData
            (
                points[pointi],
                0.0,
                pointi,                       // passive scalar
                pointNormals[patchPointi]     // surface normals
            );
        }

        // Do all calculations
        List<pointData> edgeWallDist(mesh.nEdges());
        PointEdgeWave<pointData> wallDistCalc
        (
            mesh,
            meshPoints,
            wallInfo,
            pointWallDist,
            edgeWallDist,
            mesh.globalData().nTotalPoints(),   // max iterations
            dummyTrackData
        );


        label nUnvisit = returnReduce
        (
            wallDistCalc.getUnsetPoints(),
            sumOp<label>()
        );

        if (nUnvisit > 0)
        {
            WarningInFunction
                << "Walking did not visit all points." << nl
                << "    Did not visit " << nUnvisit
                << " out of " << mesh.globalData().nTotalPoints()
                << " points. This is not necessarily a problem" << nl
                << "    and might be due to faceZones splitting of part"
                << " of the domain." << nl << endl;
        }
    }


    // Check sync of wall distance
    if (debug&meshRefinement::MESH || debug&meshRefinement::LAYERINFO)
    {
        pointField origin(pointWallDist.size());
        scalarField distSqr(pointWallDist.size());
        // NA scalarField passiveS(pointWallDist.size());
        pointField passiveV(pointWallDist.size());
        forAll(pointWallDist, pointi)
        {
            origin[pointi] = pointWallDist[pointi].origin();
            distSqr[pointi] = pointWallDist[pointi].distSqr();
            // passiveS[pointi] = pointWallDist[pointi].s();
            passiveV[pointi] = pointWallDist[pointi].v();
        }
        meshRefinement::testSyncPointList("origin", mesh, origin);
        meshRefinement::testSyncPointList("distSqr", mesh, distSqr);
        // meshRefinement::testSyncPointList("passiveS", mesh, passiveS);
        meshRefinement::testSyncPointList("passiveV", mesh, passiveV);
    }


    // 2. Find points with max distance and transport information back to
    //    wall.
    {
        List<pointData> pointMedialDist(mesh.nPoints());
        List<pointData> edgeMedialDist(mesh.nEdges());

        // Seed point data.
        DynamicList<pointData> maxInfo(meshPoints.size());
        DynamicList<label> maxPoints(meshPoints.size());

        // 1. Medial axis points

        const edgeList& edges = mesh.edges();

        forAll(edges, edgeI)
        {
            const edge& e = edges[edgeI];

            if
            (
                !pointWallDist[e[0]].valid(dummyTrackData)
             || !pointWallDist[e[1]].valid(dummyTrackData)
            )
            {
                // Unvisited point. See above about nUnvisit warning
            }
            else if (isMaxEdge(pointWallDist, edgeI, minMedialAxisAngleCos))
            {
                // Both end points of edge have very different nearest wall
                // point. Mark both points as medial axis points.

                // Approximate medial axis location on edge.
                // const point medialAxisPt = e.centre(points);
                vector eVec = e.vec(points);
                scalar eMag = mag(eVec);
                if (eMag > vSmall)
                {
                    eVec /= eMag;

                    // Calculate distance along edge
                    const point& p0 = points[e[0]];
                    const point& p1 = points[e[1]];
                    scalar dist0 = (p0-pointWallDist[e[0]].origin()) & eVec;
                    scalar dist1 = (pointWallDist[e[1]].origin()-p1) & eVec;
                    scalar s = 0.5*(dist1+eMag+dist0);

                    point medialAxisPt;
                    if (s <= dist0)
                    {
                        medialAxisPt = p0;
                    }
                    else if (s >= dist0+eMag)
                    {
                        medialAxisPt = p1;
                    }
                    else
                    {
                        medialAxisPt = p0+(s-dist0)*eVec;
                    }

                    forAll(e, ep)
                    {
                        label pointi = e[ep];

                        if (!pointMedialDist[pointi].valid(dummyTrackData))
                        {
                            maxPoints.append(pointi);
                            maxInfo.append
                            (
                                pointData
                                (
                                    medialAxisPt,   // points[pointi],
                                    magSqr(points[pointi]-medialAxisPt),//0.0,
                                    pointi,         // passive data
                                    Zero    // passive data
                                )
                            );
                            pointMedialDist[pointi] = maxInfo.last();
                        }
                    }
                }
            }
        }


        // 2. Seed non-adapt patches
        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        labelHashSet adaptPatches(meshMover.adaptPatchIDs());


        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];
            const pointPatchVectorField& pvf =
                meshMover.displacement().boundaryField()[patchi];

            if
            (
                !pp.coupled()
             && !isA<emptyPolyPatch>(pp)
             && !adaptPatches.found(patchi)
            )
            {
                const labelList& meshPoints = pp.meshPoints();

                // Check the type of the patchField. The types are
                //  - fixedValue (0 or more layers) but the >0 layers have
                //    already been handled in the adaptPatches loop
                //  - constraint (but not coupled) types, e.g. symmetryPlane,
                //    slip.
                if (pvf.fixesValue())
                {
                    // Disable all movement on fixedValue patchFields
                    Info<< "Inserting all points on patch " << pp.name()
                        << endl;

                    forAll(meshPoints, i)
                    {
                        label pointi = meshPoints[i];
                        if (!pointMedialDist[pointi].valid(dummyTrackData))
                        {
                            maxPoints.append(pointi);
                            maxInfo.append
                            (
                                pointData
                                (
                                    points[pointi],
                                    0.0,
                                    pointi,         // passive data
                                    Zero    // passive data
                                )
                            );
                            pointMedialDist[pointi] = maxInfo.last();
                        }
                    }
                }
                else
                {
                    // Based on geometry: analyse angle w.r.t. nearest moving
                    // point. In the pointWallDist we transported the
                    // normal as the passive vector. Note that this points
                    // out of the originating wall so inside of the domain
                    // on this patch.
                    Info<< "Inserting points on patch " << pp.name()
                        << " if angle to nearest layer patch > "
                        << featureAngle << " degrees." << endl;

                    scalar featureAngleCos = Foam::cos(degToRad(featureAngle));
                    pointField pointNormals(PatchTools::pointNormals(mesh, pp));

                    forAll(meshPoints, i)
                    {
                        label pointi = meshPoints[i];

                        if
                        (
                            pointWallDist[pointi].valid(dummyTrackData)
                        && !pointMedialDist[pointi].valid(dummyTrackData)
                        )
                        {
                            // Check if angle not too large.
                            scalar cosAngle =
                            (
                               -pointWallDist[pointi].v()
                              & pointNormals[i]
                            );
                            if (cosAngle > featureAngleCos)
                            {
                                // Extrusion direction practically perpendicular
                                // to the patch. Disable movement at the patch.

                                maxPoints.append(pointi);
                                maxInfo.append
                                (
                                    pointData
                                    (
                                        points[pointi],
                                        0.0,
                                        pointi,         // passive data
                                        Zero    // passive data
                                    )
                                );
                                pointMedialDist[pointi] = maxInfo.last();
                            }
                            else
                            {
                                // Extrusion direction makes angle with patch
                                // so allow sliding - don't insert zero points
                            }
                        }
                    }
                }
            }
        }

        maxInfo.shrink();
        maxPoints.shrink();

        // Do all calculations
        PointEdgeWave<pointData> medialDistCalc
        (
            mesh,
            maxPoints,
            maxInfo,

            pointMedialDist,
            edgeMedialDist,
            mesh.globalData().nTotalPoints(),   // max iterations
            dummyTrackData
        );

        // Extract medial axis distance as pointScalarField
        forAll(pointMedialDist, pointi)
        {
            medialDist[pointi] = Foam::sqrt(pointMedialDist[pointi].distSqr());
            medialVec[pointi] = pointMedialDist[pointi].origin();
        }

        // Check
        if (debug&meshRefinement::MESH || debug&meshRefinement::LAYERINFO)
        {
            meshRefinement::testSyncPointList("medialDist", mesh, medialDist);
            meshRefinement::testSyncPointList("medialVec", mesh, medialVec);
        }
    }

    // Extract transported surface normals as pointVectorField
    forAll(dispVec, i)
    {
        if (!pointWallDist[i].valid(dummyTrackData))
        {
            dispVec[i] = vector(1, 0, 0);
        }
        else
        {
            dispVec[i] = pointWallDist[i].v();
        }
    }

    // Smooth normal vectors. Do not change normals on pp.meshPoints
    smoothNormals
    (
        nSmoothNormals,
        isMasterPoint,
        isMasterEdge,
        meshPoints,
        dispVec
    );

    if (debug&meshRefinement::MESH || debug&meshRefinement::LAYERINFO)
    {
        meshRefinement::testSyncPointList("smoothed dispVec", mesh, dispVec);
    }

    // Calculate ratio point medial distance to point wall distance
    forAll(medialRatio, pointi)
    {
        if (!pointWallDist[pointi].valid(dummyTrackData))
        {
            medialRatio[pointi] = 0.0;
        }
        else
        {
            scalar wDist2 = pointWallDist[pointi].distSqr();
            scalar mDist = medialDist[pointi];

            if (wDist2 < sqr(small) && mDist < small)
            //- Note: maybe less strict:
            //(
            //    wDist2 < sqr(meshRefiner_.mergeDistance())
            // && mDist < meshRefiner_.mergeDistance()
            //)
            {
                medialRatio[pointi] = 0.0;
            }
            else
            {
                medialRatio[pointi] = mDist / (Foam::sqrt(wDist2) + mDist);
            }
        }
    }


    if (debug&meshRefinement::MESH || debug&meshRefinement::LAYERINFO)
    {
        // medialRatio
        meshRefinement::testSyncPointList("medialRatio", mesh, medialRatio);
    }


    if (debug&meshRefinement::MESH || debug&meshRefinement::LAYERINFO)
    {
        Info<< "medialAxisSmoothingInfo :"
            << " Writing:" << nl
            << "    " << dispVec.name()
            << " : normalised direction of nearest displacement" << nl
            << "    " << medialDist.name()
            << " : distance to medial axis" << nl
            << "    " << medialVec.name()
            << " : nearest point on medial axis" << nl
            << "    " << medialRatio.name()
            << " : ratio of medial distance to wall distance" << nl
            << endl;
        meshRefiner_.mesh().setInstance(meshRefiner_.timeName());
        meshRefiner_.write
        (
            meshRefinement::debugType(debug),
            meshRefinement::writeType
            (
                meshRefinement::writeLevel()
              | meshRefinement::WRITEMESH
            ),
            mesh.time().path()/meshRefiner_.timeName()
        );
        dispVec.write();
        medialDist.write();
        medialVec.write();
        medialRatio.write();
    }
}


void Foam::snappyLayerDriver::shrinkMeshMedialDistance
(
    motionSmoother& meshMover,
    const dictionary& meshQualityDict,
    const List<labelPair>& baffles,
    const label nSmoothPatchThickness,
    const label nSmoothDisplacement,
    const scalar maxThicknessToMedialRatio,
    const label nAllowableErrors,
    const label nSnap,
    const scalar minCosLayerTermination,

    const scalarField& layerThickness,
    const scalarField& minThickness,

    const pointVectorField& dispVec,
    const pointScalarField& medialRatio,
    const pointScalarField& medialDist,
    const pointVectorField& medialVec,

    List<extrudeMode>& extrudeStatus,
    pointField& patchDisp,
    labelList& patchNLayers
) const
{
    Info<< "shrinkMeshMedialDistance : Smoothing using Medial Axis ..." << endl;

    const polyMesh& mesh = meshMover.mesh();

    const indirectPrimitivePatch& pp = meshMover.patch();
    const labelList& meshPoints = pp.meshPoints();

    // Precalulate master points/edge (only relevant for shared points/edges)
    const PackedBoolList isMasterPoint(syncTools::getMasterPoints(mesh));
    const PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh));
    // Precalculate meshEdge per pp edge
    const labelList meshEdges
    (
        pp.meshEdges
        (
            mesh.edges(),
            mesh.pointEdges()
        )
    );


    scalarField thickness(layerThickness.size());

    thickness = mag(patchDisp);

    forAll(thickness, patchPointi)
    {
        if (extrudeStatus[patchPointi] == NOEXTRUDE)
        {
            thickness[patchPointi] = 0.0;
        }
    }

    label numThicknessRatioExclude = 0;

    // reduce thickness where thickness/medial axis distance large
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    autoPtr<OBJstream> str;
    if (debug)
    {
        str.reset
        (
            new OBJstream
            (
                mesh.time().path()
              / "thicknessRatioExcludePoints_"
              + meshRefiner_.timeName()
              + ".obj"
            )
        );
        Info<< "Writing points with too large an extrusion distance to "
            << str().name() << endl;
    }

    autoPtr<OBJstream> medialVecStr;
    if (debug)
    {
        medialVecStr.reset
        (
            new OBJstream
            (
                mesh.time().path()
              / "thicknessRatioExcludeMedialVec_"
              + meshRefiner_.timeName()
              + ".obj"
            )
        );
        Info<< "Writing points with too large an extrusion distance to "
            << medialVecStr().name() << endl;
    }

    forAll(meshPoints, patchPointi)
    {
        if (extrudeStatus[patchPointi] != NOEXTRUDE)
        {
            label pointi = meshPoints[patchPointi];

            //- Option 1: look only at extrusion thickness v.s. distance
            //  to nearest (medial axis or static) point.
            scalar mDist = medialDist[pointi];
            scalar thicknessRatio = thickness[patchPointi]/(mDist+vSmall);

            //- Option 2: Look at component in the direction
            //  of nearest (medial axis or static) point.
            vector n =
                patchDisp[patchPointi]
              / (mag(patchDisp[patchPointi]) + vSmall);
            vector mVec = mesh.points()[pointi]-medialVec[pointi];
            mVec /= mag(mVec)+vSmall;
            thicknessRatio *= (n&mVec);

            if (thicknessRatio > maxThicknessToMedialRatio)
            {
                // Truncate thickness.
                if (debug)
                {
                    Pout<< "truncating displacement at "
                        << mesh.points()[pointi]
                        << " from " << thickness[patchPointi]
                        << " to "
                        <<  0.5
                           *(
                                minThickness[patchPointi]
                               +thickness[patchPointi]
                            )
                        << " medial direction:" << mVec
                        << " extrusion direction:" << n
                        << " with thicknessRatio:" << thicknessRatio
                        << endl;
                }

                thickness[patchPointi] =
                    0.5*(minThickness[patchPointi]+thickness[patchPointi]);

                patchDisp[patchPointi] = thickness[patchPointi]*n;

                if (isMasterPoint[pointi])
                {
                    numThicknessRatioExclude++;
                }

                if (str.valid())
                {
                    const point& pt = mesh.points()[pointi];
                    str().write(linePointRef(pt, pt+patchDisp[patchPointi]));
                }
                if (medialVecStr.valid())
                {
                    const point& pt = mesh.points()[pointi];
                    medialVecStr().write
                    (
                        linePointRef
                        (
                            pt,
                            medialVec[pointi]
                        )
                    );
                }
            }
        }
    }

    reduce(numThicknessRatioExclude, sumOp<label>());
    Info<< "shrinkMeshMedialDistance : Reduce layer thickness at "
        << numThicknessRatioExclude
        << " nodes where thickness to medial axis distance is large " << endl;


    // find points where layer growth isolated to a lone point, edge or face

    findIsolatedRegions
    (
        minCosLayerTermination,
        isMasterPoint,
        isMasterEdge,
        pp,
        meshEdges,
        minThickness,

        extrudeStatus,
        patchDisp,
        patchNLayers
    );

    // Update thickess for changed extrusion
    forAll(thickness, patchPointi)
    {
        if (extrudeStatus[patchPointi] == NOEXTRUDE)
        {
            thickness[patchPointi] = 0.0;
        }
    }


    // smooth layer thickness on moving patch
    smoothField
    (
        meshMover,
        isMasterPoint,
        isMasterEdge,
        meshEdges,
        minThickness,
        nSmoothPatchThickness,

        thickness
    );

    // Dummy additional info for PointEdgeWave
    int dummyTrackData = 0;

    List<pointData> pointWallDist(mesh.nPoints());

    const pointField& points = mesh.points();
    // 1. Calculate distance to points where displacement is specified.
    // This wave is used to transport layer thickness
    {
        // Distance to wall and medial axis on edges.
        List<pointData> edgeWallDist(mesh.nEdges());
        labelList wallPoints(meshPoints.size());

        // Seed data.
        List<pointData> wallInfo(meshPoints.size());

        forAll(meshPoints, patchPointi)
        {
            label pointi = meshPoints[patchPointi];
            wallPoints[patchPointi] = pointi;
            wallInfo[patchPointi] = pointData
            (
                points[pointi],
                0.0,
                thickness[patchPointi],       // transport layer thickness
                Zero                  // passive vector
            );
        }

        // Do all calculations
        PointEdgeWave<pointData> wallDistCalc
        (
            mesh,
            wallPoints,
            wallInfo,
            pointWallDist,
            edgeWallDist,
            mesh.globalData().nTotalPoints(),   // max iterations
            dummyTrackData
        );
    }

    // Calculate scaled displacement vector
    pointVectorField& displacement = meshMover.displacement();

    forAll(displacement, pointi)
    {
        if (!pointWallDist[pointi].valid(dummyTrackData))
        {
            displacement[pointi] = Zero;
        }
        else
        {
            // 1) displacement on nearest wall point, scaled by medialRatio
            //    (wall distance / medial distance)
            // 2) pointWallDist[pointi].s() is layer thickness transported
            //    from closest wall point.
            // 3) shrink in opposite direction of addedPoints
            displacement[pointi] =
                -medialRatio[pointi]
                *pointWallDist[pointi].s()
                *dispVec[pointi];
        }
    }



    // Smear displacement away from fixed values (medialRatio=0 or 1)
    if (nSmoothDisplacement > 0)
    {
        scalarField invSumWeight(mesh.nPoints());
        sumWeights
        (
            isMasterEdge,
            identity(mesh.nEdges()),
            identity(mesh.nPoints()),
            mesh.edges(),
            invSumWeight
        );


        // Get smoothly varying patch field.
        Info<< "shrinkMeshDistance : Smoothing displacement ..." << endl;

        const scalar lambda = 0.33;
        const scalar mu = -0.34;

        pointField average(mesh.nPoints());
        for (label iter = 0; iter < nSmoothDisplacement; iter++)
        {
            // Calculate average of field
            averageNeighbours
            (
                mesh,
                isMasterEdge,
                identity(mesh.nEdges()),    // meshEdges,
                identity(mesh.nPoints()),   // meshPoints,
                mesh.edges(),               // edges,
                invSumWeight,
                displacement,
                average
            );

            forAll(displacement, i)
            {
                if (medialRatio[i] > small && medialRatio[i] < 1-small)
                {
                    displacement[i] =
                        (1-lambda)*displacement[i]
                       +lambda*average[i];
                }
            }


            // Calculate average of field
            averageNeighbours
            (
                mesh,
                isMasterEdge,
                identity(mesh.nEdges()),    // meshEdges,
                identity(mesh.nPoints()),   // meshPoints,
                mesh.edges(),               // edges,
                invSumWeight,
                displacement,
                average
            );

            forAll(displacement, i)
            {
                if (medialRatio[i] > small && medialRatio[i] < 1-small)
                {
                    displacement[i] = (1-mu)*displacement[i]+mu*average[i];
                }
            }


            // Do residual calculation every so often.
            if ((iter % 10) == 0)
            {
                scalar resid = meshRefinement::gAverage
                (
                    mesh,
                    syncTools::getMasterPoints(mesh),
                    mag(displacement-average)()
                );
                Info<< "    Iteration " << iter << "   residual " << resid
                    << endl;
            }
        }
    }

    // Make sure displacement boundary conditions is uptodate with
    // internal field
    meshMover.setDisplacementPatchFields();


    // Check a bit the sync of displacements
    if (debug&meshRefinement::MESH || debug&meshRefinement::LAYERINFO)
    {
        // initial mesh
        meshRefinement::testSyncPointList("mesh.points()", mesh, mesh.points());

        // pointWallDist
        scalarField pWallDist(pointWallDist.size());
        forAll(pointWallDist, pointi)
        {
            pWallDist[pointi] = pointWallDist[pointi].s();
        }
        meshRefinement::testSyncPointList("pointWallDist", mesh, pWallDist);

        // dispVec
        meshRefinement::testSyncPointList("dispVec", mesh, dispVec);

        // displacement before and after correction
        meshRefinement::testSyncPointList
        (
            "displacement BEFORE",
            mesh,
            displacement
        );

        meshMover.correctBoundaryConditions(displacement);
        meshRefinement::testSyncPointList
        (
            "displacement AFTER",
            mesh,
            displacement
        );
    }


    if (debug&meshRefinement::MESH || debug&meshRefinement::LAYERINFO)
    {
        const_cast<Time&>(mesh.time())++;
        Info<< "Writing wanted-displacement mesh (possibly illegal) to "
            << meshRefiner_.timeName() << endl;
        pointField oldPoints(mesh.points());

        meshRefiner_.mesh().movePoints(meshMover.curPoints());
        // Warn meshMover for changed geometry
        meshMover.movePoints();

        // Above move will have changed the instance only on the points (which
        // is correct).
        // However the previous mesh written will be the mesh with layers
        // (see snappyLayerDriver.C) so we now have to force writing all files
        // so we can easily step through time steps. Note that if you
        // don't write the mesh with layers this is not necessary.
        meshRefiner_.mesh().setInstance(meshRefiner_.timeName());

        meshRefiner_.write
        (
            meshRefinement::debugType(debug),
            meshRefinement::writeType
            (
                meshRefinement::writeLevel()
              | meshRefinement::WRITEMESH
            ),
            mesh.time().path()/meshRefiner_.timeName()
        );
        dispVec.write();
        medialDist.write();
        medialRatio.write();

        // Move mesh back
        meshRefiner_.mesh().movePoints(oldPoints);
        // Warn meshMover for changed geometry
        meshMover.movePoints();
    }


    // Current faces to check. Gets modified in meshMover.scaleMesh
    labelList checkFaces(identity(mesh.nFaces()));

    Info<< "shrinkMeshMedialDistance : Moving mesh ..." << endl;
    scalar oldErrorReduction = -1;

    for (label iter = 0; iter < 2*nSnap ; iter++)
    {
        Info<< "Iteration " << iter << endl;
        if (iter == nSnap)
        {
            Info<< "Displacement scaling for error reduction set to 0." << endl;
            oldErrorReduction = meshMover.setErrorReduction(0.0);
        }

        if
        (
            meshMover.scaleMesh
            (
                checkFaces,
                baffles,
                meshMover.paramDict(),
                meshQualityDict,
                true,
                nAllowableErrors
            )
        )
        {
            Info<< "shrinkMeshMedialDistance : Successfully moved mesh" << endl;
            break;
        }
    }

    if (oldErrorReduction >= 0)
    {
        meshMover.setErrorReduction(oldErrorReduction);
    }

    Info<< "shrinkMeshMedialDistance : Finished moving mesh ..." << endl;
}


// ************************************************************************* //
