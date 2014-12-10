/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

#include "medialAxisMeshMover.H"
#include "addToRunTimeSelectionTable.H"
#include "pointFields.H"
#include "valuePointPatchFields.H"
#include "PointEdgeWave.H"
#include "meshRefinement.H"
#include "unitConversion.H"
#include "PatchTools.H"
#include "OBJstream.H"
#include "pointData.H"
#include "zeroFixedValuePointPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(medialAxisMeshMover, 0);

    addToRunTimeSelectionTable
    (
        externalDisplacementMeshMover,
        medialAxisMeshMover,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::medialAxisMeshMover::getFixedValueBCs
(
    const pointVectorField& fld
)
{
    DynamicList<label> adaptPatchIDs;
    forAll(fld.boundaryField(), patchI)
    {
        const pointPatchField<vector>& patchFld =
            fld.boundaryField()[patchI];

        if (isA<valuePointPatchField<vector> >(patchFld))
        {
            if (isA<zeroFixedValuePointPatchField<vector> >(patchFld))
            {
                // Special condition of fixed boundary condition. Does not
                // get adapted
            }
            else
            {
                adaptPatchIDs.append(patchI);
            }
        }
    }
    return adaptPatchIDs;
}


Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::medialAxisMeshMover::getPatch
(
    const polyMesh& mesh,
    const labelList& patchIDs
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Count faces.
    label nFaces = 0;

    forAll(patchIDs, i)
    {
        const polyPatch& pp = patches[patchIDs[i]];

        nFaces += pp.size();
    }

    // Collect faces.
    labelList addressing(nFaces);
    nFaces = 0;

    forAll(patchIDs, i)
    {
        const polyPatch& pp = patches[patchIDs[i]];

        label meshFaceI = pp.start();

        forAll(pp, i)
        {
            addressing[nFaces++] = meshFaceI++;
        }
    }

    return autoPtr<indirectPrimitivePatch>
    (
        new indirectPrimitivePatch
        (
            IndirectList<face>(mesh.faces(), addressing),
            mesh.points()
        )
    );
}


void Foam::medialAxisMeshMover::smoothPatchNormals
(
    const label nSmoothDisp,
    const PackedBoolList& isPatchMasterPoint,
    const PackedBoolList& isPatchMasterEdge,
    pointField& normals
) const
{
    const indirectPrimitivePatch& pp = adaptPatchPtr_();
    const edgeList& edges = pp.edges();
    const labelList& meshPoints = pp.meshPoints();

    // Get smoothly varying internal normals field.
    Info<< typeName << " : Smoothing normals ..." << endl;

    scalarField edgeWeights(edges.size());
    scalarField invSumWeight(meshPoints.size());
    meshRefinement::calculateEdgeWeights
    (
        mesh(),
        isPatchMasterEdge,
        meshPoints,
        edges,
        edgeWeights,
        invSumWeight
    );


    vectorField average;
    for (label iter = 0; iter < nSmoothDisp; iter++)
    {
        meshRefinement::weightedSum
        (
            mesh(),
            isPatchMasterEdge,
            meshPoints,
            edges,
            edgeWeights,
            normals,
            average
        );
        average *= invSumWeight;

        // Do residual calculation every so often.
        if ((iter % 10) == 0)
        {
            scalar resid = meshRefinement::gAverage
            (
                isPatchMasterPoint,
                mag(normals-average)()
            );
            Info<< "    Iteration " << iter << "   residual " << resid << endl;
        }

        // Transfer to normals vector field
        forAll(average, pointI)
        {
            // full smoothing neighbours + point value
            average[pointI] = 0.5*(normals[pointI]+average[pointI]);
            normals[pointI] = average[pointI];
            normals[pointI] /= mag(normals[pointI]) + VSMALL;
        }
    }
}


// Smooth normals in interior.
void Foam::medialAxisMeshMover::smoothNormals
(
    const label nSmoothDisp,
    const PackedBoolList& isMeshMasterPoint,
    const PackedBoolList& isMeshMasterEdge,
    const labelList& fixedPoints,
    pointVectorField& normals
) const
{
    // Get smoothly varying internal normals field.
    Info<< typeName
        << " : Smoothing normals in interior ..." << endl;

    const edgeList& edges = mesh().edges();

    // Points that do not change.
    PackedBoolList isFixedPoint(mesh().nPoints());

    // Internal points that are fixed
    forAll(fixedPoints, i)
    {
        label meshPointI = fixedPoints[i];
        isFixedPoint.set(meshPointI, 1);
    }

    // Make sure that points that are coupled to meshPoints but not on a patch
    // are fixed as well
    syncTools::syncPointList(mesh(), isFixedPoint, maxEqOp<unsigned int>(), 0);


    // Correspondence between local edges/points and mesh edges/points
    const labelList meshPoints(identity(mesh().nPoints()));

    // Calculate inverse sum of weights

    scalarField edgeWeights(mesh().nEdges());
    scalarField invSumWeight(meshPoints.size());
    meshRefinement::calculateEdgeWeights
    (
        mesh(),
        isMeshMasterEdge,
        meshPoints,
        edges,
        edgeWeights,
        invSumWeight
    );

    vectorField average;
    for (label iter = 0; iter < nSmoothDisp; iter++)
    {
        meshRefinement::weightedSum
        (
            mesh(),
            isMeshMasterEdge,
            meshPoints,
            edges,
            edgeWeights,
            normals,
            average
        );
        average *= invSumWeight;

        // Do residual calculation every so often.
        if ((iter % 10) == 0)
        {
            scalar resid = meshRefinement::gAverage
            (
                isMeshMasterPoint,
                mag(normals-average)()
            );
            Info<< "    Iteration " << iter << "   residual " << resid << endl;
        }


        // Transfer to normals vector field
        forAll(average, pointI)
        {
            if (isFixedPoint.get(pointI) == 0)
            {
                //full smoothing neighbours + point value
                average[pointI] = 0.5*(normals[pointI]+average[pointI]);
                normals[pointI] = average[pointI];
                normals[pointI] /= mag(normals[pointI]) + VSMALL;
            }
        }
    }
}


// Tries and find a medial axis point. Done by comparing vectors to nearest
// wall point for both vertices of edge.
bool Foam::medialAxisMeshMover::isMaxEdge
(
    const List<pointData>& pointWallDist,
    const label edgeI,
    const scalar minCos
) const
{
    const pointField& points = mesh().points();

    // Do not mark edges with one side on moving wall.

    const edge& e = mesh().edges()[edgeI];

    vector v0(points[e[0]] - pointWallDist[e[0]].origin());
    scalar magV0(mag(v0));

    if (magV0 < SMALL)
    {
        return false;
    }

    vector v1(points[e[1]] - pointWallDist[e[1]].origin());
    scalar magV1(mag(v1));

    if (magV1 < SMALL)
    {
        return false;
    }


    //- Detect based on vector to nearest point differing for both endpoints
    //v0 /= magV0;
    //v1 /= magV1;
    //
    //// Test angle.
    //if ((v0 & v1) < minCos)
    //{
    //    return true;
    //}
    //else
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


void Foam::medialAxisMeshMover::update(const dictionary& coeffDict)
{
    Info<< typeName
        << " : Calculating distance to Medial Axis ..." << endl;

    const pointField& points = mesh().points();

    const indirectPrimitivePatch& pp = adaptPatchPtr_();
    const labelList& meshPoints = pp.meshPoints();


    // Read a few parameters
    // ~~~~~~~~~~~~~~~~~~~~~

    //- Smooth surface normals
    const label nSmoothSurfaceNormals = readLabel
    (
        coeffDict.lookup("nSmoothSurfaceNormals")
    );

    //- When is medial axis
    word angleKey = "minMedialAxisAngle";
    if (!coeffDict.found(angleKey))
    {
        // Backwards compatibility
        angleKey = "minMedianAxisAngle";
    }
    scalar minMedialAxisAngleCos = Foam::cos
    (
        degToRad(readScalar(coeffDict.lookup(angleKey)))
    );

    //- Feature angle when to stop adding layers
    const scalar featureAngle = readScalar(coeffDict.lookup("featureAngle"));

    //- When to slip along wall
    const scalar slipFeatureAngle =
    (
        coeffDict.found("slipFeatureAngle")
      ? readScalar(coeffDict.lookup("slipFeatureAngle"))
      : 0.5*featureAngle
    );

    //- Smooth internal normals
    const label nSmoothNormals = readLabel
    (
        coeffDict.lookup("nSmoothNormals")
    );

    //- Number of edges walking out
    const label nMedialAxisIter = coeffDict.lookupOrDefault<label>
    (
        "nMedialAxisIter",
        mesh().globalData().nTotalPoints()
    );


    // Predetermine mesh edges
    // ~~~~~~~~~~~~~~~~~~~~~~~

    // Precalulate (mesh) master point/edge (only relevant for shared pts/edges)
    const PackedBoolList isMeshMasterPoint(syncTools::getMasterPoints(mesh()));
    const PackedBoolList isMeshMasterEdge(syncTools::getMasterEdges(mesh()));
    // Precalculate meshEdge per pp edge
    const labelList meshEdges
    (
        pp.meshEdges
        (
            mesh().edges(),
            mesh().pointEdges()
        )
    );

    // Precalulate (patch) master point/edge
    const PackedBoolList isPatchMasterPoint
    (
        meshRefinement::getMasterPoints
        (
            mesh(),
            meshPoints
        )
    );
    const PackedBoolList isPatchMasterEdge
    (
        meshRefinement::getMasterEdges
        (
            mesh(),
            meshEdges
        )
    );

    // Determine pointNormal
    // ~~~~~~~~~~~~~~~~~~~~~

    pointField pointNormals(PatchTools::pointNormals(mesh(), pp));

    // Smooth patch normal vectors
    smoothPatchNormals
    (
        nSmoothSurfaceNormals,
        isPatchMasterPoint,
        isPatchMasterEdge,
        pointNormals
    );


    // Calculate distance to pp points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Distance to wall
    List<pointData> pointWallDist(mesh().nPoints());

    // Dummy additional info for PointEdgeWave
    int dummyTrackData = 0;


    // 1. Calculate distance to points where displacement is specified.
    {
        // Seed data.
        List<pointData> wallInfo(meshPoints.size());

        forAll(meshPoints, patchPointI)
        {
            label pointI = meshPoints[patchPointI];
            wallInfo[patchPointI] = pointData
            (
                points[pointI],
                0.0,
                pointI,                       // passive scalar
                pointNormals[patchPointI]     // surface normals
            );
        }

        // Do all calculations
        List<pointData> edgeWallDist(mesh().nEdges());
        PointEdgeWave<pointData> wallDistCalc
        (
            mesh(),
            meshPoints,
            wallInfo,
            pointWallDist,
            edgeWallDist,
            0,   // max iterations
            dummyTrackData
        );
        wallDistCalc.iterate(nMedialAxisIter);

        label nUnvisit = returnReduce
        (
            wallDistCalc.getUnsetPoints(),
            sumOp<label>()
        );

        if (nUnvisit > 0)
        {
            if (nMedialAxisIter > 0)
            {
                Info<< typeName
                    << " : Limited walk to " << nMedialAxisIter
                    << " steps. Not visited " << nUnvisit
                    << " out of " << mesh().globalData().nTotalPoints()
                    << " points" << endl;
            }
            else
            {
                WarningIn
                (
                    "medialAxisMeshMover::"
                    "medialAxisSmoothingInfo(..)"
                )   << "Walking did not visit all points." << nl
                    << "    Did not visit " << nUnvisit
                    << " out of " << mesh().globalData().nTotalPoints()
                    << " points. This is not necessarily a problem" << nl
                    << "    and might be due to faceZones splitting of part"
                    << " of the domain." << nl << endl;
            }
        }
    }


    // 2. Find points with max distance and transport information back to
    //    wall.
    {
        List<pointData> pointMedialDist(mesh().nPoints());
        List<pointData> edgeMedialDist(mesh().nEdges());

        // Seed point data.
        DynamicList<pointData> maxInfo(meshPoints.size());
        DynamicList<label> maxPoints(meshPoints.size());

        // 1. Medial axis points

        const edgeList& edges = mesh().edges();

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
                forAll(e, ep)
                {
                    label pointI = e[ep];

                    if (!pointMedialDist[pointI].valid(dummyTrackData))
                    {
                        maxPoints.append(pointI);
                        maxInfo.append
                        (
                            pointData
                            (
                                points[pointI],
                                0.0,
                                pointI,         // passive data
                                vector::zero    // passive data
                            )
                        );
                        pointMedialDist[pointI] = maxInfo.last();
                    }
                }

            }
            else if (isMaxEdge(pointWallDist, edgeI, minMedialAxisAngleCos))
            {
                // Both end points of edge have very different nearest wall
                // point. Mark both points as medial axis points.

                // Approximate medial axis location on edge.
                //const point medialAxisPt = e.centre(points);
                vector eVec = e.vec(points);
                scalar eMag = mag(eVec);
                if (eMag > VSMALL)
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
                        label pointI = e[ep];

                        if (!pointMedialDist[pointI].valid(dummyTrackData))
                        {
                            maxPoints.append(pointI);
                            maxInfo.append
                            (
                                pointData
                                (
                                    medialAxisPt,   //points[pointI],
                                    magSqr(points[pointI]-medialAxisPt),//0.0,
                                    pointI,         // passive data
                                    vector::zero    // passive data
                                )
                            );
                            pointMedialDist[pointI] = maxInfo.last();
                        }
                    }
                }
            }
        }


        // 2. Seed non-adapt patches
        const polyBoundaryMesh& patches = mesh().boundaryMesh();

        labelHashSet adaptPatches(adaptPatchIDs_);


        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];
            const pointPatchVectorField& pvf =
                pointDisplacement().boundaryField()[patchI];

            if
            (
                !pp.coupled()
             && !isA<emptyPolyPatch>(pp)
             && !adaptPatches.found(patchI)
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
                    Info<< typeName
                        << " : Inserting all points on patch " << pp.name()
                        << endl;

                    forAll(meshPoints, i)
                    {
                        label pointI = meshPoints[i];
                        if (!pointMedialDist[pointI].valid(dummyTrackData))
                        {
                            maxPoints.append(pointI);
                            maxInfo.append
                            (
                                pointData
                                (
                                    points[pointI],
                                    0.0,
                                    pointI,         // passive data
                                    vector::zero    // passive data
                                )
                            );
                            pointMedialDist[pointI] = maxInfo.last();
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
                    Info<< typeName
                        << " : Inserting points on patch " << pp.name()
                        << " if angle to nearest layer patch > "
                        << slipFeatureAngle << " degrees." << endl;

                    scalar slipFeatureAngleCos = Foam::cos
                    (
                        degToRad(slipFeatureAngle)
                    );
                    pointField pointNormals
                    (
                        PatchTools::pointNormals(mesh(), pp)
                    );

                    forAll(meshPoints, i)
                    {
                        label pointI = meshPoints[i];

                        if
                        (
                            pointWallDist[pointI].valid(dummyTrackData)
                        && !pointMedialDist[pointI].valid(dummyTrackData)
                        )
                        {
                            // Check if angle not too large.
                            scalar cosAngle =
                            (
                               -pointWallDist[pointI].v()
                              & pointNormals[i]
                            );
                            if (cosAngle > slipFeatureAngleCos)
                            {
                                // Extrusion direction practically perpendicular
                                // to the patch. Disable movement at the patch.

                                maxPoints.append(pointI);
                                maxInfo.append
                                (
                                    pointData
                                    (
                                        points[pointI],
                                        0.0,
                                        pointI,         // passive data
                                        vector::zero    // passive data
                                    )
                                );
                                pointMedialDist[pointI] = maxInfo.last();
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
            mesh(),
            maxPoints,
            maxInfo,

            pointMedialDist,
            edgeMedialDist,
            0,   // max iterations
            dummyTrackData
        );
        medialDistCalc.iterate(2*nMedialAxisIter);


        // Extract medial axis distance as pointScalarField
        forAll(pointMedialDist, pointI)
        {
            if (pointMedialDist[pointI].valid(dummyTrackData))
            {
                medialDist_[pointI] = Foam::sqrt
                (
                    pointMedialDist[pointI].distSqr()
                );
                medialVec_[pointI] = pointMedialDist[pointI].origin();
            }
            else
            {
                // Unvisited. Do as if on medial axis so unmoving
                medialDist_[pointI] = 0.0;
                medialVec_[pointI] = point(1, 0, 0);
            }
        }
    }

    // Extract transported surface normals as pointVectorField
    forAll(dispVec_, i)
    {
        if (!pointWallDist[i].valid(dummyTrackData))
        {
            dispVec_[i] = vector(1, 0, 0);
        }
        else
        {
            dispVec_[i] = pointWallDist[i].v();
        }
    }

    // Smooth normal vectors. Do not change normals on pp.meshPoints
    smoothNormals
    (
        nSmoothNormals,
        isMeshMasterPoint,
        isMeshMasterEdge,
        meshPoints,
        dispVec_
    );

    // Calculate ratio point medial distance to point wall distance
    forAll(medialRatio_, pointI)
    {
        if (!pointWallDist[pointI].valid(dummyTrackData))
        {
            medialRatio_[pointI] = 0.0;
        }
        else
        {
            scalar wDist2 = pointWallDist[pointI].distSqr();
            scalar mDist = medialDist_[pointI];

            if (wDist2 < sqr(SMALL) && mDist < SMALL)
            //- Note: maybe less strict:
            //(
            //    wDist2 < sqr(meshRefiner_.mergeDistance())
            // && mDist < meshRefiner_.mergeDistance()
            //)
            {
                medialRatio_[pointI] = 0.0;
            }
            else
            {
                medialRatio_[pointI] = mDist / (Foam::sqrt(wDist2) + mDist);
            }
        }
    }


    if (debug)
    {
        Info<< typeName
            << " : Writing medial axis fields:" << nl
            << incrIndent
            << "ratio of medial distance to wall distance : "
            << medialRatio_.name() << nl
            << "distance to nearest medial axis           : "
            << medialDist_.name() << nl
            << "nearest medial axis location              : "
            << medialVec_.name() << nl
            << "normal at nearest wall                    : "
            << dispVec_.name() << nl
            << decrIndent << nl
            << endl;

        dispVec_.write();
        medialRatio_.write();
        medialDist_.write();
        medialVec_.write();
    }
}


bool Foam::medialAxisMeshMover::unmarkExtrusion
(
    const label patchPointI,
    pointField& patchDisp,
    List<autoLayerDriver::extrudeMode>& extrudeStatus
)
{
    if (extrudeStatus[patchPointI] == autoLayerDriver::EXTRUDE)
    {
        extrudeStatus[patchPointI] = autoLayerDriver::NOEXTRUDE;
        patchDisp[patchPointI] = vector::zero;
        return true;
    }
    else if (extrudeStatus[patchPointI] == autoLayerDriver::EXTRUDEREMOVE)
    {
        extrudeStatus[patchPointI] = autoLayerDriver::NOEXTRUDE;
        patchDisp[patchPointI] = vector::zero;
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::medialAxisMeshMover::syncPatchDisplacement
(
    const scalarField& minThickness,
    pointField& patchDisp,
    List<autoLayerDriver::extrudeMode>& extrudeStatus
) const
{
    const indirectPrimitivePatch& pp = adaptPatchPtr_();
    const labelList& meshPoints = pp.meshPoints();

    label nChangedTotal = 0;

    while (true)
    {
        label nChanged = 0;

        // Sync displacement (by taking min)
        syncTools::syncPointList
        (
            mesh(),
            meshPoints,
            patchDisp,
            minMagSqrEqOp<vector>(),
            point::rootMax           // null value
        );

        // Unmark if displacement too small
        forAll(patchDisp, i)
        {
            if (mag(patchDisp[i]) < minThickness[i])
            {
                if (unmarkExtrusion(i, patchDisp, extrudeStatus))
                {
                    nChanged++;
                }
            }
        }

        //labelList syncPatchNLayers(patchNLayers);
        //
        //syncTools::syncPointList
        //(
        //    mesh(),
        //    meshPoints,
        //    syncPatchNLayers,
        //    minEqOp<label>(),
        //    labelMax            // null value
        //);
        //
        //// Reset if differs
        //// 1. take max
        //forAll(syncPatchNLayers, i)
        //{
        //    if (syncPatchNLayers[i] != patchNLayers[i])
        //    {
        //        if
        //        (
        //            unmarkExtrusion
        //            (
        //                i,
        //                patchDisp,
        //                patchNLayers,
        //                extrudeStatus
        //            )
        //        )
        //        {
        //            nChanged++;
        //        }
        //    }
        //}
        //
        //syncTools::syncPointList
        //(
        //    mesh(),
        //    meshPoints,
        //    syncPatchNLayers,
        //    maxEqOp<label>(),
        //    labelMin            // null value
        //);
        //
        //// Reset if differs
        //// 2. take min
        //forAll(syncPatchNLayers, i)
        //{
        //    if (syncPatchNLayers[i] != patchNLayers[i])
        //    {
        //        if
        //        (
        //            unmarkExtrusion
        //            (
        //                i,
        //                patchDisp,
        //                patchNLayers,
        //                extrudeStatus
        //            )
        //        )
        //        {
        //            nChanged++;
        //        }
        //    }
        //}

        nChangedTotal += nChanged;

        if (!returnReduce(nChanged, sumOp<label>()))
        {
            break;
        }
    }

    //Info<< "Prevented extrusion on "
    //    << returnReduce(nChangedTotal, sumOp<label>())
    //    << " coupled patch points during syncPatchDisplacement." << endl;
}


void Foam::medialAxisMeshMover::minSmoothField
(
    const label nSmoothDisp,
    const PackedBoolList& isPatchMasterPoint,
    const PackedBoolList& isPatchMasterEdge,
    const scalarField& fieldMin,
    scalarField& field
) const
{
    const indirectPrimitivePatch& pp = adaptPatchPtr_();
    const edgeList& edges = pp.edges();
    const labelList& meshPoints = pp.meshPoints();

    scalarField edgeWeights(edges.size());
    scalarField invSumWeight(meshPoints.size());
    meshRefinement::calculateEdgeWeights
    (
        mesh(),
        isPatchMasterEdge,
        meshPoints,
        edges,
        edgeWeights,
        invSumWeight
    );

    // Get smoothly varying patch field.
    Info<< typeName << " : Smoothing field ..." << endl;

    for (label iter = 0; iter < nSmoothDisp; iter++)
    {
        scalarField average(pp.nPoints());
        meshRefinement::weightedSum
        (
            mesh(),
            isPatchMasterEdge,
            meshPoints,
            edges,
            edgeWeights,
            field,
            average
        );
        average *= invSumWeight;

        // Transfer to field
        forAll(field, pointI)
        {
            //full smoothing neighbours + point value
            average[pointI] = 0.5*(field[pointI]+average[pointI]);

            // perform monotonic smoothing
            if
            (
                average[pointI] < field[pointI]
             && average[pointI] >= fieldMin[pointI]
            )
            {
                field[pointI] = average[pointI];
            }
        }

        // Do residual calculation every so often.
        if ((iter % 10) == 0)
        {
            scalar resid = meshRefinement::gAverage
            (
                isPatchMasterPoint,
                mag(field-average)()
            );
            Info<< "    Iteration " << iter << "   residual " << resid << endl;
        }
    }
}

// Stop layer growth where mesh wraps around edge with a
// large feature angle
void Foam::medialAxisMeshMover::
handleFeatureAngleLayerTerminations
(
    const scalar minCos,
    const PackedBoolList& isPatchMasterPoint,
    const labelList& meshEdges,
    List<autoLayerDriver::extrudeMode>& extrudeStatus,
    pointField& patchDisp,
    label& nPointCounter
) const
{
    const indirectPrimitivePatch& pp = adaptPatchPtr_();

    // Mark faces that have all points extruded
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    boolList extrudedFaces(pp.size(), true);

    forAll(pp.localFaces(), faceI)
    {
        const face& f = pp.localFaces()[faceI];

        forAll(f, fp)
        {
            if (extrudeStatus[f[fp]] == autoLayerDriver::NOEXTRUDE)
            {
                extrudedFaces[faceI] = false;
                break;
            }
        }
    }



    //label nOldPointCounter = nPointCounter;

    // Detect situation where two featureedge-neighbouring faces are partly or
    // not extruded and the edge itself is extruded. In this case unmark the
    // edge for extrusion.


    List<List<point> > edgeFaceNormals(pp.nEdges());
    List<List<bool> > edgeFaceExtrude(pp.nEdges());

    const labelListList& edgeFaces = pp.edgeFaces();
    const vectorField& faceNormals = pp.faceNormals();

    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        edgeFaceNormals[edgeI].setSize(eFaces.size());
        edgeFaceExtrude[edgeI].setSize(eFaces.size());
        forAll(eFaces, i)
        {
            label faceI = eFaces[i];
            edgeFaceNormals[edgeI][i] = faceNormals[faceI];
            edgeFaceExtrude[edgeI][i] = extrudedFaces[faceI];
        }
    }

    syncTools::syncEdgeList
    (
        mesh(),
        meshEdges,
        edgeFaceNormals,
        globalMeshData::ListPlusEqOp<List<point> >(),   // combine operator
        List<point>()               // null value
    );

    syncTools::syncEdgeList
    (
        mesh(),
        meshEdges,
        edgeFaceExtrude,
        globalMeshData::ListPlusEqOp<List<bool> >(),    // combine operator
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
                extrudeStatus[v0] != autoLayerDriver::NOEXTRUDE
             || extrudeStatus[v1] != autoLayerDriver::NOEXTRUDE
            )
            {
                if (!eFaceExtrude[0] || !eFaceExtrude[1])
                {
                    const vector& n0 = eFaceNormals[0];
                    const vector& n1 = eFaceNormals[1];

                    if ((n0 & n1) < minCos)
                    {
                        if (unmarkExtrusion(v0, patchDisp, extrudeStatus))
                        {
                            if (isPatchMasterPoint[v0])
                            {
                                nPointCounter++;
                            }
                        }
                        if (unmarkExtrusion(v1, patchDisp, extrudeStatus))
                        {
                            if (isPatchMasterPoint[v1])
                            {
                                nPointCounter++;
                            }
                        }
                    }
                }
            }
        }
    }

    //Info<< "Added "
    //    << returnReduce(nPointCounter-nOldPointCounter, sumOp<label>())
    //    << " point not to extrude." << endl;
}


// Find isolated islands (points, edges and faces and layer terminations)
// in the layer mesh and stop any layer growth at these points.
void Foam::medialAxisMeshMover::findIsolatedRegions
(
    const scalar minCosLayerTermination,
    const bool detectExtrusionIsland,
    const PackedBoolList& isPatchMasterPoint,
    const PackedBoolList& isPatchMasterEdge,
    const labelList& meshEdges,
    const scalarField& minThickness,
    List<autoLayerDriver::extrudeMode>& extrudeStatus,
    pointField& patchDisp
) const
{
    const indirectPrimitivePatch& pp = adaptPatchPtr_();
    const labelListList& pointFaces = pp.pointFaces();
    const labelList& meshPoints = pp.meshPoints();

    Info<< typeName << " : Removing isolated regions ..." << endl;

    // Keep count of number of points unextruded
    label nPointCounter = 0;


    autoPtr<OBJstream> str;
    if (debug)
    {
        str.reset
        (
            new OBJstream
            (
                mesh().time().path()
              / "islandExcludePoints_"
              + mesh().time().timeName()
              + ".obj"
            )
        );
        Info<< typeName
            << " : Writing points surrounded by non-extruded points to "
            << str().name() << endl;
    }

    while (true)
    {
        // Stop layer growth where mesh wraps around edge with a
        // large feature angle
        handleFeatureAngleLayerTerminations
        (
            minCosLayerTermination,
            isPatchMasterPoint,
            meshEdges,

            extrudeStatus,
            patchDisp,
            nPointCounter
        );

        syncPatchDisplacement(minThickness, patchDisp, extrudeStatus);



        // Detect either:
        // - point where all surrounding points are not extruded
        //   (detectExtrusionIsland)
        // or
        // - point where all the faces surrounding it are not fully
        //   extruded

        boolList keptPoints(pp.nPoints(), false);

        if (detectExtrusionIsland)
        {
            // Do not extrude from point where all neighbouring
            // points are not grown
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            labelList islandPoint(pp.size(), -1);
            forAll(pp, faceI)
            {
                const face& f = pp.localFaces()[faceI];

                forAll(f, fp)
                {
                    if (extrudeStatus[f[fp]] != autoLayerDriver::NOEXTRUDE)
                    {
                        if (islandPoint[faceI] == -1)
                        {
                            // First point to extrude
                            islandPoint[faceI] = f[fp];
                        }
                        else if (islandPoint[faceI] != -2)
                        {
                            // Second or more point to extrude
                            islandPoint[faceI] = -2;
                        }
                    }
                }
            }

            // islandPoint:
            //  -1 : no point extruded on face
            //  -2 : >= 2 points extruded on face
            //  >=0: label of point extruded

            // Check all surrounding faces that I am the islandPoint
            forAll(pointFaces, patchPointI)
            {
                if (extrudeStatus[patchPointI] != autoLayerDriver::NOEXTRUDE)
                {
                    const labelList& pFaces = pointFaces[patchPointI];

                    forAll(pFaces, i)
                    {
                        label faceI = pFaces[i];
                        if (islandPoint[faceI] != patchPointI)
                        {
                            keptPoints[patchPointI] = true;
                            break;
                        }
                    }
                }
            }
        }
        else
        {
            // Do not extrude from point where all neighbouring
            // faces are not grown
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            boolList extrudedFaces(pp.size(), true);
            forAll(pp.localFaces(), faceI)
            {
                const face& f = pp.localFaces()[faceI];
                forAll(f, fp)
                {
                    if (extrudeStatus[f[fp]] == autoLayerDriver::NOEXTRUDE)
                    {
                        extrudedFaces[faceI] = false;
                        break;
                    }
                }
            }

            const labelListList& pointFaces = pp.pointFaces();

            forAll(keptPoints, patchPointI)
            {
                const labelList& pFaces = pointFaces[patchPointI];

                forAll(pFaces, i)
                {
                    label faceI = pFaces[i];
                    if (extrudedFaces[faceI])
                    {
                        keptPoints[patchPointI] = true;
                        break;
                    }
                }
            }
        }


        syncTools::syncPointList
        (
            mesh(),
            meshPoints,
            keptPoints,
            orEqOp<bool>(),
            false               // null value
        );

        label nChanged = 0;

        forAll(keptPoints, patchPointI)
        {
            if (!keptPoints[patchPointI])
            {
                if (unmarkExtrusion(patchPointI, patchDisp, extrudeStatus))
                {
                    nPointCounter++;
                    nChanged++;

                    if (str.valid())
                    {
                        str().write(pp.points()[meshPoints[patchPointI]]);
                    }
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
        if (isPatchMasterEdge[edgeI])
        {
            const edge& e = edges[edgeI];

            label v0 = e[0];
            label v1 = e[1];

            if (extrudeStatus[v1] != autoLayerDriver::NOEXTRUDE)
            {
                isolatedPoint[v0] += 1;
            }
            if (extrudeStatus[v0] != autoLayerDriver::NOEXTRUDE)
            {
                isolatedPoint[v1] += 1;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh(),
        meshPoints,
        isolatedPoint,
        plusEqOp<label>(),
        label(0)        // null value
    );

    // stop layer growth on isolated faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    forAll(pp, faceI)
    {
        const face& f = pp.localFaces()[faceI];
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
                if (extrudeStatus[f[fp]] == autoLayerDriver::NOEXTRUDE)
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
                            extrudeStatus
                        )
                    )
                    {
                        nPointCounter++;

                        if (str.valid())
                        {
                            str().write(pp.points()[meshPoints[f[fp]]]);
                        }
                    }
                }
            }
        }
    }

    reduce(nPointCounter, sumOp<label>());
    Info<< typeName
        << " : Number of isolated points extrusion stopped : "<< nPointCounter
        << endl;
}


void Foam::medialAxisMeshMover::smoothLambdaMuDisplacement
(
    const label nSmoothDisp,
    const PackedBoolList& isMeshMasterPoint,
    const PackedBoolList& isMeshMasterEdge,
    vectorField& displacement
) const
{
    const edgeList& edges = mesh().edges();

    // Correspondence between local edges/points and mesh edges/points
    const labelList meshPoints(identity(mesh().nPoints()));

    // Calculate inverse sum of weights
    scalarField edgeWeights(mesh().nEdges());
    scalarField invSumWeight(meshPoints.size());
    meshRefinement::calculateEdgeWeights
    (
        mesh(),
        isMeshMasterEdge,
        meshPoints,
        edges,
        edgeWeights,
        invSumWeight
    );

    // Get smoothly varying patch field.
    Info<< typeName << " : Smoothing displacement ..." << endl;

    const scalar lambda = 0.33;
    const scalar mu = -0.34;

    vectorField average;

    for (label iter = 0; iter < nSmoothDisp; iter++)
    {
        meshRefinement::weightedSum
        (
            mesh(),
            isMeshMasterEdge,
            meshPoints,
            edges,
            edgeWeights,
            displacement,
            average
        );
        average *= invSumWeight;

        forAll(displacement, i)
        {
            if (medialRatio_[i] > SMALL && medialRatio_[i] < 1-SMALL)
            {
                displacement[i] = (1-lambda)*displacement[i]+lambda*average[i];
            }
        }

        meshRefinement::weightedSum
        (
            mesh(),
            isMeshMasterEdge,
            meshPoints,
            edges,
            edgeWeights,
            displacement,
            average
        );
        average *= invSumWeight;


        forAll(displacement, i)
        {
            if (medialRatio_[i] > SMALL && medialRatio_[i] < 1-SMALL)
            {
                displacement[i] = (1-mu)*displacement[i]+mu*average[i];
            }
        }


        // Do residual calculation every so often.
        if ((iter % 10) == 0)
        {
            scalar resid = meshRefinement::gAverage
            (
                isMeshMasterPoint,
                mag(displacement-average)()
            );
            Info<< "    Iteration " << iter << "   residual " << resid << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::medialAxisMeshMover::medialAxisMeshMover
(
    const dictionary& dict,
    const List<labelPair>& baffles,
    pointVectorField& pointDisplacement
)
:
    externalDisplacementMeshMover(dict, baffles, pointDisplacement),
    adaptPatchIDs_(getFixedValueBCs(pointDisplacement)),
    adaptPatchPtr_(getPatch(mesh(), adaptPatchIDs_)),
    scale_
    (
        IOobject
        (
            "scale",
            pointDisplacement.time().timeName(),
            pointDisplacement.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedScalar("scale", dimless, 1.0)
    ),
    oldPoints_(mesh().points()),
    meshMover_
    (
        const_cast<polyMesh&>(mesh()),
        const_cast<pointMesh&>(pMesh()),
        adaptPatchPtr_(),
        pointDisplacement,
        scale_,
        oldPoints_,
        adaptPatchIDs_,
        dict
    ),
    dispVec_
    (
        IOobject
        (
            "dispVec",
            pointDisplacement.time().timeName(),
            pointDisplacement.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pMesh(),
        dimensionedVector("dispVec", dimLength, vector::zero)
    ),
    medialRatio_
    (
        IOobject
        (
            "medialRatio",
            pointDisplacement.time().timeName(),
            pointDisplacement.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pMesh(),
        dimensionedScalar("medialRatio", dimless, 0.0)
    ),
    medialDist_
    (
        IOobject
        (
            "pointMedialDist",
            pointDisplacement.time().timeName(),
            pointDisplacement.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pMesh(),
        dimensionedScalar("pointMedialDist", dimLength, 0.0)
    ),
    medialVec_
    (
        IOobject
        (
            "medialVec",
            pointDisplacement.time().timeName(),
            pointDisplacement.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pMesh(),
        dimensionedVector("medialVec", dimLength, vector::zero)
    )
{
    update(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::medialAxisMeshMover::~medialAxisMeshMover()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::medialAxisMeshMover::calculateDisplacement
(
    const dictionary& coeffDict,
    const scalarField& minThickness,
    List<autoLayerDriver::extrudeMode>& extrudeStatus,
    pointField& patchDisp
)
{
    Info<< typeName << " : Smoothing using Medial Axis ..." << endl;

    const indirectPrimitivePatch& pp = adaptPatchPtr_;
    const labelList& meshPoints = pp.meshPoints();


    // Read settings
    // ~~~~~~~~~~~~~

    //- (lambda-mu) smoothing of internal displacement
    const label nSmoothDisplacement =  coeffDict.lookupOrDefault
    (
        "nSmoothDisplacement",
        0
    );

    //- Layer thickness too big
    const scalar maxThicknessToMedialRatio  = readScalar
    (
        coeffDict.lookup("maxThicknessToMedialRatio")
    );

    //- Feature angle when to stop adding layers
    const scalar featureAngle = readScalar(coeffDict.lookup("featureAngle"));

    //- Stop layer growth where mesh wraps around sharp edge
    const scalar minCosLayerTermination = Foam::cos
    (
        degToRad(0.5*featureAngle)
    );

    //- Smoothing wanted patch thickness
    const label nSmoothPatchThickness = readLabel
    (
        coeffDict.lookup("nSmoothThickness")
    );

    //- Number of edges walking out
    const label nMedialAxisIter = coeffDict.lookupOrDefault<label>
    (
        "nMedialAxisIter",
        mesh().globalData().nTotalPoints()
    );

    //- Use strict extrusionIsland detection
    const Switch detectExtrusionIsland = coeffDict.lookupOrDefault<Switch>
    (
        "detectExtrusionIsland",
        false
    );


    // Precalulate master points/edge (only relevant for shared points/edges)
    const PackedBoolList isMeshMasterPoint(syncTools::getMasterPoints(mesh()));
    const PackedBoolList isMeshMasterEdge(syncTools::getMasterEdges(mesh()));
    // Precalculate meshEdge per pp edge
    const labelList meshEdges
    (
        pp.meshEdges
        (
            mesh().edges(),
            mesh().pointEdges()
        )
    );

    // Precalulate (patch) master point/edge
    const PackedBoolList isPatchMasterPoint
    (
        meshRefinement::getMasterPoints
        (
            mesh(),
            meshPoints
        )
    );
    const PackedBoolList isPatchMasterEdge
    (
        meshRefinement::getMasterEdges
        (
            mesh(),
            meshEdges
        )
    );


    scalarField thickness(patchDisp.size());

    thickness = mag(patchDisp);

    forAll(thickness, patchPointI)
    {
        if (extrudeStatus[patchPointI] == autoLayerDriver::NOEXTRUDE)
        {
            thickness[patchPointI] = 0.0;
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
                mesh().time().path()
              / "thicknessRatioExcludePoints_"
              + mesh().time().timeName()
              + ".obj"
            )
        );
        Info<< typeName
            << " : Writing points with too large an extrusion distance to "
            << str().name() << endl;
    }

    autoPtr<OBJstream> medialVecStr;
    if (debug)
    {
        medialVecStr.reset
        (
            new OBJstream
            (
                mesh().time().path()
              / "thicknessRatioExcludeMedialVec_"
              + mesh().time().timeName()
              + ".obj"
            )
        );
        Info<< typeName
            << " : Writing medial axis vectors on points with too large"
            << " an extrusion distance to " << medialVecStr().name() << endl;
    }

    forAll(meshPoints, patchPointI)
    {
        if (extrudeStatus[patchPointI] != autoLayerDriver::NOEXTRUDE)
        {
            label pointI = meshPoints[patchPointI];

            //- Option 1: look only at extrusion thickness v.s. distance
            //  to nearest (medial axis or static) point.
            scalar mDist = medialDist_[pointI];
            scalar thicknessRatio = thickness[patchPointI]/(mDist+VSMALL);

            //- Option 2: Look at component in the direction
            //  of nearest (medial axis or static) point.
            vector n =
                patchDisp[patchPointI]
              / (mag(patchDisp[patchPointI]) + VSMALL);
            vector mVec = mesh().points()[pointI]-medialVec_[pointI];
            mVec /= mag(mVec)+VSMALL;
            thicknessRatio *= (n&mVec);

            if (thicknessRatio > maxThicknessToMedialRatio)
            {
                // Truncate thickness.
                if (debug&2)
                {
                    Pout<< "truncating displacement at "
                        << mesh().points()[pointI]
                        << " from " << thickness[patchPointI]
                        << " to "
                        <<  0.5
                           *(
                                minThickness[patchPointI]
                               +thickness[patchPointI]
                            )
                        << " medial direction:" << mVec
                        << " extrusion direction:" << n
                        << " with thicknessRatio:" << thicknessRatio
                        << endl;
                }

                thickness[patchPointI] =
                    0.5*(minThickness[patchPointI]+thickness[patchPointI]);

                patchDisp[patchPointI] = thickness[patchPointI]*n;

                if (isPatchMasterPoint[patchPointI])
                {
                    numThicknessRatioExclude++;
                }

                if (str.valid())
                {
                    const point& pt = mesh().points()[pointI];
                    str().write(linePointRef(pt, pt+patchDisp[patchPointI]));
                }
                if (medialVecStr.valid())
                {
                    const point& pt = mesh().points()[pointI];
                    medialVecStr().write
                    (
                        linePointRef
                        (
                            pt,
                            medialVec_[pointI]
                        )
                    );
                }
            }
        }
    }

    reduce(numThicknessRatioExclude, sumOp<label>());
    Info<< typeName << " : Reducing layer thickness at "
        << numThicknessRatioExclude
        << " nodes where thickness to medial axis distance is large " << endl;


    // find points where layer growth isolated to a lone point, edge or face

    findIsolatedRegions
    (
        minCosLayerTermination,
        detectExtrusionIsland,

        isPatchMasterPoint,
        isPatchMasterEdge,
        meshEdges,
        minThickness,

        extrudeStatus,
        patchDisp
    );

    // Update thickess for changed extrusion
    forAll(thickness, patchPointI)
    {
        if (extrudeStatus[patchPointI] == autoLayerDriver::NOEXTRUDE)
        {
            thickness[patchPointI] = 0.0;
        }
    }


    // smooth layer thickness on moving patch
    minSmoothField
    (
        nSmoothPatchThickness,
        isPatchMasterPoint,
        isPatchMasterEdge,
        minThickness,

        thickness
    );


    // Dummy additional info for PointEdgeWave
    int dummyTrackData = 0;

    List<pointData> pointWallDist(mesh().nPoints());

    const pointField& points = mesh().points();
    // 1. Calculate distance to points where displacement is specified.
    // This wave is used to transport layer thickness
    {
        // Distance to wall and medial axis on edges.
        List<pointData> edgeWallDist(mesh().nEdges());
        labelList wallPoints(meshPoints.size());

        // Seed data.
        List<pointData> wallInfo(meshPoints.size());

        forAll(meshPoints, patchPointI)
        {
            label pointI = meshPoints[patchPointI];
            wallPoints[patchPointI] = pointI;
            wallInfo[patchPointI] = pointData
            (
                points[pointI],
                0.0,
                thickness[patchPointI],       // transport layer thickness
                vector::zero                  // passive vector
            );
        }

        // Do all calculations
        PointEdgeWave<pointData> wallDistCalc
        (
            mesh(),
            wallPoints,
            wallInfo,
            pointWallDist,
            edgeWallDist,
            0,   // max iterations
            dummyTrackData
        );
        wallDistCalc.iterate(nMedialAxisIter);
    }


    // Calculate scaled displacement vector
    pointField& displacement = pointDisplacement_;

    forAll(displacement, pointI)
    {
        if (!pointWallDist[pointI].valid(dummyTrackData))
        {
            displacement[pointI] = vector::zero;
        }
        else
        {
            // 1) displacement on nearest wall point, scaled by medialRatio
            //    (wall distance / medial distance)
            // 2) pointWallDist[pointI].s() is layer thickness transported
            //    from closest wall point.
            // 3) shrink in opposite direction of addedPoints
            displacement[pointI] =
                -medialRatio_[pointI]
                *pointWallDist[pointI].s()
                *dispVec_[pointI];
        }
    }


    // Smear displacement away from fixed values (medialRatio=0 or 1)
    if (nSmoothDisplacement > 0)
    {
        smoothLambdaMuDisplacement
        (
            nSmoothDisplacement,
            isMeshMasterPoint,
            isMeshMasterEdge,
            displacement
        );
    }
}


bool Foam::medialAxisMeshMover::shrinkMesh
(
    const dictionary& meshQualityDict,
    const label nAllowableErrors,
    labelList& checkFaces
)
{
    //- Number of attempts shrinking the mesh
    const label nSnap  = readLabel(meshQualityDict.lookup("nRelaxIter"));




    // Make sure displacement boundary conditions is uptodate with
    // internal field
    meshMover_.setDisplacementPatchFields();

    Info<< typeName << " : Moving mesh ..." << endl;
    scalar oldErrorReduction = -1;

    bool meshOk = false;

    for (label iter = 0; iter < 2*nSnap ; iter++)
    {
        Info<< typeName
            << " : Iteration " << iter << endl;
        if (iter == nSnap)
        {
            Info<< typeName
                << " : Displacement scaling for error reduction set to 0."
                << endl;
            oldErrorReduction = meshMover_.setErrorReduction(0.0);
        }

        if
        (
            meshMover_.scaleMesh
            (
                checkFaces,
                baffles_,
                meshMover_.paramDict(),
                meshQualityDict,
                true,
                nAllowableErrors
            )
        )
        {
            Info<< typeName << " : Successfully moved mesh" << endl;
            meshOk = true;
            break;
        }
    }

    if (oldErrorReduction >= 0)
    {
        meshMover_.setErrorReduction(oldErrorReduction);
    }

    Info<< typeName << " : Finished moving mesh ..." << endl;

    return meshOk;
}


bool Foam::medialAxisMeshMover::move
(
    const dictionary& moveDict,
    const label nAllowableErrors,
    labelList& checkFaces
)
{
    // Read a few settings
    // ~~~~~~~~~~~~~~~~~~~

    //- Name of field specifying min thickness
    const word minThicknessName = word(moveDict.lookup("minThicknessName"));


    // The points have moved so before calculation update
    // the mesh and motionSolver accordingly
    movePoints(mesh().points());
    //
    //// Update any point motion bcs (e.g. timevarying)
    //pointDisplacement_.boundaryField().updateCoeffs();


    // Extract out patch-wise displacement
    const indirectPrimitivePatch& pp = adaptPatchPtr_();

    scalarField zeroMinThickness;
    if (minThicknessName == "none")
    {
        zeroMinThickness = scalarField(pp.nPoints(), 0.0);
    }
    const scalarField& minThickness =
    (
        (minThicknessName == "none")
      ? zeroMinThickness
      : mesh().lookupObject<scalarField>(minThicknessName)
    );


    pointField patchDisp
    (
        pointDisplacement_.internalField(),
        pp.meshPoints()
    );

    List<autoLayerDriver::extrudeMode> extrudeStatus
    (
        pp.nPoints(),
        autoLayerDriver::EXTRUDE
    );
    forAll(extrudeStatus, pointI)
    {
        if (mag(patchDisp[pointI]) <= minThickness[pointI]+SMALL)
        {
            extrudeStatus[pointI] = autoLayerDriver::NOEXTRUDE;
        }
    }


    // Solve displacement
    calculateDisplacement(moveDict, minThickness, extrudeStatus, patchDisp);

    //- Move mesh according to calculated displacement
    return shrinkMesh
    (
        moveDict,           // meshQualityDict,
        nAllowableErrors,   // nAllowableErrors
        checkFaces
    );
}


void Foam::medialAxisMeshMover::movePoints(const pointField& p)
{
    externalDisplacementMeshMover::movePoints(p);

    // Update local data for new geometry
    adaptPatchPtr_().movePoints(p);

    // Update motionSmoother for new geometry
    meshMover_.movePoints();

    // Assume corrent mesh location is correct
    meshMover_.correct();
}


// ************************************************************************* //
