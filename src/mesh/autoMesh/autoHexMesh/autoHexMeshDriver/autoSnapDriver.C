/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    All to do with snapping to the surface

\*----------------------------------------------------------------------------*/

#include "autoSnapDriver.H"
#include "motionSmoother.H"
#include "polyTopoChange.H"
#include "syncTools.H"
#include "fvMesh.H"
#include "Time.H"
#include "OFstream.H"
#include "OBJstream.H"
#include "mapPolyMesh.H"
#include "pointEdgePoint.H"
#include "PointEdgeWave.H"
#include "mergePoints.H"
#include "snapParameters.H"
#include "refinementSurfaces.H"
#include "unitConversion.H"
#include "localPointRegion.H"
#include "PatchTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(autoSnapDriver, 0);

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Calculate geometrically collocated points, Requires PackedList to be
// sized and initalised!
Foam::label Foam::autoSnapDriver::getCollocatedPoints
(
    const scalar tol,
    const pointField& points,
    PackedBoolList& isCollocatedPoint
)
{
    labelList pointMap;
    label nUnique = mergePoints
    (
        points,                         // points
        tol,                            // mergeTol
        false,                          // verbose
        pointMap
    );
    bool hasMerged = (nUnique < points.size());

    if (!returnReduce(hasMerged, orOp<bool>()))
    {
        return 0;
    }

    // Determine which merged points are referenced more than once
    label nCollocated = 0;

    // Per old point the newPoint. Or -1 (not set yet) or -2 (already seen
    // twice)
    labelList firstOldPoint(nUnique, -1);
    forAll(pointMap, oldPointI)
    {
        label newPointI = pointMap[oldPointI];

        if (firstOldPoint[newPointI] == -1)
        {
            // First use of oldPointI. Store.
            firstOldPoint[newPointI] = oldPointI;
        }
        else if (firstOldPoint[newPointI] == -2)
        {
            // Third or more reference of oldPointI -> non-manifold
            isCollocatedPoint.set(oldPointI, 1u);
            nCollocated++;
        }
        else
        {
            // Second reference of oldPointI -> non-manifold
            isCollocatedPoint.set(firstOldPoint[newPointI], 1u);
            nCollocated++;

            isCollocatedPoint.set(oldPointI, 1u);
            nCollocated++;

            // Mark with special value to save checking next time round
            firstOldPoint[newPointI] = -2;
        }
    }
    return returnReduce(nCollocated, sumOp<label>());
}


// Calculate displacement as average of patch points.
Foam::pointField Foam::autoSnapDriver::smoothPatchDisplacement
(
    const motionSmoother& meshMover,
    const List<labelPair>& baffles
)
{
    const indirectPrimitivePatch& pp = meshMover.patch();

    // Calculate geometrically non-manifold points on the patch to be moved.
    PackedBoolList nonManifoldPoint(pp.nPoints());
    label nNonManifoldPoints = getCollocatedPoints
    (
        SMALL,
        pp.localPoints(),
        nonManifoldPoint
    );
    Info<< "Found " << nNonManifoldPoints << " non-manifold point(s)."
        << endl;


    // Average points
    // ~~~~~~~~~~~~~~

    // We determine three points:
    // - average of (centres of) connected patch faces
    // - average of (centres of) connected internal mesh faces
    // - as fallback: centre of any connected cell
    // so we can do something moderately sensible for non/manifold points.

    // Note: the averages are calculated properly parallel. This is
    // necessary to get the points shared by processors correct.


    const labelListList& pointFaces = pp.pointFaces();
    const labelList& meshPoints = pp.meshPoints();
    const pointField& points = pp.points();
    const polyMesh& mesh = meshMover.mesh();

    // Get labels of faces to count (master of coupled faces and baffle pairs)
    PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh));

    {
        forAll(baffles, i)
        {
            label f0 = baffles[i].first();
            label f1 = baffles[i].second();

            if (isMasterFace.get(f0))
            {
                // Make f1 a slave
                isMasterFace.unset(f1);
            }
            else if (isMasterFace.get(f1))
            {
                isMasterFace.unset(f0);
            }
            else
            {
                FatalErrorIn("autoSnapDriver::smoothPatchDisplacement(..)")
                    << "Both sides of baffle consisting of faces " << f0
                    << " and " << f1 << " are already slave faces."
                    << abort(FatalError);
            }
        }
    }


    // Get average position of boundary face centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    vectorField avgBoundary(pointFaces.size(), vector::zero);
    labelList nBoundary(pointFaces.size(), 0);

    forAll(pointFaces, patchPointI)
    {
        const labelList& pFaces = pointFaces[patchPointI];

        forAll(pFaces, pfI)
        {
            label faceI = pFaces[pfI];

            if (isMasterFace.get(pp.addressing()[faceI]))
            {
                avgBoundary[patchPointI] += pp[faceI].centre(points);
                nBoundary[patchPointI]++;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        avgBoundary,
        plusEqOp<point>(),  // combine op
        vector::zero        // null value
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        nBoundary,
        plusEqOp<label>(),  // combine op
        label(0)            // null value
    );

    forAll(avgBoundary, i)
    {
        avgBoundary[i] /= nBoundary[i];
    }


    // Get average position of internal face centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    vectorField avgInternal;
    labelList nInternal;
    {
        vectorField globalSum(mesh.nPoints(), vector::zero);
        labelList globalNum(mesh.nPoints(), 0);

        // Note: no use of pointFaces
        const faceList& faces = mesh.faces();

        for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
        {
            const face& f = faces[faceI];
            const point& fc = mesh.faceCentres()[faceI];

            forAll(f, fp)
            {
                globalSum[f[fp]] += fc;
                globalNum[f[fp]]++;
            }
        }

        // Count coupled faces as internal ones (but only once)
        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        forAll(patches, patchI)
        {
            if
            (
                patches[patchI].coupled()
             && refCast<const coupledPolyPatch>(patches[patchI]).owner()
            )
            {
                const coupledPolyPatch& pp =
                    refCast<const coupledPolyPatch>(patches[patchI]);

                const vectorField::subField faceCentres = pp.faceCentres();

                forAll(pp, i)
                {
                    const face& f = pp[i];
                    const point& fc = faceCentres[i];

                    forAll(f, fp)
                    {
                        globalSum[f[fp]] += fc;
                        globalNum[f[fp]]++;
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            globalSum,
            plusEqOp<vector>(), // combine op
            vector::zero        // null value
        );
        syncTools::syncPointList
        (
            mesh,
            globalNum,
            plusEqOp<label>(),  // combine op
            label(0)            // null value
        );

        avgInternal.setSize(meshPoints.size());
        nInternal.setSize(meshPoints.size());

        forAll(avgInternal, patchPointI)
        {
            label meshPointI = meshPoints[patchPointI];

            nInternal[patchPointI] = globalNum[meshPointI];

            if (nInternal[patchPointI] == 0)
            {
                avgInternal[patchPointI] = globalSum[meshPointI];
            }
            else
            {
                avgInternal[patchPointI] =
                    globalSum[meshPointI]
                  / nInternal[patchPointI];
            }
        }
    }


    // Precalculate any cell using mesh point (replacement of pointCells()[])
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList anyCell(mesh.nPoints(), -1);
    forAll(mesh.faceNeighbour(), faceI)
    {
        label own = mesh.faceOwner()[faceI];
        const face& f = mesh.faces()[faceI];

        forAll(f, fp)
        {
            anyCell[f[fp]] = own;
        }
    }
    for (label faceI = mesh.nInternalFaces(); faceI < mesh.nFaces(); faceI++)
    {
        label own = mesh.faceOwner()[faceI];

        const face& f = mesh.faces()[faceI];

        forAll(f, fp)
        {
            anyCell[f[fp]] = own;
        }
    }


    // Displacement to calculate.
    pointField patchDisp(meshPoints.size(), vector::zero);

    forAll(pointFaces, i)
    {
        label meshPointI = meshPoints[i];
        const point& currentPos = pp.points()[meshPointI];

        // Now we have the two average points: avgBoundary and avgInternal
        // and how many boundary/internal faces connect to the point
        // (nBoundary, nInternal)
        // Do some blending between the two.
        // Note: the following section has some reasoning behind it but the
        // blending factors can be experimented with.

        point newPos;

        if (!nonManifoldPoint.get(i))
        {
            // Points that are manifold. Weight the internal and boundary
            // by their number of faces and blend with
            scalar internalBlend = 0.1;
            scalar blend = 0.1;

            point avgPos =
                (
                   internalBlend*nInternal[i]*avgInternal[i]
                  +(1-internalBlend)*nBoundary[i]*avgBoundary[i]
                )
              / (internalBlend*nInternal[i]+(1-internalBlend)*nBoundary[i]);

            newPos = (1-blend)*avgPos + blend*currentPos;
        }
        else if (nInternal[i] == 0)
        {
            // Non-manifold without internal faces. Use any connected cell
            // as internal point instead. Use precalculated any cell to avoid
            // e.g. pointCells()[meshPointI][0]

            const point& cc = mesh.cellCentres()[anyCell[meshPointI]];

            scalar cellCBlend = 0.8;
            scalar blend = 0.1;

            point avgPos = (1-cellCBlend)*avgBoundary[i] + cellCBlend*cc;

            newPos = (1-blend)*avgPos + blend*currentPos;
        }
        else
        {
            // Non-manifold point with internal faces connected to them
            scalar internalBlend = 0.9;
            scalar blend = 0.1;

            point avgPos =
                internalBlend*avgInternal[i]
              + (1-internalBlend)*avgBoundary[i];

            newPos = (1-blend)*avgPos + blend*currentPos;
        }

        patchDisp[i] = newPos - currentPos;
    }

    return patchDisp;
}
//XXXXXXX
//Foam::tmp<Foam::pointField> Foam::autoSnapDriver::avg
//(
//    const indirectPrimitivePatch& pp,
//    const pointField& localPoints
//)
//{
//    const labelListList& pointEdges = pp.pointEdges();
//    const edgeList& edges = pp.edges();
//
//    tmp<pointField> tavg(new pointField(pointEdges.size(), vector::zero));
//    pointField& avg = tavg();
//
//    forAll(pointEdges, vertI)
//    {
//        vector& avgPos = avg[vertI];
//
//        const labelList& pEdges = pointEdges[vertI];
//
//        forAll(pEdges, myEdgeI)
//        {
//            const edge& e = edges[pEdges[myEdgeI]];
//
//            label otherVertI = e.otherVertex(vertI);
//
//            avgPos += localPoints[otherVertI];
//        }
//
//        avgPos /= pEdges.size();
//    }
//    return tavg;
//}
//Foam::tmp<Foam::pointField>
//Foam::autoSnapDriver::smoothLambdaMuPatchDisplacement
//(
//    const motionSmoother& meshMover,
//    const List<labelPair>& baffles
//)
//{
//    const indirectPrimitivePatch& pp = meshMover.patch();
//    pointField newLocalPoints(pp.localPoints());
//
//    const label iters = 90;
//    const scalar lambda = 0.33;
//    const scalar mu = 0.34;
//
//    for (label iter = 0; iter < iters; iter++)
//    {
//        // Lambda
//        newLocalPoints =
//            (1 - lambda)*newLocalPoints
//          + lambda*avg(pp, newLocalPoints);
//
//        // Mu
//        newLocalPoints =
//            (1 + mu)*newLocalPoints
//          - mu*avg(pp, newLocalPoints);
//    }
//    return newLocalPoints-pp.localPoints();
//}
//XXXXXXX


Foam::tmp<Foam::scalarField> Foam::autoSnapDriver::edgePatchDist
(
    const pointMesh& pMesh,
    const indirectPrimitivePatch& pp
)
{
    const polyMesh& mesh = pMesh();

    // Set initial changed points to all the patch points
    List<pointEdgePoint> wallInfo(pp.nPoints());

    forAll(pp.localPoints(), ppI)
    {
        wallInfo[ppI] = pointEdgePoint(pp.localPoints()[ppI], 0.0);
    }

    // Current info on points
    List<pointEdgePoint> allPointInfo(mesh.nPoints());

    // Current info on edges
    List<pointEdgePoint> allEdgeInfo(mesh.nEdges());

    PointEdgeWave<pointEdgePoint> wallCalc
    (
        mesh,
        pp.meshPoints(),
        wallInfo,

        allPointInfo,
        allEdgeInfo,
        mesh.globalData().nTotalPoints()  // max iterations
    );

    // Copy edge values into scalarField
    tmp<scalarField> tedgeDist(new scalarField(mesh.nEdges()));
    scalarField& edgeDist = tedgeDist();

    forAll(allEdgeInfo, edgeI)
    {
        edgeDist[edgeI] = Foam::sqrt(allEdgeInfo[edgeI].distSqr());
    }


    //{
    //    // For debugging: dump to file
    //    pointScalarField pointDist
    //    (
    //        IOobject
    //        (
    //            "pointDist",
    //            meshRefiner_.timeName(),
    //            mesh.DB(),
    //            IOobject::NO_READ,
    //            IOobject::AUTO_WRITE
    //        ),
    //        pMesh,
    //        dimensionedScalar("pointDist", dimless, 0.0)
    //    );
    //
    //    forAll(allEdgeInfo, edgeI)
    //    {
    //        scalar d = Foam::sqrt(allEdgeInfo[edgeI].distSqr());
    //
    //        const edge& e = mesh.edges()[edgeI];
    //
    //        pointDist[e[0]] += d;
    //        pointDist[e[1]] += d;
    //    }
    //    forAll(pointDist, pointI)
    //    {
    //        pointDist[pointI] /= mesh.pointEdges()[pointI].size();
    //    }
    //    Info<< "Writing patch distance to " << pointDist.name()
    //        << " at time " << meshRefiner_.timeName() << endl;
    //
    //    pointDist.write();
    //}

    return tedgeDist;
}


void Foam::autoSnapDriver::dumpMove
(
    const fileName& fName,
    const pointField& meshPts,
    const pointField& surfPts
)
{
    // Dump direction of growth into file
    Info<< "Dumping move direction to " << fName << endl;

    OFstream nearestStream(fName);

    label vertI = 0;

    forAll(meshPts, ptI)
    {
        meshTools::writeOBJ(nearestStream, meshPts[ptI]);
        vertI++;

        meshTools::writeOBJ(nearestStream, surfPts[ptI]);
        vertI++;

        nearestStream<< "l " << vertI-1 << ' ' << vertI << nl;
    }
}


// Check whether all displacement vectors point outwards of patch. Return true
// if so.
bool Foam::autoSnapDriver::outwardsDisplacement
(
    const indirectPrimitivePatch& pp,
    const vectorField& patchDisp
)
{
    const vectorField& faceNormals = pp.faceNormals();
    const labelListList& pointFaces = pp.pointFaces();

    forAll(pointFaces, pointI)
    {
        const labelList& pFaces = pointFaces[pointI];

        vector disp(patchDisp[pointI]);

        scalar magDisp = mag(disp);

        if (magDisp > SMALL)
        {
            disp /= magDisp;

            bool outwards = meshTools::visNormal(disp, faceNormals, pFaces);

            if (!outwards)
            {
                Warning<< "Displacement " << patchDisp[pointI]
                    << " at mesh point " << pp.meshPoints()[pointI]
                    << " coord " << pp.points()[pp.meshPoints()[pointI]]
                    << " points through the surrounding patch faces" << endl;
                return false;
            }
        }
        else
        {
            //? Displacement small but in wrong direction. Would probably be ok.
        }
    }
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoSnapDriver::autoSnapDriver
(
    meshRefinement& meshRefiner,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch
)
:
    meshRefiner_(meshRefiner),
    globalToMasterPatch_(globalToMasterPatch),
    globalToSlavePatch_(globalToSlavePatch)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mapPolyMesh> Foam::autoSnapDriver::mergeZoneBaffles
(
    const List<labelPair>& baffles
)
{
    labelList zonedSurfaces =
        surfaceZonesInfo::getNamedSurfaces(meshRefiner_.surfaces().surfZones());

    autoPtr<mapPolyMesh> map;

    // No need to sync; all processors will have all same zonedSurfaces.
    label nBaffles = returnReduce(baffles.size(), sumOp<label>());
    if (zonedSurfaces.size() && nBaffles > 0)
    {
        // Merge any baffles
        Info<< "Converting " << nBaffles << " baffles back into zoned faces ..."
            << endl;

        map = meshRefiner_.mergeBaffles(baffles);

        Info<< "Converted baffles in = "
            << meshRefiner_.mesh().time().cpuTimeIncrement()
            << " s\n" << nl << endl;
    }

    return map;
}


Foam::scalarField Foam::autoSnapDriver::calcSnapDistance
(
    const fvMesh& mesh,
    const snapParameters& snapParams,
    const indirectPrimitivePatch& pp
)
{
    const edgeList& edges = pp.edges();
    const labelListList& pointEdges = pp.pointEdges();
    const pointField& localPoints = pp.localPoints();

    scalarField maxEdgeLen(localPoints.size(), -GREAT);

    forAll(pointEdges, pointI)
    {
        const labelList& pEdges = pointEdges[pointI];

        forAll(pEdges, pEdgeI)
        {
            const edge& e = edges[pEdges[pEdgeI]];

            scalar len = e.mag(localPoints);

            maxEdgeLen[pointI] = max(maxEdgeLen[pointI], len);
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        maxEdgeLen,
        maxEqOp<scalar>(),  // combine op
        -GREAT              // null value
    );

    return scalarField(snapParams.snapTol()*maxEdgeLen);
}


void Foam::autoSnapDriver::preSmoothPatch
(
    const meshRefinement& meshRefiner,
    const snapParameters& snapParams,
    const label nInitErrors,
    const List<labelPair>& baffles,
    motionSmoother& meshMover
)
{
    const fvMesh& mesh = meshRefiner.mesh();

    labelList checkFaces;

    Info<< "Smoothing patch points ..." << endl;
    for
    (
        label smoothIter = 0;
        smoothIter < snapParams.nSmoothPatch();
        smoothIter++
    )
    {
        Info<< "Smoothing iteration " << smoothIter << endl;
        checkFaces.setSize(mesh.nFaces());
        forAll(checkFaces, faceI)
        {
            checkFaces[faceI] = faceI;
        }

        pointField patchDisp(smoothPatchDisplacement(meshMover, baffles));
        //pointField patchDisp
        //(
        //  smoothLambdaMuPatchDisplacement(meshMover, baffles)
        //);

        // The current mesh is the starting mesh to smooth from.
        meshMover.setDisplacement(patchDisp);

        meshMover.correct();

        scalar oldErrorReduction = -1;

        for (label snapIter = 0; snapIter < 2*snapParams.nSnap(); snapIter++)
        {
            Info<< nl << "Scaling iteration " << snapIter << endl;

            if (snapIter == snapParams.nSnap())
            {
                Info<< "Displacement scaling for error reduction set to 0."
                    << endl;
                oldErrorReduction = meshMover.setErrorReduction(0.0);
            }

            // Try to adapt mesh to obtain displacement by smoothly
            // decreasing displacement at error locations.
            if (meshMover.scaleMesh(checkFaces, baffles, true, nInitErrors))
            {
                Info<< "Successfully moved mesh" << endl;
                break;
            }
        }

        if (oldErrorReduction >= 0)
        {
            meshMover.setErrorReduction(oldErrorReduction);
        }
        Info<< endl;
    }


    // The current mesh is the starting mesh to smooth from.
    meshMover.correct();

    if (debug&meshRefinement::MESH)
    {
        const_cast<Time&>(mesh.time())++;
        Info<< "Writing patch smoothed mesh to time "
            << meshRefiner.timeName() << '.' << endl;
        meshRefiner.write
        (
            meshRefinement::debugType(debug),
            meshRefinement::writeType
            (
                meshRefinement::writeLevel()
              | meshRefinement::WRITEMESH
            ),
            mesh.time().path()/meshRefiner.timeName()
        );
        Info<< "Dumped mesh in = "
            << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;
    }

    Info<< "Patch points smoothed in = "
        << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;
}


// Get (pp-local) indices of points that are both on zone and on patched surface
Foam::labelList Foam::autoSnapDriver::getZoneSurfacePoints
(
    const fvMesh& mesh,
    const indirectPrimitivePatch& pp,
    const word& zoneName
)
{
    label zoneI = mesh.faceZones().findZoneID(zoneName);

    if (zoneI == -1)
    {
        FatalErrorIn
        (
            "autoSnapDriver::getZoneSurfacePoints"
            "(const fvMesh&, const indirectPrimitivePatch&, const word&)"
        )   << "Cannot find zone " << zoneName
            << exit(FatalError);
    }

    const faceZone& fZone = mesh.faceZones()[zoneI];


    // Could use PrimitivePatch & localFaces to extract points but might just
    // as well do it ourselves.

    boolList pointOnZone(pp.nPoints(), false);

    forAll(fZone, i)
    {
        const face& f = mesh.faces()[fZone[i]];

        forAll(f, fp)
        {
            label meshPointI = f[fp];

            Map<label>::const_iterator iter =
                pp.meshPointMap().find(meshPointI);

            if (iter != pp.meshPointMap().end())
            {
                label pointI = iter();
                pointOnZone[pointI] = true;
            }
        }
    }

    return findIndices(pointOnZone, true);
}


Foam::tmp<Foam::pointField> Foam::autoSnapDriver::avgCellCentres
(
    const fvMesh& mesh,
    const indirectPrimitivePatch& pp
)
{
    const labelListList& pointFaces = pp.pointFaces();


    tmp<pointField> tavgBoundary
    (
        new pointField(pointFaces.size(), vector::zero)
    );
    pointField& avgBoundary = tavgBoundary();
    labelList nBoundary(pointFaces.size(), 0);

    forAll(pointFaces, pointI)
    {
        const labelList& pFaces = pointFaces[pointI];

        forAll(pFaces, pfI)
        {
            label faceI = pFaces[pfI];
            label meshFaceI = pp.addressing()[faceI];

            label own = mesh.faceOwner()[meshFaceI];
            avgBoundary[pointI] += mesh.cellCentres()[own];
            nBoundary[pointI]++;
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        avgBoundary,
        plusEqOp<point>(),  // combine op
        vector::zero        // null value
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        nBoundary,
        plusEqOp<label>(),  // combine op
        label(0)            // null value
    );

    forAll(avgBoundary, i)
    {
        avgBoundary[i] /= nBoundary[i];
    }
    return tavgBoundary;
}


//Foam::tmp<Foam::scalarField> Foam::autoSnapDriver::calcEdgeLen
//(
//    const indirectPrimitivePatch& pp
//) const
//{
//    // Get local edge length based on refinement level
//    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//    // (Ripped from autoLayerDriver)
//
//    tmp<scalarField> tedgeLen(new scalarField(pp.nPoints()));
//    scalarField& edgeLen = tedgeLen();
//    {
//        const fvMesh& mesh = meshRefiner_.mesh();
//        const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();
//        const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
//
//        labelList maxPointLevel(pp.nPoints(), labelMin);
//
//        forAll(pp, i)
//        {
//            label ownLevel = cellLevel[mesh.faceOwner()[pp.addressing()[i]]];
//            const face& f = pp.localFaces()[i];
//            forAll(f, fp)
//            {
//                maxPointLevel[f[fp]] = max(maxPointLevel[f[fp]], ownLevel);
//            }
//        }
//
//        syncTools::syncPointList
//        (
//            mesh,
//            pp.meshPoints(),
//            maxPointLevel,
//            maxEqOp<label>(),
//            labelMin            // null value
//        );
//
//
//        forAll(maxPointLevel, pointI)
//        {
//            // Find undistorted edge size for this level.
//            edgeLen[pointI] = edge0Len/(1<<maxPointLevel[pointI]);
//        }
//    }
//    return tedgeLen;
//}


void Foam::autoSnapDriver::detectNearSurfaces
(
    const scalar planarCos,
    const indirectPrimitivePatch& pp,
    const pointField& nearestPoint,
    const vectorField& nearestNormal,

    vectorField& disp
) const
{
    Info<< "Detecting near surfaces ..." << endl;

    const pointField& localPoints = pp.localPoints();
    const labelList& meshPoints = pp.meshPoints();
    const refinementSurfaces& surfaces = meshRefiner_.surfaces();
    const fvMesh& mesh = meshRefiner_.mesh();

    //// Get local edge length based on refinement level
    //const scalarField edgeLen(calcEdgeLen(pp));
    //
    //// Generate rays for every surface point
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    //{
    //    const scalar cos45 = Foam::cos(degToRad(45));
    //    vector n(cos45, cos45, cos45);
    //    n /= mag(n);
    //
    //    pointField start(14*pp.nPoints());
    //    pointField end(start.size());
    //
    //    label rayI = 0;
    //    forAll(localPoints, pointI)
    //    {
    //        const point& pt = localPoints[pointI];
    //
    //        // Along coordinate axes
    //
    //        {
    //            start[rayI] = pt;
    //            point& endPt = end[rayI++];
    //            endPt = pt;
    //            endPt.x() -= edgeLen[pointI];
    //        }
    //        {
    //            start[rayI] = pt;
    //            point& endPt = end[rayI++];
    //            endPt = pt;
    //            endPt.x() += edgeLen[pointI];
    //        }
    //        {
    //            start[rayI] = pt;
    //            point& endPt = end[rayI++];
    //            endPt = pt;
    //            endPt.y() -= edgeLen[pointI];
    //        }
    //        {
    //            start[rayI] = pt;
    //            point& endPt = end[rayI++];
    //            endPt = pt;
    //            endPt.y() += edgeLen[pointI];
    //        }
    //        {
    //            start[rayI] = pt;
    //            point& endPt = end[rayI++];
    //            endPt = pt;
    //            endPt.z() -= edgeLen[pointI];
    //        }
    //        {
    //            start[rayI] = pt;
    //            point& endPt = end[rayI++];
    //            endPt = pt;
    //            endPt.z() += edgeLen[pointI];
    //        }
    //
    //        // At 45 degrees
    //
    //        const vector vec(edgeLen[pointI]*n);
    //
    //        {
    //            start[rayI] = pt;
    //            point& endPt = end[rayI++];
    //            endPt = pt;
    //            endPt.x() += vec.x();
    //            endPt.y() += vec.y();
    //            endPt.z() += vec.z();
    //        }
    //        {
    //            start[rayI] = pt;
    //            point& endPt = end[rayI++];
    //            endPt = pt;
    //            endPt.x() -= vec.x();
    //            endPt.y() += vec.y();
    //            endPt.z() += vec.z();
    //        }
    //        {
    //            start[rayI] = pt;
    //            point& endPt = end[rayI++];
    //            endPt = pt;
    //            endPt.x() += vec.x();
    //            endPt.y() -= vec.y();
    //            endPt.z() += vec.z();
    //        }
    //        {
    //            start[rayI] = pt;
    //            point& endPt = end[rayI++];
    //            endPt = pt;
    //            endPt.x() -= vec.x();
    //            endPt.y() -= vec.y();
    //            endPt.z() += vec.z();
    //        }
    //        {
    //            start[rayI] = pt;
    //            point& endPt = end[rayI++];
    //            endPt = pt;
    //            endPt.x() += vec.x();
    //            endPt.y() += vec.y();
    //            endPt.z() -= vec.z();
    //        }
    //        {
    //            start[rayI] = pt;
    //            point& endPt = end[rayI++];
    //            endPt = pt;
    //            endPt.x() -= vec.x();
    //            endPt.y() += vec.y();
    //            endPt.z() -= vec.z();
    //        }
    //        {
    //            start[rayI] = pt;
    //            point& endPt = end[rayI++];
    //            endPt = pt;
    //            endPt.x() += vec.x();
    //            endPt.y() -= vec.y();
    //            endPt.z() -= vec.z();
    //        }
    //        {
    //            start[rayI] = pt;
    //            point& endPt = end[rayI++];
    //            endPt = pt;
    //            endPt.x() -= vec.x();
    //            endPt.y() -= vec.y();
    //            endPt.z() -= vec.z();
    //        }
    //    }
    //
    //    labelList surface1;
    //    List<pointIndexHit> hit1;
    //    labelList region1;
    //    vectorField normal1;
    //
    //    labelList surface2;
    //    List<pointIndexHit> hit2;
    //    labelList region2;
    //    vectorField normal2;
    //    surfaces.findNearestIntersection
    //    (
    //        unzonedSurfaces,    // surfacesToTest,
    //        start,
    //        end,
    //
    //        surface1,
    //        hit1,
    //        region1,
    //        normal1,
    //
    //        surface2,
    //        hit2,
    //        region2,
    //        normal2
    //    );
    //
    //    // All intersections
    //    {
    //        OBJstream str
    //        (
    //            mesh.time().path()
    //          / "surfaceHits_" + meshRefiner_.timeName() + ".obj"
    //        );
    //
    //        Info<< "Dumping intersections with rays to " << str.name()
    //            << endl;
    //
    //        forAll(hit1, i)
    //        {
    //            if (hit1[i].hit())
    //            {
    //                str.write(linePointRef(start[i], hit1[i].hitPoint()));
    //            }
    //            if (hit2[i].hit())
    //            {
    //                str.write(linePointRef(start[i], hit2[i].hitPoint()));
    //            }
    //        }
    //    }
    //
    //    // Co-planar intersections
    //    {
    //        OBJstream str
    //        (
    //            mesh.time().path()
    //          / "coplanarHits_" + meshRefiner_.timeName() + ".obj"
    //        );
    //
    //        Info<< "Dumping intersections with co-planar surfaces to "
    //            << str.name() << endl;
    //
    //        forAll(localPoints, pointI)
    //        {
    //            bool hasNormal = false;
    //            point surfPointA;
    //            vector surfNormalA;
    //            point surfPointB;
    //            vector surfNormalB;
    //
    //            bool isCoplanar = false;
    //
    //            label rayI = 14*pointI;
    //            for (label i = 0; i < 14; i++)
    //            {
    //                if (hit1[rayI].hit())
    //                {
    //                    const point& pt = hit1[rayI].hitPoint();
    //                    const vector& n = normal1[rayI];
    //
    //                    if (!hasNormal)
    //                    {
    //                        hasNormal = true;
    //                        surfPointA = pt;
    //                        surfNormalA = n;
    //                    }
    //                    else
    //                    {
    //                        if
    //                        (
    //                            meshRefiner_.isGap
    //                            (
    //                                planarCos,
    //                                surfPointA,
    //                                surfNormalA,
    //                                pt,
    //                                n
    //                            )
    //                        )
    //                        {
    //                            isCoplanar = true;
    //                            surfPointB = pt;
    //                            surfNormalB = n;
    //                            break;
    //                        }
    //                    }
    //                }
    //                if (hit2[rayI].hit())
    //                {
    //                    const point& pt = hit2[rayI].hitPoint();
    //                    const vector& n = normal2[rayI];
    //
    //                    if (!hasNormal)
    //                    {
    //                        hasNormal = true;
    //                        surfPointA = pt;
    //                        surfNormalA = n;
    //                    }
    //                    else
    //                    {
    //                        if
    //                        (
    //                            meshRefiner_.isGap
    //                            (
    //                                planarCos,
    //                                surfPointA,
    //                                surfNormalA,
    //                                pt,
    //                                n
    //                            )
    //                        )
    //                        {
    //                            isCoplanar = true;
    //                            surfPointB = pt;
    //                            surfNormalB = n;
    //                            break;
    //                        }
    //                    }
    //                }
    //
    //                rayI++;
    //            }
    //
    //            if (isCoplanar)
    //            {
    //                str.write(linePointRef(surfPointA, surfPointB));
    //            }
    //        }
    //    }
    //}


    const pointField avgCc(avgCellCentres(mesh, pp));

    // Construct rays through localPoints to beyond cell centre
    pointField start(pp.nPoints());
    pointField end(pp.nPoints());
    forAll(localPoints, pointI)
    {
        const point& pt = localPoints[pointI];
        const vector d = 2*(avgCc[pointI]-pt);
        start[pointI] = pt - d;
        end[pointI] = pt + d;
    }


    autoPtr<OBJstream> gapStr;
    if (debug&meshRefinement::ATTRACTION)
    {
        gapStr.reset
        (
            new OBJstream
            (
                mesh.time().path()
              / "detectNearSurfaces_" + meshRefiner_.timeName() + ".obj"
            )
        );
    }


    const PackedBoolList isPatchMasterPoint
    (
        meshRefinement::getMasterPoints
        (
            mesh,
            meshPoints
        )
    );

    label nOverride = 0;

    // 1. All points to non-interface surfaces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    {
        const labelList unzonedSurfaces =
            surfaceZonesInfo::getUnnamedSurfaces
            (
                meshRefiner_.surfaces().surfZones()
            );

        // Do intersection test
        labelList surface1;
        List<pointIndexHit> hit1;
        labelList region1;
        vectorField normal1;

        labelList surface2;
        List<pointIndexHit> hit2;
        labelList region2;
        vectorField normal2;
        surfaces.findNearestIntersection
        (
            unzonedSurfaces,
            start,
            end,

            surface1,
            hit1,
            region1,
            normal1,

            surface2,
            hit2,
            region2,
            normal2
        );


        forAll(localPoints, pointI)
        {
            // Current location
            const point& pt = localPoints[pointI];

            bool override = false;

            //if (hit1[pointI].hit())
            //{
            //    if
            //    (
            //        meshRefiner_.isGap
            //        (
            //            planarCos,
            //            nearestPoint[pointI],
            //            nearestNormal[pointI],
            //            hit1[pointI].hitPoint(),
            //            normal1[pointI]
            //        )
            //    )
            //    {
            //        disp[pointI] = hit1[pointI].hitPoint()-pt;
            //        override = true;
            //    }
            //}
            //if (hit2[pointI].hit())
            //{
            //    if
            //    (
            //        meshRefiner_.isGap
            //        (
            //            planarCos,
            //            nearestPoint[pointI],
            //            nearestNormal[pointI],
            //            hit2[pointI].hitPoint(),
            //            normal2[pointI]
            //        )
            //    )
            //    {
            //        disp[pointI] = hit2[pointI].hitPoint()-pt;
            //        override = true;
            //    }
            //}

            if (hit1[pointI].hit() && hit2[pointI].hit())
            {
                if
                (
                    meshRefiner_.isGap
                    (
                        planarCos,
                        hit1[pointI].hitPoint(),
                        normal1[pointI],
                        hit2[pointI].hitPoint(),
                        normal2[pointI]
                    )
                )
                {
                    // TBD: check if the attraction (to nearest) would attract
                    // good enough and not override attraction

                    if (gapStr.valid())
                    {
                        const point& intPt = hit2[pointI].hitPoint();
                        gapStr().write(linePointRef(pt, intPt));
                    }

                    // Choose hit2 : nearest to end point (so inside the domain)
                    disp[pointI] = hit2[pointI].hitPoint()-pt;
                    override = true;
                }
            }

            if (override && isPatchMasterPoint[pointI])
            {
                nOverride++;
            }
        }
    }


    // 2. All points on zones to their respective surface
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        // Surfaces with zone information
        const PtrList<surfaceZonesInfo>& surfZones = surfaces.surfZones();

        const labelList zonedSurfaces = surfaceZonesInfo::getNamedSurfaces
        (
            surfZones
        );

        forAll(zonedSurfaces, i)
        {
            label zoneSurfI = zonedSurfaces[i];

            const word& faceZoneName = surfZones[zoneSurfI].faceZoneName();

            const labelList surfacesToTest(1, zoneSurfI);

            // Get indices of points both on faceZone and on pp.
            labelList zonePointIndices
            (
                getZoneSurfacePoints
                (
                    mesh,
                    pp,
                    faceZoneName
                )
            );

            // Do intersection test
            labelList surface1;
            List<pointIndexHit> hit1;
            labelList region1;
            vectorField normal1;

            labelList surface2;
            List<pointIndexHit> hit2;
            labelList region2;
            vectorField normal2;
            surfaces.findNearestIntersection
            (
                surfacesToTest,
                pointField(start, zonePointIndices),
                pointField(end, zonePointIndices),

                surface1,
                hit1,
                region1,
                normal1,

                surface2,
                hit2,
                region2,
                normal2
            );


            forAll(hit1, i)
            {
                label pointI = zonePointIndices[i];

                // Current location
                const point& pt = localPoints[pointI];

                bool override = false;

                //if (hit1[i].hit())
                //{
                //    if
                //    (
                //        meshRefiner_.isGap
                //        (
                //            planarCos,
                //            nearestPoint[pointI],
                //            nearestNormal[pointI],
                //            hit1[i].hitPoint(),
                //            normal1[i]
                //        )
                //    )
                //    {
                //        disp[pointI] = hit1[i].hitPoint()-pt;
                //        override = true;
                //    }
                //}
                //if (hit2[i].hit())
                //{
                //    if
                //    (
                //        meshRefiner_.isGap
                //        (
                //            planarCos,
                //            nearestPoint[pointI],
                //            nearestNormal[pointI],
                //            hit2[i].hitPoint(),
                //            normal2[i]
                //        )
                //    )
                //    {
                //        disp[pointI] = hit2[i].hitPoint()-pt;
                //        override = true;
                //    }
                //}

                if (hit1[i].hit() && hit2[i].hit())
                {
                    if
                    (
                        meshRefiner_.isGap
                        (
                            planarCos,
                            hit1[i].hitPoint(),
                            normal1[i],
                            hit2[i].hitPoint(),
                            normal2[i]
                        )
                    )
                    {
                        if (gapStr.valid())
                        {
                            const point& intPt = hit2[i].hitPoint();
                            gapStr().write(linePointRef(pt, intPt));
                        }

                        disp[pointI] = hit2[i].hitPoint()-pt;
                        override = true;
                    }
                }

                if (override && isPatchMasterPoint[pointI])
                {
                    nOverride++;
                }
            }
        }
    }

    Info<< "Overriding nearest with intersection of close gaps at "
        << returnReduce(nOverride, sumOp<label>())
        << " out of " << returnReduce(pp.nPoints(), sumOp<label>())
        << " points." << endl;
}


Foam::vectorField Foam::autoSnapDriver::calcNearestSurface
(
    const meshRefinement& meshRefiner,
    const scalarField& snapDist,
    const indirectPrimitivePatch& pp,
    pointField& nearestPoint,
    vectorField& nearestNormal
)
{
    Info<< "Calculating patchDisplacement as distance to nearest surface"
        << " point ..." << endl;

    const pointField& localPoints = pp.localPoints();
    const refinementSurfaces& surfaces = meshRefiner.surfaces();
    const fvMesh& mesh = meshRefiner.mesh();

    // Displacement per patch point
    vectorField patchDisp(localPoints.size(), vector::zero);

    if (returnReduce(localPoints.size(), sumOp<label>()) > 0)
    {
        // Current surface snapped to
        labelList snapSurf(localPoints.size(), -1);

        // Divide surfaces into zoned and unzoned
        const labelList zonedSurfaces =
            surfaceZonesInfo::getNamedSurfaces
            (
                meshRefiner.surfaces().surfZones()
            );
        const labelList unzonedSurfaces =
            surfaceZonesInfo::getUnnamedSurfaces
            (
                meshRefiner.surfaces().surfZones()
            );


        // 1. All points to non-interface surfaces
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        {
            List<pointIndexHit> hitInfo;
            labelList hitSurface;

            if (nearestNormal.size() == localPoints.size())
            {
                labelList hitRegion;
                vectorField hitNormal;
                surfaces.findNearestRegion
                (
                    unzonedSurfaces,
                    localPoints,
                    sqr(snapDist),
                    hitSurface,
                    hitInfo,
                    hitRegion,
                    hitNormal
                );

                forAll(hitInfo, pointI)
                {
                    if (hitInfo[pointI].hit())
                    {
                        nearestPoint[pointI] = hitInfo[pointI].hitPoint();
                        nearestNormal[pointI] = hitNormal[pointI];
                    }
                }
            }
            else
            {
                surfaces.findNearest
                (
                    unzonedSurfaces,
                    localPoints,
                    sqr(snapDist),        // sqr of attract distance
                    hitSurface,
                    hitInfo
                );
            }

            forAll(hitInfo, pointI)
            {
                if (hitInfo[pointI].hit())
                {
                    patchDisp[pointI] =
                        hitInfo[pointI].hitPoint()
                      - localPoints[pointI];

                    snapSurf[pointI] = hitSurface[pointI];
                }
            }
        }



        // 2. All points on zones to their respective surface
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Surfaces with zone information
        const PtrList<surfaceZonesInfo>& surfZones = surfaces.surfZones();

        // Current best snap distance
        scalarField minSnapDist(snapDist);

        forAll(zonedSurfaces, i)
        {
            label zoneSurfI = zonedSurfaces[i];

            const word& faceZoneName = surfZones[zoneSurfI].faceZoneName();

            const labelList surfacesToTest(1, zoneSurfI);

            // Get indices of points both on faceZone and on pp.
            labelList zonePointIndices
            (
                getZoneSurfacePoints
                (
                    mesh,
                    pp,
                    faceZoneName
                )
            );

            // Find nearest for points both on faceZone and pp.
            List<pointIndexHit> hitInfo;
            labelList hitSurface;

            if (nearestNormal.size() == localPoints.size())
            {
                labelList hitRegion;
                vectorField hitNormal;
                surfaces.findNearestRegion
                (
                    surfacesToTest,
                    pointField(localPoints, zonePointIndices),
                    sqr(scalarField(minSnapDist, zonePointIndices)),
                    hitSurface,
                    hitInfo,
                    hitRegion,
                    hitNormal
                );

                forAll(hitInfo, i)
                {
                    if (hitInfo[i].hit())
                    {
                        label pointI = zonePointIndices[i];
                        nearestPoint[pointI] = hitInfo[i].hitPoint();
                        nearestNormal[pointI] = hitNormal[i];
                    }
                }
            }
            else
            {
                surfaces.findNearest
                (
                    surfacesToTest,
                    pointField(localPoints, zonePointIndices),
                    sqr(scalarField(minSnapDist, zonePointIndices)),
                    hitSurface,
                    hitInfo
                );
            }

            forAll(hitInfo, i)
            {
                label pointI = zonePointIndices[i];

                if (hitInfo[i].hit())
                {
                    patchDisp[pointI] =
                        hitInfo[i].hitPoint()
                      - localPoints[pointI];

                    minSnapDist[pointI] = min
                    (
                        minSnapDist[pointI],
                        mag(patchDisp[pointI])
                    );

                    snapSurf[pointI] = zoneSurfI;
                }
            }
        }

        // Check if all points are being snapped
        forAll(snapSurf, pointI)
        {
            if (snapSurf[pointI] == -1)
            {
                WarningIn("autoSnapDriver::calcNearestSurface(..)")
                    << "For point:" << pointI
                    << " coordinate:" << localPoints[pointI]
                    << " did not find any surface within:"
                    << minSnapDist[pointI]
                    << " metre." << endl;
            }
        }

        {
            const PackedBoolList isPatchMasterPoint
            (
                meshRefinement::getMasterPoints
                (
                    mesh,
                    pp.meshPoints()
                )
            );

            scalarField magDisp(mag(patchDisp));

            Info<< "Wanted displacement : average:"
                <<  meshRefinement::gAverage(isPatchMasterPoint, magDisp)
                << " min:" << gMin(magDisp)
                << " max:" << gMax(magDisp) << endl;
        }
    }

    Info<< "Calculated surface displacement in = "
        << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;


    // Limit amount of movement.
    forAll(patchDisp, patchPointI)
    {
        scalar magDisp = mag(patchDisp[patchPointI]);

        if (magDisp > snapDist[patchPointI])
        {
            patchDisp[patchPointI] *= snapDist[patchPointI] / magDisp;

            Pout<< "Limiting displacement for " << patchPointI
                << " from " << magDisp << " to " << snapDist[patchPointI]
                << endl;
        }
    }

    // Points on zones in one domain but only present as point on other
    // will not do condition 2 on all. Sync explicitly.
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        patchDisp,
        minMagSqrEqOp<point>(),         // combine op
        vector(GREAT, GREAT, GREAT)     // null value (note: cannot use VGREAT)
    );

    return patchDisp;
}


////XXXXXXXXX
//// Get (pp-local) indices of points that are on both patches
//Foam::labelList Foam::autoSnapDriver::getPatchSurfacePoints
//(
//    const fvMesh& mesh,
//    const indirectPrimitivePatch& allPp,
//    const polyPatch& pp
//)
//{
//    // Could use PrimitivePatch & localFaces to extract points but might just
//    // as well do it ourselves.
//
//    boolList pointOnZone(allPp.nPoints(), false);
//
//    forAll(pp, i)
//    {
//        const face& f = pp[i];
//
//        forAll(f, fp)
//        {
//            label meshPointI = f[fp];
//
//            Map<label>::const_iterator iter =
//                allPp.meshPointMap().find(meshPointI);
//
//            if (iter != allPp.meshPointMap().end())
//            {
//                label pointI = iter();
//                pointOnZone[pointI] = true;
//            }
//        }
//    }
//
//    return findIndices(pointOnZone, true);
//}
//Foam::vectorField Foam::autoSnapDriver::calcNearestLocalSurface
//(
//    const meshRefinement& meshRefiner,
//    const scalarField& snapDist,
//    const indirectPrimitivePatch& pp
//)
//{
//    Info<< "Calculating patchDisplacement as distance to nearest"
//        << " local surface point ..." << endl;
//
//    const pointField& localPoints = pp.localPoints();
//    const refinementSurfaces& surfaces = meshRefiner.surfaces();
//    const fvMesh& mesh = meshRefiner.mesh();
//    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
//
//
////    // Assume that all patch-internal points get attracted to their surface
////    // only. So we want to know if point is on multiple regions
////
////    labelList minPatch(mesh.nPoints(), labelMax);
////    labelList maxPatch(mesh.nPoints(), labelMin);
////
////    forAll(meshMover.adaptPatchIDs(), i)
////    {
////        label patchI = meshMover.adaptPatchIDs()[i];
////        const labelList& meshPoints = pbm[patchI].meshPoints();
////
////        forAll(meshPoints, meshPointI)
////        {
////            label meshPointI = meshPoints[meshPointI];
////            minPatch[meshPointI] = min(minPatch[meshPointI], patchI);
////            maxPatch[meshPointI] = max(maxPatch[meshPointI], patchI);
////        }
////    }
////
////    syncTools::syncPointList
////    (
////        mesh,
////        minPatch,
////        minEqOp<label>(),   // combine op
////        labelMax            // null value
////    );
////    syncTools::syncPointList
////    (
////        mesh,
////        maxPatch,
////        maxEqOp<label>(),   // combine op
////        labelMin            // null value
////    );
//
//    // Now all points with minPatch != maxPatch will be on the outside of
//    // the patch.
//
//    // Displacement per patch point
//    vectorField patchDisp(localPoints.size(), vector::zero);
//    // Current best snap distance
//    scalarField minSnapDist(snapDist);
//    // Current surface snapped to
//    labelList snapSurf(localPoints.size(), -1);
//
//    const labelList& surfaceGeometry = surfaces.surfaces();
//    forAll(surfaceGeometry, surfI)
//    {
//        label geomI = surfaceGeometry[surfI];
//        const wordList& regNames = allGeometry.regionNames()[geomI];
//        forAll(regNames, regionI)
//        {
//            label globalRegionI = surfaces.globalRegion(surfI, regionI);
//            // Collect master patch points
//            label masterPatchI = globalToMasterPatch_[globalRegionI];
//            label slavePatchI = globalToSlavePatch_[globalRegionI];
//
//            labelList patchPointIndices
//            (
//                getPatchSurfacePoints
//                (
//                    mesh,
//                    pp,
//                    pbm[masterPatchI]
//                )
//            );
//
//            // Find nearest for points both on faceZone and pp.
//            List<pointIndexHit> hitInfo;
//            surfaces.findNearest
//            (
//                surfI,
//                regionI,
//                pointField(localPoints, patchPointIndices),
//                sqr(scalarField(minSnapDist, patchPointIndices)),
//                hitInfo
//            );
//
//            forAll(hitInfo, i)
//            {
//                label pointI = patchPointIndices[i];
//
//                if (hitInfo[i].hit())
//                {
//                    const point& pt = hitInfo[i].hitPoint();
//                    patchDisp[pointI] = pt-localPoints[pointI];
//                    minSnapDist[pointI] = min
//                    (
//                        minSnapDist[pointI],
//                        mag(patchDisp[pointI])
//                    );
//                    snapSurf[pointI] = surfI;
//                }
//            }
//
//            // Slave patch
//            if (slavePatchI != masterPatchI)
//            {
//                labelList patchPointIndices
//                (
//                    getPatchSurfacePoints
//                    (
//                        mesh,
//                        pp,
//                        pbm[slavePatchI]
//                    )
//                );
//
//                // Find nearest for points both on faceZone and pp.
//                List<pointIndexHit> hitInfo;
//                surfaces.findNearest
//                (
//                    surfI,
//                    regionI,
//                    pointField(localPoints, patchPointIndices),
//                    sqr(scalarField(minSnapDist, patchPointIndices)),
//                    hitInfo
//                );
//
//                forAll(hitInfo, i)
//                {
//                    label pointI = patchPointIndices[i];
//
//                    if (hitInfo[i].hit())
//                    {
//                        const point& pt = hitInfo[i].hitPoint();
//                        patchDisp[pointI] = pt-localPoints[pointI];
//                        minSnapDist[pointI] = min
//                        (
//                            minSnapDist[pointI],
//                            mag(patchDisp[pointI])
//                        );
//                        snapSurf[pointI] = surfI;
//                    }
//                }
//            }
//        }
//    }
//
//
//    // Check if all points are being snapped
//    forAll(snapSurf, pointI)
//    {
//        if (snapSurf[pointI] == -1)
//        {
//            WarningIn("autoSnapDriver::calcNearestLocalSurface(..)")
//                << "For point:" << pointI
//                << " coordinate:" << localPoints[pointI]
//                << " did not find any surface within:"
//                << minSnapDist[pointI]
//                << " metre." << endl;
//        }
//    }
//
//    {
//        scalarField magDisp(mag(patchDisp));
//
//        Info<< "Wanted displacement : average:"
//            << gSum(magDisp)/returnReduce(patchDisp.size(), sumOp<label>())
//            << " min:" << gMin(magDisp)
//            << " max:" << gMax(magDisp) << endl;
//    }
//
//
//    Info<< "Calculated surface displacement in = "
//        << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;
//
//
//    // Limit amount of movement.
//    forAll(patchDisp, patchPointI)
//    {
//        scalar magDisp = mag(patchDisp[patchPointI]);
//
//        if (magDisp > snapDist[patchPointI])
//        {
//            patchDisp[patchPointI] *= snapDist[patchPointI] / magDisp;
//
//            Pout<< "Limiting displacement for " << patchPointI
//                << " from " << magDisp << " to " << snapDist[patchPointI]
//                << endl;
//        }
//    }
//
//    // Points on zones in one domain but only present as point on other
//    // will not do condition 2 on all. Sync explicitly.
//    syncTools::syncPointList
//    (
//        mesh,
//        pp.meshPoints(),
//        patchDisp,
//        minMagSqrEqOp<point>(),    // combine op
//        vector(GREAT, GREAT, GREAT)// null value (note: cannot use VGREAT)
//    );
//
//    return patchDisp;
//}
////XXXXXXXXX

void Foam::autoSnapDriver::smoothDisplacement
(
    const snapParameters& snapParams,
    motionSmoother& meshMover
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const indirectPrimitivePatch& pp = meshMover.patch();

    Info<< "Smoothing displacement ..." << endl;

    // Set edge diffusivity as inverse of distance to patch
    scalarField edgeGamma(1.0/(edgePatchDist(meshMover.pMesh(), pp) + SMALL));
    //scalarField edgeGamma(mesh.nEdges(), 1.0);
    //scalarField edgeGamma(wallGamma(mesh, pp, 10, 1));

    // Get displacement field
    pointVectorField& disp = meshMover.displacement();

    for (label iter = 0; iter < snapParams.nSmoothDispl(); iter++)
    {
        if ((iter % 10) == 0)
        {
            Info<< "Iteration " << iter << endl;
        }
        pointVectorField oldDisp(disp);
        meshMover.smooth(oldDisp, edgeGamma, disp);
    }
    Info<< "Displacement smoothed in = "
        << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;

    if (debug&meshRefinement::MESH)
    {
        const_cast<Time&>(mesh.time())++;
        Info<< "Writing smoothed mesh to time " << meshRefiner_.timeName()
            << endl;

        // Moving mesh creates meshPhi. Can be cleared out by a mesh.clearOut
        // but this will also delete all pointMesh but not pointFields which
        // gives an illegal situation.

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
        Info<< "Writing displacement field ..." << endl;
        disp.write();
        tmp<pointScalarField> magDisp(mag(disp));
        magDisp().write();

        Info<< "Writing actual patch displacement ..." << endl;
        vectorField actualPatchDisp(disp, pp.meshPoints());
        dumpMove
        (
            mesh.time().path()
          / "actualPatchDisplacement_" + meshRefiner_.timeName() + ".obj",
            pp.localPoints(),
            pp.localPoints() + actualPatchDisp
        );
    }
}


bool Foam::autoSnapDriver::scaleMesh
(
    const snapParameters& snapParams,
    const label nInitErrors,
    const List<labelPair>& baffles,
    motionSmoother& meshMover
)
{
    const fvMesh& mesh = meshRefiner_.mesh();

    // Relax displacement until correct mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    labelList checkFaces(identity(mesh.nFaces()));

    scalar oldErrorReduction = -1;

    bool meshOk = false;

    Info<< "Moving mesh ..." << endl;
    for (label iter = 0; iter < 2*snapParams.nSnap(); iter++)
    {
        Info<< nl << "Iteration " << iter << endl;

        if (iter == snapParams.nSnap())
        {
            Info<< "Displacement scaling for error reduction set to 0." << endl;
            oldErrorReduction = meshMover.setErrorReduction(0.0);
        }

        meshOk = meshMover.scaleMesh(checkFaces, baffles, true, nInitErrors);

        if (meshOk)
        {
            Info<< "Successfully moved mesh" << endl;
            break;
        }
        if (debug&meshRefinement::MESH)
        {
            const_cast<Time&>(mesh.time())++;
            Info<< "Writing scaled mesh to time " << meshRefiner_.timeName()
                << endl;
            mesh.write();

            Info<< "Writing displacement field ..." << endl;
            meshMover.displacement().write();
            tmp<pointScalarField> magDisp(mag(meshMover.displacement()));
            magDisp().write();
        }
    }

    if (oldErrorReduction >= 0)
    {
        meshMover.setErrorReduction(oldErrorReduction);
    }
    Info<< "Moved mesh in = "
        << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;

    return meshOk;
}


// After snapping: correct patching according to nearest surface.
// Code is very similar to calcNearestSurface.
// - calculate face-wise snap distance as max of point-wise
// - calculate face-wise nearest surface point
// - repatch face according to patch for surface point.
Foam::autoPtr<Foam::mapPolyMesh> Foam::autoSnapDriver::repatchToSurface
(
    const snapParameters& snapParams,
    const labelList& adaptPatchIDs,
    const labelList& preserveFaces
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const refinementSurfaces& surfaces = meshRefiner_.surfaces();

    Info<< "Repatching faces according to nearest surface ..." << endl;

    // Get the labels of added patches.
    autoPtr<indirectPrimitivePatch> ppPtr
    (
        meshRefinement::makePatch
        (
            mesh,
            adaptPatchIDs
        )
    );
    indirectPrimitivePatch& pp = ppPtr();

    // Divide surfaces into zoned and unzoned
    labelList zonedSurfaces =
        surfaceZonesInfo::getNamedSurfaces(surfaces.surfZones());
    labelList unzonedSurfaces =
        surfaceZonesInfo::getUnnamedSurfaces(surfaces.surfZones());


    // Faces that do not move
    PackedBoolList isZonedFace(mesh.nFaces());
    {
        // 1. Preserve faces in preserveFaces list
        forAll(preserveFaces, faceI)
        {
            if (preserveFaces[faceI] != -1)
            {
                isZonedFace.set(faceI, 1);
            }
        }

        // 2. All faces on zoned surfaces
        const PtrList<surfaceZonesInfo>& surfZones = surfaces.surfZones();
        const faceZoneMesh& fZones = mesh.faceZones();

        forAll(zonedSurfaces, i)
        {
            const label zoneSurfI = zonedSurfaces[i];
            const faceZone& fZone = fZones[surfZones[zoneSurfI].faceZoneName()];

            forAll(fZone, i)
            {
                isZonedFace.set(fZone[i], 1);
            }
        }
    }


    // Determine per pp face which patch it should be in
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Patch that face should be in
    labelList closestPatch(pp.size(), -1);
    {
        // face snap distance as max of point snap distance
        scalarField faceSnapDist(pp.size(), -GREAT);
        {
            // Distance to attract to nearest feature on surface
            const scalarField snapDist
            (
                calcSnapDistance
                (
                    mesh,
                    snapParams,
                    pp
                )
            );

            const faceList& localFaces = pp.localFaces();

            forAll(localFaces, faceI)
            {
                const face& f = localFaces[faceI];

                forAll(f, fp)
                {
                    faceSnapDist[faceI] = max
                    (
                        faceSnapDist[faceI],
                        snapDist[f[fp]]
                    );
                }
            }
        }

        pointField localFaceCentres(mesh.faceCentres(), pp.addressing());

        // Get nearest surface and region
        labelList hitSurface;
        labelList hitRegion;
        surfaces.findNearestRegion
        (
            unzonedSurfaces,
            localFaceCentres,
            sqr(faceSnapDist),    // sqr of attract distance
            hitSurface,
            hitRegion
        );

        // Get patch
        forAll(pp, i)
        {
            label faceI = pp.addressing()[i];

            if (hitSurface[i] != -1 && !isZonedFace.get(faceI))
            {
                closestPatch[i] = globalToMasterPatch_
                [
                    surfaces.globalRegion
                    (
                        hitSurface[i],
                        hitRegion[i]
                    )
                ];
            }
        }
    }


    // Change those faces for which there is a different closest patch
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList ownPatch(mesh.nFaces(), -1);
    labelList neiPatch(mesh.nFaces(), -1);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        forAll(pp, i)
        {
            ownPatch[pp.start()+i] = patchI;
            neiPatch[pp.start()+i] = patchI;
        }
    }

    label nChanged = 0;
    forAll(closestPatch, i)
    {
        label faceI = pp.addressing()[i];

        if (closestPatch[i] != -1 && closestPatch[i] != ownPatch[faceI])
        {
            ownPatch[faceI] = closestPatch[i];
            neiPatch[faceI] = closestPatch[i];
            nChanged++;
        }
    }

    Info<< "Repatched " << returnReduce(nChanged, sumOp<label>())
        << " faces in = " << mesh.time().cpuTimeIncrement() << " s\n" << nl
        << endl;

    return meshRefiner_.createBaffles(ownPatch, neiPatch);
}


void Foam::autoSnapDriver::detectWarpedFaces
(
    const scalar featureCos,
    const indirectPrimitivePatch& pp,

    DynamicList<label>& splitFaces,
    DynamicList<labelPair>& splits
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const faceList& localFaces = pp.localFaces();
    const pointField& localPoints = pp.localPoints();
    const labelList& bFaces = pp.addressing();

    splitFaces.clear();
    splitFaces.setCapacity(bFaces.size());
    splits.clear();
    splits.setCapacity(bFaces.size());

    // Determine parallel consistent normals on points
    const vectorField pointNormals(PatchTools::pointNormals(mesh, pp));

    face f0(4);
    face f1(4);

    forAll(localFaces, faceI)
    {
        const face& f = localFaces[faceI];

        if (f.size() >= 4)
        {
            // See if splitting face across diagonal would make two faces with
            // biggish normal angle

            labelPair minDiag(-1, -1);
            scalar minCos(GREAT);

            for (label startFp = 0; startFp < f.size()-2; startFp++)
            {
                label minFp = f.rcIndex(startFp);

                for
                (
                    label endFp = f.fcIndex(f.fcIndex(startFp));
                    endFp < f.size() && endFp != minFp;
                    endFp++
                )
                {
                    // Form two faces
                    f0.setSize(endFp-startFp+1);
                    label i0 = 0;
                    for (label fp = startFp; fp <= endFp; fp++)
                    {
                        f0[i0++] = f[fp];
                    }
                    f1.setSize(f.size()+2-f0.size());
                    label i1 = 0;
                    for (label fp = endFp; fp != startFp; fp = f.fcIndex(fp))
                    {
                        f1[i1++] = f[fp];
                    }
                    f1[i1++] = f[startFp];

                    //Info<< "Splitting face:" << f << " into f0:" << f0
                    //    << " f1:" << f1 << endl;

                    vector n0 = f0.normal(localPoints);
                    scalar n0Mag = mag(n0);
                    vector n1 = f1.normal(localPoints);
                    scalar n1Mag = mag(n1);

                    if (n0Mag > ROOTVSMALL && n1Mag > ROOTVSMALL)
                    {
                        scalar cosAngle = (n0/n0Mag) & (n1/n1Mag);
                        if (cosAngle < minCos)
                        {
                            minCos = cosAngle;
                            minDiag = labelPair(startFp, endFp);
                        }
                    }
                }
            }


            if (minCos < featureCos)
            {
                splitFaces.append(bFaces[faceI]);
                splits.append(minDiag);
            }
        }
    }
}


void Foam::autoSnapDriver::doSnap
(
    const dictionary& snapDict,
    const dictionary& motionDict,
    const scalar featureCos,
    const scalar planarAngle,
    const snapParameters& snapParams
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl
        << "Morphing phase" << nl
        << "--------------" << nl
        << endl;

    // Get the labels of added patches.
    labelList adaptPatchIDs(meshRefiner_.meshedPatches());


    // faceZone handling
    // ~~~~~~~~~~~~~~~~~
    //
    // We convert all faceZones into baffles during snapping so we can use
    // a standard mesh motion (except for the mesh checking which for baffles
    // created from internal faces should check across the baffles). The state
    // is stored in two variables:
    //      baffles : pairs of boundary faces
    //      duplicateFace : from mesh face to its baffle colleague (or -1 for
    //                      normal faces)
    // There are three types of faceZones according to the faceType property:
    //
    // internal
    // --------
    // - baffles: contains all faces on faceZone so
    //      - mesh checks check across baffles
    //      - they get back merged into internal faces
    // - duplicateFace: from face to duplicate face. Contains
    //   all faces on faceZone to prevents merging patch faces.
    //
    // baffle
    // ------
    // - baffles: contains no faces on faceZone since need not be merged/checked
    //   across
    // - duplicateFace: contains all faces on faceZone to prevent
    //   merging patch faces.
    //
    // boundary
    // --------
    // - baffles: contains no faces on faceZone since need not be merged/checked
    //   across
    // - duplicateFace: contains no faces on faceZone since both sides can
    //   merge faces independently.


    // Create baffles (pairs of faces that share the same points)
    // Baffles stored as owner and neighbour face that have been created.
    List<labelPair> baffles;
    meshRefiner_.createZoneBaffles
    (
        globalToMasterPatch_,
        globalToSlavePatch_,
        baffles
    );

    // Maintain map from face to baffle face (-1 for non-baffle faces). Used
    // later on to prevent patchface merging if faceType=baffle
    labelList duplicateFace(mesh.nFaces(), -1);

    forAll(baffles, i)
    {
        const labelPair& baffle = baffles[i];
        duplicateFace[baffle.first()] = baffle.second();
        duplicateFace[baffle.second()] = baffle.first();
    }

    // Selectively 'forget' about the baffles, i.e. not check across them
    // or merge across them.
    {
        const faceZoneMesh& fZones = mesh.faceZones();
        const refinementSurfaces& surfaces = meshRefiner_.surfaces();
        const PtrList<surfaceZonesInfo>& surfZones = surfaces.surfZones();

        // Determine which
        //  - faces to remove from list of baffles (so not merge)
        //  - points to duplicate

        // Per face if is on faceType 'baffle' or 'boundary'
        labelList filterFace(mesh.nFaces(), -1);
        label nFilterFaces = 0;
        // Per point whether it need to be duplicated
        PackedBoolList duplicatePoint(mesh.nPoints());
        label nDuplicatePoints = 0;
        forAll(surfZones, surfI)
        {
            const word& faceZoneName = surfZones[surfI].faceZoneName();

            if (faceZoneName.size())
            {
                const surfaceZonesInfo::faceZoneType& faceType =
                    surfZones[surfI].faceType();

                if
                (
                    faceType == surfaceZonesInfo::BAFFLE
                 || faceType == surfaceZonesInfo::BOUNDARY
                )
                {
                    // Filter out all faces for this zone.
                    label zoneI = fZones.findZoneID(faceZoneName);
                    const faceZone& fZone = fZones[zoneI];
                    forAll(fZone, i)
                    {
                        label faceI = fZone[i];
                        filterFace[faceI] = zoneI;
                        nFilterFaces++;
                    }

                    if (faceType == surfaceZonesInfo::BOUNDARY)
                    {
                        forAll(fZone, i)
                        {
                            label faceI = fZone[i];

                            // Allow combining patch faces across this face
                            duplicateFace[faceI] = -1;

                            const face& f = mesh.faces()[faceI];
                            forAll(f, fp)
                            {
                                if (!duplicatePoint[f[fp]])
                                {
                                    duplicatePoint[f[fp]] = 1;
                                    nDuplicatePoints++;
                                }
                            }
                        }
                    }

                    Info<< "Surface : " << surfaces.names()[surfI] << nl
                        << "    faces to become baffle : "
                        << returnReduce(nFilterFaces, sumOp<label>()) << nl
                        << "    points to duplicate    : "
                        << returnReduce(nDuplicatePoints, sumOp<label>())
                        << endl;
                }
            }
        }

        // Duplicate points only if all points agree
        syncTools::syncPointList
        (
            mesh,
            duplicatePoint,
            andEqOp<unsigned int>(),    // combine op
            0u                          // null value
        );
        // Mark as duplicate (avoids combining patch faces) if one or both
        syncTools::syncFaceList(mesh, duplicateFace, maxEqOp<label>());
        // Mark as resulting from baffle/boundary face zone only if both agree
        syncTools::syncFaceList(mesh, filterFace, minEqOp<label>());

        // Duplicate points
        if (returnReduce(nDuplicatePoints, sumOp<label>()) > 0)
        {
            // Collect all points (recount since syncPointList might have
            // increased set)
            nDuplicatePoints = 0;
            forAll(duplicatePoint, pointI)
            {
                if (duplicatePoint[pointI])
                {
                    nDuplicatePoints++;
                }
            }
            labelList candidatePoints(nDuplicatePoints);
            nDuplicatePoints = 0;
            forAll(duplicatePoint, pointI)
            {
                if (duplicatePoint[pointI])
                {
                    candidatePoints[nDuplicatePoints++] = pointI;
                }
            }


            localPointRegion regionSide(mesh, candidatePoints);
            autoPtr<mapPolyMesh> mapPtr = meshRefiner_.dupNonManifoldPoints
            (
                regionSide
            );
            meshRefinement::updateList(mapPtr().faceMap(), -1, filterFace);
            meshRefinement::updateList(mapPtr().faceMap(), -1, duplicateFace);

            // Update baffles and baffle-to-baffle addressing

            const labelList& reverseFaceMap = mapPtr().reverseFaceMap();

            forAll(baffles, i)
            {
                labelPair& baffle = baffles[i];
                baffle.first() = reverseFaceMap[baffle.first()];
                baffle.second() = reverseFaceMap[baffle.second()];
            }

            if (debug&meshRefinement::MESH)
            {
                const_cast<Time&>(mesh.time())++;
                Pout<< "Writing duplicatedPoints mesh to time "
                    << meshRefiner_.timeName()
                    << endl;
                meshRefiner_.write
                (
                    meshRefinement::debugType(debug),
                    meshRefinement::writeType
                    (
                        meshRefinement::writeLevel()
                      | meshRefinement::WRITEMESH
                    ),
                    mesh.time().path()/"duplicatedPoints"
                );
            }
        }


        // Forget about baffles in a BAFFLE/BOUNDARY type zone
        DynamicList<labelPair> newBaffles(baffles.size());
        forAll(baffles, i)
        {
            const labelPair& baffle = baffles[i];
            if
            (
                filterFace[baffle.first()] == -1
             && filterFace[baffles[i].second()] == -1
            )
            {
                newBaffles.append(baffle);
            }
        }

        if (newBaffles.size() < baffles.size())
        {
            //Info<< "Splitting baffles into" << nl
            //    << "    internal : " << newBaffles.size() << nl
            //    << "    baffle   : " << baffles.size()-newBaffles.size()
            //    << nl << endl;
            baffles.transfer(newBaffles);
        }
        Info<< endl;
    }


    bool doFeatures = false;
    label nFeatIter = 1;
    if (snapParams.nFeatureSnap() > 0)
    {
        doFeatures = true;
        nFeatIter = snapParams.nFeatureSnap();

        Info<< "Snapping to features in " << nFeatIter
            << " iterations ..." << endl;
    }


    bool meshOk = false;

    {
        autoPtr<indirectPrimitivePatch> ppPtr
        (
            meshRefinement::makePatch
            (
                mesh,
                adaptPatchIDs
            )
        );


        // Distance to attract to nearest feature on surface
        const scalarField snapDist(calcSnapDistance(mesh, snapParams, ppPtr()));


        // Construct iterative mesh mover.
        Info<< "Constructing mesh displacer ..." << endl;
        Info<< "Using mesh parameters " << motionDict << nl << endl;

        autoPtr<motionSmoother> meshMoverPtr
        (
            new motionSmoother
            (
                mesh,
                ppPtr(),
                adaptPatchIDs,
                meshRefinement::makeDisplacementField
                (
                    pointMesh::New(mesh),
                    adaptPatchIDs
                ),
                motionDict
            )
        );


        // Check initial mesh
        Info<< "Checking initial mesh ..." << endl;
        labelHashSet wrongFaces(mesh.nFaces()/100);
        motionSmoother::checkMesh(false, mesh, motionDict, wrongFaces);
        const label nInitErrors = returnReduce
        (
            wrongFaces.size(),
            sumOp<label>()
        );

        Info<< "Detected " << nInitErrors << " illegal faces"
            << " (concave, zero area or negative cell pyramid volume)"
            << endl;


        Info<< "Checked initial mesh in = "
            << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;

        // Pre-smooth patch vertices (so before determining nearest)
        preSmoothPatch
        (
            meshRefiner_,
            snapParams,
            nInitErrors,
            baffles,
            meshMoverPtr()
        );



        //- Only if in feature attraction mode:
        // Nearest feature
        vectorField patchAttraction;
        // Constraints at feature
        List<pointConstraint> patchConstraints;


        for (label iter = 0; iter < nFeatIter; iter++)
        {
            //if (doFeatures && (iter == 0 || iter == nFeatIter/2))
            //{
            //    Info<< "Splitting diagonal attractions" << endl;
            //
            //    indirectPrimitivePatch& pp = ppPtr();
            //    motionSmoother& meshMover = meshMoverPtr();
            //
            //    // Calculate displacement at every patch point. Insert into
            //    // meshMover.
            //    // Calculate displacement at every patch point
            //    pointField nearestPoint;
            //    vectorField nearestNormal;
            //
            //    if (snapParams.detectNearSurfacesSnap())
            //    {
            //        nearestPoint.setSize(pp.nPoints(), vector::max);
            //        nearestNormal.setSize(pp.nPoints(), vector::zero);
            //    }
            //
            //    vectorField disp = calcNearestSurface
            //    (
            //        meshRefiner_,
            //        snapDist,
            //        pp,
            //        nearestPoint,
            //        nearestNormal
            //    );
            //
            //
            //    // Override displacement at thin gaps
            //    if (snapParams.detectNearSurfacesSnap())
            //    {
            //        detectNearSurfaces
            //        (
            //            Foam::cos(degToRad(planarAngle)),// planar gaps
            //            pp,
            //            nearestPoint,   // surfacepoint from nearest test
            //            nearestNormal,  // surfacenormal from nearest test
            //
            //            disp
            //        );
            //    }
            //
            //    // Override displacement with feature edge attempt
            //    const label iter = 0;
            //    calcNearestSurfaceFeature
            //    (
            //        snapParams,
            //        false,   // avoidSnapProblems
            //        iter,
            //        featureCos,
            //        scalar(iter+1)/nFeatIter,
            //        snapDist,
            //        disp,
            //        meshMover,
            //        patchAttraction,
            //        patchConstraints
            //    );
            //
            //
            //    const labelList& bFaces = ppPtr().addressing();
            //    DynamicList<label> splitFaces(bFaces.size());
            //    DynamicList<labelPair> splits(bFaces.size());
            //
            //    forAll(bFaces, faceI)
            //    {
            //        const labelPair split
            //        (
            //            findDiagonalAttraction
            //            (
            //                ppPtr(),
            //                patchAttraction,
            //                patchConstraints,
            //                faceI
            //            )
            //        );
            //
            //        if (split != labelPair(-1, -1))
            //        {
            //            splitFaces.append(bFaces[faceI]);
            //            splits.append(split);
            //        }
            //    }
            //
            //    Info<< "Splitting "
            //        << returnReduce(splitFaces.size(), sumOp<label>())
            //        << " faces along diagonal attractions" << endl;
            //
            //    autoPtr<mapPolyMesh> mapPtr = meshRefiner_.splitFaces
            //    (
            //        splitFaces,
            //        splits
            //    );
            //
            //    const labelList& faceMap = mapPtr().faceMap();
            //    meshRefinement::updateList(faceMap, -1, duplicateFace);
            //    const labelList& reverseFaceMap = mapPtr().reverseFaceMap();
            //    forAll(baffles, i)
            //    {
            //        labelPair& baffle = baffles[i];
            //        baffle.first() = reverseFaceMap[baffle.first()];
            //        baffle.second() = reverseFaceMap[baffle.second()];
            //    }
            //
            //    meshMoverPtr.clear();
            //    ppPtr.clear();
            //
            //    ppPtr = meshRefinement::makePatch(mesh, adaptPatchIDs);
            //    meshMoverPtr.reset
            //    (
            //        new motionSmoother
            //        (
            //            mesh,
            //            ppPtr(),
            //            adaptPatchIDs,
            //            meshRefinement::makeDisplacementField
            //            (
            //                pointMesh::New(mesh),
            //                adaptPatchIDs
            //            ),
            //            motionDict
            //        )
            //    );
            //
            //    if (debug&meshRefinement::MESH)
            //    {
            //        const_cast<Time&>(mesh.time())++;
            //        Info<< "Writing split diagonal mesh to time "
            //            << meshRefiner_.timeName() << endl;
            //        meshRefiner_.write
            //        (
            //            meshRefinement::debugType(debug),
            //            meshRefinement::writeType
            //            (
            //                meshRefinement::writeLevel()
            //              | meshRefinement::WRITEMESH
            //            ),
            //            mesh.time().path()/meshRefiner_.timeName()
            //        );
            //    }
            //}
            //else
            //if
            //(
            //    doFeatures
            // && (iter == 1 || iter == nFeatIter/2+1 || iter == nFeatIter-1)
            //)
            //{
            //    Info<< "Splitting warped faces" << endl;
            //
            //    const labelList& bFaces = ppPtr().addressing();
            //    DynamicList<label> splitFaces(bFaces.size());
            //    DynamicList<labelPair> splits(bFaces.size());
            //
            //    detectWarpedFaces
            //    (
            //        featureCos,
            //        ppPtr(),
            //
            //        splitFaces,
            //        splits
            //    );
            //
            //    Info<< "Splitting "
            //        << returnReduce(splitFaces.size(), sumOp<label>())
            //        << " faces along diagonal to avoid warpage" << endl;
            //
            //    autoPtr<mapPolyMesh> mapPtr = meshRefiner_.splitFaces
            //    (
            //        splitFaces,
            //        splits
            //    );
            //
            //    const labelList& faceMap = mapPtr().faceMap();
            //    meshRefinement::updateList(faceMap, -1, duplicateFace);
            //    const labelList& reverseFaceMap = mapPtr().reverseFaceMap();
            //    forAll(baffles, i)
            //    {
            //        labelPair& baffle = baffles[i];
            //        baffle.first() = reverseFaceMap[baffle.first()];
            //        baffle.second() = reverseFaceMap[baffle.second()];
            //    }
            //
            //    meshMoverPtr.clear();
            //    ppPtr.clear();
            //
            //    ppPtr = meshRefinement::makePatch(mesh, adaptPatchIDs);
            //    meshMoverPtr.reset
            //    (
            //        new motionSmoother
            //        (
            //            mesh,
            //            ppPtr(),
            //            adaptPatchIDs,
            //            meshRefinement::makeDisplacementField
            //            (
            //                pointMesh::New(mesh),
            //                adaptPatchIDs
            //            ),
            //            motionDict
            //        )
            //    );
            //
            //    if (debug&meshRefinement::MESH)
            //    {
            //        const_cast<Time&>(mesh.time())++;
            //        Info<< "Writing split warped mesh to time "
            //            << meshRefiner_.timeName() << endl;
            //        meshRefiner_.write
            //        (
            //            meshRefinement::debugType(debug),
            //            meshRefinement::writeType
            //            (
            //                meshRefinement::writeLevel()
            //              | meshRefinement::WRITEMESH
            //            ),
            //            mesh.time().path()/meshRefiner_.timeName()
            //        );
            //    }
            //}



            Info<< nl
                << "Morph iteration " << iter << nl
                << "-----------------" << endl;

            indirectPrimitivePatch& pp = ppPtr();
            motionSmoother& meshMover = meshMoverPtr();


            // Calculate displacement at every patch point. Insert into
            // meshMover.
            // Calculate displacement at every patch point
            pointField nearestPoint;
            vectorField nearestNormal;

            if (snapParams.detectNearSurfacesSnap())
            {
                nearestPoint.setSize(pp.nPoints(), vector::max);
                nearestNormal.setSize(pp.nPoints(), vector::zero);
            }

            vectorField disp = calcNearestSurface
            (
                meshRefiner_,
                snapDist,
                pp,
                nearestPoint,
                nearestNormal
            );


            // Override displacement at thin gaps
            if (snapParams.detectNearSurfacesSnap())
            {
                detectNearSurfaces
                (
                    Foam::cos(degToRad(planarAngle)),// planar cos for gaps
                    pp,
                    nearestPoint,   // surfacepoint from nearest test
                    nearestNormal,  // surfacenormal from nearest test

                    disp
                );
            }

            // Override displacement with feature edge attempt
            if (doFeatures)
            {
                disp = calcNearestSurfaceFeature
                (
                    snapParams,
                    true,   // avoidSnapProblems
                    iter,
                    featureCos,
                    scalar(iter+1)/nFeatIter,
                    snapDist,
                    disp,
                    meshMover,
                    patchAttraction,
                    patchConstraints
                );
            }

            // Check for displacement being outwards.
            outwardsDisplacement(pp, disp);

            // Set initial distribution of displacement field (on patches)
            // from patchDisp and make displacement consistent with b.c.
            // on displacement pointVectorField.
            meshMover.setDisplacement(disp);


            if (debug&meshRefinement::ATTRACTION)
            {
                dumpMove
                (
                    mesh.time().path()
                  / "patchDisplacement_" + name(iter) + ".obj",
                    pp.localPoints(),
                    pp.localPoints() + disp
                );
            }

            // Get smoothly varying internal displacement field.
            smoothDisplacement(snapParams, meshMover);

            // Apply internal displacement to mesh.
            meshOk = scaleMesh
            (
                snapParams,
                nInitErrors,
                baffles,
                meshMover
            );

            if (!meshOk)
            {
                WarningIn("autoSnapDriver::doSnap(..)")
                    << "Did not succesfully snap mesh."
                    << " Continuing to snap to resolve easy" << nl
                    << "    surfaces but the"
                    << " resulting mesh will not satisfy your quality"
                    << " constraints" << nl << endl;
            }

            if (debug&meshRefinement::MESH)
            {
                const_cast<Time&>(mesh.time())++;
                Info<< "Writing scaled mesh to time "
                    << meshRefiner_.timeName() << endl;
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
                Info<< "Writing displacement field ..." << endl;
                meshMover.displacement().write();
                tmp<pointScalarField> magDisp
                (
                    mag(meshMover.displacement())
                );
                magDisp().write();
            }

            // Use current mesh as base mesh
            meshMover.correct();
        }
    }


    // Merge any introduced baffles (from faceZones of faceType 'internal')
    {
        autoPtr<mapPolyMesh> mapPtr = mergeZoneBaffles(baffles);

        if (mapPtr.valid())
        {
            forAll(duplicateFace, faceI)
            {
                if (duplicateFace[faceI] != -1)
                {
                    duplicateFace[faceI] = mapPtr().reverseFaceMap()[faceI];
                }
            }
        }
    }

    // Repatch faces according to nearest. Do not repatch baffle faces.
    {
        autoPtr<mapPolyMesh> mapPtr = repatchToSurface
        (
            snapParams,
            adaptPatchIDs,
            duplicateFace
        );
        meshRefinement::updateList(mapPtr().faceMap(), -1, duplicateFace);
    }

    // Repatching might have caused faces to be on same patch and hence
    // mergeable so try again to merge coplanar faces. Do not merge baffle
    // faces to ensure they both stay the same.
    label nChanged = meshRefiner_.mergePatchFacesUndo
    (
        featureCos,     // minCos
        featureCos,     // concaveCos
        meshRefiner_.meshedPatches(),
        motionDict,
        duplicateFace   // faces not to merge
    );

    nChanged += meshRefiner_.mergeEdgesUndo(featureCos, motionDict);

    if (nChanged > 0 && debug & meshRefinement::MESH)
    {
        const_cast<Time&>(mesh.time())++;
        Info<< "Writing patchFace merged mesh to time "
            << meshRefiner_.timeName() << endl;
        meshRefiner_.write
        (
            meshRefinement::debugType(debug),
            meshRefinement::writeType
            (
                meshRefinement::writeLevel()
              | meshRefinement::WRITEMESH
            ),
            meshRefiner_.timeName()
        );
    }
}


// ************************************************************************* //
