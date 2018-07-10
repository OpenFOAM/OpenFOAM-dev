/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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
    All to do with adding cell layers

\*---------------------------------------------------------------------------*/

#include "snappyLayerDriver.H"
#include "fvMesh.H"
#include "Time.H"
#include "meshRefinement.H"
#include "removePoints.H"
#include "pointFields.H"
#include "motionSmoother.H"
#include "unitConversion.H"
#include "pointSet.H"
#include "faceSet.H"
#include "cellSet.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "addPatchCellLayer.H"
#include "mapDistributePolyMesh.H"
#include "OBJstream.H"
#include "layerParameters.H"
#include "combineFaces.H"
#include "IOmanip.H"
#include "globalIndex.H"
#include "DynamicField.H"
#include "PatchTools.H"
#include "slipPointPatchFields.H"
#include "fixedValuePointPatchFields.H"
#include "zeroFixedValuePointPatchFields.H"
#include "calculatedPointPatchFields.H"
#include "cyclicSlipPointPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "localPointRegion.H"

#include "externalDisplacementMeshMover.H"
#include "medialAxisMeshMover.H"
#include "scalarIOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(snappyLayerDriver, 0);

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// For debugging: Dump displacement to .obj files
void Foam::snappyLayerDriver::dumpDisplacement
(
    const fileName& prefix,
    const indirectPrimitivePatch& pp,
    const vectorField& patchDisp,
    const List<extrudeMode>& extrudeStatus
)
{
    OBJstream dispStr(prefix + "_disp.obj");
    Info<< "Writing all displacements to " << dispStr.name() << endl;

    forAll(patchDisp, patchPointi)
    {
        const point& pt = pp.localPoints()[patchPointi];
        dispStr.write(linePointRef(pt, pt + patchDisp[patchPointi]));
    }


    OBJstream illStr(prefix + "_illegal.obj");
    Info<< "Writing invalid displacements to " << illStr.name() << endl;

    forAll(patchDisp, patchPointi)
    {
        if (extrudeStatus[patchPointi] != EXTRUDE)
        {
            const point& pt = pp.localPoints()[patchPointi];
            illStr.write(linePointRef(pt, pt + patchDisp[patchPointi]));
        }
    }
}


Foam::tmp<Foam::scalarField> Foam::snappyLayerDriver::avgPointData
(
    const indirectPrimitivePatch& pp,
    const scalarField& pointFld
)
{
    tmp<scalarField> tfaceFld(new scalarField(pp.size(), 0.0));
    scalarField& faceFld = tfaceFld.ref();

    forAll(pp.localFaces(), facei)
    {
        const face& f = pp.localFaces()[facei];
        if (f.size())
        {
            forAll(f, fp)
            {
                faceFld[facei] += pointFld[f[fp]];
            }
            faceFld[facei] /= f.size();
        }
    }
    return tfaceFld;
}


// Check that primitivePatch is not multiply connected. Collect non-manifold
// points in pointSet.
void Foam::snappyLayerDriver::checkManifold
(
    const indirectPrimitivePatch& fp,
    pointSet& nonManifoldPoints
)
{
    // Check for non-manifold points (surface pinched at point)
    fp.checkPointManifold(false, &nonManifoldPoints);

    // Check for edge-faces (surface pinched at edge)
    const labelListList& edgeFaces = fp.edgeFaces();

    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        if (eFaces.size() > 2)
        {
            const edge& e = fp.edges()[edgeI];

            nonManifoldPoints.insert(fp.meshPoints()[e[0]]);
            nonManifoldPoints.insert(fp.meshPoints()[e[1]]);
        }
    }
}


void Foam::snappyLayerDriver::checkMeshManifold() const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl << "Checking mesh manifoldness ..." << endl;

    // Get all outside faces
    labelList outsideFaces(mesh.nFaces() - mesh.nInternalFaces());

    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
         outsideFaces[facei - mesh.nInternalFaces()] = facei;
    }

    pointSet nonManifoldPoints
    (
        mesh,
        "nonManifoldPoints",
        mesh.nPoints() / 100
    );

    // Build primitivePatch out of faces and check it for problems.
    checkManifold
    (
        indirectPrimitivePatch
        (
            IndirectList<face>(mesh.faces(), outsideFaces),
            mesh.points()
        ),
        nonManifoldPoints
    );

    label nNonManif = returnReduce(nonManifoldPoints.size(), sumOp<label>());

    if (nNonManif > 0)
    {
        Info<< "Outside of mesh is multiply connected across edges or"
            << " points." << nl
            << "This is not a fatal error but might cause some unexpected"
            << " behaviour." << nl
            //<< "Writing " << nNonManif
            //<< " points where this happens to pointSet "
            //<< nonManifoldPoints.name()
            << endl;

        // nonManifoldPoints.instance() = meshRefiner_.timeName();
        // nonManifoldPoints.write();
    }
    Info<< endl;
}



// Unset extrusion on point. Returns true if anything unset.
bool Foam::snappyLayerDriver::unmarkExtrusion
(
    const label patchPointi,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
)
{
    if (extrudeStatus[patchPointi] == EXTRUDE)
    {
        extrudeStatus[patchPointi] = NOEXTRUDE;
        patchNLayers[patchPointi] = 0;
        patchDisp[patchPointi] = Zero;
        return true;
    }
    else if (extrudeStatus[patchPointi] == EXTRUDEREMOVE)
    {
        extrudeStatus[patchPointi] = NOEXTRUDE;
        patchNLayers[patchPointi] = 0;
        patchDisp[patchPointi] = Zero;
        return true;
    }
    else
    {
        return false;
    }
}


// Unset extrusion on face. Returns true if anything unset.
bool Foam::snappyLayerDriver::unmarkExtrusion
(
    const face& localFace,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
)
{
    bool unextruded = false;

    forAll(localFace, fp)
    {
        if
        (
            unmarkExtrusion
            (
                localFace[fp],
                patchDisp,
                patchNLayers,
                extrudeStatus
            )
        )
        {
            unextruded = true;
        }
    }
    return unextruded;
}


// No extrusion at non-manifold points.
void Foam::snappyLayerDriver::handleNonManifolds
(
    const indirectPrimitivePatch& pp,
    const labelList& meshEdges,
    const labelListList& edgeGlobalFaces,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl << "Handling non-manifold points ..." << endl;

    // Detect non-manifold points
    Info<< nl << "Checking patch manifoldness ..." << endl;

    pointSet nonManifoldPoints(mesh, "nonManifoldPoints", pp.nPoints());

    // 1. Local check
    checkManifold(pp, nonManifoldPoints);

    // 2. Remote check for boundary edges on coupled boundaries
    forAll(edgeGlobalFaces, edgeI)
    {
        if
        (
            pp.edgeFaces()[edgeI].size() == 1
         && edgeGlobalFaces[edgeI].size() > 2
        )
        {
            // So boundary edges that are connected to more than 2 processors
            // i.e. a non-manifold edge which is exactly on a processor
            // boundary.
            const edge& e = pp.edges()[edgeI];
            nonManifoldPoints.insert(pp.meshPoints()[e[0]]);
            nonManifoldPoints.insert(pp.meshPoints()[e[1]]);
        }
    }

    // 3. Remote check for end of layer across coupled boundaries
    {
        PackedBoolList isCoupledEdge(mesh.nEdges());

        const labelList& cpEdges = mesh.globalData().coupledPatchMeshEdges();
        forAll(cpEdges, i)
        {
            isCoupledEdge[cpEdges[i]] = true;
        }
        syncTools::syncEdgeList
        (
            mesh,
            isCoupledEdge,
            orEqOp<unsigned int>(),
            0
        );

        forAll(edgeGlobalFaces, edgeI)
        {
            label meshEdgeI = meshEdges[edgeI];

            if
            (
                pp.edgeFaces()[edgeI].size() == 1
             && edgeGlobalFaces[edgeI].size() == 1
             && isCoupledEdge[meshEdgeI]
            )
            {
                // Edge of patch but no continuation across processor.
                const edge& e = pp.edges()[edgeI];
                // Pout<< "** Stopping extrusion on edge "
                //    << pp.localPoints()[e[0]]
                //    << pp.localPoints()[e[1]] << endl;
                nonManifoldPoints.insert(pp.meshPoints()[e[0]]);
                nonManifoldPoints.insert(pp.meshPoints()[e[1]]);
            }
        }
    }



    label nNonManif = returnReduce(nonManifoldPoints.size(), sumOp<label>());

    Info<< "Outside of local patch is multiply connected across edges or"
        << " points at " << nNonManif << " points." << endl;

    if (nNonManif > 0)
    {
        const labelList& meshPoints = pp.meshPoints();

        forAll(meshPoints, patchPointi)
        {
            if (nonManifoldPoints.found(meshPoints[patchPointi]))
            {
                unmarkExtrusion
                (
                    patchPointi,
                    patchDisp,
                    patchNLayers,
                    extrudeStatus
                );
            }
        }
    }

    Info<< "Set displacement to zero for all " << nNonManif
        << " non-manifold points" << endl;
}


// Parallel feature edge detection. Assumes non-manifold edges already handled.
void Foam::snappyLayerDriver::handleFeatureAngle
(
    const indirectPrimitivePatch& pp,
    const labelList& meshEdges,
    const scalar minCos,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl << "Handling feature edges ..." << endl;

    if (minCos < 1-small)
    {
        // Normal component of normals of connected faces.
        vectorField edgeNormal(mesh.nEdges(), point::max);

        const labelListList& edgeFaces = pp.edgeFaces();

        forAll(edgeFaces, edgeI)
        {
            const labelList& eFaces = pp.edgeFaces()[edgeI];

            label meshEdgeI = meshEdges[edgeI];

            forAll(eFaces, i)
            {
                nomalsCombine()
                (
                    edgeNormal[meshEdgeI],
                    pp.faceNormals()[eFaces[i]]
                );
            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            edgeNormal,
            nomalsCombine(),
            point::max          // null value
        );

        autoPtr<OBJstream> str;
        if (debug&meshRefinement::MESH)
        {
            str.reset
            (
                new OBJstream
                (
                    mesh.time().path()
                  / "featureEdges_"
                  + meshRefiner_.timeName()
                  + ".obj"
                )
            );
            Info<< "Writing feature edges to " << str().name() << endl;
        }

        label nFeats = 0;

        // Now on coupled edges the edgeNormal will have been truncated and
        // only be still be the old value where two faces have the same normal
        forAll(edgeFaces, edgeI)
        {
            const labelList& eFaces = pp.edgeFaces()[edgeI];

            label meshEdgeI = meshEdges[edgeI];

            const vector& n = edgeNormal[meshEdgeI];

            if (n != point::max)
            {
                scalar cos = n & pp.faceNormals()[eFaces[0]];

                if (cos < minCos)
                {
                    const edge& e = pp.edges()[edgeI];

                    unmarkExtrusion
                    (
                        e[0],
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    );
                    unmarkExtrusion
                    (
                        e[1],
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    );

                    nFeats++;

                    if (str.valid())
                    {
                        const point& p0 = pp.localPoints()[e[0]];
                        const point& p1 = pp.localPoints()[e[1]];
                        str().write(linePointRef(p0, p1));
                    }
                }
            }
        }

        Info<< "Set displacement to zero for points on "
            << returnReduce(nFeats, sumOp<label>())
            << " feature edges" << endl;
    }
}


// No extrusion on cells with warped faces. Calculates the thickness of the
// layer and compares it to the space the warped face takes up. Disables
// extrusion if layer thickness is more than faceRatio of the thickness of
// the face.
void Foam::snappyLayerDriver::handleWarpedFaces
(
    const indirectPrimitivePatch& pp,
    const scalar faceRatio,
    const scalar edge0Len,
    const labelList& cellLevel,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl << "Handling cells with warped patch faces ..." << nl;

    const pointField& points = mesh.points();

    label nWarpedFaces = 0;

    forAll(pp, i)
    {
        const face& f = pp[i];

        if (f.size() > 3)
        {
            label facei = pp.addressing()[i];

            label ownLevel = cellLevel[mesh.faceOwner()[facei]];
            scalar edgeLen = edge0Len/(1<<ownLevel);

            // Normal distance to face centre plane
            const point& fc = mesh.faceCentres()[facei];
            const vector& fn = pp.faceNormals()[i];

            scalarField vProj(f.size());

            forAll(f, fp)
            {
                vector n = points[f[fp]] - fc;
                vProj[fp] = (n & fn);
            }

            // Get normal 'span' of face
            scalar minVal = min(vProj);
            scalar maxVal = max(vProj);

            if ((maxVal - minVal) > faceRatio * edgeLen)
            {
                if
                (
                    unmarkExtrusion
                    (
                        pp.localFaces()[i],
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    )
                )
                {
                    nWarpedFaces++;
                }
            }
        }
    }

    Info<< "Set displacement to zero on "
        << returnReduce(nWarpedFaces, sumOp<label>())
        << " warped faces since layer would be > " << faceRatio
        << " of the size of the bounding box." << endl;
}


//// No extrusion on cells with multiple patch faces. There usually is a reason
//// why combinePatchFaces hasn't succeeded.
//void Foam::snappyLayerDriver::handleMultiplePatchFaces
//(
//    const indirectPrimitivePatch& pp,
//    pointField& patchDisp,
//    labelList& patchNLayers,
//    List<extrudeMode>& extrudeStatus
//) const
//{
//    const fvMesh& mesh = meshRefiner_.mesh();
//
//    Info<< nl << "Handling cells with multiple patch faces ..." << nl;
//
//    const labelListList& pointFaces = pp.pointFaces();
//
//    // Cells that should not get an extrusion layer
//    cellSet multiPatchCells(mesh, "multiPatchCells", pp.size());
//
//    // Detect points that use multiple faces on same cell.
//    forAll(pointFaces, patchPointi)
//    {
//        const labelList& pFaces = pointFaces[patchPointi];
//
//        labelHashSet pointCells(pFaces.size());
//
//        forAll(pFaces, i)
//        {
//            label celli = mesh.faceOwner()[pp.addressing()[pFaces[i]]];
//
//            if (!pointCells.insert(celli))
//            {
//                // Second or more occurrence of cell so cell has two or more
//                // pp faces connected to this point.
//                multiPatchCells.insert(celli);
//            }
//        }
//    }
//
//    label nMultiPatchCells = returnReduce
//    (
//        multiPatchCells.size(),
//        sumOp<label>()
//    );
//
//    Info<< "Detected " << nMultiPatchCells
//        << " cells with multiple (connected) patch faces." << endl;
//
//    label nChanged = 0;
//
//    if (nMultiPatchCells > 0)
//    {
//        multiPatchCells.instance() = meshRefiner_.timeName();
//        Info<< "Writing " << nMultiPatchCells
//            << " cells with multiple (connected) patch faces to cellSet "
//            << multiPatchCells.objectPath() << endl;
//        multiPatchCells.write();
//
//
//        // Go through all points and remove extrusion on any cell in
//        // multiPatchCells
//        // (has to be done in separate loop since having one point on
//        // multipatches has to reset extrusion on all points of cell)
//
//        forAll(pointFaces, patchPointi)
//        {
//            if (extrudeStatus[patchPointi] != NOEXTRUDE)
//            {
//                const labelList& pFaces = pointFaces[patchPointi];
//
//                forAll(pFaces, i)
//                {
//                    label celli =
//                        mesh.faceOwner()[pp.addressing()[pFaces[i]]];
//
//                    if (multiPatchCells.found(celli))
//                    {
//                        if
//                        (
//                            unmarkExtrusion
//                            (
//                                patchPointi,
//                                patchDisp,
//                                patchNLayers,
//                                extrudeStatus
//                            )
//                        )
//                        {
//                            nChanged++;
//                        }
//                    }
//                }
//            }
//        }
//
//        reduce(nChanged, sumOp<label>());
//    }
//
//    Info<< "Prevented extrusion on " << nChanged
//        << " points due to multiple patch faces." << nl << endl;
//}


void Foam::snappyLayerDriver::setNumLayers
(
    const labelList& patchToNLayers,
    const labelList& patchIDs,
    const indirectPrimitivePatch& pp,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus,
    label& nAddedCells
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl << "Handling points with inconsistent layer specification ..."
        << endl;

    // Get for every point (really only necessary on patch external points)
    // the max and min of any patch faces using it.
    labelList maxLayers(patchNLayers.size(), labelMin);
    labelList minLayers(patchNLayers.size(), labelMax);

    forAll(patchIDs, i)
    {
        label patchi = patchIDs[i];

        const labelList& meshPoints = mesh.boundaryMesh()[patchi].meshPoints();

        label wantedLayers = patchToNLayers[patchi];

        forAll(meshPoints, patchPointi)
        {
            label ppPointi = pp.meshPointMap()[meshPoints[patchPointi]];

            maxLayers[ppPointi] = max(wantedLayers, maxLayers[ppPointi]);
            minLayers[ppPointi] = min(wantedLayers, minLayers[ppPointi]);
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        maxLayers,
        maxEqOp<label>(),
        labelMin            // null value
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        minLayers,
        minEqOp<label>(),
        labelMax            // null value
    );

    // Unmark any point with different min and max
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // label nConflicts = 0;

    forAll(maxLayers, i)
    {
        if (maxLayers[i] == labelMin || minLayers[i] == labelMax)
        {
            FatalErrorInFunction
                << "Patchpoint:" << i << " coord:" << pp.localPoints()[i]
                << " maxLayers:" << maxLayers
                << " minLayers:" << minLayers
                << abort(FatalError);
        }
        else if (maxLayers[i] == minLayers[i])
        {
            // Ok setting.
            patchNLayers[i] = maxLayers[i];
        }
        else
        {
            // Inconsistent num layers between patch faces using point
            // if
            //(
            //    unmarkExtrusion
            //    (
            //        i,
            //        patchDisp,
            //        patchNLayers,
            //        extrudeStatus
            //    )
            //)
            //{
            //    nConflicts++;
            //}
            patchNLayers[i] = maxLayers[i];
        }
    }


    // Calculate number of cells to create
    nAddedCells = 0;
    forAll(pp.localFaces(), facei)
    {
        const face& f = pp.localFaces()[facei];

        // Get max of extrusion per point
        label nCells = 0;
        forAll(f, fp)
        {
            nCells = max(nCells, patchNLayers[f[fp]]);
        }

        nAddedCells += nCells;
    }
    reduce(nAddedCells, sumOp<label>());

    // reduce(nConflicts, sumOp<label>());
    //
    // Info<< "Set displacement to zero for " << nConflicts
    //    << " points due to points being on multiple regions"
    //    << " with inconsistent nLayers specification." << endl;
}


// Construct pointVectorField with correct boundary conditions for adding
// layers
Foam::tmp<Foam::pointVectorField>
Foam::snappyLayerDriver::makeLayerDisplacementField
(
    const pointMesh& pMesh,
    const labelList& numLayers
)
{
    // Construct displacement field.
    const pointBoundaryMesh& pointPatches = pMesh.boundary();

    wordList patchFieldTypes
    (
        pointPatches.size(),
        slipPointPatchVectorField::typeName
    );
    wordList actualPatchTypes(patchFieldTypes.size());
    forAll(pointPatches, patchi)
    {
        actualPatchTypes[patchi] = pointPatches[patchi].type();
    }

    forAll(numLayers, patchi)
    {
        //  0 layers: do not allow slip so fixedValue 0
        // >0 layers: fixedValue which gets adapted
        if (numLayers[patchi] == 0)
        {
            patchFieldTypes[patchi] =
                zeroFixedValuePointPatchVectorField::typeName;
        }
        else if (numLayers[patchi] > 0)
        {
            patchFieldTypes[patchi] = fixedValuePointPatchVectorField::typeName;
        }
    }

    forAll(pointPatches, patchi)
    {
        if (isA<processorPointPatch>(pointPatches[patchi]))
        {
            patchFieldTypes[patchi] = calculatedPointPatchVectorField::typeName;
        }
        else if (isA<cyclicPointPatch>(pointPatches[patchi]))
        {
            patchFieldTypes[patchi] = cyclicSlipPointPatchVectorField::typeName;
        }
    }


    const polyMesh& mesh = pMesh();

    // Note: time().timeName() instead of meshRefinement::timeName() since
    // postprocessable field.
    tmp<pointVectorField> tfld
    (
        new pointVectorField
        (
            IOobject
            (
                "pointDisplacement",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            dimensionedVector("displacement", dimLength, Zero),
            patchFieldTypes,
            actualPatchTypes
        )
    );
    return tfld;
}


void Foam::snappyLayerDriver::growNoExtrusion
(
    const indirectPrimitivePatch& pp,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    Info<< nl << "Growing non-extrusion points by one layer ..." << endl;

    List<extrudeMode> grownExtrudeStatus(extrudeStatus);

    const faceList& localFaces = pp.localFaces();

    label nGrown = 0;

    forAll(localFaces, facei)
    {
        const face& f = localFaces[facei];

        bool hasSqueeze = false;
        forAll(f, fp)
        {
            if (extrudeStatus[f[fp]] == NOEXTRUDE)
            {
                hasSqueeze = true;
                break;
            }
        }

        if (hasSqueeze)
        {
            // Squeeze all points of face
            forAll(f, fp)
            {
                if
                (
                    extrudeStatus[f[fp]] == EXTRUDE
                 && grownExtrudeStatus[f[fp]] != NOEXTRUDE
                )
                {
                    grownExtrudeStatus[f[fp]] = NOEXTRUDE;
                    nGrown++;
                }
            }
        }
    }

    extrudeStatus.transfer(grownExtrudeStatus);


    // Synchronise since might get called multiple times.
    // Use the fact that NOEXTRUDE is the minimum value.
    {
        labelList status(extrudeStatus.size());
        forAll(status, i)
        {
            status[i] = extrudeStatus[i];
        }
        syncTools::syncPointList
        (
            meshRefiner_.mesh(),
            pp.meshPoints(),
            status,
            minEqOp<label>(),
            labelMax            // null value
        );
        forAll(status, i)
        {
            extrudeStatus[i] = extrudeMode(status[i]);
        }
    }


    forAll(extrudeStatus, patchPointi)
    {
        if (extrudeStatus[patchPointi] == NOEXTRUDE)
        {
            patchDisp[patchPointi] = Zero;
            patchNLayers[patchPointi] = 0;
        }
    }

    reduce(nGrown, sumOp<label>());

    Info<< "Set displacement to zero for an additional " << nGrown
        << " points." << endl;
}


void Foam::snappyLayerDriver::determineSidePatches
(
    const globalIndex& globalFaces,
    const labelListList& edgeGlobalFaces,
    const indirectPrimitivePatch& pp,

    labelList& sidePatchID
)
{
    // Sometimes edges-to-be-extruded are on more than 2 processors.
    // Work out which 2 hold the faces to be extruded and thus which procpatch
    // the side-face should be in. As an additional complication this might
    // mean that 2 procesors that were only edge-connected now suddenly need
    // to become face-connected i.e. have a processor patch between them.

    fvMesh& mesh = meshRefiner_.mesh();

    // Determine sidePatchID. Any additional processor boundary gets added to
    // patchToNbrProc,nbrProcToPatch and nPatches gets set to the new number
    // of patches.
    label nPatches;
    Map<label> nbrProcToPatch;
    Map<label> patchToNbrProc;
    addPatchCellLayer::calcSidePatch
    (
        mesh,
        globalFaces,
        edgeGlobalFaces,
        pp,

        sidePatchID,
        nPatches,
        nbrProcToPatch,
        patchToNbrProc
    );

    label nOldPatches = mesh.boundaryMesh().size();
    label nAdded = returnReduce(nPatches-nOldPatches, sumOp<label>());
    Info<< nl << "Adding in total " << nAdded/2 << " inter-processor patches to"
        << " handle extrusion of non-manifold processor boundaries."
        << endl;

    if (nAdded > 0)
    {
        // We might not add patches in same order as in patchToNbrProc
        // so prepare to renumber sidePatchID
        Map<label> wantedToAddedPatch;

        for (label patchi = nOldPatches; patchi < nPatches; patchi++)
        {
            label nbrProci = patchToNbrProc[patchi];
            word name
            (
                processorPolyPatch::newName(Pstream::myProcNo(), nbrProci)
            );

            dictionary patchDict;
            patchDict.add("type", processorPolyPatch::typeName);
            patchDict.add("myProcNo", Pstream::myProcNo());
            patchDict.add("neighbProcNo", nbrProci);
            patchDict.add("nFaces", 0);
            patchDict.add("startFace", mesh.nFaces());

            label procPatchi = meshRefiner_.appendPatch
            (
                mesh,
                mesh.boundaryMesh().size(), // new patch index
                name,
                patchDict
            );
            wantedToAddedPatch.insert(patchi, procPatchi);
        }

        // Renumber sidePatchID
        forAll(sidePatchID, i)
        {
            label patchi = sidePatchID[i];
            Map<label>::const_iterator fnd = wantedToAddedPatch.find(patchi);
            if (fnd != wantedToAddedPatch.end())
            {
                sidePatchID[i] = fnd();
            }
        }

        mesh.clearOut();
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh()).updateMesh();
    }
}


void Foam::snappyLayerDriver::calculateLayerThickness
(
    const indirectPrimitivePatch& pp,
    const labelList& patchIDs,
    const layerParameters& layerParams,
    const labelList& cellLevel,
    const labelList& patchNLayers,
    const scalar edge0Len,

    scalarField& thickness,
    scalarField& minThickness,
    scalarField& expansionRatio
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();


    // Rework patch-wise layer parameters into minimum per point
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Note: only layer parameters consistent with layer specification
    // method (see layerParameters) will be correct.
    scalarField firstLayerThickness(pp.nPoints(), great);
    scalarField finalLayerThickness(pp.nPoints(), great);
    scalarField totalThickness(pp.nPoints(), great);
    scalarField expRatio(pp.nPoints(), great);

    minThickness.setSize(pp.nPoints());
    minThickness = great;

    forAll(patchIDs, i)
    {
        label patchi = patchIDs[i];

        const labelList& meshPoints = patches[patchi].meshPoints();

        forAll(meshPoints, patchPointi)
        {
            label ppPointi = pp.meshPointMap()[meshPoints[patchPointi]];

            firstLayerThickness[ppPointi] = min
            (
                firstLayerThickness[ppPointi],
                layerParams.firstLayerThickness()[patchi]
            );
            finalLayerThickness[ppPointi] = min
            (
                finalLayerThickness[ppPointi],
                layerParams.finalLayerThickness()[patchi]
            );
            totalThickness[ppPointi] = min
            (
                totalThickness[ppPointi],
                layerParams.thickness()[patchi]
            );
            expRatio[ppPointi] = min
            (
                expRatio[ppPointi],
                layerParams.expansionRatio()[patchi]
            );
            minThickness[ppPointi] = min
            (
                minThickness[ppPointi],
                layerParams.minThickness()[patchi]
            );
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        firstLayerThickness,
        minEqOp<scalar>(),
        great               // null value
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        finalLayerThickness,
        minEqOp<scalar>(),
        great               // null value
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        totalThickness,
        minEqOp<scalar>(),
        great               // null value
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        expRatio,
        minEqOp<scalar>(),
        great               // null value
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        minThickness,
        minEqOp<scalar>(),
        great               // null value
    );


    // Now the thicknesses are set according to the minimum of connected
    // patches.


    // Rework relative thickness into absolute
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // by multiplying with the internal cell size.

    if (layerParams.relativeSizes())
    {
        if
        (
            min(layerParams.minThickness()) < 0
         || max(layerParams.minThickness()) > 2
        )
        {
            FatalErrorInFunction
                << "Thickness should be factor of local undistorted cell size."
                << " Valid values are [0..2]." << nl
                << " minThickness:" << layerParams.minThickness()
                << exit(FatalError);
        }


        // Determine per point the max cell level of connected cells
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        labelList maxPointLevel(pp.nPoints(), labelMin);

        forAll(pp, i)
        {
            label ownLevel = cellLevel[mesh.faceOwner()[pp.addressing()[i]]];

            const face& f = pp.localFaces()[i];

            forAll(f, fp)
            {
                maxPointLevel[f[fp]] = max(maxPointLevel[f[fp]], ownLevel);
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            maxPointLevel,
            maxEqOp<label>(),
            labelMin            // null value
        );


        forAll(maxPointLevel, pointi)
        {
            // Find undistorted edge size for this level.
            scalar edgeLen = edge0Len/(1<<maxPointLevel[pointi]);
            firstLayerThickness[pointi] *= edgeLen;
            finalLayerThickness[pointi] *= edgeLen;
            totalThickness[pointi] *= edgeLen;
            minThickness[pointi] *= edgeLen;
        }
    }



    // Rework thickness parameters into overall thickness
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(firstLayerThickness, pointi)
    {
        thickness[pointi] = layerParams.layerThickness
        (
            patchNLayers[pointi],
            firstLayerThickness[pointi],
            finalLayerThickness[pointi],
            totalThickness[pointi],
            expRatio[pointi]
        );

        expansionRatio[pointi] = layerParams.layerExpansionRatio
        (
            patchNLayers[pointi],
            firstLayerThickness[pointi],
            finalLayerThickness[pointi],
            totalThickness[pointi],
            expRatio[pointi]
        );
    }

    // Info<< "calculateLayerThickness : min:" << gMin(thickness)
    //    << " max:" << gMax(thickness) << endl;

    // Print a bit
    {
        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        int oldPrecision = Info().precision();

        // Find maximum length of a patch name, for a nicer output
        label maxPatchNameLen = 0;
        forAll(patchIDs, i)
        {
            label patchi = patchIDs[i];
            word patchName = patches[patchi].name();
            maxPatchNameLen = max(maxPatchNameLen, label(patchName.size()));
        }

        Info<< nl
            << setf(ios_base::left) << setw(maxPatchNameLen) << "patch"
            << setw(0) << " faces    layers avg thickness[m]" << nl
            << setf(ios_base::left) << setw(maxPatchNameLen) << " "
            << setw(0) << "                 near-wall overall" << nl
            << setf(ios_base::left) << setw(maxPatchNameLen) << "-----"
            << setw(0) << " -----    ------ --------- -------" << endl;


        const PackedBoolList isMasterPoint(syncTools::getMasterPoints(mesh));

        forAll(patchIDs, i)
        {
            label patchi = patchIDs[i];

            const labelList& meshPoints = patches[patchi].meshPoints();

            scalar sumThickness = 0;
            scalar sumNearWallThickness = 0;
            label nMasterPoints = 0;

            forAll(meshPoints, patchPointi)
            {
                label meshPointi = meshPoints[patchPointi];
                if (isMasterPoint[meshPointi])
                {
                    label ppPointi = pp.meshPointMap()[meshPointi];

                    sumThickness += thickness[ppPointi];
                    sumNearWallThickness += layerParams.firstLayerThickness
                    (
                        patchNLayers[ppPointi],
                        firstLayerThickness[ppPointi],
                        finalLayerThickness[ppPointi],
                        thickness[ppPointi],
                        expansionRatio[ppPointi]
                    );
                    nMasterPoints++;
                }
            }

            label totNPoints = returnReduce(nMasterPoints, sumOp<label>());

            // For empty patches, totNPoints is 0.
            scalar avgThickness = 0;
            scalar avgNearWallThickness = 0;

            if (totNPoints > 0)
            {
                avgThickness =
                    returnReduce(sumThickness, sumOp<scalar>())
                  / totNPoints;
                avgNearWallThickness =
                    returnReduce(sumNearWallThickness, sumOp<scalar>())
                  / totNPoints;
            }

            Info<< setf(ios_base::left) << setw(maxPatchNameLen)
                << patches[patchi].name() << setprecision(3)
                << " " << setw(8)
                << returnReduce(patches[patchi].size(), sumOp<scalar>())
                << " " << setw(6) << layerParams.numLayers()[patchi]
                << " " << setw(8) << avgNearWallThickness
                << "  " << setw(8) << avgThickness
                << endl;
        }
        Info<< setprecision(oldPrecision) << endl;
    }
}


// Synchronize displacement among coupled patches.
void Foam::snappyLayerDriver::syncPatchDisplacement
(
    const indirectPrimitivePatch& pp,
    const scalarField& minThickness,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const labelList& meshPoints = pp.meshPoints();

    label nChangedTotal = 0;

    while (true)
    {
        label nChanged = 0;

        // Sync displacement (by taking min)
        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            patchDisp,
            minMagSqrEqOp<vector>(),
            point::rootMax      // null value
        );

        // Unmark if displacement too small
        forAll(patchDisp, i)
        {
            if (mag(patchDisp[i]) < minThickness[i])
            {
                if
                (
                    unmarkExtrusion
                    (
                        i,
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    )
                )
                {
                    nChanged++;
                }
            }
        }

        labelList syncPatchNLayers(patchNLayers);

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            syncPatchNLayers,
            minEqOp<label>(),
            labelMax            // null value
        );

        // Reset if differs
        // 1. take max
        forAll(syncPatchNLayers, i)
        {
            if (syncPatchNLayers[i] != patchNLayers[i])
            {
                if
                (
                    unmarkExtrusion
                    (
                        i,
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    )
                )
                {
                    nChanged++;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            syncPatchNLayers,
            maxEqOp<label>(),
            labelMin            // null value
        );

        // Reset if differs
        // 2. take min
        forAll(syncPatchNLayers, i)
        {
            if (syncPatchNLayers[i] != patchNLayers[i])
            {
                if
                (
                    unmarkExtrusion
                    (
                        i,
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    )
                )
                {
                    nChanged++;
                }
            }
        }
        nChangedTotal += nChanged;

        if (!returnReduce(nChanged, sumOp<label>()))
        {
            break;
        }
    }

    // Info<< "Prevented extrusion on "
    //    << returnReduce(nChangedTotal, sumOp<label>())
    //    << " coupled patch points during syncPatchDisplacement." << endl;
}


// Calculate displacement vector for all patch points. Uses pointNormal.
// Checks that displaced patch point would be visible from all centres
// of the faces using it.
// extrudeStatus is both input and output and gives the status of each
// patch point.
void Foam::snappyLayerDriver::getPatchDisplacement
(
    const indirectPrimitivePatch& pp,
    const scalarField& thickness,
    const scalarField& minThickness,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    Info<< nl << "Determining displacement for added points"
        << " according to pointNormal ..." << endl;

    const fvMesh& mesh = meshRefiner_.mesh();
    const vectorField& faceNormals = pp.faceNormals();
    const labelListList& pointFaces = pp.pointFaces();
    const pointField& localPoints = pp.localPoints();

    // Determine pointNormal
    // ~~~~~~~~~~~~~~~~~~~~~

    pointField pointNormals(PatchTools::pointNormals(mesh, pp));


    // Determine local length scale on patch
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Start off from same thickness everywhere (except where no extrusion)
    patchDisp = thickness*pointNormals;


    label nNoVisNormal = 0;
    label nExtrudeRemove = 0;


    // Check if no extrude possible.
    forAll(pointNormals, patchPointi)
    {
        label meshPointi = pp.meshPoints()[patchPointi];

        if (extrudeStatus[patchPointi] == NOEXTRUDE)
        {
            // Do not use unmarkExtrusion; forcibly set to zero extrusion.
            patchNLayers[patchPointi] = 0;
            patchDisp[patchPointi] = Zero;
        }
        else
        {
            // Get normal
            const vector& n = pointNormals[patchPointi];

            if (!meshTools::visNormal(n, faceNormals, pointFaces[patchPointi]))
            {
                if (debug&meshRefinement::ATTRACTION)
                {
                    Pout<< "No valid normal for point " << meshPointi
                        << ' ' << pp.points()[meshPointi]
                        << "; setting displacement to "
                        << patchDisp[patchPointi]
                        << endl;
                }

                extrudeStatus[patchPointi] = EXTRUDEREMOVE;
                nNoVisNormal++;
            }
        }
    }

    // At illegal points make displacement average of new neighbour positions
    forAll(extrudeStatus, patchPointi)
    {
        if (extrudeStatus[patchPointi] == EXTRUDEREMOVE)
        {
            point avg(Zero);
            label nPoints = 0;

            const labelList& pEdges = pp.pointEdges()[patchPointi];

            forAll(pEdges, i)
            {
                label edgeI = pEdges[i];

                label otherPointi = pp.edges()[edgeI].otherVertex(patchPointi);

                if (extrudeStatus[otherPointi] != NOEXTRUDE)
                {
                    avg += localPoints[otherPointi] + patchDisp[otherPointi];
                    nPoints++;
                }
            }

            if (nPoints > 0)
            {
                if (debug&meshRefinement::ATTRACTION)
                {
                    Pout<< "Displacement at illegal point "
                        << localPoints[patchPointi]
                        << " set to "
                        << (avg / nPoints - localPoints[patchPointi])
                        << endl;
                }

                patchDisp[patchPointi] =
                    avg / nPoints
                  - localPoints[patchPointi];

                nExtrudeRemove++;
            }
            else
            {
                // All surrounding points are not extruded. Leave patchDisp
                // intact.
            }
        }
    }

    Info<< "Detected " << returnReduce(nNoVisNormal, sumOp<label>())
        << " points with point normal pointing through faces." << nl
        << "Reset displacement at "
        << returnReduce(nExtrudeRemove, sumOp<label>())
        << " points to average of surrounding points." << endl;

    // Make sure displacement is equal on both sides of coupled patches.
    syncPatchDisplacement
    (
        pp,
        minThickness,
        patchDisp,
        patchNLayers,
        extrudeStatus
    );

    Info<< endl;
}


bool Foam::snappyLayerDriver::sameEdgeNeighbour
(
    const labelListList& globalEdgeFaces,
    const label myGlobalFacei,
    const label nbrGlobFacei,
    const label edgeI
) const
{
    const labelList& eFaces = globalEdgeFaces[edgeI];
    if (eFaces.size() == 2)
    {
        return edge(myGlobalFacei, nbrGlobFacei) == edge(eFaces[0], eFaces[1]);
    }
    else
    {
        return false;
    }
}


void Foam::snappyLayerDriver::getVertexString
(
    const indirectPrimitivePatch& pp,
    const labelListList& globalEdgeFaces,
    const label facei,
    const label edgeI,
    const label myGlobFacei,
    const label nbrGlobFacei,
    DynamicList<label>& vertices
) const
{
    const labelList& fEdges = pp.faceEdges()[facei];
    label fp = findIndex(fEdges, edgeI);

    if (fp == -1)
    {
        FatalErrorInFunction
            << "problem." << abort(FatalError);
    }

    // Search back
    label startFp = fp;

    forAll(fEdges, i)
    {
        label prevFp = fEdges.rcIndex(startFp);
        if
        (
           !sameEdgeNeighbour
            (
                globalEdgeFaces,
                myGlobFacei,
                nbrGlobFacei,
                fEdges[prevFp]
            )
        )
        {
            break;
        }
        startFp = prevFp;
    }

    label endFp = fp;
    forAll(fEdges, i)
    {
        label nextFp = fEdges.fcIndex(endFp);
        if
        (
           !sameEdgeNeighbour
            (
                globalEdgeFaces,
                myGlobFacei,
                nbrGlobFacei,
                fEdges[nextFp]
            )
        )
        {
            break;
        }
        endFp = nextFp;
    }

    const face& f = pp.localFaces()[facei];
    vertices.clear();
    fp = startFp;
    while (fp != endFp)
    {
        vertices.append(f[fp]);
        fp = f.fcIndex(fp);
    }
    vertices.append(f[fp]);
    fp = f.fcIndex(fp);
    vertices.append(f[fp]);
}


// Truncates displacement
// - for all patchFaces in the faceset displacement gets set to zero
// - all displacement < minThickness gets set to zero
Foam::label Foam::snappyLayerDriver::truncateDisplacement
(
    const globalIndex& globalFaces,
    const labelListList& edgeGlobalFaces,
    const indirectPrimitivePatch& pp,
    const scalarField& minThickness,
    const faceSet& illegalPatchFaces,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    label nChanged = 0;

    const Map<label>& meshPointMap = pp.meshPointMap();

    forAllConstIter(faceSet, illegalPatchFaces, iter)
    {
        label facei = iter.key();

        if (mesh.isInternalFace(facei))
        {
            FatalErrorInFunction
                << "Faceset " << illegalPatchFaces.name()
                << " contains internal face " << facei << nl
                << "It should only contain patch faces" << abort(FatalError);
        }

        const face& f = mesh.faces()[facei];


        forAll(f, fp)
        {
            if (meshPointMap.found(f[fp]))
            {
                label patchPointi = meshPointMap[f[fp]];

                if (extrudeStatus[patchPointi] != NOEXTRUDE)
                {
                    unmarkExtrusion
                    (
                        patchPointi,
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    );
                    nChanged++;
                }
            }
        }
    }

    forAll(patchDisp, patchPointi)
    {
        if (mag(patchDisp[patchPointi]) < minThickness[patchPointi])
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
                nChanged++;
            }
        }
        else if (extrudeStatus[patchPointi] == NOEXTRUDE)
        {
            // Make sure displacement is 0. Should already be so but ...
            patchDisp[patchPointi] = Zero;
            patchNLayers[patchPointi] = 0;
        }
    }


    const faceList& localFaces = pp.localFaces();

    while (true)
    {
        syncPatchDisplacement
        (
            pp,
            minThickness,
            patchDisp,
            patchNLayers,
            extrudeStatus
        );


        // Pinch
        // ~~~~~

        // Make sure that a face doesn't have two non-consecutive areas
        // not extruded (e.g. quad where vertex 0 and 2 are not extruded
        // but 1 and 3 are) since this gives topological errors.

        label nPinched = 0;

        forAll(localFaces, i)
        {
            const face& localF = localFaces[i];

            // Count number of transitions from unsnapped to snapped.
            label nTrans = 0;

            extrudeMode prevMode = extrudeStatus[localF.prevLabel(0)];

            forAll(localF, fp)
            {
                extrudeMode fpMode = extrudeStatus[localF[fp]];

                if (prevMode == NOEXTRUDE && fpMode != NOEXTRUDE)
                {
                    nTrans++;
                }
                prevMode = fpMode;
            }

            if (nTrans > 1)
            {
                // Multiple pinches. Reset whole face as unextruded.
                if
                (
                    unmarkExtrusion
                    (
                        localF,
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    )
                )
                {
                    nPinched++;
                    nChanged++;
                }
            }
        }

        reduce(nPinched, sumOp<label>());

        Info<< "truncateDisplacement : Unextruded " << nPinched
            << " faces due to non-consecutive vertices being extruded." << endl;


        // Butterfly
        // ~~~~~~~~~

        // Make sure that a string of edges becomes a single face so
        // not a butterfly. Occasionally an 'edge' will have a single dangling
        // vertex due to face combining. These get extruded as a single face
        // (with a dangling vertex) so make sure this extrusion forms a single
        // shape.
        //  - continuous i.e. no butterfly:
        //      +     +
        //      |\   /|
        //      | \ / |
        //      +--+--+
        //  - extrudes from all but the endpoints i.e. no partial
        //    extrude
        //            +
        //           /|
        //          / |
        //      +--+--+
        // The common error topology is a pinch somewhere in the middle
        label nButterFly = 0;
        {
            DynamicList<label> stringedVerts;
            forAll(pp.edges(), edgeI)
            {
                const labelList& globFaces = edgeGlobalFaces[edgeI];

                if (globFaces.size() == 2)
                {
                    label myFacei = pp.edgeFaces()[edgeI][0];
                    label myGlobalFacei = globalFaces.toGlobal
                    (
                        pp.addressing()[myFacei]
                    );
                    label nbrGlobalFacei =
                    (
                        globFaces[0] != myGlobalFacei
                      ? globFaces[0]
                      : globFaces[1]
                    );
                    getVertexString
                    (
                        pp,
                        edgeGlobalFaces,
                        myFacei,
                        edgeI,
                        myGlobalFacei,
                        nbrGlobalFacei,
                        stringedVerts
                    );

                    if
                    (
                        extrudeStatus[stringedVerts[0]] != NOEXTRUDE
                     || extrudeStatus[stringedVerts.last()] != NOEXTRUDE
                    )
                    {
                        // Any pinch in the middle
                        bool pinch = false;
                        for (label i = 1; i < stringedVerts.size()-1; i++)
                        {
                            if
                            (
                                extrudeStatus[stringedVerts[i]] == NOEXTRUDE
                            )
                            {
                                pinch = true;
                                break;
                            }
                        }
                        if (pinch)
                        {
                            forAll(stringedVerts, i)
                            {
                                if
                                (
                                    unmarkExtrusion
                                    (
                                        stringedVerts[i],
                                        patchDisp,
                                        patchNLayers,
                                        extrudeStatus
                                    )
                                )
                                {
                                    nButterFly++;
                                    nChanged++;
                                }
                            }
                        }
                    }
                }
            }
        }

        reduce(nButterFly, sumOp<label>());

        Info<< "truncateDisplacement : Unextruded " << nButterFly
            << " faces due to stringed edges with inconsistent extrusion."
            << endl;



        // Consistent number of layers
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Make sure that a face has consistent number of layers for all
        // its vertices.

        label nDiffering = 0;

        // forAll(localFaces, i)
        //{
        //    const face& localF = localFaces[i];
        //
        //    label numLayers = -1;
        //
        //    forAll(localF, fp)
        //    {
        //        if (patchNLayers[localF[fp]] > 0)
        //        {
        //            if (numLayers == -1)
        //            {
        //                numLayers = patchNLayers[localF[fp]];
        //            }
        //            else if (numLayers != patchNLayers[localF[fp]])
        //            {
        //                // Differing number of layers
        //                if
        //                (
        //                    unmarkExtrusion
        //                    (
        //                        localF,
        //                        patchDisp,
        //                        patchNLayers,
        //                        extrudeStatus
        //                    )
        //                )
        //                {
        //                    nDiffering++;
        //                    nChanged++;
        //                }
        //                break;
        //            }
        //        }
        //    }
        //}
        //
        // reduce(nDiffering, sumOp<label>());
        //
        // Info<< "truncateDisplacement : Unextruded " << nDiffering
        //    << " faces due to having differing number of layers." << endl;

        if (nPinched+nButterFly+nDiffering == 0)
        {
            break;
        }
    }

    return nChanged;
}


// Setup layer information (at points and faces) to modify mesh topology in
// regions where layer mesh terminates.
void Foam::snappyLayerDriver::setupLayerInfoTruncation
(
    const indirectPrimitivePatch& pp,
    const labelList& patchNLayers,
    const List<extrudeMode>& extrudeStatus,
    const label nBufferCellsNoExtrude,
    labelList& nPatchPointLayers,
    labelList& nPatchFaceLayers
) const
{
    Info<< nl << "Setting up information for layer truncation ..." << endl;

    const fvMesh& mesh = meshRefiner_.mesh();

    if (nBufferCellsNoExtrude < 0)
    {
        Info<< nl << "Performing no layer truncation."
            << " nBufferCellsNoExtrude set to less than 0  ..." << endl;

        // Face layers if any point gets extruded
        forAll(pp.localFaces(), patchFacei)
        {
            const face& f = pp.localFaces()[patchFacei];

            forAll(f, fp)
            {
                if (patchNLayers[f[fp]] > 0)
                {
                    nPatchFaceLayers[patchFacei] = patchNLayers[f[fp]];
                    break;
                }
            }
        }
        nPatchPointLayers = patchNLayers;

        // Set any unset patch face layers
        forAll(nPatchFaceLayers, patchFacei)
        {
            if (nPatchFaceLayers[patchFacei] == -1)
            {
                nPatchFaceLayers[patchFacei] = 0;
            }
        }
    }
    else
    {
        // Determine max point layers per face.
        labelList maxLevel(pp.size(), 0);

        forAll(pp.localFaces(), patchFacei)
        {
            const face& f = pp.localFaces()[patchFacei];

            // find patch faces where layer terminates (i.e contains extrude
            // and noextrude points).

            bool noExtrude = false;
            label mLevel = 0;

            forAll(f, fp)
            {
                if (extrudeStatus[f[fp]] == NOEXTRUDE)
                {
                    noExtrude = true;
                }
                mLevel = max(mLevel, patchNLayers[f[fp]]);
            }

            if (mLevel > 0)
            {
                // So one of the points is extruded. Check if all are extruded
                // or is a mix.

                if (noExtrude)
                {
                    nPatchFaceLayers[patchFacei] = 1;
                    maxLevel[patchFacei] = mLevel;
                }
                else
                {
                    maxLevel[patchFacei] = mLevel;
                }
            }
        }

        // We have the seed faces (faces with nPatchFaceLayers != maxLevel)
        // Now do a meshwave across the patch where we pick up neighbours
        // of seed faces.
        // Note: quite inefficient. Could probably be coded better.

        const labelListList& pointFaces = pp.pointFaces();

        label nLevels = gMax(patchNLayers);

        // flag neighbouring patch faces with number of layers to grow
        for (label ilevel = 1; ilevel < nLevels; ilevel++)
        {
            label nBuffer;

            if (ilevel == 1)
            {
                nBuffer = nBufferCellsNoExtrude - 1;
            }
            else
            {
                nBuffer = nBufferCellsNoExtrude;
            }

            for (label ibuffer = 0; ibuffer < nBuffer + 1; ibuffer++)
            {
                labelList tempCounter(nPatchFaceLayers);

                boolList foundNeighbour(pp.nPoints(), false);

                forAll(pp.meshPoints(), patchPointi)
                {
                    forAll(pointFaces[patchPointi], pointFacei)
                    {
                        label facei = pointFaces[patchPointi][pointFacei];

                        if
                        (
                            nPatchFaceLayers[facei] != -1
                         && maxLevel[facei] > 0
                        )
                        {
                            foundNeighbour[patchPointi] = true;
                            break;
                        }
                    }
                }

                syncTools::syncPointList
                (
                    mesh,
                    pp.meshPoints(),
                    foundNeighbour,
                    orEqOp<bool>(),
                    false               // null value
                );

                forAll(pp.meshPoints(), patchPointi)
                {
                    if (foundNeighbour[patchPointi])
                    {
                        forAll(pointFaces[patchPointi], pointFacei)
                        {
                            label facei = pointFaces[patchPointi][pointFacei];
                            if
                            (
                                nPatchFaceLayers[facei] == -1
                             && maxLevel[facei] > 0
                             && ilevel < maxLevel[facei]
                            )
                            {
                                tempCounter[facei] = ilevel;
                            }
                        }
                    }
                }
                nPatchFaceLayers = tempCounter;
            }
        }

        forAll(pp.localFaces(), patchFacei)
        {
            if (nPatchFaceLayers[patchFacei] == -1)
            {
                nPatchFaceLayers[patchFacei] = maxLevel[patchFacei];
            }
        }

        forAll(pp.meshPoints(), patchPointi)
        {
            if (extrudeStatus[patchPointi] != NOEXTRUDE)
            {
                forAll(pointFaces[patchPointi], pointFacei)
                {
                    label face = pointFaces[patchPointi][pointFacei];
                    nPatchPointLayers[patchPointi] = max
                    (
                        nPatchPointLayers[patchPointi],
                        nPatchFaceLayers[face]
                    );
                }
            }
            else
            {
                nPatchPointLayers[patchPointi] = 0;
            }
        }
        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            nPatchPointLayers,
            maxEqOp<label>(),
            label(0)        // null value
        );
    }
}


// Does any of the cells use a face from faces?
bool Foam::snappyLayerDriver::cellsUseFace
(
    const polyMesh& mesh,
    const labelList& cellLabels,
    const labelHashSet& faces
)
{
    forAll(cellLabels, i)
    {
        const cell& cFaces = mesh.cells()[cellLabels[i]];

        forAll(cFaces, cFacei)
        {
            if (faces.found(cFaces[cFacei]))
            {
                return true;
            }
        }
    }
    return false;
}


// Checks the newly added cells and locally unmarks points so they
// will not get extruded next time round. Returns global number of unmarked
// points (0 if all was fine)
Foam::label Foam::snappyLayerDriver::checkAndUnmark
(
    const addPatchCellLayer& addLayer,
    const dictionary& meshQualityDict,
    const bool additionalReporting,
    const List<labelPair>& baffles,
    const indirectPrimitivePatch& pp,
    const fvMesh& newMesh,

    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
)
{
    // Check the resulting mesh for errors
    Info<< nl << "Checking mesh with layer ..." << endl;
    faceSet wrongFaces(newMesh, "wrongFaces", newMesh.nFaces()/1000);
    motionSmoother::checkMesh
    (
        false,
        newMesh,
        meshQualityDict,
        identity(newMesh.nFaces()),
        baffles,
        wrongFaces
    );
    Info<< "Detected " << returnReduce(wrongFaces.size(), sumOp<label>())
        << " illegal faces"
        << " (concave, zero area or negative cell pyramid volume)"
        << endl;

    // Undo local extrusion if
    // - any of the added cells in error

    label nChanged = 0;

    // Get all cells in the layer.
    labelListList addedCells
    (
        addPatchCellLayer::addedCells
        (
            newMesh,
            addLayer.layerFaces()
        )
    );

    // Check if any of the faces in error uses any face of an added cell
    // - if additionalReporting print the few remaining areas for ease of
    //   finding out where the problems are.

    const label nReportMax = 10;
    DynamicField<point> disabledFaceCentres(nReportMax);

    forAll(addedCells, oldPatchFacei)
    {
        // Get the cells (in newMesh labels) per old patch face (in mesh
        // labels)
        const labelList& fCells = addedCells[oldPatchFacei];

        if (cellsUseFace(newMesh, fCells, wrongFaces))
        {
            // Unmark points on old mesh
            if
            (
                unmarkExtrusion
                (
                    pp.localFaces()[oldPatchFacei],
                    patchDisp,
                    patchNLayers,
                    extrudeStatus
                )
            )
            {
                if (additionalReporting && (nChanged < nReportMax))
                {
                    disabledFaceCentres.append
                    (
                        pp.faceCentres()[oldPatchFacei]
                    );
                }

                nChanged++;
            }
        }
    }


    label nChangedTotal = returnReduce(nChanged, sumOp<label>());

    if (additionalReporting)
    {
        // Limit the number of points to be printed so that
        // not too many points are reported when running in parallel
        // Not accurate, i.e. not always nReportMax points are written,
        // but this estimation avoid some communication here.
        // The important thing, however, is that when only a few faces
        // are disabled, their coordinates are printed, and this should be
        // the case
        label nReportLocal = nChanged;
        if (nChangedTotal > nReportMax)
        {
            nReportLocal = min
            (
                max(nChangedTotal / Pstream::nProcs(), 1),
                min
                (
                    nChanged,
                    max(nReportMax / Pstream::nProcs(), 1)
                )
            );
        }

        if (nReportLocal)
        {
            Pout<< "Checked mesh with layers. Disabled extrusion at " << endl;
            for (label i=0; i < nReportLocal; i++)
            {
                Pout<< "    " << disabledFaceCentres[i] << endl;
            }
        }

        label nReportTotal = returnReduce(nReportLocal, sumOp<label>());

        if (nReportTotal < nChangedTotal)
        {
            Info<< "Suppressed disabled extrusion message for other "
                << nChangedTotal - nReportTotal << " faces." << endl;
        }
    }

    return nChangedTotal;
}


Foam::label Foam::snappyLayerDriver::countExtrusion
(
    const indirectPrimitivePatch& pp,
    const List<extrudeMode>& extrudeStatus
)
{
    // Count number of extruded patch faces
    label nExtruded = 0;
    {
        const faceList& localFaces = pp.localFaces();

        forAll(localFaces, i)
        {
            const face& localFace = localFaces[i];

            forAll(localFace, fp)
            {
                if (extrudeStatus[localFace[fp]] != NOEXTRUDE)
                {
                    nExtruded++;
                    break;
                }
            }
        }
    }

    return returnReduce(nExtruded, sumOp<label>());
}


// Collect layer faces and layer cells into mesh fields for ease of handling
void Foam::snappyLayerDriver::getLayerCellsFaces
(
    const polyMesh& mesh,
    const addPatchCellLayer& addLayer,
    const scalarField& oldRealThickness,

    labelList& cellNLayers,
    scalarField& faceRealThickness
)
{
    cellNLayers.setSize(mesh.nCells());
    cellNLayers = 0;
    faceRealThickness.setSize(mesh.nFaces());
    faceRealThickness = 0;

    // Mark all faces in the layer
    const labelListList& layerFaces = addLayer.layerFaces();

    // Mark all cells in the layer.
    labelListList addedCells(addPatchCellLayer::addedCells(mesh, layerFaces));

    forAll(addedCells, oldPatchFacei)
    {
        const labelList& added = addedCells[oldPatchFacei];

        const labelList& layer = layerFaces[oldPatchFacei];

        if (layer.size())
        {
            forAll(added, i)
            {
                cellNLayers[added[i]] = layer.size()-1;
            }
        }
    }

    forAll(layerFaces, oldPatchFacei)
    {
        const labelList& layer = layerFaces[oldPatchFacei];
        const scalar realThickness = oldRealThickness[oldPatchFacei];

        if (layer.size())
        {
            // Layer contains both original boundary face and new boundary
            // face so is nLayers+1. Leave out old internal face.
            for (label i = 1; i < layer.size(); i++)
            {
                faceRealThickness[layer[i]] = realThickness;
            }
        }
    }
}


void Foam::snappyLayerDriver::printLayerData
(
    const fvMesh& mesh,
    const labelList& patchIDs,
    const labelList& cellNLayers,
    const scalarField& faceWantedThickness,
    const scalarField& faceRealThickness
) const
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    int oldPrecision = Info().precision();

    // Find maximum length of a patch name, for a nicer output
    label maxPatchNameLen = 0;
    forAll(patchIDs, i)
    {
        label patchi = patchIDs[i];
        word patchName = pbm[patchi].name();
        maxPatchNameLen = max(maxPatchNameLen, label(patchName.size()));
    }

    Info<< nl
        << setf(ios_base::left) << setw(maxPatchNameLen) << "patch"
        << setw(0) << " faces    layers   overall thickness" << nl
        << setf(ios_base::left) << setw(maxPatchNameLen) << " "
        << setw(0) << "                   [m]       [%]" << nl
        << setf(ios_base::left) << setw(maxPatchNameLen) << "-----"
        << setw(0) << " -----    ------   ---       ---" << endl;


    forAll(patchIDs, i)
    {
        label patchi = patchIDs[i];
        const polyPatch& pp = pbm[patchi];

        label sumSize = pp.size();

        // Number of layers
        const labelList& faceCells = pp.faceCells();
        label sumNLayers = 0;
        forAll(faceCells, i)
        {
            sumNLayers += cellNLayers[faceCells[i]];
        }

        // Thickness
        scalarField::subField patchWanted = pbm[patchi].patchSlice
        (
            faceWantedThickness
        );
        scalarField::subField patchReal = pbm[patchi].patchSlice
        (
            faceRealThickness
        );

        scalar sumRealThickness = sum(patchReal);
        scalar sumFraction = 0;
        forAll(patchReal, i)
        {
            if (patchWanted[i] > vSmall)
            {
                sumFraction += (patchReal[i]/patchWanted[i]);
            }
        }


        reduce(sumSize, sumOp<label>());
        reduce(sumNLayers, sumOp<label>());
        reduce(sumRealThickness, sumOp<scalar>());
        reduce(sumFraction, sumOp<scalar>());


        scalar avgLayers = 0;
        scalar avgReal = 0;
        scalar avgFraction = 0;
        if (sumSize > 0)
        {
            avgLayers = scalar(sumNLayers)/sumSize;
            avgReal = sumRealThickness/sumSize;
            avgFraction = sumFraction/sumSize;
        }

        Info<< setf(ios_base::left) << setw(maxPatchNameLen)
            << pbm[patchi].name() << setprecision(3)
            << " " << setw(8) << sumSize
            << " " << setw(8) << avgLayers
            << " " << setw(8) << avgReal
            << "  " << setw(8) << 100*avgFraction
            << endl;
    }
    Info<< setprecision(oldPrecision) << endl;
}


bool Foam::snappyLayerDriver::writeLayerData
(
    const fvMesh& mesh,
    const labelList& patchIDs,
    const labelList& cellNLayers,
    const scalarField& faceWantedThickness,
    const scalarField& faceRealThickness
) const
{
    bool allOk = true;

    if (meshRefinement::writeLevel() & meshRefinement::WRITELAYERSETS)
    {
        {
            label nAdded = 0;
            forAll(cellNLayers, celli)
            {
                if (cellNLayers[celli] > 0)
                {
                    nAdded++;
                }
            }
            cellSet addedCellSet(mesh, "addedCells", nAdded);
            forAll(cellNLayers, celli)
            {
                if (cellNLayers[celli] > 0)
                {
                    addedCellSet.insert(celli);
                }
            }
            addedCellSet.instance() = meshRefiner_.timeName();
            Info<< "Writing "
                << returnReduce(addedCellSet.size(), sumOp<label>())
                << " added cells to cellSet "
                << addedCellSet.name() << endl;
            bool ok = addedCellSet.write();
            allOk = allOk & ok;
        }
        {
            label nAdded = 0;
            for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
            {
                if (faceRealThickness[facei] > 0)
                {
                    nAdded++;
                }
            }

            faceSet layerFacesSet(mesh, "layerFaces", nAdded);
            for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
            {
                if (faceRealThickness[facei] > 0)
                {
                    layerFacesSet.insert(facei);
                }
            }
            layerFacesSet.instance() = meshRefiner_.timeName();
            Info<< "Writing "
                << returnReduce(layerFacesSet.size(), sumOp<label>())
                << " faces inside added layer to faceSet "
                << layerFacesSet.name() << endl;
            bool ok = layerFacesSet.write();
            allOk = allOk & ok;
        }
    }

    if (meshRefinement::writeLevel() & meshRefinement::WRITELAYERFIELDS)
    {
        Info<< nl << "Writing fields with layer information:" << incrIndent
            << endl;
        {
            volScalarField fld
            (
                IOobject
                (
                    "nSurfaceLayers",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false
                ),
                mesh,
                dimensionedScalar("zero", dimless, 0),
                fixedValueFvPatchScalarField::typeName
            );
            const polyBoundaryMesh& pbm = mesh.boundaryMesh();

            volScalarField::Boundary& fldBf =
                fld.boundaryFieldRef();

            forAll(patchIDs, i)
            {
                label patchi = patchIDs[i];
                const polyPatch& pp = pbm[patchi];
                const labelList& faceCells = pp.faceCells();
                scalarField pfld(faceCells.size());
                forAll(faceCells, i)
                {
                    pfld[i] = cellNLayers[faceCells[i]];
                }
                fldBf[patchi] == pfld;
            }
            Info<< indent << fld.name() << "    : actual number of layers"
                << endl;
            bool ok = fld.write();
            allOk = allOk & ok;
        }
        {
            volScalarField fld
            (
                IOobject
                (
                    "thickness",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false
                ),
                mesh,
                dimensionedScalar("zero", dimless, 0),
                fixedValueFvPatchScalarField::typeName
            );

            const polyBoundaryMesh& pbm = mesh.boundaryMesh();

            volScalarField::Boundary& fldBf =
                fld.boundaryFieldRef();

            forAll(patchIDs, i)
            {
                label patchi = patchIDs[i];
                fldBf[patchi] == pbm[patchi].patchSlice
                (
                    faceRealThickness
                );
            }
            Info<< indent << fld.name() << "         : overall layer thickness"
                << endl;
            bool ok = fld.write();
            allOk = allOk & ok;
        }
        {
            volScalarField fld
            (
                IOobject
                (
                    "thicknessFraction",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false
                ),
                mesh,
                dimensionedScalar("zero", dimless, 0),
                fixedValueFvPatchScalarField::typeName
            );

            const polyBoundaryMesh& pbm = mesh.boundaryMesh();

            volScalarField::Boundary& fldBf =
                fld.boundaryFieldRef();

            forAll(patchIDs, i)
            {
                label patchi = patchIDs[i];

                scalarField::subField patchWanted = pbm[patchi].patchSlice
                (
                    faceWantedThickness
                );
                scalarField::subField patchReal = pbm[patchi].patchSlice
                (
                    faceRealThickness
                );

                // Convert patchReal to relavtive thickness
                scalarField pfld(patchReal.size(), 0.0);
                forAll(patchReal, i)
                {
                    if (patchWanted[i] > vSmall)
                    {
                        pfld[i] = patchReal[i]/patchWanted[i];
                    }
                }

                fldBf[patchi] == pfld;
            }
            Info<< indent << fld.name()
                << " : overall layer thickness (fraction"
                << " of desired thickness)" << endl;
            bool ok = fld.write();
            allOk = allOk & ok;
        }
        Info<< decrIndent<< endl;
    }

    // if (meshRefinement::outputLevel() & meshRefinement::OUTPUTLAYERINFO)
    {
        printLayerData
        (
            mesh,
            patchIDs,
            cellNLayers,
            faceWantedThickness,
            faceRealThickness
        );
    }

    return allOk;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::snappyLayerDriver::snappyLayerDriver
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

void Foam::snappyLayerDriver::mergePatchFacesUndo
(
    const layerParameters& layerParams,
    const dictionary& motionDict
)
{
    // Clip to 30 degrees. Not helpful!
    // scalar planarAngle = min(30.0, layerParams.featureAngle());
    scalar planarAngle = layerParams.featureAngle();
    scalar minCos = Foam::cos(degToRad(planarAngle));

    scalar concaveCos = Foam::cos(degToRad(layerParams.concaveAngle()));

    Info<< nl
        << "Merging all faces of a cell" << nl
        << "---------------------------" << nl
        << "    - which are on the same patch" << nl
        << "    - which make an angle < " << planarAngle
        << " degrees"
        << nl
        << "      (cos:" << minCos << ')' << nl
        << "    - as long as the resulting face doesn't become concave"
        << " by more than "
        << layerParams.concaveAngle() << " degrees" << nl
        << "      (0=straight, 180=fully concave)" << nl
        << endl;

    const fvMesh& mesh = meshRefiner_.mesh();

    List<labelPair> couples(localPointRegion::findDuplicateFacePairs(mesh));

    labelList duplicateFace(mesh.nFaces(), -1);
    forAll(couples, i)
    {
        const labelPair& cpl = couples[i];
        duplicateFace[cpl[0]] = cpl[1];
        duplicateFace[cpl[1]] = cpl[0];
    }

    label nChanged = meshRefiner_.mergePatchFacesUndo
    (
        minCos,
        concaveCos,
        meshRefiner_.meshedPatches(),
        motionDict,
        duplicateFace
    );

    nChanged += meshRefiner_.mergeEdgesUndo(minCos, motionDict);
}


void Foam::snappyLayerDriver::addLayers
(
    const layerParameters& layerParams,
    const dictionary& motionDict,
    const labelList& patchIDs,
    const label nAllowableErrors,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    // Create baffles (pairs of faces that share the same points)
    // Baffles stored as owner and neighbour face that have been created.
    List<labelPair> baffles;
    meshRefiner_.createZoneBaffles
    (
        globalToMasterPatch_,
        globalToSlavePatch_,
        baffles
    );

    if (debug&meshRefinement::MESH)
    {
        const_cast<Time&>(mesh.time())++;
        Info<< "Writing baffled mesh to time "
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
    }


    autoPtr<indirectPrimitivePatch> pp
    (
        meshRefinement::makePatch
        (
            mesh,
            patchIDs
        )
    );


    // Global face indices engine
    const globalIndex globalFaces(mesh.nFaces());

    // Determine extrudePatch.edgeFaces in global numbering (so across
    // coupled patches). This is used only to string up edges between coupled
    // faces (all edges between same (global)face indices get extruded).
    labelListList edgeGlobalFaces
    (
        addPatchCellLayer::globalEdgeFaces
        (
            mesh,
            globalFaces,
            pp
        )
    );

    // Determine patches for extruded boundary edges. Calculates if any
    // additional processor patches need to be constructed.

    labelList sidePatchID;
    determineSidePatches
    (
        globalFaces,
        edgeGlobalFaces,
        pp,

        sidePatchID
    );


    // Point-wise extrusion data
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    // Displacement for all pp.localPoints.
    vectorField patchDisp(pp().nPoints(), vector::one);

    // Number of layers for all pp.localPoints. Note: only valid if
    // extrudeStatus = EXTRUDE.
    labelList patchNLayers(pp().nPoints(), 0);

    // Ideal number of cells added
    label nIdealTotAddedCells = 0;

    // Whether to add edge for all pp.localPoints.
    List<extrudeMode> extrudeStatus(pp().nPoints(), EXTRUDE);


    {
        // Get number of layer per point from number of layers per patch
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        setNumLayers
        (
            layerParams.numLayers(),    // per patch the num layers
            patchIDs,                   // patches that are being moved
            pp,                         // indirectpatch for all faces moving

            patchDisp,
            patchNLayers,
            extrudeStatus,
            nIdealTotAddedCells
        );

        // Precalculate mesh edge labels for patch edges
        labelList meshEdges(pp().meshEdges(mesh.edges(), mesh.pointEdges()));

        // Disable extrusion on non-manifold points
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        handleNonManifolds
        (
            pp,
            meshEdges,
            edgeGlobalFaces,

            patchDisp,
            patchNLayers,
            extrudeStatus
        );

        // Disable extrusion on feature angles
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        handleFeatureAngle
        (
            pp,
            meshEdges,
            degToRad(layerParams.featureAngle()),

            patchDisp,
            patchNLayers,
            extrudeStatus
        );

        // Disable extrusion on warped faces
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // It is hard to calculate some length scale if not in relative
        // mode so disable this check.
        if (layerParams.relativeSizes())
        {
            // Undistorted edge length
            const scalar edge0Len =
                meshRefiner_.meshCutter().level0EdgeLength();
            const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();

            handleWarpedFaces
            (
                pp,
                layerParams.maxFaceThicknessRatio(),
                edge0Len,
                cellLevel,

                patchDisp,
                patchNLayers,
                extrudeStatus
            );
        }

        //// Disable extrusion on cells with multiple patch faces
        //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //
        // handleMultiplePatchFaces
        //(
        //    pp,
        //
        //    patchDisp,
        //    patchNLayers,
        //    extrudeStatus
        //);


        // Grow out region of non-extrusion
        for (label i = 0; i < layerParams.nGrow(); i++)
        {
            growNoExtrusion
            (
                pp,
                patchDisp,
                patchNLayers,
                extrudeStatus
            );
        }
    }


    // Undistorted edge length
    const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();

    // Determine (wanted) point-wise overall layer thickness and expansion
    // ratio
    scalarField thickness(pp().nPoints());
    scalarIOField minThickness
    (
        IOobject
        (
            "minThickness",
            meshRefiner_.timeName(),
            mesh,
            IOobject::NO_READ
        ),
        pp().nPoints()
    );
    scalarField expansionRatio(pp().nPoints());
    calculateLayerThickness
    (
        pp,
        patchIDs,
        layerParams,
        cellLevel,
        patchNLayers,
        edge0Len,

        thickness,
        minThickness,
        expansionRatio
    );



    // Overall displacement field
    pointVectorField displacement
    (
        makeLayerDisplacementField
        (
            pointMesh::New(mesh),
            layerParams.numLayers()
        )
    );

    // Allocate run-time selectable mesh mover
    autoPtr<externalDisplacementMeshMover> medialAxisMoverPtr;
    {
        // Set up controls for meshMover
        dictionary combinedDict(layerParams.dict());
        // Add mesh quality constraints
        combinedDict.merge(motionDict);
        // Where to get minThickness from
        combinedDict.add("minThicknessName", minThickness.name());

        // Take over patchDisp as boundary conditions on displacement
        // pointVectorField
        medialAxisMoverPtr = externalDisplacementMeshMover::New
        (
            layerParams.meshShrinker(),
            combinedDict,
            baffles,
            displacement
        );
    }


    // Saved old points
    pointField oldPoints(mesh.points());

    // Current set of topology changes. (changing mesh clears out
    // polyTopoChange)
    polyTopoChange savedMeshMod(mesh.boundaryMesh().size());
    // Per cell 0 or number of layers in the cell column it is part of
    labelList cellNLayers;
    // Per face actual overall layer thickness
    scalarField faceRealThickness;
    // Per face wanted overall layer thickness
    scalarField faceWantedThickness(mesh.nFaces(), 0.0);
    {
        UIndirectList<scalar>(faceWantedThickness, pp().addressing()) =
            avgPointData(pp, thickness);
    }


    for (label iteration = 0; iteration < layerParams.nLayerIter(); iteration++)
    {
        Info<< nl
            << "Layer addition iteration " << iteration << nl
            << "--------------------------" << endl;


        // Unset the extrusion at the pp.
        const dictionary& meshQualityDict =
        (
            iteration < layerParams.nRelaxedIter()
          ? motionDict
          : motionDict.subDict("relaxed")
        );

        if (iteration >= layerParams.nRelaxedIter())
        {
            Info<< "Switched to relaxed meshQuality constraints." << endl;
        }



        // Make sure displacement is equal on both sides of coupled patches.
        syncPatchDisplacement
        (
            pp,
            minThickness,
            patchDisp,
            patchNLayers,
            extrudeStatus
        );

        // Displacement acc. to pointnormals
        getPatchDisplacement
        (
            pp,
            thickness,
            minThickness,
            patchDisp,
            patchNLayers,
            extrudeStatus
        );

        // Shrink mesh by displacement value first.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        {
            pointField oldPatchPos(pp().localPoints());

            // Take over patchDisp into pointDisplacement field and
            // adjust both for multi-patch constraints
            motionSmootherAlgo::setDisplacement
            (
                patchIDs,
                pp,
                patchDisp,
                displacement
            );


            // Move mesh
            // ~~~~~~~~~

            // Set up controls for meshMover
            dictionary combinedDict(layerParams.dict());
            // Add standard quality constraints
            combinedDict.merge(motionDict);
            // Add relaxed constraints (overrides standard ones)
            combinedDict.merge(meshQualityDict);
            // Where to get minThickness from
            combinedDict.add("minThicknessName", minThickness.name());

            labelList checkFaces(identity(mesh.nFaces()));
            medialAxisMoverPtr().move
            (
                combinedDict,
                nAllowableErrors,
                checkFaces
            );

            pp().movePoints(mesh.points());

            // Update patchDisp (since not all might have been honoured)
            patchDisp = oldPatchPos - pp().localPoints();
        }

        // Truncate displacements that are too small (this will do internal
        // ones, coupled ones have already been truncated by
        // syncPatchDisplacement)
        faceSet dummySet(mesh, "wrongPatchFaces", 0);
        truncateDisplacement
        (
            globalFaces,
            edgeGlobalFaces,
            pp,
            minThickness,
            dummySet,
            patchDisp,
            patchNLayers,
            extrudeStatus
        );


        // Dump to .obj file for debugging.
        if (debug&meshRefinement::MESH || debug&meshRefinement::LAYERINFO)
        {
            dumpDisplacement
            (
                mesh.time().path()/"layer_" + meshRefiner_.timeName(),
                pp(),
                patchDisp,
                extrudeStatus
            );

            const_cast<Time&>(mesh.time())++;
            Info<< "Writing shrunk mesh to time "
                << meshRefiner_.timeName() << endl;

            // See comment in snappySnapDriver why we should not remove meshPhi
            // using mesh.clearOut().

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
        }


        // Mesh topo change engine
        polyTopoChange meshMod(mesh);

        // Grow layer of cells on to patch. Handles zero sized displacement.
        addPatchCellLayer addLayer(mesh);

        // Determine per point/per face number of layers to extrude. Also
        // handles the slow termination of layers when going switching layers

        labelList nPatchPointLayers(pp().nPoints(), -1);
        labelList nPatchFaceLayers(pp().size(), -1);
        setupLayerInfoTruncation
        (
            pp,
            patchNLayers,
            extrudeStatus,
            layerParams.nBufferCellsNoExtrude(),
            nPatchPointLayers,
            nPatchFaceLayers
        );

        // Calculate displacement for final layer for addPatchLayer.
        // (layer of cells next to the original mesh)
        vectorField finalDisp(patchNLayers.size(), Zero);

        forAll(nPatchPointLayers, i)
        {
            scalar ratio = layerParams.finalLayerThicknessRatio
            (
                nPatchPointLayers[i],
                expansionRatio[i]
            );
            finalDisp[i] = ratio*patchDisp[i];
        }


        const scalarField invExpansionRatio(1.0/expansionRatio);

        // Add topo regardless of whether extrudeStatus is extruderemove.
        // Not add layer if patchDisp is zero.
        addLayer.setRefinement
        (
            globalFaces,
            edgeGlobalFaces,

            invExpansionRatio,
            pp(),
            sidePatchID,        // boundary patch for extruded boundary edges
            labelList(0),       // exposed patchIDs, not used for adding layers
            nPatchFaceLayers,   // layers per face
            nPatchPointLayers,  // layers per point
            finalDisp,          // thickness of layer nearest internal mesh
            meshMod
        );

        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }

        // Store mesh changes for if mesh is correct.
        savedMeshMod = meshMod;


        // With the stored topo changes we create a new mesh so we can
        // undo if necessary.

        autoPtr<fvMesh> newMeshPtr;
        autoPtr<mapPolyMesh> map = meshMod.makeMesh
        (
            newMeshPtr,
            IOobject
            (
                // mesh.name()+"_layer",
                mesh.name(),
                static_cast<polyMesh&>(mesh).instance(),
                mesh.time(),  // register with runTime
                IOobject::NO_READ,
                static_cast<polyMesh&>(mesh).writeOpt()
            ),              // io params from original mesh but new name
            mesh,           // original mesh
            true            // parallel sync
        );
        fvMesh& newMesh = newMeshPtr();

        //?necessary? Update fields
        newMesh.updateMesh(map);

        newMesh.setInstance(meshRefiner_.timeName());

        // Update numbering on addLayer:
        // - cell/point labels to be newMesh.
        // - patchFaces to remain in oldMesh order.
        addLayer.updateMesh
        (
            map,
            identity(pp().size()),
            identity(pp().nPoints())
        );

        // Update numbering of baffles
        List<labelPair> newMeshBaffles(baffles.size());
        forAll(baffles, i)
        {
            const labelPair& p = baffles[i];
            newMeshBaffles[i][0] = map().reverseFaceMap()[p[0]];
            newMeshBaffles[i][1] = map().reverseFaceMap()[p[1]];
        }

        // Collect layer faces and cells for outside loop.
        getLayerCellsFaces
        (
            newMesh,
            addLayer,
            avgPointData(pp, mag(patchDisp))(), // current thickness

            cellNLayers,
            faceRealThickness
        );


        // Count number of added cells
        label nAddedCells = 0;
        forAll(cellNLayers, celli)
        {
            if (cellNLayers[celli] > 0)
            {
                nAddedCells++;
            }
        }


        if (debug&meshRefinement::MESH)
        {
            Info<< "Writing layer mesh to time " << meshRefiner_.timeName()
                << endl;
            newMesh.write();

            cellSet addedCellSet(newMesh, "addedCells", nAddedCells);
            forAll(cellNLayers, celli)
            {
                if (cellNLayers[celli] > 0)
                {
                    addedCellSet.insert(celli);
                }
            }
            addedCellSet.instance() = meshRefiner_.timeName();
            Info<< "Writing "
                << returnReduce(addedCellSet.size(), sumOp<label>())
                << " added cells to cellSet " << addedCellSet.name() << endl;
            addedCellSet.write();

            faceSet layerFacesSet(newMesh, "layerFaces", newMesh.nFaces()/100);
            for (label facei = 0; facei < newMesh.nInternalFaces(); facei++)
            {
                if (faceRealThickness[facei] > 0)
                {
                    layerFacesSet.insert(facei);
                }
            }
            layerFacesSet.instance() = meshRefiner_.timeName();
            Info<< "Writing "
                << returnReduce(layerFacesSet.size(), sumOp<label>())
                << " faces inside added layer to faceSet "
                << layerFacesSet.name() << endl;
            layerFacesSet.write();
        }


        label nTotChanged = checkAndUnmark
        (
            addLayer,
            meshQualityDict,
            layerParams.additionalReporting(),
            newMeshBaffles,
            pp(),
            newMesh,

            patchDisp,
            patchNLayers,
            extrudeStatus
        );

        label nTotExtruded = countExtrusion(pp, extrudeStatus);
        label nTotFaces = returnReduce(pp().size(), sumOp<label>());
        label nTotAddedCells = returnReduce(nAddedCells, sumOp<label>());

        Info<< "Extruding " << nTotExtruded
            << " out of " << nTotFaces
            << " faces (" << 100.0*nTotExtruded/nTotFaces << "%)."
            << " Removed extrusion at " << nTotChanged << " faces."
            << endl
            << "Added " << nTotAddedCells << " out of " << nIdealTotAddedCells
            << " cells (" << 100.0*nTotAddedCells/nIdealTotAddedCells << "%)."
            << endl;

        if (nTotChanged == 0)
        {
            break;
        }

        // Reset mesh points and start again
        mesh.movePoints(oldPoints);
        pp().movePoints(mesh.points());

        // Grow out region of non-extrusion
        for (label i = 0; i < layerParams.nGrow(); i++)
        {
            growNoExtrusion
            (
                pp,
                patchDisp,
                patchNLayers,
                extrudeStatus
            );
        }

        Info<< endl;
    }


    // At this point we have a (shrunk) mesh and a set of topology changes
    // which will make a valid mesh with layer. Apply these changes to the
    // current mesh.

    // Apply the stored topo changes to the current mesh.
    autoPtr<mapPolyMesh> map = savedMeshMod.changeMesh(mesh, false);

    // Hack to remove meshPhi - mapped incorrectly. TBD.
    mesh.clearOut();

    // Update fields
    mesh.updateMesh(map);

    // Move mesh (since morphing does not do this)
    if (map().hasMotionPoints())
    {
        mesh.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes.
        mesh.clearOut();
    }

    // Reset the instance for if in overwrite mode
    mesh.setInstance(meshRefiner_.timeName());

    meshRefiner_.updateMesh(map, labelList(0));

    // Update numbering of faceWantedThickness
    meshRefinement::updateList(map().faceMap(), scalar(0), faceWantedThickness);

    // Update numbering on baffles
    forAll(baffles, i)
    {
        labelPair& p = baffles[i];
        p[0] = map().reverseFaceMap()[p[0]];
        p[1] = map().reverseFaceMap()[p[1]];
    }


    label nBaffles = returnReduce(baffles.size(), sumOp<label>());
    if (nBaffles > 0)
    {
        // Merge any baffles
        Info<< "Converting " << nBaffles
            << " baffles back into zoned faces ..."
            << endl;

        autoPtr<mapPolyMesh> map = meshRefiner_.mergeBaffles(baffles);

        inplaceReorder(map().reverseCellMap(), cellNLayers);
        inplaceReorder(map().reverseFaceMap(), faceWantedThickness);
        inplaceReorder(map().reverseFaceMap(), faceRealThickness);

        Info<< "Converted baffles in = "
            << meshRefiner_.mesh().time().cpuTimeIncrement()
            << " s\n" << nl << endl;
    }


    // Do final balancing
    // ~~~~~~~~~~~~~~~~~~

    if (Pstream::parRun())
    {
        Info<< nl
            << "Doing final balancing" << nl
            << "---------------------" << nl
            << endl;

        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }

        // Balance. No restriction on face zones and baffles.
        autoPtr<mapDistributePolyMesh> map = meshRefiner_.balance
        (
            false,
            false,
            scalarField(mesh.nCells(), 1.0),
            decomposer,
            distributor
        );

        // Re-distribute flag of layer faces and cells
        map().distributeCellData(cellNLayers);
        map().distributeFaceData(faceWantedThickness);
        map().distributeFaceData(faceRealThickness);
    }


    // Write mesh data
    // ~~~~~~~~~~~~~~~

    writeLayerData
    (
        mesh,
        patchIDs,
        cellNLayers,
        faceWantedThickness,
        faceRealThickness
    );
}


void Foam::snappyLayerDriver::doLayers
(
    const dictionary& shrinkDict,
    const dictionary& motionDict,
    const layerParameters& layerParams,
    const bool preBalance,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor
)
{
    const fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl
        << "Shrinking and layer addition phase" << nl
        << "----------------------------------" << nl
        << endl;

    Info<< "Using mesh parameters " << motionDict << nl << endl;

    // Merge coplanar boundary faces
    mergePatchFacesUndo(layerParams, motionDict);

    // Per patch the number of layers (-1 or 0 if no layer)
    const labelList& numLayers = layerParams.numLayers();

    // Patches that need to get a layer
    DynamicList<label> patchIDs(numLayers.size());
    label nFacesWithLayers = 0;
    forAll(numLayers, patchi)
    {
        if (numLayers[patchi] > 0)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchi];

            if (!pp.coupled())
            {
                patchIDs.append(patchi);
                nFacesWithLayers += mesh.boundaryMesh()[patchi].size();
            }
            else
            {
                WarningInFunction
                    << "Ignoring layers on coupled patch " << pp.name()
                    << endl;
            }
        }
    }
    patchIDs.shrink();

    if (returnReduce(nFacesWithLayers, sumOp<label>()) == 0)
    {
        Info<< nl << "No layers to generate ..." << endl;
    }
    else
    {
        // Check that outside of mesh is not multiply connected.
        checkMeshManifold();

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


        // Balance
        if (Pstream::parRun() && preBalance)
        {
            Info<< nl
                << "Doing initial balancing" << nl
                << "-----------------------" << nl
                << endl;

            scalarField cellWeights(mesh.nCells(), 1);
            forAll(numLayers, patchi)
            {
                if (numLayers[patchi] > 0)
                {
                    const polyPatch& pp = mesh.boundaryMesh()[patchi];
                    forAll(pp.faceCells(), i)
                    {
                        cellWeights[pp.faceCells()[i]] += numLayers[patchi];
                    }
                }
            }

            // Balance mesh (and meshRefinement). Restrict faceZones to
            // be on internal faces only since they will be converted into
            // baffles.
            autoPtr<mapDistributePolyMesh> map = meshRefiner_.balance
            (
                true,   // false,    // keepZoneFaces
                false,
                cellWeights,
                decomposer,
                distributor
            );
        }


        // Do all topo changes
        addLayers
        (
            layerParams,
            motionDict,
            patchIDs,
            nInitErrors,
            decomposer,
            distributor
        );
    }
}


// ************************************************************************* //
