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

\*---------------------------------------------------------------------------*/

#include "autoSnapDriver.H"
#include "polyTopoChange.H"
#include "syncTools.H"
#include "fvMesh.H"
#include "OBJstream.H"
#include "motionSmoother.H"
#include "refinementSurfaces.H"
#include "refinementFeatures.H"
#include "unitConversion.H"
#include "plane.H"
#include "featureEdgeMesh.H"
#include "treeDataPoint.H"
#include "indexedOctree.H"
#include "snapParameters.H"
#include "PatchTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<class T>
    class listPlusEqOp
    {
    public:

        void operator()(List<T>& x, const List<T>& y) const
        {
            label sz = x.size();
            x.setSize(sz+y.size());
            forAll(y, i)
            {
                x[sz++] = y[i];
            }
        }
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::autoSnapDriver::isFeaturePoint
(
    const scalar featureCos,
    const indirectPrimitivePatch& pp,
    const PackedBoolList& isFeatureEdge,
    const label pointI
) const
{
    const pointField& points = pp.localPoints();
    const edgeList& edges = pp.edges();
    const labelList& pEdges = pp.pointEdges()[pointI];

    label nFeatEdges = 0;

    forAll(pEdges, i)
    {
        if (isFeatureEdge[pEdges[i]])
        {
            nFeatEdges++;

            for (label j = i+1; j < pEdges.size(); j++)
            {
                if (isFeatureEdge[pEdges[j]])
                {
                    const edge& eI = edges[pEdges[i]];
                    const edge& eJ = edges[pEdges[j]];

                    const point& p = points[pointI];
                    const point& pI = points[eI.otherVertex(pointI)];
                    const point& pJ = points[eJ.otherVertex(pointI)];

                    vector vI = p-pI;
                    scalar vIMag = mag(vI);

                    vector vJ = pJ-p;
                    scalar vJMag = mag(vJ);

                    if
                    (
                        vIMag > SMALL
                     && vJMag > SMALL
                     && ((vI/vIMag & vJ/vJMag) < featureCos)
                    )
                    {
                        return true;
                    }
                }
            }
        }
    }

    if (nFeatEdges == 1)
    {
        // End of feature-edge string
        return true;
    }

    return false;
}


void Foam::autoSnapDriver::smoothAndConstrain
(
    const PackedBoolList& isPatchMasterEdge,
    const indirectPrimitivePatch& pp,
    const labelList& meshEdges,
    const List<pointConstraint>& constraints,
    vectorField& disp
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    for (label avgIter = 0; avgIter < 20; avgIter++)
    {
        // Calculate average displacement of neighbours
        // - unconstrained (i.e. surface) points use average of all
        //   neighbouring points
        // - from testing it has been observed that it is not beneficial
        //   to have edge constrained points use average of all edge or point
        //   constrained neighbours since they're already attracted to
        //   the nearest point on the feature.
        //   Having them attract to point-constrained neighbours does not
        //   make sense either since there is usually just one of them so
        //   it severely distorts it.
        // - same for feature points. They are already attracted to the
        //   nearest feature point.

        vectorField dispSum(pp.nPoints(), vector::zero);
        labelList dispCount(pp.nPoints(), 0);

        const labelListList& pointEdges = pp.pointEdges();
        const edgeList& edges = pp.edges();

        forAll(pointEdges, pointI)
        {
            const labelList& pEdges = pointEdges[pointI];

            label nConstraints = constraints[pointI].first();

            if (nConstraints <= 1)
            {
                forAll(pEdges, i)
                {
                    label edgeI = pEdges[i];

                    if (isPatchMasterEdge[edgeI])
                    {
                        label nbrPointI = edges[edgeI].otherVertex(pointI);
                        if (constraints[nbrPointI].first() >= nConstraints)
                        {
                            dispSum[pointI] += disp[nbrPointI];
                            dispCount[pointI]++;
                        }
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            dispSum,
            plusEqOp<point>(),
            vector::zero,
            mapDistribute::transform()
        );
        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            dispCount,
            plusEqOp<label>(),
            0,
            mapDistribute::transform()
        );

        // Constraints
        forAll(constraints, pointI)
        {
            if (dispCount[pointI] > 0)
            {
                // Mix my displacement with neighbours' displacement
                disp[pointI] =
                    0.5
                   *(disp[pointI] + dispSum[pointI]/dispCount[pointI]);
            }
        }
    }
}


void Foam::autoSnapDriver::calcNearestFace
(
    const label iter,
    const indirectPrimitivePatch& pp,
    const scalarField& faceSnapDist,
    vectorField& faceDisp,
    vectorField& faceSurfaceNormal,
    labelList& faceSurfaceGlobalRegion,
    vectorField& faceRotation
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const refinementSurfaces& surfaces = meshRefiner_.surfaces();

    // Displacement and orientation per pp face.
    faceDisp.setSize(pp.size());
    faceDisp = vector::zero;
    faceSurfaceNormal.setSize(pp.size());
    faceSurfaceNormal = vector::zero;
    faceSurfaceGlobalRegion.setSize(pp.size());
    faceSurfaceGlobalRegion = -1;

    // Divide surfaces into zoned and unzoned
    const labelList zonedSurfaces =
        surfaceZonesInfo::getNamedSurfaces(surfaces.surfZones());
    const labelList unzonedSurfaces =
        surfaceZonesInfo::getUnnamedSurfaces(surfaces.surfZones());

    // Per pp face the current surface snapped to
    labelList snapSurf(pp.size(), -1);


    // Do zoned surfaces
    // ~~~~~~~~~~~~~~~~~
    // Zoned faces only attract to corresponding surface

    // Extract faces per zone
    const PtrList<surfaceZonesInfo>& surfZones = surfaces.surfZones();

    forAll(zonedSurfaces, i)
    {
        label zoneSurfI = zonedSurfaces[i];

        const word& faceZoneName = surfZones[zoneSurfI].faceZoneName();

        // Get indices of faces on pp that are also in zone
        label zoneI = mesh.faceZones().findZoneID(faceZoneName);
        if (zoneI == -1)
        {
            FatalErrorIn
            (
                "autoSnapDriver::calcNearestFace(..)"
            )   << "Problem. Cannot find zone " << faceZoneName
                << exit(FatalError);
        }
        const faceZone& fZone = mesh.faceZones()[zoneI];
        PackedBoolList isZonedFace(mesh.nFaces());
        forAll(fZone, i)
        {
            isZonedFace[fZone[i]] = 1;
        }

        DynamicList<label> ppFaces(fZone.size());
        DynamicList<label> meshFaces(fZone.size());
        forAll(pp.addressing(), i)
        {
            if (isZonedFace[pp.addressing()[i]])
            {
                snapSurf[i] = zoneSurfI;
                ppFaces.append(i);
                meshFaces.append(pp.addressing()[i]);
            }
        }

        //Pout<< "For faceZone " << fZone.name()
        //    << " found " << ppFaces.size() << " out of " << pp.size()
        //    << endl;

        pointField fc
        (
            indirectPrimitivePatch
            (
                IndirectList<face>(mesh.faces(), meshFaces),
                mesh.points()
            ).faceCentres()
        );

        List<pointIndexHit> hitInfo;
        labelList hitSurface;
        labelList hitRegion;
        vectorField hitNormal;
        surfaces.findNearestRegion
        (
            labelList(1, zoneSurfI),
            fc,
            sqr(faceSnapDist),// sqr of attract dist
            hitSurface,
            hitInfo,
            hitRegion,
            hitNormal
        );

        forAll(hitInfo, hitI)
        {
            if (hitInfo[hitI].hit())
            {
                label faceI = ppFaces[hitI];
                faceDisp[faceI] = hitInfo[hitI].hitPoint() - fc[hitI];
                faceSurfaceNormal[faceI] = hitNormal[hitI];
                faceSurfaceGlobalRegion[faceI] = surfaces.globalRegion
                (
                    hitSurface[hitI],
                    hitRegion[hitI]
                );
            }
            else
            {
                WarningIn
                (
                    "autoSnapDriver::calcNearestFace(..)"
                )   << "Did not find surface near face centre " << fc[hitI]
                    << endl;
            }
        }
    }


    // Do unzoned surfaces
    // ~~~~~~~~~~~~~~~~~~~
    // Unzoned faces attract to any unzoned surface

    DynamicList<label> ppFaces(pp.size());
    DynamicList<label> meshFaces(pp.size());
    forAll(pp.addressing(), i)
    {
        if (snapSurf[i] == -1)
        {
            ppFaces.append(i);
            meshFaces.append(pp.addressing()[i]);
        }
    }
    //Pout<< "Found " << ppFaces.size() << " unzoned faces out of "
    //   << pp.size() << endl;

    pointField fc
    (
        indirectPrimitivePatch
        (
            IndirectList<face>(mesh.faces(), meshFaces),
            mesh.points()
        ).faceCentres()
    );

    List<pointIndexHit> hitInfo;
    labelList hitSurface;
    labelList hitRegion;
    vectorField hitNormal;
    surfaces.findNearestRegion
    (
        unzonedSurfaces,
        fc,
        sqr(faceSnapDist),// sqr of attract dist
        hitSurface,
        hitInfo,
        hitRegion,
        hitNormal
    );

    forAll(hitInfo, hitI)
    {
        if (hitInfo[hitI].hit())
        {
            label faceI = ppFaces[hitI];
            faceDisp[faceI] = hitInfo[hitI].hitPoint() - fc[hitI];
            faceSurfaceNormal[faceI] = hitNormal[hitI];
            faceSurfaceGlobalRegion[faceI] = surfaces.globalRegion
            (
                hitSurface[hitI],
                hitRegion[hitI]
            );
        }
        else
        {
            WarningIn
            (
                "autoSnapDriver::calcNearestFace(..)"
            )   << "Did not find surface near face centre " << fc[hitI]
                << endl;
        }
    }


    // Determine rotation
    // ~~~~~~~~~~~~~~~~~~

    // Determine rotation axis
    faceRotation.setSize(pp.size());
    faceRotation = vector::zero;

    forAll(faceRotation, faceI)
    {
        // Note: extend to >180 degrees checking
        faceRotation[faceI] =
            pp.faceNormals()[faceI]
          ^ faceSurfaceNormal[faceI];
    }

    if (debug&meshRefinement::ATTRACTION)
    {
        dumpMove
        (
            mesh.time().path()
          / "faceDisp_" + name(iter) + ".obj",
            pp.faceCentres(),
            pp.faceCentres() + faceDisp
        );
        dumpMove
        (
            mesh.time().path()
          / "faceRotation_" + name(iter) + ".obj",
            pp.faceCentres(),
            pp.faceCentres() + faceRotation
        );
    }
}


// Collect (possibly remote) per point data of all surrounding faces
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// - faceSurfaceNormal
// - faceDisp
// - faceCentres&faceNormal
void Foam::autoSnapDriver::calcNearestFacePointProperties
(
    const label iter,
    const indirectPrimitivePatch& pp,

    const vectorField& faceDisp,
    const vectorField& faceSurfaceNormal,
    const labelList& faceSurfaceGlobalRegion,

    List<List<point> >& pointFaceSurfNormals,
    List<List<point> >& pointFaceDisp,
    List<List<point> >& pointFaceCentres,
    List<labelList>&    pointFacePatchID
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    const PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh));


    // For now just get all surrounding face data. Expensive - should just
    // store and sync data on coupled points only
    // (see e.g PatchToolsNormals.C)

    pointFaceSurfNormals.setSize(pp.nPoints());
    pointFaceDisp.setSize(pp.nPoints());
    pointFaceCentres.setSize(pp.nPoints());
    pointFacePatchID.setSize(pp.nPoints());

    // Fill local data
    forAll(pp.pointFaces(), pointI)
    {
        const labelList& pFaces = pp.pointFaces()[pointI];

        // Count valid face normals
        label nFaces = 0;
        forAll(pFaces, i)
        {
            label faceI = pFaces[i];
            if (isMasterFace[faceI] && faceSurfaceGlobalRegion[faceI] != -1)
            {
                nFaces++;
            }
        }


        List<point>& pNormals = pointFaceSurfNormals[pointI];
        pNormals.setSize(nFaces);
        List<point>& pDisp = pointFaceDisp[pointI];
        pDisp.setSize(nFaces);
        List<point>& pFc = pointFaceCentres[pointI];
        pFc.setSize(nFaces);
        labelList& pFid = pointFacePatchID[pointI];
        pFid.setSize(nFaces);

        nFaces = 0;
        forAll(pFaces, i)
        {
            label faceI = pFaces[i];
            label globalRegionI = faceSurfaceGlobalRegion[faceI];

            if (isMasterFace[faceI] && globalRegionI != -1)
            {
                pNormals[nFaces] = faceSurfaceNormal[faceI];
                pDisp[nFaces] = faceDisp[faceI];
                pFc[nFaces] = pp.faceCentres()[faceI];
                pFid[nFaces] = globalToMasterPatch_[globalRegionI];
                nFaces++;
            }
        }
    }


    // Collect additionally 'normal' boundary faces for boundaryPoints of pp
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // points on the boundary of pp should pick up non-pp normals
    // as well for the feature-reconstruction to behave correctly.
    // (the movement is already constrained outside correctly so it
    //  is only that the unconstrained attraction vector is calculated
    //  correctly)
    {
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();
        labelList patchID(pbm.patchID());

        // Unmark all non-coupled boundary faces
        forAll(pbm, patchI)
        {
            const polyPatch& pp = pbm[patchI];

            if (pp.coupled() || isA<emptyPolyPatch>(pp))
            {
                forAll(pp, i)
                {
                    label meshFaceI = pp.start()+i;
                    patchID[meshFaceI-mesh.nInternalFaces()] = -1;
                }
            }
        }

        // Remove any meshed faces
        forAll(pp.addressing(), i)
        {
            label meshFaceI = pp.addressing()[i];
            patchID[meshFaceI-mesh.nInternalFaces()] = -1;
        }

        // See if pp point uses any non-meshed boundary faces

        const labelList& boundaryPoints = pp.boundaryPoints();
        forAll(boundaryPoints, i)
        {
            label pointI = boundaryPoints[i];
            label meshPointI = pp.meshPoints()[pointI];
            const point& pt = mesh.points()[meshPointI];
            const labelList& pFaces = mesh.pointFaces()[meshPointI];

            List<point>& pNormals = pointFaceSurfNormals[pointI];
            List<point>& pDisp = pointFaceDisp[pointI];
            List<point>& pFc = pointFaceCentres[pointI];
            labelList& pFid = pointFacePatchID[pointI];

            forAll(pFaces, i)
            {
                label meshFaceI = pFaces[i];
                if (!mesh.isInternalFace(meshFaceI))
                {
                    label patchI = patchID[meshFaceI-mesh.nInternalFaces()];

                    if (patchI != -1)
                    {
                        vector fn = mesh.faceAreas()[meshFaceI];
                        pNormals.append(fn/mag(fn));
                        pDisp.append(mesh.faceCentres()[meshFaceI]-pt);
                        pFc.append(mesh.faceCentres()[meshFaceI]);
                        pFid.append(patchI);
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        pointFaceSurfNormals,
        listPlusEqOp<point>(),
        List<point>(),
        mapDistribute::transform()
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        pointFaceDisp,
        listPlusEqOp<point>(),
        List<point>(),
        mapDistribute::transform()
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        pointFaceCentres,
        listPlusEqOp<point>(),
        List<point>(),
        mapDistribute::transformPosition()
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        pointFacePatchID,
        listPlusEqOp<label>(),
        List<label>()
    );


    // Sort the data according to the face centres. This is only so we get
    // consistent behaviour serial and parallel.
    labelList visitOrder;
    forAll(pointFaceDisp, pointI)
    {
        List<point>& pNormals = pointFaceSurfNormals[pointI];
        List<point>& pDisp = pointFaceDisp[pointI];
        List<point>& pFc = pointFaceCentres[pointI];
        labelList& pFid = pointFacePatchID[pointI];

        sortedOrder(mag(pFc)(), visitOrder);

        pNormals = List<point>(pNormals, visitOrder);
        pDisp = List<point>(pDisp, visitOrder);
        pFc = List<point>(pFc, visitOrder);
        pFid = UIndirectList<label>(pFid, visitOrder)();
    }
}


// Gets passed in offset to nearest point on feature edge. Calculates
// if the point has a different number of faces on either side of the feature
// and if so attracts the point to that non-dominant plane.
void Foam::autoSnapDriver::correctAttraction
(
    const DynamicList<point>& surfacePoints,
    const DynamicList<label>& surfaceCounts,
    const point& edgePt,
    const vector& edgeNormal,       // normalised normal
    const point& pt,

    vector& edgeOffset              // offset from pt to point on edge
) const
{
    // Tangential component along edge
    scalar tang = ((pt-edgePt)&edgeNormal);

    labelList order;
    Foam::sortedOrder(surfaceCounts, order);

    if (order[0] < order[1])
    {
        // There is a non-dominant plane. Use the point on the plane to
        // attract to.
        vector attractD = surfacePoints[order[0]]-edgePt;
        // Tangential component along edge
        scalar tang2 = (attractD&edgeNormal);
        // Normal component
        attractD -= tang2*edgeNormal;
        // Calculate fraction of normal distances
        scalar magAttractD = mag(attractD);
        scalar fraction = magAttractD/(magAttractD+mag(edgeOffset));

        point linePt =
            edgePt
          + ((1.0-fraction)*tang2 + fraction*tang)*edgeNormal;
        edgeOffset = linePt-pt;
    }
}


Foam::pointIndexHit Foam::autoSnapDriver::findMultiPatchPoint
(
    const point& pt,
    const labelList& patchIDs,
    const List<point>& faceCentres
) const
{
    // Determine if multiple patchIDs
    if (patchIDs.size())
    {
        label patch0 = patchIDs[0];

        for (label i = 1; i < patchIDs.size(); i++)
        {
            if (patchIDs[i] != patch0)
            {
                return pointIndexHit(true, pt, labelMax);
            }
        }
    }
    return pointIndexHit(false, vector::zero, labelMax);
}


Foam::label Foam::autoSnapDriver::findNormal
(
    const scalar featureCos,
    const vector& n,
    const DynamicList<vector>& surfaceNormals
) const
{
    label index = -1;

    forAll(surfaceNormals, j)
    {
        scalar cosAngle = (n&surfaceNormals[j]);

        if
        (
            (cosAngle >= featureCos)
         || (cosAngle < (-1+0.001)) // triangle baffles
        )
        {
            index = j;
            break;
        }
    }
    return index;
}


// Detect multiple patches. Returns pointIndexHit:
// - false, index=-1 : single patch
// - true , index=0  : multiple patches but on different normals planes
//                     (so geometric feature edge is also a region edge)
// - true , index=1  : multiple patches on same normals plane i.e. flat region
//                     edge
Foam::pointIndexHit Foam::autoSnapDriver::findMultiPatchPoint
(
    const point& pt,
    const labelList& patchIDs,
    const DynamicList<vector>& surfaceNormals,
    const labelList& faceToNormalBin
) const
{
    if (patchIDs.empty())
    {
        return pointIndexHit(false, pt, -1);
    }

    // Detect single patch situation (to avoid allocation)
    label patch0 = patchIDs[0];

    for (label i = 1; i < patchIDs.size(); i++)
    {
        if (patchIDs[i] != patch0)
        {
            patch0 = -1;
            break;
        }
    }

    if (patch0 >= 0)
    {
        // Single patch
        return pointIndexHit(false, pt, -1);
    }
    else
    {
        if (surfaceNormals.size() == 1)
        {
            // Same normals plane, flat region edge.
            return pointIndexHit(true, pt, 1);
        }
        else
        {
            // Detect per normals bin
            labelList normalToPatch(surfaceNormals.size(), -1);
            forAll(faceToNormalBin, i)
            {
                if (faceToNormalBin[i] != -1)
                {
                    label& patch = normalToPatch[faceToNormalBin[i]];
                    if (patch == -1)
                    {
                        // First occurence
                        patch = patchIDs[i];
                    }
                    else if (patch == -2)
                    {
                        // Already marked as being on multiple patches
                    }
                    else if (patch != patchIDs[i])
                    {
                        // Mark as being on multiple patches
                        patch = -2;
                    }
                }
            }

            forAll(normalToPatch, normalI)
            {
                if (normalToPatch[normalI] == -2)
                {
                    // Multiple patches on same normals plane, flat region
                    // edge
                    return pointIndexHit(true, pt, 1);
                }
            }

            // All patches on either side of geometric feature anyway
            return pointIndexHit(true, pt, 0);
        }
    }
}


void Foam::autoSnapDriver::featureAttractionUsingReconstruction
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,
    const vectorField& nearestDisp,
    const label pointI,

    const List<List<point> >& pointFaceSurfNormals,
    const List<List<point> >& pointFaceDisp,
    const List<List<point> >& pointFaceCentres,
    const labelListList& pointFacePatchID,

    DynamicList<point>& surfacePoints,
    DynamicList<vector>& surfaceNormals,
    labelList& faceToNormalBin,

    vector& patchAttraction,
    pointConstraint& patchConstraint
) const
{
    patchAttraction = vector::zero;
    patchConstraint = pointConstraint();

    const List<point>& pfSurfNormals = pointFaceSurfNormals[pointI];
    const List<point>& pfDisp = pointFaceDisp[pointI];
    const List<point>& pfCentres = pointFaceCentres[pointI];

    // Bin according to surface normal
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //- bins of differing normals:
    //  - one normal   : flat(tish) surface
    //  - two normals  : geometric feature edge
    //  - three normals: geometric feature point
    //  - four normals : too complex a feature
    surfacePoints.clear();
    surfaceNormals.clear();

    //- from face to above normals bin
    faceToNormalBin.setSize(pfDisp.size());
    faceToNormalBin = -1;

    forAll(pfSurfNormals, i)
    {
        const point& fc = pfCentres[i];
        const vector& fSNormal = pfSurfNormals[i];
        const vector& fDisp = pfDisp[i];

        // What to do with very far attraction? For now just ignore the face
        if (magSqr(fDisp) < sqr(snapDist[pointI]) && mag(fSNormal) > VSMALL)
        {
            const point pt = fc + fDisp;

            // Do we already have surface normal?
            faceToNormalBin[i] = findNormal
            (
                featureCos,
                fSNormal,
                surfaceNormals
            );

            if (faceToNormalBin[i] != -1)
            {
                // Same normal
            }
            else
            {
                // Now check if the planes go through the same edge or point

                if (surfacePoints.size() <= 1)
                {
                    surfacePoints.append(pt);
                    faceToNormalBin[i] = surfaceNormals.size();
                    surfaceNormals.append(fSNormal);
                }
                else if (surfacePoints.size() == 2)
                {
                    plane pl0(surfacePoints[0], surfaceNormals[0]);
                    plane pl1(surfacePoints[1], surfaceNormals[1]);
                    plane::ray r(pl0.planeIntersect(pl1));
                    vector featureNormal = r.dir() / mag(r.dir());

                    if (mag(fSNormal&featureNormal) >= 0.001)
                    {
                        // Definitely makes a feature point
                        surfacePoints.append(pt);
                        faceToNormalBin[i] = surfaceNormals.size();
                        surfaceNormals.append(fSNormal);
                    }
                }
                else if (surfacePoints.size() == 3)
                {
                    // Have already feature point. See if this new plane is
                    // the same point or not.
                    plane pl0(surfacePoints[0], surfaceNormals[0]);
                    plane pl1(surfacePoints[1], surfaceNormals[1]);
                    plane pl2(surfacePoints[2], surfaceNormals[2]);
                    point p012(pl0.planePlaneIntersect(pl1, pl2));

                    plane::ray r(pl0.planeIntersect(pl1));
                    vector featureNormal = r.dir() / mag(r.dir());
                    if (mag(fSNormal&featureNormal) >= 0.001)
                    {
                        plane pl3(pt, fSNormal);
                        point p013(pl0.planePlaneIntersect(pl1, pl3));

                        if (mag(p012-p013) > snapDist[pointI])
                        {
                            // Different feature point
                            surfacePoints.append(pt);
                            faceToNormalBin[i] = surfaceNormals.size();
                            surfaceNormals.append(fSNormal);
                        }
                    }
                }
            }
        }
    }


    const point& pt = pp.localPoints()[pointI];

    // Check the number of directions
    if (surfaceNormals.size() == 1)
    {
        // Normal distance to plane
        vector d =
            ((surfacePoints[0]-pt) & surfaceNormals[0])
           *surfaceNormals[0];

        // Trim to snap distance
        if (magSqr(d) > sqr(snapDist[pointI]))
        {
            d *= Foam::sqrt(sqr(snapDist[pointI])/magSqr(d));
        }

        patchAttraction = d;

        // Store constraints
        patchConstraint.applyConstraint(surfaceNormals[0]);
    }
    else if (surfaceNormals.size() == 2)
    {
        plane pl0(surfacePoints[0], surfaceNormals[0]);
        plane pl1(surfacePoints[1], surfaceNormals[1]);
        plane::ray r(pl0.planeIntersect(pl1));
        vector n = r.dir() / mag(r.dir());

        // Get nearest point on infinite ray
        vector d = r.refPoint()-pt;
        d -= (d&n)*n;

        // Trim to snap distance
        if (magSqr(d) > sqr(snapDist[pointI]))
        {
            d *= Foam::sqrt(sqr(snapDist[pointI])/magSqr(d));
        }

        patchAttraction = d;

        // Store constraints
        patchConstraint.applyConstraint(surfaceNormals[0]);
        patchConstraint.applyConstraint(surfaceNormals[1]);
    }
    else if (surfaceNormals.size() == 3)
    {
        // Calculate point from the faces.
        plane pl0(surfacePoints[0], surfaceNormals[0]);
        plane pl1(surfacePoints[1], surfaceNormals[1]);
        plane pl2(surfacePoints[2], surfaceNormals[2]);
        point cornerPt(pl0.planePlaneIntersect(pl1, pl2));
        vector d = cornerPt - pt;

        // Trim to snap distance
        if (magSqr(d) > sqr(snapDist[pointI]))
        {
            d *= Foam::sqrt(sqr(snapDist[pointI])/magSqr(d));
        }

        patchAttraction = d;

        // Store constraints
        patchConstraint.applyConstraint(surfaceNormals[0]);
        patchConstraint.applyConstraint(surfaceNormals[1]);
        patchConstraint.applyConstraint(surfaceNormals[2]);
    }
}


// Special version that calculates attraction in one go
void Foam::autoSnapDriver::featureAttractionUsingReconstruction
(
    const label iter,
    const bool avoidSnapProblems,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,
    const vectorField& nearestDisp,

    const List<List<point> >& pointFaceSurfNormals,
    const List<List<point> >& pointFaceDisp,
    const List<List<point> >& pointFaceCentres,
    const labelListList& pointFacePatchID,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    autoPtr<OBJstream> feStr;
    autoPtr<OBJstream> fpStr;
    if (debug&meshRefinement::ATTRACTION)
    {
        feStr.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "implicitFeatureEdge_" + name(iter) + ".obj"
            )
        );
        Info<< "Dumping implicit feature-edge direction to "
            << feStr().name() << endl;

        fpStr.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "implicitFeaturePoint_" + name(iter) + ".obj"
            )
        );
        Info<< "Dumping implicit feature-point direction to "
            << fpStr().name() << endl;
    }


    DynamicList<point> surfacePoints(4);
    DynamicList<vector> surfaceNormals(4);
    labelList faceToNormalBin;

    forAll(pp.localPoints(), pointI)
    {
        vector attraction = vector::zero;
        pointConstraint constraint;

        featureAttractionUsingReconstruction
        (
            iter,
            featureCos,

            pp,
            snapDist,
            nearestDisp,

            pointI,

            pointFaceSurfNormals,
            pointFaceDisp,
            pointFaceCentres,
            pointFacePatchID,

            surfacePoints,
            surfaceNormals,
            faceToNormalBin,

            attraction,
            constraint
        );

        if
        (
            (constraint.first() > patchConstraints[pointI].first())
         || (
                (constraint.first() == patchConstraints[pointI].first())
             && (magSqr(attraction) < magSqr(patchAttraction[pointI]))
            )
        )
        {
            patchAttraction[pointI] = attraction;
            patchConstraints[pointI] = constraint;

            const point& pt = pp.localPoints()[pointI];

            if (patchConstraints[pointI].first() == 2 && feStr.valid())
            {
                feStr().write(linePointRef(pt, pt+patchAttraction[pointI]));
            }
            else if (patchConstraints[pointI].first() == 3 && fpStr.valid())
            {
                fpStr().write(linePointRef(pt, pt+patchAttraction[pointI]));
            }
        }
    }
}


void Foam::autoSnapDriver::stringFeatureEdges
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,

    const vectorField& rawPatchAttraction,
    const List<pointConstraint>& rawPatchConstraints,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    // Snap edges to feature edges
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Walk existing edges and snap remaining ones (that are marked as
    // feature edges in rawPatchConstraints)

    // What this does is fill in any faces where not all points
    // on the face are being attracted:
    /*
           +
          / \
         /   \
      ---+    +---
         \   /
          \ /
           +
    */
    // so the top and bottom will never get attracted since the nearest
    // back from the feature edge will always be one of the left or right
    // points since the face is diamond like. So here we walk the feature edges
    // and add any non-attracted points.


    while (true)
    {
        label nChanged = 0;

        const labelListList& pointEdges = pp.pointEdges();
        forAll(pointEdges, pointI)
        {
            if (patchConstraints[pointI].first() == 2)
            {
                const point& pt = pp.localPoints()[pointI];
                const labelList& pEdges = pointEdges[pointI];
                const vector& featVec = patchConstraints[pointI].second();

                // Detect whether there are edges in both directions.
                // (direction along the feature edge that is)
                bool hasPos = false;
                bool hasNeg = false;

                forAll(pEdges, pEdgeI)
                {
                    const edge& e = pp.edges()[pEdges[pEdgeI]];
                    label nbrPointI = e.otherVertex(pointI);

                    if (patchConstraints[nbrPointI].first() > 1)
                    {
                        const point& nbrPt = pp.localPoints()[nbrPointI];
                        const point featPt =
                            nbrPt + patchAttraction[nbrPointI];
                        const scalar cosAngle = (featVec & (featPt-pt));

                        if (cosAngle > 0)
                        {
                            hasPos = true;
                        }
                        else
                        {
                            hasNeg = true;
                        }
                    }
                }

                if (!hasPos || !hasNeg)
                {
                    //Pout<< "**Detected feature string end at  "
                    //    << pp.localPoints()[pointI] << endl;

                    // No string. Assign best choice on either side
                    label bestPosPointI = -1;
                    scalar minPosDistSqr = GREAT;
                    label bestNegPointI = -1;
                    scalar minNegDistSqr = GREAT;

                    forAll(pEdges, pEdgeI)
                    {
                        const edge& e = pp.edges()[pEdges[pEdgeI]];
                        label nbrPointI = e.otherVertex(pointI);

                        if
                        (
                            patchConstraints[nbrPointI].first() <= 1
                         && rawPatchConstraints[nbrPointI].first() > 1
                        )
                        {
                            const vector& nbrFeatVec =
                                rawPatchConstraints[pointI].second();

                            if (mag(featVec&nbrFeatVec) > featureCos)
                            {
                                // nbrPointI attracted to sameish feature
                                // Note: also check on position.

                                scalar d2 = magSqr
                                (
                                    rawPatchAttraction[nbrPointI]
                                );

                                const point featPt =
                                    pp.localPoints()[nbrPointI]
                                  + rawPatchAttraction[nbrPointI];
                                const scalar cosAngle =
                                    (featVec & (featPt-pt));

                                if (cosAngle > 0)
                                {
                                    if (!hasPos && d2 < minPosDistSqr)
                                    {
                                        minPosDistSqr = d2;
                                        bestPosPointI = nbrPointI;
                                    }
                                }
                                else
                                {
                                    if (!hasNeg && d2 < minNegDistSqr)
                                    {
                                        minNegDistSqr = d2;
                                        bestNegPointI = nbrPointI;
                                    }
                                }
                            }
                        }
                    }

                    if (bestPosPointI != -1)
                    {
                        // Use reconstructed-feature attraction. Use only
                        // part of it since not sure...
                        //const point& bestPt =
                        //    pp.localPoints()[bestPosPointI];
                        //Pout<< "**Overriding point " << bestPt
                        //    << " on reconstructed feature edge at "
                        //    << rawPatchAttraction[bestPosPointI]+bestPt
                        //    << " to attracted-to-feature-edge." << endl;
                        patchAttraction[bestPosPointI] =
                            0.5*rawPatchAttraction[bestPosPointI];
                        patchConstraints[bestPosPointI] =
                            rawPatchConstraints[bestPosPointI];

                        nChanged++;
                    }
                    if (bestNegPointI != -1)
                    {
                        // Use reconstructed-feature attraction. Use only
                        // part of it since not sure...
                        //const point& bestPt =
                        //    pp.localPoints()[bestNegPointI];
                        //Pout<< "**Overriding point " << bestPt
                        //    << " on reconstructed feature edge at "
                        //    << rawPatchAttraction[bestNegPointI]+bestPt
                        //    << " to attracted-to-feature-edge." << endl;
                        patchAttraction[bestNegPointI] =
                            0.5*rawPatchAttraction[bestNegPointI];
                        patchConstraints[bestNegPointI] =
                            rawPatchConstraints[bestNegPointI];

                        nChanged++;
                    }
                }
            }
        }


        reduce(nChanged, sumOp<label>());
        Info<< "Stringing feature edges : changed " << nChanged << " points"
            << endl;
        if (nChanged == 0)
        {
            break;
        }
    }
}


void Foam::autoSnapDriver::releasePointsNextToMultiPatch
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,

    const List<List<point> >& pointFaceCentres,
    const labelListList& pointFacePatchID,

    const vectorField& rawPatchAttraction,
    const List<pointConstraint>& rawPatchConstraints,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    autoPtr<OBJstream> multiPatchStr;
    if (debug&meshRefinement::ATTRACTION)
    {
        multiPatchStr.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "multiPatch_" + name(iter) + ".obj"
            )
        );
        Info<< "Dumping removed constraints due to same-face"
            << " multi-patch points to "
            << multiPatchStr().name() << endl;
    }


    // 1. Mark points on multiple patches
    PackedBoolList isMultiPatchPoint(pp.size());

    forAll(pointFacePatchID, pointI)
    {
        pointIndexHit multiPatchPt = findMultiPatchPoint
        (
            pp.localPoints()[pointI],
            pointFacePatchID[pointI],
            pointFaceCentres[pointI]
        );
        isMultiPatchPoint[pointI] = multiPatchPt.hit();
    }

    // 2. Make sure multi-patch points are also attracted
    forAll(isMultiPatchPoint, pointI)
    {
        if (isMultiPatchPoint[pointI])
        {
            if
            (
                patchConstraints[pointI].first() <= 1
             && rawPatchConstraints[pointI].first() > 1
            )
            {
                patchAttraction[pointI] = rawPatchAttraction[pointI];
                patchConstraints[pointI] = rawPatchConstraints[pointI];

                //if (multiPatchStr.valid())
                //{
                //    Pout<< "Adding constraint on multiPatchPoint:"
                //        << pp.localPoints()[pointI]
                //        << " constraint:" << patchConstraints[pointI]
                //        << " attraction:" << patchAttraction[pointI]
                //        << endl;
                //}
            }
        }
    }

    // Up to here it is all parallel ok.


    // 3. Knock out any attraction on faces with multi-patch points
    label nChanged = 0;
    forAll(pp.localFaces(), faceI)
    {
        const face& f = pp.localFaces()[faceI];

        label nMultiPatchPoints = 0;
        forAll(f, fp)
        {
            label pointI = f[fp];
            if
            (
                isMultiPatchPoint[pointI]
             && patchConstraints[pointI].first() > 1
            )
            {
                nMultiPatchPoints++;
            }
        }

        if (nMultiPatchPoints > 0)
        {
            forAll(f, fp)
            {
                label pointI = f[fp];
                if
                (
                   !isMultiPatchPoint[pointI]
                 && patchConstraints[pointI].first() > 1
                )
                {
                    //Pout<< "Knocking out constraint"
                    //    << " on non-multiPatchPoint:"
                    //    << pp.localPoints()[pointI] << endl;
                    patchAttraction[pointI] = vector::zero;
                    patchConstraints[pointI] = pointConstraint();
                    nChanged++;

                    if (multiPatchStr.valid())
                    {
                        multiPatchStr().write(pp.localPoints()[pointI]);
                    }
                }
            }
        }
    }

    reduce(nChanged, sumOp<label>());
    Info<< "Removing constraints near multi-patch points : changed "
        << nChanged << " points" << endl;
}


// If only two attractions and across face return the face indices
Foam::labelPair Foam::autoSnapDriver::findDiagonalAttraction
(
    const indirectPrimitivePatch& pp,
    const vectorField& patchAttraction,
    const List<pointConstraint>& patchConstraints,
    const label faceI
) const
{
    const face& f = pp.localFaces()[faceI];
    // For now just detect any attraction. Improve this to look at
    // actual attraction position and orientation

    labelPair attractIndices(-1, -1);

    if (f.size() >= 4)
    {
        forAll(f, fp)
        {
            label pointI = f[fp];
            if (patchConstraints[pointI].first() >= 2)
            {
                // Attract to feature edge or point
                if (attractIndices[0] == -1)
                {
                    // First attraction. Store
                    attractIndices[0] = fp;
                }
                else if (attractIndices[1] == -1)
                {
                    // Second attraction. Check if not consecutive to first
                    // attraction
                    label fp0 = attractIndices[0];
                    if (f.fcIndex(fp0) == fp || f.fcIndex(fp) == fp0)
                    {
                        // Consecutive. Skip.
                        attractIndices = labelPair(-1, -1);
                        break;
                    }
                    else
                    {
                        attractIndices[1] = fp;
                    }
                }
                else
                {
                    // More than two attractions. Skip.
                    attractIndices = labelPair(-1, -1);
                    break;
                }
            }
        }


        if (attractIndices[1] == -1)
        {
            // Found only one attraction. Skip.
            attractIndices = labelPair(-1, -1);
        }
    }
    return attractIndices;
}


void Foam::autoSnapDriver::avoidDiagonalAttraction
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    forAll(pp.localFaces(), faceI)
    {
        const face& f = pp.localFaces()[faceI];

        labelPair diag = findDiagonalAttraction
        (
            pp,
            patchAttraction,
            patchConstraints,
            faceI
        );

        if (diag[0] != -1 && diag[1] != -1)
        {
            // Found two diagonal points that being attracted.
            // For now just attract my one to the average of those.
            const label i0 = f[diag[0]];
            const point pt0 =
                pp.localPoints()[i0]+patchAttraction[i0];
            const label i1 = f[diag[1]];
            const point pt1 =
                pp.localPoints()[i1]+patchAttraction[i1];
            const point mid = 0.5*(pt0+pt1);

            const scalar cosAngle = mag
            (
                patchConstraints[i0].second()
              & patchConstraints[i1].second()
            );

            //Pout<< "Found diagonal attraction at indices:"
            //    << diag[0]
            //    << " and " << diag[1]
            //    << " with cosAngle:" << cosAngle
            //    << " mid:" << mid << endl;

            if (cosAngle > featureCos)
            {
                // Add the nearest of the other two points as
                // attractor
                label minFp = -1;
                scalar minDistSqr = GREAT;
                forAll(f, fp)
                {
                    label pointI = f[fp];
                    if (patchConstraints[pointI].first() <= 1)
                    {
                        const point& pt = pp.localPoints()[pointI];
                        scalar distSqr = magSqr(mid-pt);
                        if (distSqr < minDistSqr)
                        {
                            distSqr = minDistSqr;
                            minFp = fp;
                        }
                    }
                }
                if (minFp != -1)
                {
                    label minPointI = f[minFp];
                    patchAttraction[minPointI] =
                        mid-pp.localPoints()[minPointI];
                    patchConstraints[minPointI] =
                        patchConstraints[f[diag[0]]];
                }
            }
            else
            {
                //Pout<< "Diagonal attractors at" << nl
                //    << "    pt0:" << pt0
                //    << "    constraint:"
                //    << patchConstraints[i0].second() << nl
                //    << "    pt1:" << pt1
                //    << "    constraint:"
                //    << patchConstraints[i1].second() << nl
                //    << "    make too large an angle:"
                //    <<  mag
                //        (
                //            patchConstraints[i0].second()
                //          & patchConstraints[i1].second()
                //        )
                //    << endl;
            }
        }
    }
}


Foam::Tuple2<Foam::label, Foam::pointIndexHit>
Foam::autoSnapDriver::findNearFeatureEdge
(
    const bool isRegionEdge,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,
    const label pointI,
    const point& estimatedPt,

    List<List<DynamicList<point> > >& edgeAttractors,
    List<List<DynamicList<pointConstraint> > >& edgeConstraints,
    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    const refinementFeatures& features = meshRefiner_.features();

    labelList nearEdgeFeat;
    List<pointIndexHit> nearEdgeInfo;
    vectorField nearNormal;

    if (isRegionEdge)
    {
        features.findNearestRegionEdge
        (
            pointField(1, estimatedPt),
            scalarField(1, sqr(snapDist[pointI])),
            nearEdgeFeat,
            nearEdgeInfo,
            nearNormal
        );
    }
    else
    {
        features.findNearestEdge
        (
            pointField(1, estimatedPt),
            scalarField(1, sqr(snapDist[pointI])),
            nearEdgeFeat,
            nearEdgeInfo,
            nearNormal
        );
    }

    const pointIndexHit& nearInfo = nearEdgeInfo[0];
    label featI = nearEdgeFeat[0];

    if (nearInfo.hit())
    {
        // So we have a point on the feature edge. Use this
        // instead of our estimate from planes.
        edgeAttractors[featI][nearInfo.index()].append
        (
            nearInfo.hitPoint()
        );
        pointConstraint c(Tuple2<label, vector>(2, nearNormal[0]));
        edgeConstraints[featI][nearInfo.index()].append(c);

        // Store for later use
        patchAttraction[pointI] =
            nearInfo.hitPoint()-pp.localPoints()[pointI];
        patchConstraints[pointI] = c;
    }
    return Tuple2<label, pointIndexHit>(featI, nearInfo);
}


Foam::Tuple2<Foam::label, Foam::pointIndexHit>
Foam::autoSnapDriver::findNearFeaturePoint
(
    const bool isRegionPoint,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,
    const label pointI,
    const point& estimatedPt,

    // Feature-point to pp point
    List<labelList>& pointAttractor,
    List<List<pointConstraint> >& pointConstraints,
    // Feature-edge to pp point
    List<List<DynamicList<point> > >& edgeAttractors,
    List<List<DynamicList<pointConstraint> > >& edgeConstraints,
    // pp point to nearest feature
    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    const refinementFeatures& features = meshRefiner_.features();

    labelList nearFeat;
    List<pointIndexHit> nearInfo;
    features.findNearestPoint
    (
        pointField(1, estimatedPt),
        scalarField(1, sqr(snapDist[pointI])),
        nearFeat,
        nearInfo
    );

    label featI = nearFeat[0];

    if (featI != -1)
    {
        const point& pt = pp.localPoints()[pointI];

        label featPointI = nearInfo[0].index();
        const point& featPt = nearInfo[0].hitPoint();
        scalar distSqr = magSqr(featPt-pt);

        // Check if already attracted
        label oldPointI = pointAttractor[featI][featPointI];

        if (oldPointI != -1)
        {
            // Check distance
            if (distSqr >= magSqr(featPt-pp.localPoints()[oldPointI]))
            {
                // oldPointI nearest. Keep.
                featI = -1;
                featPointI = -1;
            }
            else
            {
                // Current pointI nearer.
                pointAttractor[featI][featPointI] = pointI;
                pointConstraints[featI][featPointI].first() = 3;
                pointConstraints[featI][featPointI].second() = vector::zero;

                // Store for later use
                patchAttraction[pointI] = featPt-pt;
                patchConstraints[pointI] =
                    pointConstraints[featI][featPointI];

                // Reset oldPointI to nearest on feature edge
                patchAttraction[oldPointI] = vector::zero;
                patchConstraints[oldPointI] = pointConstraint();

                findNearFeatureEdge
                (
                    isRegionPoint,      // search region edges only

                    pp,
                    snapDist,
                    oldPointI,
                    pp.localPoints()[oldPointI],

                    edgeAttractors,
                    edgeConstraints,
                    patchAttraction,
                    patchConstraints
                );
            }
        }
        else
        {
            // Current pointI nearer.
            pointAttractor[featI][featPointI] = pointI;
            pointConstraints[featI][featPointI].first() = 3;
            pointConstraints[featI][featPointI].second() = vector::zero;

            // Store for later use
            patchAttraction[pointI] = featPt-pt;
            patchConstraints[pointI] = pointConstraints[featI][featPointI];
        }
    }

    return Tuple2<label, pointIndexHit>(featI, nearInfo[0]);
}


// Determines for every pp point - that is on multiple faces that form
// a feature - the nearest feature edge/point.
void Foam::autoSnapDriver::determineFeatures
(
    const label iter,
    const scalar featureCos,
    const bool multiRegionFeatureSnap,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,
    const vectorField& nearestDisp,

    const List<List<point> >& pointFaceSurfNormals,
    const List<List<point> >& pointFaceDisp,
    const List<List<point> >& pointFaceCentres,
    const labelListList& pointFacePatchID,

    // Feature-point to pp point
    List<labelList>& pointAttractor,
    List<List<pointConstraint> >& pointConstraints,
    // Feature-edge to pp point
    List<List<DynamicList<point> > >& edgeAttractors,
    List<List<DynamicList<pointConstraint> > >& edgeConstraints,
    // pp point to nearest feature
    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    autoPtr<OBJstream> featureEdgeStr;
    autoPtr<OBJstream> missedEdgeStr;
    autoPtr<OBJstream> featurePointStr;

    if (debug&meshRefinement::ATTRACTION)
    {
        featureEdgeStr.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "featureEdge_" + name(iter) + ".obj"
            )
        );
        Info<< "Dumping feature-edge sampling to "
            << featureEdgeStr().name() << endl;

        missedEdgeStr.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "missedFeatureEdge_" + name(iter) + ".obj"
            )
        );
        Info<< "Dumping feature-edges that are too far away to "
            << missedEdgeStr().name() << endl;

        featurePointStr.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "featurePoint_" + name(iter) + ".obj"
            )
        );
        Info<< "Dumping feature-point sampling to "
            << featurePointStr().name() << endl;
    }


    DynamicList<point> surfacePoints(4);
    DynamicList<vector> surfaceNormals(4);
    labelList faceToNormalBin;

    forAll(pp.localPoints(), pointI)
    {
        const point& pt = pp.localPoints()[pointI];

        vector attraction = vector::zero;
        pointConstraint constraint;

        featureAttractionUsingReconstruction
        (
            iter,
            featureCos,

            pp,
            snapDist,
            nearestDisp,

            pointI,

            pointFaceSurfNormals,
            pointFaceDisp,
            pointFaceCentres,
            pointFacePatchID,

            surfacePoints,
            surfaceNormals,
            faceToNormalBin,

            attraction,
            constraint
        );


        if
        (
            (constraint.first() > patchConstraints[pointI].first())
         || (
                (constraint.first() == patchConstraints[pointI].first())
             && (magSqr(attraction) < magSqr(patchAttraction[pointI]))
            )
        )
        {
            patchAttraction[pointI] = attraction;
            patchConstraints[pointI] = constraint;

            // Check the number of directions
            if (patchConstraints[pointI].first() == 1)
            {
                // Flat surface. Check for different patchIDs
                if (multiRegionFeatureSnap)
                {
                    const point estimatedPt(pt + nearestDisp[pointI]);
                    pointIndexHit multiPatchPt
                    (
                        findMultiPatchPoint
                        (
                            estimatedPt,
                            pointFacePatchID[pointI],
                            surfaceNormals,
                            faceToNormalBin
                        )
                    );

                    if (multiPatchPt.hit())
                    {
                        // Behave like when having two surface normals so
                        // attract to nearest feature edge (with a guess for
                        // the multipatch point as starting point)
                        Tuple2<label, pointIndexHit> nearInfo =
                        findNearFeatureEdge
                        (
                            true,                       // isRegionPoint
                            pp,
                            snapDist,
                            pointI,
                            multiPatchPt.hitPoint(),    // estimatedPt

                            edgeAttractors,
                            edgeConstraints,

                            patchAttraction,
                            patchConstraints
                        );

                        const pointIndexHit& info = nearInfo.second();
                        if (info.hit())
                        {
                            // Dump
                            if (featureEdgeStr.valid())
                            {
                                featureEdgeStr().write
                                (
                                    linePointRef(pt, info.hitPoint())
                                );
                            }
                        }
                        else
                        {
                            if (missedEdgeStr.valid())
                            {
                                missedEdgeStr().write
                                (
                                    linePointRef(pt, info.missPoint())
                                );
                            }
                        }
                    }
                }
            }
            else if (patchConstraints[pointI].first() == 2)
            {
                // Mark point on the nearest feature edge. Note that we
                // only search within the surrounding since the plane
                // reconstruction might find a feature where there isn't one.
                const point estimatedPt(pt + patchAttraction[pointI]);

                Tuple2<label, pointIndexHit> nearInfo(-1, pointIndexHit());

                // Geometric feature edge. Check for different patchIDs
                bool hasSnapped = false;
                if (multiRegionFeatureSnap)
                {
                    pointIndexHit multiPatchPt
                    (
                        findMultiPatchPoint
                        (
                            estimatedPt,
                            pointFacePatchID[pointI],
                            surfaceNormals,
                            faceToNormalBin
                        )
                    );
                    if (multiPatchPt.hit())
                    {
                        if (multiPatchPt.index() == 0)
                        {
                            // Geometric feature edge is also region edge
                            nearInfo = findNearFeatureEdge
                            (
                                true,               // isRegionPoint
                                pp,
                                snapDist,
                                pointI,
                                estimatedPt,

                                edgeAttractors,
                                edgeConstraints,

                                patchAttraction,
                                patchConstraints
                            );
                            hasSnapped = true;
                        }
                        else
                        {
                            // One of planes of feature contains multiple
                            // regions
                            nearInfo = findNearFeaturePoint
                            (
                                true,           // isRegionPoint
                                pp,
                                snapDist,
                                pointI,
                                estimatedPt,

                                // Feature-point to pp point
                                pointAttractor,
                                pointConstraints,
                                // Feature-edge to pp point
                                edgeAttractors,
                                edgeConstraints,
                                // pp point to nearest feature
                                patchAttraction,
                                patchConstraints
                            );
                            hasSnapped = true;
                        }
                    }
                }

                if (!hasSnapped)
                {
                    // Determine nearest point on feature edge. Store
                    // constraint
                    // (calculated from feature edge, alternative would be to
                    //  use constraint calculated from both surfaceNormals)
                    nearInfo = findNearFeatureEdge
                    (
                        false,      // isRegionPoint
                        pp,
                        snapDist,
                        pointI,
                        estimatedPt,

                        edgeAttractors,
                        edgeConstraints,

                        patchAttraction,
                        patchConstraints
                    );
                    hasSnapped = true;
                }

                // Dump to obj
                const pointIndexHit& info = nearInfo.second();
                if (info.hit())
                {
                    if
                    (
                        patchConstraints[pointI].first() == 3
                     && featurePointStr.valid()
                    )
                    {
                        featurePointStr().write
                        (
                            linePointRef(pt, info.hitPoint())
                        );
                    }
                    else if
                    (
                        patchConstraints[pointI].first() == 2
                     && featureEdgeStr.valid()
                    )
                    {
                        featureEdgeStr().write
                        (
                            linePointRef(pt, info.hitPoint())
                        );
                    }
                }
                else
                {
                    if (missedEdgeStr.valid())
                    {
                        missedEdgeStr().write
                        (
                            linePointRef(pt, info.missPoint())
                        );
                    }
                }
            }
            else if (patchConstraints[pointI].first() == 3)
            {
                // Mark point on the nearest feature point.
                const point estimatedPt(pt + patchAttraction[pointI]);

                Tuple2<label, pointIndexHit> nearInfo(-1, pointIndexHit());

                if (multiRegionFeatureSnap)
                {
                    pointIndexHit multiPatchPt
                    (
                        findMultiPatchPoint
                        (
                            estimatedPt,
                            pointFacePatchID[pointI],
                            surfaceNormals,
                            faceToNormalBin
                        )
                    );
                    if (multiPatchPt.hit())
                    {
                        // Multiple regions
                        nearInfo = findNearFeaturePoint
                        (
                            true,           // isRegionPoint
                            pp,
                            snapDist,
                            pointI,
                            estimatedPt,

                            // Feature-point to pp point
                            pointAttractor,
                            pointConstraints,
                            // Feature-edge to pp point
                            edgeAttractors,
                            edgeConstraints,
                            // pp point to nearest feature
                            patchAttraction,
                            patchConstraints
                        );
                    }
                    else
                    {
                        nearInfo = findNearFeaturePoint
                        (
                            false,              // isRegionPoint
                            pp,
                            snapDist,
                            pointI,
                            estimatedPt,

                            // Feature-point to pp point
                            pointAttractor,
                            pointConstraints,
                            // Feature-edge to pp point
                            edgeAttractors,
                            edgeConstraints,
                            // pp point to nearest feature
                            patchAttraction,
                            patchConstraints
                        );
                    }
                }
                else
                {
                    // No multi-patch snapping
                    nearInfo = findNearFeaturePoint
                    (
                        false,              // isRegionPoint
                        pp,
                        snapDist,
                        pointI,
                        estimatedPt,

                        // Feature-point to pp point
                        pointAttractor,
                        pointConstraints,
                        // Feature-edge to pp point
                        edgeAttractors,
                        edgeConstraints,
                        // pp point to nearest feature
                        patchAttraction,
                        patchConstraints
                    );
                }

                const pointIndexHit& info = nearInfo.second();
                if (info.hit() && featurePointStr.valid())
                {
                    featurePointStr().write
                    (
                        linePointRef(pt, info.hitPoint())
                    );
                }
            }
        }
    }
}


// Baffle handling
// ~~~~~~~~~~~~~~~
// Override pointAttractor, edgeAttractor, patchAttration etc. to
// implement 'baffle' handling.
// Baffle: the mesh pp point originates from a loose standing
// baffle.
// Sampling the surface with the surrounding face-centres only picks up
// a single triangle normal so above determineFeatures will not have
// detected anything. So explicitly pick up feature edges on the pp
// (after duplicating points & smoothing so will already have been
// expanded) and match these to the features.
void Foam::autoSnapDriver::determineBaffleFeatures
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,

    // Feature-point to pp point
    List<labelList>& pointAttractor,
    List<List<pointConstraint> >& pointConstraints,
    // Feature-edge to pp point
    List<List<DynamicList<point> > >& edgeAttractors,
    List<List<DynamicList<pointConstraint> > >& edgeConstraints,
    // pp point to nearest feature
    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const refinementFeatures& features = meshRefiner_.features();

    // Calculate edge-faces
    List<List<point> > edgeFaceNormals(pp.nEdges());

    // Fill local data
    forAll(pp.edgeFaces(), edgeI)
    {
        const labelList& eFaces = pp.edgeFaces()[edgeI];
        List<point>& eFc = edgeFaceNormals[edgeI];
        eFc.setSize(eFaces.size());
        forAll(eFaces, i)
        {
            label faceI = eFaces[i];
            eFc[i] = pp.faceNormals()[faceI];
        }
    }

    {
        // Precalculate mesh edges for pp.edges.
        const labelList meshEdges
        (
            pp.meshEdges(mesh.edges(), mesh.pointEdges())
        );
        syncTools::syncEdgeList
        (
            mesh,
            meshEdges,
            edgeFaceNormals,
            listPlusEqOp<point>(),
            List<point>(),
            mapDistribute::transform()
        );
    }

    // Detect baffle edges. Assume initial mesh will have 0,90 or 180
    // (baffle) degree angles so smoothing should make 0,90
    // to be less than 90.
    const scalar baffleFeatureCos = Foam::cos(degToRad(91));


    autoPtr<OBJstream> baffleEdgeStr;
    if (debug&meshRefinement::ATTRACTION)
    {
        baffleEdgeStr.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "baffleEdge_" + name(iter) + ".obj"
            )
        );
        Info<< nl << "Dumping baffle-edges to "
            << baffleEdgeStr().name() << endl;
    }


    // Is edge on baffle
    PackedBoolList isBaffleEdge(pp.nEdges());
    label nBaffleEdges = 0;
    // Is point on
    //  0 : baffle-edge (0)
    //  1 : baffle-feature-point (1)
    // -1 : rest
    labelList pointStatus(pp.nPoints(), -1);

    forAll(edgeFaceNormals, edgeI)
    {
        const List<point>& efn = edgeFaceNormals[edgeI];

        if (efn.size() == 2 && (efn[0]&efn[1]) < baffleFeatureCos)
        {
            isBaffleEdge[edgeI] = true;
            nBaffleEdges++;
            const edge& e = pp.edges()[edgeI];
            pointStatus[e[0]] = 0;
            pointStatus[e[1]] = 0;

            if (baffleEdgeStr.valid())
            {
                const point& p0 = pp.localPoints()[e[0]];
                const point& p1 = pp.localPoints()[e[1]];
                baffleEdgeStr().write(linePointRef(p0, p1));
            }
        }
    }

    Info<< "Detected "
        << returnReduce(nBaffleEdges, sumOp<label>())
        << " baffle edges out of "
        << returnReduce(pp.nEdges(), sumOp<label>())
        << " edges." << endl;


    forAll(pp.pointEdges(), pointI)
    {
        if
        (
            isFeaturePoint
            (
                featureCos,
                pp,
                isBaffleEdge,
                pointI
            )
        )
        {
            //Pout<< "Detected feature point:" << pp.localPoints()[pointI]
            //    << endl;
            //-TEMPORARILY DISABLED:
            //pointStatus[pointI] = 1;
        }
    }


    forAll(pointStatus, pointI)
    {
        const point& pt = pp.localPoints()[pointI];

        if (pointStatus[pointI] == 0)   // baffle edge
        {
            Tuple2<label, pointIndexHit> nearInfo = findNearFeatureEdge
            (
                false,          // isRegionPoint?
                pp,
                snapDist,
                pointI,
                pt,

                edgeAttractors,
                edgeConstraints,
                patchAttraction,
                patchConstraints
            );

            if (!nearInfo.second().hit())
            {
                //Pout<< "*** Failed to find close edge to point " << pt
                //    << endl;
            }
        }
        else if (pointStatus[pointI] == 1)   // baffle point
        {
            labelList nearFeat;
            List<pointIndexHit> nearInfo;
            features.findNearestPoint
            (
                pointField(1, pt),
                scalarField(1, sqr(snapDist[pointI])),
                nearFeat,
                nearInfo
            );

            label featI = nearFeat[0];

            if (featI != -1)
            {
                label featPointI = nearInfo[0].index();
                const point& featPt = nearInfo[0].hitPoint();
                scalar distSqr = magSqr(featPt-pt);

                // Check if already attracted
                label oldPointI = pointAttractor[featI][featPointI];

                if
                (
                    oldPointI == -1
                 || (
                        distSqr
                      < magSqr(featPt-pp.localPoints()[oldPointI])
                    )
                )
                {
                    pointAttractor[featI][featPointI] = pointI;
                    pointConstraints[featI][featPointI].first() = 3;
                    pointConstraints[featI][featPointI].second() =
                        vector::zero;

                    // Store for later use
                    patchAttraction[pointI] = featPt-pt;
                    patchConstraints[pointI] =
                        pointConstraints[featI][featPointI];

                    if (oldPointI != -1)
                    {
                        // The current point is closer so wins. Reset
                        // the old point to attract to nearest edge
                        // instead.
                        findNearFeatureEdge
                        (
                            false,              // isRegionPoint
                            pp,
                            snapDist,
                            oldPointI,
                            pp.localPoints()[oldPointI],

                            edgeAttractors,
                            edgeConstraints,
                            patchAttraction,
                            patchConstraints
                        );
                    }
                }
                else
                {
                    // Make it fall through to check below
                    featI = -1;
                }
            }

            // Not found a feature point or another point is already
            // closer to that feature
            if (featI == -1)
            {
                //Pout<< "*** Falling back to finding nearest feature"
                //    << " edge"
                //    << " for baffle-feature-point " << pt
                //    << endl;

                findNearFeatureEdge
                (
                    false,                  // isRegionPoint
                    pp,
                    snapDist,
                    pointI,
                    pt,                     // starting point

                    edgeAttractors,
                    edgeConstraints,
                    patchAttraction,
                    patchConstraints
                );
            }
        }
    }
}


void Foam::autoSnapDriver::reverseAttractMeshPoints
(
    const label iter,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,

    // Feature-point to pp point
    const List<labelList>& pointAttractor,
    const List<List<pointConstraint> >& pointConstraints,
    // Feature-edge to pp point
    const List<List<DynamicList<point> > >& edgeAttractors,
    const List<List<DynamicList<pointConstraint> > >& edgeConstraints,

    const vectorField& rawPatchAttraction,
    const List<pointConstraint>& rawPatchConstraints,

    // pp point to nearest feature
    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    const refinementFeatures& features = meshRefiner_.features();

    // Find nearest mesh point to feature edge
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Reverse lookup : go through all edgeAttractors and find the
    // nearest point on pp

    // Get search domain and extend it a bit
    treeBoundBox bb(pp.localPoints());
    {
        // Random number generator. Bit dodgy since not exactly random ;-)
        Random rndGen(65431);

        // Slightly extended bb. Slightly off-centred just so on symmetric
        // geometry there are less face/edge aligned items.
        bb = bb.extend(rndGen, 1e-4);
        bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    }

    // Collect candidate points for attraction
    DynamicList<label> attractPoints(pp.nPoints());
    {
        const fvMesh& mesh = meshRefiner_.mesh();

        boolList isFeatureEdgeOrPoint(pp.nPoints(), false);
        label nFeats = 0;
        forAll(rawPatchConstraints, pointI)
        {
            if (rawPatchConstraints[pointI].first() >= 2)
            {
                isFeatureEdgeOrPoint[pointI] = true;
                nFeats++;
            }
        }

        Info<< "Initially selected " << returnReduce(nFeats, sumOp<label>())
            << " points out of " << returnReduce(pp.nPoints(), sumOp<label>())
            << " for reverse attraction." << endl;

        // Make sure is synchronised (note: check if constraint is already
        // synced in which case this is not needed here)
        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            isFeatureEdgeOrPoint,
            orEqOp<bool>(),         // combine op
            false
        );

        for (label nGrow = 0; nGrow < 1; nGrow++)
        {
            boolList newIsFeatureEdgeOrPoint(isFeatureEdgeOrPoint);

            forAll(pp.localFaces(), faceI)
            {
                const face& f = pp.localFaces()[faceI];

                forAll(f, fp)
                {
                    if (isFeatureEdgeOrPoint[f[fp]])
                    {
                        // Mark all points on face
                        forAll(f, fp)
                        {
                            newIsFeatureEdgeOrPoint[f[fp]] = true;
                        }
                        break;
                    }
                }
            }

            isFeatureEdgeOrPoint = newIsFeatureEdgeOrPoint;

            syncTools::syncPointList
            (
                mesh,
                pp.meshPoints(),
                isFeatureEdgeOrPoint,
                orEqOp<bool>(),         // combine op
                false
            );
        }


        // Collect attractPoints
        forAll(isFeatureEdgeOrPoint, pointI)
        {
            if (isFeatureEdgeOrPoint[pointI])
            {
                attractPoints.append(pointI);
            }
        }

        Info<< "Selected "
            << returnReduce(attractPoints.size(), sumOp<label>())
            << " points out of " << returnReduce(pp.nPoints(), sumOp<label>())
            << " for reverse attraction." << endl;
    }


    indexedOctree<treeDataPoint> ppTree
    (
        treeDataPoint(pp.localPoints(), attractPoints),
        bb,                             // overall search domain
        8,                              // maxLevel
        10,                             // leafsize
        3.0                             // duplicity
    );

    // Per mesh point the point on nearest feature edge.
    patchAttraction.setSize(pp.nPoints());
    patchAttraction = vector::zero;
    patchConstraints.setSize(pp.nPoints());
    patchConstraints = pointConstraint();

    forAll(edgeAttractors, featI)
    {
        const List<DynamicList<point> >& edgeAttr = edgeAttractors[featI];
        const List<DynamicList<pointConstraint> >& edgeConstr =
            edgeConstraints[featI];

        forAll(edgeAttr, featEdgeI)
        {
            const DynamicList<point>& attr = edgeAttr[featEdgeI];
            forAll(attr, i)
            {
                // Find nearest pp point
                const point& featPt = attr[i];
                pointIndexHit nearInfo = ppTree.findNearest
                (
                    featPt,
                    sqr(GREAT)
                );

                if (nearInfo.hit())
                {
                    label pointI =
                        ppTree.shapes().pointLabels()[nearInfo.index()];
                    const point attraction = featPt-pp.localPoints()[pointI];

                    // Check if this point is already being attracted. If so
                    // override it only if nearer.
                    if
                    (
                        patchConstraints[pointI].first() <= 1
                     || magSqr(attraction) < magSqr(patchAttraction[pointI])
                    )
                    {
                        patchAttraction[pointI] = attraction;
                        patchConstraints[pointI] = edgeConstr[featEdgeI][i];
                    }
                }
                else
                {
                    WarningIn
                    (
                        "autoSnapDriver::featureAttractionUsingFeatureEdges"
                        "(..)"
                    )   << "Did not find pp point near " << featPt
                        << endl;
                }
            }
        }
    }


    // Different procs might have different patchAttraction,patchConstraints
    // however these only contain geometric information, no topology
    // so as long as we synchronise after overriding with feature points
    // there is no problem, just possibly a small error.


    // Find nearest mesh point to feature point
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // (overrides attraction to feature edge)
    forAll(pointAttractor, featI)
    {
        const labelList& pointAttr = pointAttractor[featI];
        const List<pointConstraint>& pointConstr = pointConstraints[featI];

        forAll(pointAttr, featPointI)
        {
            if (pointAttr[featPointI] != -1)
            {
                const point& featPt = features[featI].points()
                [
                    featPointI
                ];

                // Find nearest pp point
                pointIndexHit nearInfo = ppTree.findNearest
                (
                    featPt,
                    sqr(GREAT)
                );

                if (nearInfo.hit())
                {
                    label pointI =
                        ppTree.shapes().pointLabels()[nearInfo.index()];

                    const point& pt = pp.localPoints()[pointI];
                    const point attraction = featPt-pt;

                    // - already attracted to feature edge : point always wins
                    // - already attracted to feature point: nearest wins

                    if (patchConstraints[pointI].first() <= 1)
                    {
                        patchAttraction[pointI] = attraction;
                        patchConstraints[pointI] = pointConstr[featPointI];
                    }
                    else if (patchConstraints[pointI].first() == 2)
                    {
                        patchAttraction[pointI] = attraction;
                        patchConstraints[pointI] = pointConstr[featPointI];
                    }
                    else if (patchConstraints[pointI].first() == 3)
                    {
                        // Only if nearer
                        if
                        (
                            magSqr(attraction)
                          < magSqr(patchAttraction[pointI])
                        )
                        {
                            patchAttraction[pointI] = attraction;
                            patchConstraints[pointI] =
                                pointConstr[featPointI];
                        }
                    }
                }
            }
        }
    }
}


void Foam::autoSnapDriver::featureAttractionUsingFeatureEdges
(
    const label iter,
    const bool avoidSnapProblems,
    const scalar featureCos,
    const bool multiRegionFeatureSnap,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,
    const vectorField& nearestDisp,

    const List<List<point> >& pointFaceSurfNormals,
    const List<List<point> >& pointFaceDisp,
    const List<List<point> >& pointFaceCentres,
    const labelListList& pointFacePatchID,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    const refinementFeatures& features = meshRefiner_.features();

    // Collect ordered attractions on feature edges
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Per feature, per feature-edge a list of attraction points and their
    // originating vertex.
    List<List<DynamicList<point> > > edgeAttractors(features.size());
    List<List<DynamicList<pointConstraint> > > edgeConstraints
    (
        features.size()
    );
    forAll(features, featI)
    {
        label nFeatEdges = features[featI].edges().size();
        edgeAttractors[featI].setSize(nFeatEdges);
        edgeConstraints[featI].setSize(nFeatEdges);
    }

    // Per feature, per feature-point the pp point that is attracted to it.
    // This list is only used to subset the feature-points that are actually
    // used.
    List<labelList> pointAttractor(features.size());
    List<List<pointConstraint> > pointConstraints(features.size());
    forAll(features, featI)
    {
        label nFeatPoints = features[featI].points().size();
        pointAttractor[featI].setSize(nFeatPoints, -1);
        pointConstraints[featI].setSize(nFeatPoints);
    }

    // Reverse: from pp point to nearest feature
    vectorField rawPatchAttraction(pp.nPoints(), vector::zero);
    List<pointConstraint> rawPatchConstraints(pp.nPoints());

    determineFeatures
    (
        iter,
        featureCos,
        multiRegionFeatureSnap,

        pp,
        snapDist,               // per point max distance and nearest surface
        nearestDisp,

        pointFaceSurfNormals,   // per face nearest surface
        pointFaceDisp,
        pointFaceCentres,
        pointFacePatchID,

        // Feature-point to pp point
        pointAttractor,
        pointConstraints,
        // Feature-edge to pp point
        edgeAttractors,
        edgeConstraints,
        // pp point to nearest feature
        rawPatchAttraction,
        rawPatchConstraints
    );



    // Baffle handling
    // ~~~~~~~~~~~~~~~
    // Override pointAttractor, edgeAttractor, rawPatchAttration etc. to
    // implement 'baffle' handling.
    // Baffle: the mesh pp point originates from a loose standing
    // baffle.
    // Sampling the surface with the surrounding face-centres only picks up
    // a single triangle normal so above determineFeatures will not have
    // detected anything. So explicitly pick up feature edges on the pp
    // (after duplicating points & smoothing so will already have been
    // expanded) and match these to the features.
    determineBaffleFeatures
    (
        iter,
        featureCos,

        pp,
        snapDist,

        // Feature-point to pp point
        pointAttractor,
        pointConstraints,
        // Feature-edge to pp point
        edgeAttractors,
        edgeConstraints,
        // pp point to nearest feature
        rawPatchAttraction,
        rawPatchConstraints
    );



    // Reverse lookup: Find nearest mesh point to feature edge
    // ~~~~~~~~~~~~~~~~----------------~~~~~~~~~~~~~~~~~~~~~~~
    // go through all edgeAttractors and find the nearest point on pp

    reverseAttractMeshPoints
    (
        iter,

        pp,
        snapDist,

        // Feature-point to pp point
        pointAttractor,
        pointConstraints,
        // Feature-edge to pp point
        edgeAttractors,
        edgeConstraints,

        // Estimated feature point
        rawPatchAttraction,
        rawPatchConstraints,

        // pp point to nearest feature
        patchAttraction,
        patchConstraints
    );


    // Dump
    if (debug&meshRefinement::ATTRACTION)
    {
        OBJstream featureEdgeStr
        (
            meshRefiner_.mesh().time().path()
          / "edgeAttractors_" + name(iter) + ".obj"
        );
        Info<< "Dumping feature-edge attraction to "
            << featureEdgeStr.name() << endl;

        OBJstream featurePointStr
        (
            meshRefiner_.mesh().time().path()
          / "pointAttractors_" + name(iter) + ".obj"
        );
        Info<< "Dumping feature-point attraction to "
            << featurePointStr.name() << endl;

        forAll(patchConstraints, pointI)
        {
            const point& pt = pp.localPoints()[pointI];
            const vector& attr = patchAttraction[pointI];

            if (patchConstraints[pointI].first() == 2)
            {
                featureEdgeStr.write(linePointRef(pt, pt+attr));
            }
            else if (patchConstraints[pointI].first() == 3)
            {
                featurePointStr.write(linePointRef(pt, pt+attr));
            }
        }
    }


    if (avoidSnapProblems)
    {
        //MEJ: any faces that have multi-patch points only keep the multi-patch
        //     points. The other points on the face will be dragged along
        //     (hopefully)
        if (multiRegionFeatureSnap)
        {
            releasePointsNextToMultiPatch
            (
                iter,
                featureCos,

                pp,
                snapDist,

                pointFaceCentres,
                pointFacePatchID,

                rawPatchAttraction,
                rawPatchConstraints,

                patchAttraction,
                patchConstraints
            );
        }


        // Snap edges to feature edges
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Walk existing edges and snap remaining ones (that are marked as
        // feature edges in rawPatchConstraints)

        stringFeatureEdges
        (
            iter,
            featureCos,

            pp,
            snapDist,

            rawPatchAttraction,
            rawPatchConstraints,

            patchAttraction,
            patchConstraints
        );


        // Avoid diagonal attraction
        // ~~~~~~~~~~~~~~~~~~~~~~~~~
        // Attract one of the non-diagonal points.

        avoidDiagonalAttraction
        (
            iter,
            featureCos,
            pp,
            patchAttraction,
            patchConstraints
        );
    }

    if (debug&meshRefinement::ATTRACTION)
    {
        dumpMove
        (
            meshRefiner_.mesh().time().path()
          / "patchAttraction_" + name(iter) + ".obj",
            pp.localPoints(),
            pp.localPoints() + patchAttraction
        );
    }
}


// Correct for squeezing of face
void Foam::autoSnapDriver::preventFaceSqueeze
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    pointField points;
    face singleF;
    forAll(pp.localFaces(), faceI)
    {
        const face& f = pp.localFaces()[faceI];

        if (f.size() != points.size())
        {
            points.setSize(f.size());
            singleF.setSize(f.size());
            for (label i = 0; i < f.size(); i++)
            {
                singleF[i] = i;
            }
        }
        label nConstraints = 0;
        forAll(f, fp)
        {
            label pointI = f[fp];
            const point& pt = pp.localPoints()[pointI];

            if (patchConstraints[pointI].first() > 1)
            {
                points[fp] = pt + patchAttraction[pointI];
                nConstraints++;
            }
            else
            {
                points[fp] = pt;
            }
        }

        if (nConstraints == f.size())
        {
            scalar oldArea = f.mag(pp.localPoints());
            scalar newArea = singleF.mag(points);
            if (newArea < 0.1*oldArea)
            {
                // For now remove the point with largest distance
                label maxFp = -1;
                scalar maxS = -1;
                forAll(f, fp)
                {
                    scalar s = magSqr(patchAttraction[f[fp]]);
                    if (s > maxS)
                    {
                        maxS = s;
                        maxFp = fp;
                    }
                }
                if (maxFp != -1)
                {
                    label pointI = f[maxFp];
                    // Lower attraction on pointI
                    patchAttraction[pointI] *= 0.5;
                }
            }
        }
    }
}


Foam::vectorField Foam::autoSnapDriver::calcNearestSurfaceFeature
(
    const snapParameters& snapParams,
    const bool avoidSnapProblems,
    const label iter,
    const scalar featureCos,
    const scalar featureAttract,
    const scalarField& snapDist,
    const vectorField& nearestDisp,
    motionSmoother& meshMover,
    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    const Switch implicitFeatureAttraction = snapParams.implicitFeatureSnap();
    const Switch explicitFeatureAttraction = snapParams.explicitFeatureSnap();
    const Switch multiRegionFeatureSnap = snapParams.multiRegionFeatureSnap();

    Info<< "Overriding displacement on features :" << nl
        << "   implicit features    : " << implicitFeatureAttraction << nl
        << "   explicit features    : " << explicitFeatureAttraction << nl
        << "   multi-patch features : " << multiRegionFeatureSnap << nl
        << endl;


    const indirectPrimitivePatch& pp = meshMover.patch();
    const pointField& localPoints = pp.localPoints();
    const fvMesh& mesh = meshRefiner_.mesh();



    // Per point, per surrounding face:
    // - faceSurfaceNormal
    // - faceDisp
    // - faceCentres
    List<List<point> > pointFaceSurfNormals;
    List<List<point> > pointFaceDisp;
    List<List<point> > pointFaceCentres;
    List<labelList> pointFacePatchID;

    {
        // Calculate attraction distance per face (from the attraction distance
        // per point)
        scalarField faceSnapDist(pp.size(), -GREAT);
        forAll(pp.localFaces(), faceI)
        {
            const face& f = pp.localFaces()[faceI];
            forAll(f, fp)
            {
                faceSnapDist[faceI] = max(faceSnapDist[faceI], snapDist[f[fp]]);
            }
        }


        // Displacement and orientation per pp face
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // vector from point on surface back to face centre
        vectorField faceDisp(pp.size(), vector::zero);
        // normal of surface at point on surface
        vectorField faceSurfaceNormal(pp.size(), vector::zero);
        labelList faceSurfaceGlobalRegion(pp.size(), -1);
        vectorField faceRotation(pp.size(), vector::zero);

        calcNearestFace
        (
            iter,
            pp,
            faceSnapDist,
            faceDisp,
            faceSurfaceNormal,
            faceSurfaceGlobalRegion,
            faceRotation
        );


        // Collect (possibly remote) per point data of all surrounding faces
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // - faceSurfaceNormal
        // - faceDisp
        // - faceCentres
        calcNearestFacePointProperties
        (
            iter,
            pp,

            faceDisp,
            faceSurfaceNormal,
            faceSurfaceGlobalRegion,

            pointFaceSurfNormals,
            pointFaceDisp,
            pointFaceCentres,
            pointFacePatchID
        );
    }


    // Start off with nearest point on surface
    vectorField patchDisp = nearestDisp;


    // Main calculation
    // ~~~~~~~~~~~~~~~~
    // This is the main intelligence which calculates per point the vector to
    // attract it to the nearest surface. There are lots of possibilities
    // here.

    // Nearest feature
    patchAttraction.setSize(localPoints.size());
    patchAttraction = vector::zero;
    // Constraints at feature
    patchConstraints.setSize(localPoints.size());
    patchConstraints = pointConstraint();

    if (implicitFeatureAttraction)
    {
        // Sample faces around each point and see if nearest surface normal
        // differs. Reconstruct a feature edge/point if possible and snap to
        // it.
        featureAttractionUsingReconstruction
        (
            iter,
            avoidSnapProblems,
            featureCos,

            pp,
            snapDist,
            nearestDisp,

            pointFaceSurfNormals,
            pointFaceDisp,
            pointFaceCentres,
            pointFacePatchID,

            patchAttraction,
            patchConstraints
        );
    }

    if (explicitFeatureAttraction)
    {
        // Sample faces around each point and see if nearest surface normal
        // differs. For those find the nearest real feature edge/point and
        // store the correspondence. Then loop over feature edge/point
        // and attract those nearest mesh point. (the first phase just is
        // a subsetting of candidate points, the second makes sure that only
        // one mesh point gets attracted per feature)
        featureAttractionUsingFeatureEdges
        (
            iter,
            avoidSnapProblems,
            featureCos,
            multiRegionFeatureSnap,

            pp,
            snapDist,
            nearestDisp,

            pointFaceSurfNormals,
            pointFaceDisp,
            pointFaceCentres,
            pointFacePatchID,

            patchAttraction,
            patchConstraints
        );
    }

    preventFaceSqueeze
    (
        iter,
        featureCos,

        pp,
        snapDist,

        patchAttraction,
        patchConstraints
    );

    //const PackedBoolList isMasterPoint(syncTools::getMasterPoints(mesh));
    const PackedBoolList isPatchMasterPoint
    (
        meshRefinement::getMasterPoints
        (
            mesh,
            pp.meshPoints()
        )
    );
    {
        vector avgPatchDisp = meshRefinement::gAverage
        (
            isPatchMasterPoint,
            patchDisp
        );
        vector avgPatchAttr = meshRefinement::gAverage
        (
            isPatchMasterPoint,
            patchAttraction
        );

        Info<< "Attraction:" << endl
            << "     linear   : max:" << gMaxMagSqr(patchDisp)
            << " avg:" << avgPatchDisp << endl
            << "     feature  : max:" << gMaxMagSqr(patchAttraction)
            << " avg:" << avgPatchAttr << endl;
    }

    // So now we have:
    // - patchDisp          : point movement to go to nearest point on surface
    //                       (either direct or through interpolation of
    //                        face nearest)
    // - patchAttraction    : direct attraction to features
    // - patchConstraints   : type of features

    // Use any combination of patchDisp and direct feature attraction.


    // Mix with direct feature attraction
    forAll(patchConstraints, pointI)
    {
        if (patchConstraints[pointI].first() > 1)
        {
            patchDisp[pointI] =
                (1.0-featureAttract)*patchDisp[pointI]
              + featureAttract*patchAttraction[pointI];
        }
    }



    // Count
    {
        label nMasterPoints = 0;
        label nPlanar = 0;
        label nEdge = 0;
        label nPoint = 0;

        forAll(patchConstraints, pointI)
        {
            if (isPatchMasterPoint[pointI])
            {
                nMasterPoints++;

                if (patchConstraints[pointI].first() == 1)
                {
                    nPlanar++;
                }
                else if (patchConstraints[pointI].first() == 2)
                {
                    nEdge++;
                }
                else if (patchConstraints[pointI].first() == 3)
                {
                    nPoint++;
                }
            }
        }

        reduce(nMasterPoints, sumOp<label>());
        reduce(nPlanar, sumOp<label>());
        reduce(nEdge, sumOp<label>());
        reduce(nPoint, sumOp<label>());
        Info<< "Feature analysis : total master points:"
            << nMasterPoints
            << " attraction to :" << nl
            << "    feature point   : " << nPoint << nl
            << "    feature edge    : " << nEdge << nl
            << "    nearest surface : " << nPlanar << nl
            << "    rest            : " << nMasterPoints-nPoint-nEdge-nPlanar
            << nl
            << endl;
    }


    // Now we have the displacement per patch point to move onto the surface
    // Split into tangential and normal direction.
    // - start off with all non-constrained points following the constrained
    //   ones since point normals not relevant.
    // - finish with only tangential component smoothed.
    // Note: tangential is most
    // likely to come purely from face-centre snapping, not face rotation.
    // Note: could use the constraints here (constraintTransformation())
    //       but this is not necessarily accurate and we're smoothing to
    //       get out of problems.

    if (featureAttract < 1-0.001)
    {
        //const PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh));
        const labelList meshEdges
        (
            pp.meshEdges(mesh.edges(), mesh.pointEdges())
        );
        const PackedBoolList isPatchMasterEdge
        (
            meshRefinement::getMasterEdges
            (
                mesh,
                meshEdges
            )
        );

        const vectorField pointNormals
        (
            PatchTools::pointNormals
            (
                mesh,
                pp
            )
        );


        // 1. Smoothed all displacement
        vectorField smoothedPatchDisp = patchDisp;
        smoothAndConstrain
        (
            isPatchMasterEdge,
            pp,
            meshEdges,
            patchConstraints,
            smoothedPatchDisp
        );


        // 2. Smoothed tangential component
        vectorField tangPatchDisp = patchDisp;
        tangPatchDisp -= (pointNormals & patchDisp) * pointNormals;
        smoothAndConstrain
        (
            isPatchMasterEdge,
            pp,
            meshEdges,
            patchConstraints,
            tangPatchDisp
        );

        // Re-add normal component
        tangPatchDisp += (pointNormals & patchDisp) * pointNormals;

        if (debug&meshRefinement::ATTRACTION)
        {
            dumpMove
            (
                mesh.time().path()
              / "tangPatchDispConstrained_" + name(iter) + ".obj",
                pp.localPoints(),
                pp.localPoints() + tangPatchDisp
            );
        }

        patchDisp =
             (1.0-featureAttract)*smoothedPatchDisp
           + featureAttract*tangPatchDisp;
    }


    const scalar relax = featureAttract;
    patchDisp *= relax;


    // Points on zones in one domain but only present as point on other
    // will not do condition 2 on all. Sync explicitly.
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        patchDisp,
        minMagSqrEqOp<point>(),         // combine op
        vector(GREAT, GREAT, GREAT)     // null value (note: cant use VGREAT)
    );

    return patchDisp;
}


// ************************************************************************* //
