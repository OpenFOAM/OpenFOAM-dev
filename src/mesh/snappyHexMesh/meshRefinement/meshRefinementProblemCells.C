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

#include "meshRefinement.H"
#include "fvMesh.H"
#include "syncTools.H"
#include "Time.H"
#include "refinementSurfaces.H"
#include "pointSet.H"
#include "faceSet.H"
#include "indirectPrimitivePatch.H"
#include "cellSet.H"
#include "searchableSurfaces.H"
#include "polyMeshGeometry.H"
#include "IOmanip.H"
#include "unitConversion.H"
#include "snappySnapDriver.H"

#include "snapParameters.H"
#include "motionSmoother.H"
#include "topoDistanceData.H"
#include "FaceCellWave.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::meshRefinement::markBoundaryFace
(
    const label facei,
    boolList& isBoundaryFace,
    boolList& isBoundaryEdge,
    boolList& isBoundaryPoint
) const
{
    isBoundaryFace[facei] = true;

    const labelList& fEdges = mesh_.faceEdges(facei);

    forAll(fEdges, fp)
    {
        isBoundaryEdge[fEdges[fp]] = true;
    }

    const face& f = mesh_.faces()[facei];

    forAll(f, fp)
    {
        isBoundaryPoint[f[fp]] = true;
    }
}


void Foam::meshRefinement::findNearest
(
    const labelList& meshFaces,
    List<pointIndexHit>& nearestInfo,
    labelList& nearestSurface,
    labelList& nearestRegion,
    vectorField& nearestNormal
) const
{
    pointField fc(meshFaces.size());
    forAll(meshFaces, i)
    {
        fc[i] = mesh_.faceCentres()[meshFaces[i]];
    }

    const labelList allSurfaces(identity(surfaces_.surfaces().size()));

    surfaces_.findNearest
    (
        allSurfaces,
        fc,
        scalarField(fc.size(), sqr(great)),    // sqr of attraction
        nearestSurface,
        nearestInfo
    );

    // Do normal testing per surface.
    nearestNormal.setSize(nearestInfo.size());
    nearestRegion.setSize(nearestInfo.size());

    forAll(allSurfaces, surfI)
    {
        DynamicList<pointIndexHit> localHits;

        forAll(nearestSurface, i)
        {
            if (nearestSurface[i] == surfI)
            {
                localHits.append(nearestInfo[i]);
            }
        }

        label geomI = surfaces_.surfaces()[surfI];

        pointField localNormals;
        surfaces_.geometry()[geomI].getNormal(localHits, localNormals);

        labelList localRegion;
        surfaces_.geometry()[geomI].getRegion(localHits, localRegion);

        label localI = 0;
        forAll(nearestSurface, i)
        {
            if (nearestSurface[i] == surfI)
            {
                nearestNormal[i] = localNormals[localI];
                nearestRegion[i] = localRegion[localI];
                localI++;
            }
        }
    }
}


Foam::Map<Foam::label> Foam::meshRefinement::findEdgeConnectedProblemCells
(
    const scalarField& perpendicularAngle,
    const labelList& globalToPatch
) const
{
    // Construct addressing engine from all patches added for meshing.
    autoPtr<indirectPrimitivePatch> ppPtr
    (
        meshRefinement::makePatch
        (
            mesh_,
            meshedPatches()
        )
    );
    const indirectPrimitivePatch& pp = ppPtr();


    // 1. Collect faces to test
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    DynamicList<label> candidateFaces(pp.size()/20);

    const labelListList& edgeFaces = pp.edgeFaces();

    const labelList& cellLevel = meshCutter_.cellLevel();

    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        if (eFaces.size() == 2)
        {
            label face0 = pp.addressing()[eFaces[0]];
            label face1 = pp.addressing()[eFaces[1]];

            label cell0 = mesh_.faceOwner()[face0];
            label cell1 = mesh_.faceOwner()[face1];

            if (cellLevel[cell0] > cellLevel[cell1])
            {
                // cell0 smaller.
                const vector& n0 = pp.faceNormals()[eFaces[0]];
                const vector& n1 = pp.faceNormals()[eFaces[1]];

                if (mag(n0 & n1) < 0.1)
                {
                    candidateFaces.append(face0);
                }
            }
            else if (cellLevel[cell1] > cellLevel[cell0])
            {
                // cell1 smaller.
                const vector& n0 = pp.faceNormals()[eFaces[0]];
                const vector& n1 = pp.faceNormals()[eFaces[1]];

                if (mag(n0 & n1) < 0.1)
                {
                    candidateFaces.append(face1);
                }
            }
        }
    }
    candidateFaces.shrink();

    Info<< "Testing " << returnReduce(candidateFaces.size(), sumOp<label>())
        << " faces on edge-connected cells of differing level."
        << endl;

    if (debug&meshRefinement::MESH)
    {
        faceSet fSet(mesh_, "edgeConnectedFaces", candidateFaces);
        fSet.instance() = timeName();
        Pout<< "Writing " << fSet.size()
            << " with problematic topology to faceSet "
            << fSet.objectPath() << endl;
        fSet.write();
    }


    // 2. Find nearest surface on candidate faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    List<pointIndexHit> nearestInfo;
    labelList nearestSurface;
    labelList nearestRegion;
    vectorField nearestNormal;
    findNearest
    (
        candidateFaces,
        nearestInfo,
        nearestSurface,
        nearestRegion,
        nearestNormal
    );


    // 3. Test angle to surface
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    Map<label> candidateCells(candidateFaces.size());

    faceSet perpFaces(mesh_, "perpendicularFaces", pp.size()/100);

    forAll(candidateFaces, i)
    {
        label facei = candidateFaces[i];

        vector n = mesh_.faceAreas()[facei];
        n /= mag(n);

        label region = surfaces_.globalRegion
        (
            nearestSurface[i],
            nearestRegion[i]
        );

        scalar angle = degToRad(perpendicularAngle[region]);

        if (angle >= 0)
        {
            if (mag(n & nearestNormal[i]) < Foam::sin(angle))
            {
                perpFaces.insert(facei);
                candidateCells.insert
                (
                    mesh_.faceOwner()[facei],
                    globalToPatch[region]
                );
            }
        }
    }

    if (debug&meshRefinement::MESH)
    {
        perpFaces.instance() = timeName();
        Pout<< "Writing " << perpFaces.size()
            << " faces that are perpendicular to the surface to set "
            << perpFaces.objectPath() << endl;
        perpFaces.write();
    }
    return candidateCells;
}


// Check if moving face to new points causes a 'collapsed' face.
// Uses new point position only for the face, not the neighbouring
// cell centres
bool Foam::meshRefinement::isCollapsedFace
(
    const pointField& points,
    const pointField& neiCc,
    const scalar minFaceArea,
    const scalar maxNonOrtho,
    const label facei
) const
{
    // Severe nonorthogonality threshold
    const scalar severeNonorthogonalityThreshold =
        ::cos(degToRad(maxNonOrtho));

    const vector s = mesh_.faces()[facei].area(points);
    const scalar magS = mag(s);

    // Check face area
    if (magS < minFaceArea)
    {
        return true;
    }

    // Check orthogonality
    const point& ownCc = mesh_.cellCentres()[mesh_.faceOwner()[facei]];

    if (mesh_.isInternalFace(facei))
    {
        label nei = mesh_.faceNeighbour()[facei];
        vector d = mesh_.cellCentres()[nei] - ownCc;

        scalar dDotS = (d & s)/(mag(d)*magS + vSmall);

        if (dDotS < severeNonorthogonalityThreshold)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        label patchi = mesh_.boundaryMesh().whichPatch(facei);

        if (mesh_.boundaryMesh()[patchi].coupled())
        {
            vector d = neiCc[facei-mesh_.nInternalFaces()] - ownCc;

            scalar dDotS = (d & s)/(mag(d)*magS + vSmall);

            if (dDotS < severeNonorthogonalityThreshold)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            // Collapsing normal boundary face does not cause problems with
            // non-orthogonality
            return false;
        }
    }
}


// Check if moving cell to new points causes it to collapse.
bool Foam::meshRefinement::isCollapsedCell
(
    const pointField& points,
    const scalar volFraction,
    const label celli
) const
{
    scalar vol = mesh_.cells()[celli].mag(points, mesh_.faces());

    if (vol/mesh_.cellVolumes()[celli] < volFraction)
    {
        return true;
    }
    else
    {
        return false;
    }
}


Foam::labelList Foam::meshRefinement::nearestPatch
(
    const labelList& adaptPatchIDs
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    labelList nearestAdaptPatch;

    if (adaptPatchIDs.size())
    {
        nearestAdaptPatch.setSize(mesh_.nFaces(), adaptPatchIDs[0]);


        // Count number of faces in adaptPatchIDs
        label nFaces = 0;
        forAll(adaptPatchIDs, i)
        {
            const polyPatch& pp = patches[adaptPatchIDs[i]];
            nFaces += pp.size();
        }

        // Field on cells and faces.
        List<topoDistanceData> cellData(mesh_.nCells());
        List<topoDistanceData> faceData(mesh_.nFaces());

        // Start of changes
        labelList patchFaces(nFaces);
        List<topoDistanceData> patchData(nFaces);
        nFaces = 0;
        forAll(adaptPatchIDs, i)
        {
            label patchi = adaptPatchIDs[i];
            const polyPatch& pp = patches[patchi];

            forAll(pp, i)
            {
                patchFaces[nFaces] = pp.start()+i;
                patchData[nFaces] = topoDistanceData(patchi, 0);
                nFaces++;
            }
        }

        // Propagate information inwards
        FaceCellWave<topoDistanceData> deltaCalc
        (
            mesh_,
            patchFaces,
            patchData,
            faceData,
            cellData,
            mesh_.globalData().nTotalCells()+1
        );

        // And extract

        bool haveWarned = false;
        forAll(faceData, facei)
        {
            if (!faceData[facei].valid(deltaCalc.data()))
            {
                if (!haveWarned)
                {
                    WarningInFunction
                        << "Did not visit some faces, e.g. face " << facei
                        << " at " << mesh_.faceCentres()[facei] << endl
                        << "Assigning  these cells to patch "
                        << adaptPatchIDs[0]
                        << endl;
                    haveWarned = true;
                }
            }
            else
            {
                nearestAdaptPatch[facei] = faceData[facei].data();
            }
        }
    }
    else
    {
        // Use patch 0
        nearestAdaptPatch.setSize(mesh_.nFaces(), 0);
    }

    return nearestAdaptPatch;
}


// Returns list with for every internal face -1 or the patch they should
// be baffled into. Gets run after createBaffles so all the unzoned surface
// intersections have already been turned into baffles. (Note: zoned surfaces
// are not baffled at this stage)
// Used to remove cells by baffling all their faces and have the
// splitMeshRegions chuck away non used regions.
Foam::labelList Foam::meshRefinement::markFacesOnProblemCells
(
    const dictionary& motionDict,
    const bool removeEdgeConnectedCells,
    const scalarField& perpendicularAngle,
    const labelList& globalToPatch
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const labelList& pointLevel = meshCutter_.pointLevel();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();


    // Mark all points and edges on baffle patches (so not on any inlets,
    // outlets etc.)
    boolList isBoundaryPoint(mesh_.nPoints(), false);
    boolList isBoundaryEdge(mesh_.nEdges(), false);
    boolList isBoundaryFace(mesh_.nFaces(), false);

    // Fill boundary data. All elements on meshed patches get marked.
    // Get the labels of added patches.
    labelList adaptPatchIDs(meshedPatches());

    forAll(adaptPatchIDs, i)
    {
        const polyPatch& pp = patches[adaptPatchIDs[i]];

        label facei = pp.start();

        forAll(pp, j)
        {
            markBoundaryFace
            (
                facei,
                isBoundaryFace,
                isBoundaryEdge,
                isBoundaryPoint
            );

            facei++;
        }
    }


    // Per face the nearest adaptPatch
    const labelList nearestAdaptPatch(nearestPatch(adaptPatchIDs));


    // Per internal face (boundary faces not used) the patch that the
    // baffle should get (or -1)
    labelList facePatch(mesh_.nFaces(), -1);

    // Swap neighbouring cell centres and cell level
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
    calcNeighbourData(neiLevel, neiCc);


    // Count of faces marked for baffling
    label nBaffleFaces = 0;
    PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_));

    // Count of faces not baffled since would not cause a collapse
    label nPrevented = 0;

    if (removeEdgeConnectedCells && max(perpendicularAngle) >= 0)
    {
        Info<< "markFacesOnProblemCells :"
            << " Checking for edge-connected cells of highly differing sizes."
            << endl;

        // Pick up the cells that need to be removed and (a guess for)
        // the patch they should be patched with.
        Map<label> problemCells
        (
            findEdgeConnectedProblemCells
            (
                perpendicularAngle,
                globalToPatch
            )
        );

        // Baffle all faces of cells that need to be removed
        forAllConstIter(Map<label>, problemCells, iter)
        {
            const cell& cFaces = mesh_.cells()[iter.key()];

            forAll(cFaces, i)
            {
                label facei = cFaces[i];

                if (facePatch[facei] == -1 && mesh_.isInternalFace(facei))
                {
                    facePatch[facei] = nearestAdaptPatch[facei];
                    nBaffleFaces++;

                    // Mark face as a 'boundary'
                    markBoundaryFace
                    (
                        facei,
                        isBoundaryFace,
                        isBoundaryEdge,
                        isBoundaryPoint
                    );
                }
            }
        }
        Info<< "markFacesOnProblemCells : Marked "
            << returnReduce(nBaffleFaces, sumOp<label>())
            << " additional internal faces to be converted into baffles"
            << " due to "
            << returnReduce(problemCells.size(), sumOp<label>())
            << " cells edge-connected to lower level cells." << endl;

        if (debug&meshRefinement::MESH)
        {
            cellSet problemCellSet(mesh_, "problemCells", problemCells.toc());
            problemCellSet.instance() = timeName();
            Pout<< "Writing " << problemCellSet.size()
                << " cells that are edge connected to coarser cell to set "
                << problemCellSet.objectPath() << endl;
            problemCellSet.write();
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        isBoundaryPoint,
        orEqOp<bool>(),
        false               // null value
    );

    syncTools::syncEdgeList
    (
        mesh_,
        isBoundaryEdge,
        orEqOp<bool>(),
        false               // null value
    );

    syncTools::syncFaceList
    (
        mesh_,
        isBoundaryFace,
        orEqOp<bool>()
    );


    // See if checking for collapse
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Collapse checking parameters
    const scalar volFraction =
        motionDict.lookupOrDefault<scalar>("minVolCollapseRatio", -1);

    const bool checkCollapse = (volFraction > 0);
    scalar minArea = -1;
    scalar maxNonOrtho = -1;


    // Find nearest (non-baffle) surface
    pointField newPoints;

    if (checkCollapse)
    {
        minArea = readScalar(motionDict.lookup("minArea"));
        maxNonOrtho = readScalar(motionDict.lookup("maxNonOrtho"));

        Info<< "markFacesOnProblemCells :"
            << " Deleting all-anchor surface cells only if"
            << " snapping them violates mesh quality constraints:" << nl
            << "    snapped/original cell volume < " << volFraction << nl
            << "    face area                    < " << minArea << nl
            << "    non-orthogonality            > " << maxNonOrtho << nl
            << endl;

        // Construct addressing engine.
        autoPtr<indirectPrimitivePatch> ppPtr
        (
            meshRefinement::makePatch
            (
                mesh_,
                adaptPatchIDs
            )
        );
        const indirectPrimitivePatch& pp = ppPtr();
        const pointField& localPoints = pp.localPoints();
        const labelList& meshPoints = pp.meshPoints();

        List<pointIndexHit> hitInfo;
        labelList hitSurface;
        surfaces_.findNearest
        (
            surfaceZonesInfo::getUnnamedSurfaces(surfaces_.surfZones()),
            localPoints,
            scalarField(localPoints.size(), sqr(great)),    // sqr of attraction
            hitSurface,
            hitInfo
        );

        // Start off from current points
        newPoints = mesh_.points();

        forAll(hitInfo, i)
        {
            if (hitInfo[i].hit())
            {
                newPoints[meshPoints[i]] = hitInfo[i].hitPoint();
            }
        }

        if (debug&meshRefinement::MESH)
        {
            const_cast<Time&>(mesh_.time())++;
            pointField oldPoints(mesh_.points());
            mesh_.movePoints(newPoints);
            Pout<< "Writing newPoints mesh to time " << timeName()
                << endl;
            write
            (
                debugType(debug),
                writeType(writeLevel() | WRITEMESH),
                mesh_.time().path()/"newPoints"
            );
            mesh_.movePoints(oldPoints);
        }
    }



    // For each cell count the number of anchor points that are on
    // the boundary:
    // 8 : check the number of (baffle) boundary faces. If 3 or more block
    //     off the cell since the cell would get squeezed down to a diamond
    //     (probably; if the 3 or more faces are unrefined (only use the
    //      anchor points))
    // 7 : store. Used to check later on whether there are points with
    //     3 or more of these cells. (note that on a flat surface a boundary
    //     point will only have 4 cells connected to it)


    // Does cell have exactly 7 of its 8 anchor points on the boundary?
    PackedBoolList hasSevenBoundaryAnchorPoints(mesh_.nCells());
    // If so what is the remaining non-boundary anchor point?
    labelHashSet nonBoundaryAnchors(mesh_.nCells()/10000);

    // On-the-fly addressing storage.
    DynamicList<label> dynFEdges;
    DynamicList<label> dynCPoints;

    forAll(cellLevel, celli)
    {
        const labelList& cPoints = mesh_.cellPoints(celli, dynCPoints);

        // Get number of anchor points (pointLevel <= cellLevel)

        label nBoundaryAnchors = 0;
        label nNonAnchorBoundary = 0;
        label nonBoundaryAnchor = -1;

        forAll(cPoints, i)
        {
            label pointi = cPoints[i];

            if (pointLevel[pointi] <= cellLevel[celli])
            {
                // Anchor point
                if (isBoundaryPoint[pointi])
                {
                    nBoundaryAnchors++;
                }
                else
                {
                    // Anchor point which is not on the surface
                    nonBoundaryAnchor = pointi;
                }
            }
            else if (isBoundaryPoint[pointi])
            {
                nNonAnchorBoundary++;
            }
        }

        if (nBoundaryAnchors == 8)
        {
            const cell& cFaces = mesh_.cells()[celli];

            // Count boundary faces.
            label nBfaces = 0;

            forAll(cFaces, cFacei)
            {
                if (isBoundaryFace[cFaces[cFacei]])
                {
                    nBfaces++;
                }
            }

            // If nBfaces > 1 make all non-boundary non-baffle faces baffles.
            // We assume that this situation is where there is a single
            // cell sticking out which would get flattened.

            // Eugene: delete cell no matter what.
            // if (nBfaces > 1)
            {
                if
                (
                    checkCollapse
                && !isCollapsedCell(newPoints, volFraction, celli)
                )
                {
                    nPrevented++;
                    // Pout<< "Preventing baffling/removal of 8 anchor point"
                    //    << " cell "
                    //    << celli << " at " << mesh_.cellCentres()[celli]
                    //    << " since new volume "
                    //    << mesh_.cells()[celli].mag(newPoints, mesh_.faces())
                    //    << " old volume " << mesh_.cellVolumes()[celli]
                    //    << endl;
                }
                else
                {
                    // Block all faces of cell
                    forAll(cFaces, cf)
                    {
                        label facei = cFaces[cf];

                        if
                        (
                            facePatch[facei] == -1
                         && mesh_.isInternalFace(facei)
                        )
                        {
                            facePatch[facei] = nearestAdaptPatch[facei];
                            nBaffleFaces++;

                            // Mark face as a 'boundary'
                            markBoundaryFace
                            (
                                facei,
                                isBoundaryFace,
                                isBoundaryEdge,
                                isBoundaryPoint
                            );
                        }
                    }
                }
            }
        }
        else if (nBoundaryAnchors == 7)
        {
            // Mark the cell. Store the (single!) non-boundary anchor point.
            hasSevenBoundaryAnchorPoints.set(celli, 1u);
            nonBoundaryAnchors.insert(nonBoundaryAnchor);
        }
    }


    // Loop over all points. If a point is connected to 4 or more cells
    // with 7 anchor points on the boundary set those cell's non-boundary faces
    // to baffles

    DynamicList<label> dynPCells;

    forAllConstIter(labelHashSet, nonBoundaryAnchors, iter)
    {
        label pointi = iter.key();

        const labelList& pCells = mesh_.pointCells(pointi, dynPCells);

        // Count number of 'hasSevenBoundaryAnchorPoints' cells.
        label n = 0;

        forAll(pCells, i)
        {
            if (hasSevenBoundaryAnchorPoints.get(pCells[i]) == 1u)
            {
                n++;
            }
        }

        if (n > 3)
        {
            // Point in danger of being what? Remove all 7-cells.
            forAll(pCells, i)
            {
                label celli = pCells[i];

                if (hasSevenBoundaryAnchorPoints.get(celli) == 1u)
                {
                    if
                    (
                        checkCollapse
                    && !isCollapsedCell(newPoints, volFraction, celli)
                    )
                    {
                        nPrevented++;
                        // Pout<< "Preventing baffling of 7 anchor cell "
                        //    << celli
                        //    << " at " << mesh_.cellCentres()[celli]
                        //    << " since new volume "
                        //    << mesh_.cells()[celli].mag
                        //        (newPoints, mesh_.faces())
                        //    << " old volume " << mesh_.cellVolumes()[celli]
                        //    << endl;
                    }
                    else
                    {
                        const cell& cFaces = mesh_.cells()[celli];

                        forAll(cFaces, cf)
                        {
                            label facei = cFaces[cf];

                            if
                            (
                                facePatch[facei] == -1
                             && mesh_.isInternalFace(facei)
                            )
                            {
                                facePatch[facei] = nearestAdaptPatch[facei];
                                nBaffleFaces++;

                                // Mark face as a 'boundary'
                                markBoundaryFace
                                (
                                    facei,
                                    isBoundaryFace,
                                    isBoundaryEdge,
                                    isBoundaryPoint
                                );
                            }
                        }
                    }
                }
            }
        }
    }


    // Sync all. (note that pointdata and facedata not used anymore but sync
    // anyway)

    syncTools::syncPointList
    (
        mesh_,
        isBoundaryPoint,
        orEqOp<bool>(),
        false               // null value
    );

    syncTools::syncEdgeList
    (
        mesh_,
        isBoundaryEdge,
        orEqOp<bool>(),
        false               // null value
    );

    syncTools::syncFaceList
    (
        mesh_,
        isBoundaryFace,
        orEqOp<bool>()
    );


    // Find faces with all edges on the boundary and make them baffles
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        if (facePatch[facei] == -1)
        {
            const labelList& fEdges = mesh_.faceEdges(facei, dynFEdges);
            label nFaceBoundaryEdges = 0;

            forAll(fEdges, fe)
            {
                if (isBoundaryEdge[fEdges[fe]])
                {
                    nFaceBoundaryEdges++;
                }
            }

            if (nFaceBoundaryEdges == fEdges.size())
            {
                if
                (
                    checkCollapse
                && !isCollapsedFace
                    (
                        newPoints,
                        neiCc,
                        minArea,
                        maxNonOrtho,
                        facei
                    )
                )
                {
                    nPrevented++;
                    // Pout<< "Preventing baffling (to avoid collapse) of face "
                    //    << facei
                    //    << " with all boundary edges "
                    //    << " at " << mesh_.faceCentres()[facei]
                    //    << endl;
                }
                else
                {
                    facePatch[facei] = nearestAdaptPatch[facei];
                    nBaffleFaces++;

                    // Do NOT update boundary data since this would grow blocked
                    // faces across gaps.
                }
            }
        }
    }

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled())
        {
            label facei = pp.start();

            forAll(pp, i)
            {
                if (facePatch[facei] == -1)
                {
                    const labelList& fEdges = mesh_.faceEdges(facei, dynFEdges);
                    label nFaceBoundaryEdges = 0;

                    forAll(fEdges, fe)
                    {
                        if (isBoundaryEdge[fEdges[fe]])
                        {
                            nFaceBoundaryEdges++;
                        }
                    }

                    if (nFaceBoundaryEdges == fEdges.size())
                    {
                        if
                        (
                            checkCollapse
                        && !isCollapsedFace
                            (
                                newPoints,
                                neiCc,
                                minArea,
                                maxNonOrtho,
                                facei
                            )
                        )
                        {
                            nPrevented++;
                            // Pout<< "Preventing baffling of coupled face "
                            //    << facei
                            //    << " with all boundary edges "
                            //    << " at " << mesh_.faceCentres()[facei]
                            //    << endl;
                        }
                        else
                        {
                            facePatch[facei] = nearestAdaptPatch[facei];
                            if (isMasterFace[facei])
                            {
                                nBaffleFaces++;
                            }

                            // Do NOT update boundary data since this would grow
                            // blocked faces across gaps.
                        }
                    }
                }

                facei++;
            }
        }
    }


    // Because of isCollapsedFace one side can decide not to baffle whereas
    // the other side does so sync. Baffling is preferred over not baffling.
    if (checkCollapse)  // Or always?
    {
        syncTools::syncFaceList
        (
            mesh_,
            facePatch,
            maxEqOp<label>()
        );
    }

    Info<< "markFacesOnProblemCells : marked "
        << returnReduce(nBaffleFaces, sumOp<label>())
        << " additional internal faces to be converted into baffles."
        << endl;

    if (checkCollapse)
    {
        Info<< "markFacesOnProblemCells : prevented "
            << returnReduce(nPrevented, sumOp<label>())
            << " internal faces from getting converted into baffles."
            << endl;
    }

    return facePatch;
}


// Mark faces to be baffled to prevent snapping problems. Does
// test to find nearest surface and checks which faces would get squashed.
Foam::labelList Foam::meshRefinement::markFacesOnProblemCellsGeometric
(
    const snapParameters& snapParams,
    const dictionary& motionDict
) const
{
    pointField oldPoints(mesh_.points());

    // Repeat (most of) snappySnapDriver::doSnap
    {
        labelList adaptPatchIDs(meshedPatches());

        // Construct addressing engine.
        autoPtr<indirectPrimitivePatch> ppPtr
        (
            meshRefinement::makePatch
            (
                mesh_,
                adaptPatchIDs
            )
        );
        indirectPrimitivePatch& pp = ppPtr();

        // Distance to attract to nearest feature on surface
        const scalarField snapDist
        (
            snappySnapDriver::calcSnapDistance(mesh_, snapParams, pp)
        );


        // Construct iterative mesh mover.
        Info<< "Constructing mesh displacer ..." << endl;
        Info<< "Using mesh parameters " << motionDict << nl << endl;

        const pointMesh& pMesh = pointMesh::New(mesh_);

        motionSmoother meshMover
        (
            mesh_,
            pp,
            adaptPatchIDs,
            meshRefinement::makeDisplacementField(pMesh, adaptPatchIDs)(),
            motionDict
        );


        // Check initial mesh
        Info<< "Checking initial mesh ..." << endl;
        labelHashSet wrongFaces(mesh_.nFaces()/100);
        motionSmoother::checkMesh(false, mesh_, motionDict, wrongFaces);
        const label nInitErrors = returnReduce
        (
            wrongFaces.size(),
            sumOp<label>()
        );

        Info<< "Detected " << nInitErrors << " illegal faces"
            << " (concave, zero area or negative cell pyramid volume)"
            << endl;


        Info<< "Checked initial mesh in = "
            << mesh_.time().cpuTimeIncrement() << " s\n" << nl << endl;

        // Pre-smooth patch vertices (so before determining nearest)
        snappySnapDriver::preSmoothPatch
        (
            *this,
            snapParams,
            nInitErrors,
            List<labelPair>(0), // baffles
            meshMover
        );

        pointField nearestPoint;
        vectorField nearestNormal;
        const vectorField disp
        (
            snappySnapDriver::calcNearestSurface
            (
                *this,
                snapDist,   // attraction
                pp,
                nearestPoint,
                nearestNormal
            )
        );

        const labelList& meshPoints = pp.meshPoints();

        pointField newPoints(mesh_.points());
        forAll(meshPoints, i)
        {
            newPoints[meshPoints[i]] += disp[i];
        }

        syncTools::syncPointList
        (
            mesh_,
            newPoints,
            minMagSqrEqOp<point>(),     // combine op
            vector(great, great, great) // null value (note: cannot use vGreat)
        );

        mesh_.movePoints(newPoints);
    }


    // Per face the nearest adaptPatch
    const labelList nearestAdaptPatch(nearestPatch(meshedPatches()));

    // Per face (internal or coupled!) the patch that the
    // baffle should get (or -1).
    labelList facePatch(mesh_.nFaces(), -1);
    // Count of baffled faces
    label nBaffleFaces = 0;

    {
        faceSet wrongFaces(mesh_, "wrongFaces", 100);
        {
            // motionSmoother::checkMesh(false, mesh_, motionDict, wrongFaces);

            // Just check the errors from squashing
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            const labelList allFaces(identity(mesh_.nFaces()));
            label nWrongFaces = 0;

            // const scalar minV(readScalar(motionDict.lookup("minVol", true)));
            // if (minV > -great)
            //{
            //    polyMeshGeometry::checkFacePyramids
            //    (
            //        false,
            //        minV,
            //        mesh_,
            //        mesh_.cellCentres(),
            //        mesh_.points(),
            //        allFaces,
            //        List<labelPair>(0),
            //        &wrongFaces
            //    );
            //
            //    label nNewWrongFaces = returnReduce
            //    (
            //        wrongFaces.size(),
            //        sumOp<label>()
            //    );
            //
            //    Info<< "    faces with pyramid volume < "
            //        << setw(5) << minV
            //        << " m^3                  : "
            //        << nNewWrongFaces-nWrongFaces << endl;
            //
            //    nWrongFaces = nNewWrongFaces;
            //}

            scalar minArea(readScalar(motionDict.lookup("minArea")));
            if (minArea > -small)
            {
                polyMeshGeometry::checkFaceArea
                (
                    false,
                    minArea,
                    mesh_,
                    mesh_.faceAreas(),
                    allFaces,
                    &wrongFaces
                );

                label nNewWrongFaces = returnReduce
                (
                    wrongFaces.size(),
                    sumOp<label>()
                );

                Info<< "    faces with area < "
                    << setw(5) << minArea
                    << " m^2                            : "
                    << nNewWrongFaces-nWrongFaces << endl;

                nWrongFaces = nNewWrongFaces;
            }

            scalar minDet(readScalar(motionDict.lookup("minDeterminant")));
            if (minDet > -1)
            {
                polyMeshGeometry::checkCellDeterminant
                (
                    false,
                    minDet,
                    mesh_,
                    mesh_.faceAreas(),
                    allFaces,
                    polyMeshGeometry::affectedCells(mesh_, allFaces),
                    &wrongFaces
                );

                label nNewWrongFaces = returnReduce
                (
                    wrongFaces.size(),
                    sumOp<label>()
                );

                Info<< "    faces on cells with determinant < "
                    << setw(5) << minDet << "                : "
                    << nNewWrongFaces-nWrongFaces << endl;

                nWrongFaces = nNewWrongFaces;
            }
        }


        forAllConstIter(faceSet, wrongFaces, iter)
        {
            label patchi = mesh_.boundaryMesh().whichPatch(iter.key());

            if (patchi == -1 || mesh_.boundaryMesh()[patchi].coupled())
            {
                facePatch[iter.key()] = nearestAdaptPatch[iter.key()];
                nBaffleFaces++;

                // Pout<< "    " << iter.key()
                //    //<< " on patch " << mesh_.boundaryMesh()[patchi].name()
                //    << " is destined for patch " << facePatch[iter.key()]
                //    << endl;
            }
        }
    }


    // Restore points.
    mesh_.movePoints(oldPoints);


    Info<< "markFacesOnProblemCellsGeometric : marked "
        << returnReduce(nBaffleFaces, sumOp<label>())
        << " additional internal and coupled faces"
        << " to be converted into baffles." << endl;

    syncTools::syncFaceList
    (
        mesh_,
        facePatch,
        maxEqOp<label>()
    );

    return facePatch;
}


// ************************************************************************* //
