/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "conformalVoronoiMesh.H"
#include "polyModifyFace.H"
#include "polyModifyCell.H"
#include "syncTools.H"
#include "regionSplit.H"
#include "surfaceZonesInfo.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::conformalVoronoiMesh::calcNeighbourCellCentres
(
    const polyMesh& mesh,
    const pointField& cellCentres,
    pointField& neiCc
) const
{
    label nBoundaryFaces = mesh.nFaces() - mesh.nInternalFaces();

    if (neiCc.size() != nBoundaryFaces)
    {
        FatalErrorIn("conformalVoronoiMesh::calcNeighbourCellCentres(..)")
            << "nBoundaries:" << nBoundaryFaces
            << " neiCc:" << neiCc.size()
            << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        const labelUList& faceCells = pp.faceCells();

        label bFaceI = pp.start() - mesh.nInternalFaces();

        if (pp.coupled())
        {
            forAll(faceCells, i)
            {
                neiCc[bFaceI] = cellCentres[faceCells[i]];
                bFaceI++;
            }
        }
    }

    // Swap coupled boundaries. Apply separation to cc since is coordinate.
    syncTools::swapBoundaryFacePositions(mesh, neiCc);
}


void Foam::conformalVoronoiMesh::selectSeparatedCoupledFaces
(
    const polyMesh& mesh,
    boolList& selected
) const
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchI)
    {
        // Check all coupled. Avoid using .coupled() so we also pick up AMI.
        if (isA<coupledPolyPatch>(patches[patchI]))
        {
            const coupledPolyPatch& cpp = refCast<const coupledPolyPatch>
            (
                patches[patchI]
            );

            if (cpp.separated() || !cpp.parallel())
            {
                forAll(cpp, i)
                {
                    selected[cpp.start()+i] = true;
                }
            }
        }
    }
}


void Foam::conformalVoronoiMesh::findCellZoneInsideWalk
(
    const polyMesh& mesh,
    const labelList& locationSurfaces,  // indices of surfaces with inside point
    const labelList& faceToSurface, // per face index of named surface
    labelList& cellToSurface
) const
{
    // Analyse regions. Reuse regionsplit
    boolList blockedFace(mesh.nFaces());
    selectSeparatedCoupledFaces(mesh, blockedFace);

    forAll(faceToSurface, faceI)
    {
        if (faceToSurface[faceI] == -1)
        {
            blockedFace[faceI] = false;
        }
        else
        {
            blockedFace[faceI] = true;
        }
    }
    // No need to sync since namedSurfaceIndex already is synced

    // Set region per cell based on walking
    regionSplit cellRegion(mesh, blockedFace);
    blockedFace.clear();


    // Force calculation of face decomposition (used in findCell)
    (void)mesh.tetBasePtIs();

    const PtrList<surfaceZonesInfo>& surfZones =
        geometryToConformTo().surfZones();

    // For all locationSurface find the cell
    forAll(locationSurfaces, i)
    {
        label surfI = locationSurfaces[i];

        const Foam::point& insidePoint = surfZones[surfI].zoneInsidePoint();

        const word& surfName = geometryToConformTo().geometry().names()[surfI];

        Info<< "    For surface " << surfName
            << " finding inside point " << insidePoint
            << endl;

        // Find the region containing the insidePoint
        label keepRegionI = -1;

        label cellI = mesh.findCell(insidePoint);

        if (cellI != -1)
        {
            keepRegionI = cellRegion[cellI];
        }
        reduce(keepRegionI, maxOp<label>());

        Info<< "    For surface " << surfName
            << " found point " << insidePoint << " in cell " << cellI
            << " in global region " << keepRegionI
            << " out of " << cellRegion.nRegions() << " regions." << endl;

        if (keepRegionI == -1)
        {
            FatalErrorIn
            (
                "conformalVoronoiMesh::findCellZoneInsideWalk"
                "(const polyMesh&, const labelList&"
                ", const labelList&, labelList&)"
            )   << "Point " << insidePoint
                << " is not inside the mesh." << nl
                << "Bounding box of the mesh:" << mesh.bounds()
                << exit(FatalError);
        }

        // Set all cells with this region
        forAll(cellRegion, cellI)
        {
            if (cellRegion[cellI] == keepRegionI)
            {
                if (cellToSurface[cellI] == -2)
                {
                    cellToSurface[cellI] = surfI;
                }
                else if (cellToSurface[cellI] != surfI)
                {
                    WarningIn
                    (
                        "conformalVoronoiMesh::findCellZoneInsideWalk"
                        "(const labelList&, const labelList&"
                        ", const labelList&, const labelList&)"
                    )   << "Cell " << cellI
                        << " at " << mesh.cellCentres()[cellI]
                        << " is inside surface " << surfName
                        << " but already marked as being in zone "
                        << cellToSurface[cellI] << endl
                        << "This can happen if your surfaces are not"
                        << " (sufficiently) closed."
                        << endl;
                }
            }
        }
    }
}


Foam::labelList Foam::conformalVoronoiMesh::calcCellZones
(
    const pointField& cellCentres
) const
{
    labelList cellToSurface(cellCentres.size(), -1);

    const PtrList<surfaceZonesInfo>& surfZones =
        geometryToConformTo().surfZones();

    // Get list of closed surfaces
    labelList closedNamedSurfaces
    (
        surfaceZonesInfo::getAllClosedNamedSurfaces
        (
            surfZones,
            geometryToConformTo().geometry(),
            geometryToConformTo().surfaces()
        )
    );

    forAll(closedNamedSurfaces, i)
    {
        label surfI = closedNamedSurfaces[i];

        const searchableSurface& surface =
            allGeometry()[geometryToConformTo().surfaces()[surfI]];

        const surfaceZonesInfo::areaSelectionAlgo selectionMethod =
            surfZones[surfI].zoneInside();

        if
        (
            selectionMethod != surfaceZonesInfo::INSIDE
         && selectionMethod != surfaceZonesInfo::OUTSIDE
         && selectionMethod != surfaceZonesInfo::INSIDEPOINT
        )
        {
            FatalErrorIn("conformalVoronoiMesh::calcCellZones(..)")
                << "Trying to use surface "
                << surface.name()
                << " which has non-geometric inside selection method "
                << surfaceZonesInfo::areaSelectionAlgoNames[selectionMethod]
                << exit(FatalError);
        }

        if (surface.hasVolumeType())
        {
            List<volumeType> volType;
            surface.getVolumeType(cellCentres, volType);

            bool selectInside = true;
            if (selectionMethod == surfaceZonesInfo::INSIDEPOINT)
            {
                List<volumeType> volTypeInsidePoint;
                surface.getVolumeType
                (
                    pointField(1, surfZones[surfI].zoneInsidePoint()),
                    volTypeInsidePoint
                );

                if (volTypeInsidePoint[0] == volumeType::OUTSIDE)
                {
                    selectInside = false;
                }
            }
            else if (selectionMethod == surfaceZonesInfo::OUTSIDE)
            {
                selectInside = false;
            }

            forAll(volType, pointI)
            {
                if (cellToSurface[pointI] == -1)
                {
                    if
                    (
                        (
                            volType[pointI] == volumeType::INSIDE
                         && selectInside
                        )
                     || (
                            volType[pointI] == volumeType::OUTSIDE
                         && !selectInside
                        )
                    )
                    {
                        cellToSurface[pointI] = surfI;
                    }
                }
            }
        }
    }

    return cellToSurface;
}


void Foam::conformalVoronoiMesh::calcFaceZones
(
    const polyMesh& mesh,
    const pointField& cellCentres,
    const labelList& cellToSurface,
    labelList& faceToSurface,
    boolList& flipMap
) const
{
    faceToSurface.setSize(mesh.nFaces(), -1);
    flipMap.setSize(mesh.nFaces(), false);

    const faceList& faces = mesh.faces();
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeighbour = mesh.faceNeighbour();

    labelList neiFaceOwner(mesh.nFaces() - mesh.nInternalFaces(), -1);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        const labelUList& faceCells = pp.faceCells();

        label bFaceI = pp.start() - mesh.nInternalFaces();

        if (pp.coupled())
        {
            forAll(faceCells, i)
            {
                neiFaceOwner[bFaceI] = cellToSurface[faceCells[i]];
                bFaceI++;
            }
        }
    }

    syncTools::swapBoundaryFaceList(mesh, neiFaceOwner);

    forAll(faces, faceI)
    {
        const label ownerSurfaceI = cellToSurface[faceOwner[faceI]];

        if (faceToSurface[faceI] >= 0)
        {
            continue;
        }

        if (mesh.isInternalFace(faceI))
        {
            const label neiSurfaceI = cellToSurface[faceNeighbour[faceI]];

            if
            (
                (ownerSurfaceI >= 0 || neiSurfaceI >= 0)
             && ownerSurfaceI != neiSurfaceI
            )
            {
                flipMap[faceI] =
                    (
                        ownerSurfaceI == max(ownerSurfaceI, neiSurfaceI)
                      ? false
                      : true
                    );

                faceToSurface[faceI] = max(ownerSurfaceI, neiSurfaceI);
            }
        }
        else
        {
            label patchID = mesh.boundaryMesh().whichPatch(faceI);

            if (mesh.boundaryMesh()[patchID].coupled())
            {
                const label neiSurfaceI =
                    neiFaceOwner[faceI - mesh.nInternalFaces()];

                if
                (
                    (ownerSurfaceI >= 0 || neiSurfaceI >= 0)
                 && ownerSurfaceI != neiSurfaceI
                )
                {
                    flipMap[faceI] =
                        (
                            ownerSurfaceI == max(ownerSurfaceI, neiSurfaceI)
                          ? false
                          : true
                        );

                    faceToSurface[faceI] = max(ownerSurfaceI, neiSurfaceI);
                }
            }
            else
            {
                if (ownerSurfaceI >= 0)
                {
                    faceToSurface[faceI] = ownerSurfaceI;
                }
            }
        }
    }


    const PtrList<surfaceZonesInfo>& surfZones =
        geometryToConformTo().surfZones();

    labelList unclosedSurfaces
    (
        surfaceZonesInfo::getUnclosedNamedSurfaces
        (
            surfZones,
            geometryToConformTo().geometry(),
            geometryToConformTo().surfaces()
        )
    );

    pointField neiCc(mesh.nFaces() - mesh.nInternalFaces());
    calcNeighbourCellCentres
    (
        mesh,
        cellCentres,
        neiCc
    );

    // Use intersection of cellCentre connections
    forAll(faces, faceI)
    {
        if (faceToSurface[faceI] >= 0)
        {
            continue;
        }

        label patchID = mesh.boundaryMesh().whichPatch(faceI);

        const label own = faceOwner[faceI];

        List<pointIndexHit> surfHit;
        labelList hitSurface;

        if (mesh.isInternalFace(faceI))
        {
            const label nei = faceNeighbour[faceI];

            geometryToConformTo().findSurfaceAllIntersections
            (
                cellCentres[own],
                cellCentres[nei],
                surfHit,
                hitSurface
            );
        }
        else if (patchID != -1 && mesh.boundaryMesh()[patchID].coupled())
        {
            geometryToConformTo().findSurfaceAllIntersections
            (
                cellCentres[own],
                neiCc[faceI - mesh.nInternalFaces()],
                surfHit,
                hitSurface
            );
        }

        // If there are multiple intersections then do not add to
        // a faceZone
        if (surfHit.size() == 1 && surfHit[0].hit())
        {
            if (findIndex(unclosedSurfaces, hitSurface[0]) != -1)
            {
                vectorField norm;
                geometryToConformTo().getNormal
                (
                    hitSurface[0],
                    List<pointIndexHit>(1, surfHit[0]),
                    norm
                );

                vector fN = faces[faceI].normal(mesh.points());
                fN /= mag(fN) + SMALL;

                if ((norm[0] & fN) < 0)
                {
                    flipMap[faceI] = true;
                }
                else
                {
                    flipMap[faceI] = false;
                }

                faceToSurface[faceI] = hitSurface[0];
            }
        }
    }


//    labelList neiCellSurface(mesh.nFaces()-mesh.nInternalFaces());
//
//    forAll(patches, patchI)
//    {
//        const polyPatch& pp = patches[patchI];
//
//        if (pp.coupled())
//        {
//            forAll(pp, i)
//            {
//                label faceI = pp.start()+i;
//                label ownSurface = cellToSurface[faceOwner[faceI]];
//                neiCellSurface[faceI - mesh.nInternalFaces()] = ownSurface;
//            }
//        }
//    }
//    syncTools::swapBoundaryFaceList(mesh, neiCellSurface);
//
//    forAll(patches, patchI)
//    {
//        const polyPatch& pp = patches[patchI];
//
//        if (pp.coupled())
//        {
//            forAll(pp, i)
//            {
//                label faceI = pp.start()+i;
//                label ownSurface = cellToSurface[faceOwner[faceI]];
//                label neiSurface =
//                    neiCellSurface[faceI-mesh.nInternalFaces()];
//
//                if (faceToSurface[faceI] == -1 && (ownSurface != neiSurface))
//                {
//                    // Give face the max cell zone
//                    faceToSurface[faceI] =  max(ownSurface, neiSurface);
//                }
//            }
//        }
//    }

    // Sync
    syncTools::syncFaceList(mesh, faceToSurface, maxEqOp<label>());
}


void Foam::conformalVoronoiMesh::addZones
(
    polyMesh& mesh,
    const pointField& cellCentres
) const
{
    Info<< "    Adding zones to mesh" << endl;

    const PtrList<surfaceZonesInfo>& surfZones =
        geometryToConformTo().surfZones();

    labelList cellToSurface(calcCellZones(cellCentres));

    labelList faceToSurface;
    boolList flipMap;

    calcFaceZones
    (
        mesh,
        cellCentres,
        cellToSurface,
        faceToSurface,
        flipMap
    );

    labelList insidePointNamedSurfaces
    (
        surfaceZonesInfo::getInsidePointNamedSurfaces(surfZones)
    );

    findCellZoneInsideWalk
    (
        mesh,
        insidePointNamedSurfaces,
        faceToSurface,
        cellToSurface
    );

    labelList namedSurfaces(surfaceZonesInfo::getNamedSurfaces(surfZones));

    forAll(namedSurfaces, i)
    {
        label surfI = namedSurfaces[i];

        Info<< incrIndent << indent << "Surface : "
            << geometryToConformTo().geometry().names()[surfI] << nl
            << indent << "    faceZone : "
            << surfZones[surfI].faceZoneName() << nl
            << indent << "    cellZone : "
            << surfZones[surfI].cellZoneName()
            << decrIndent << endl;
    }

    // Add zones to mesh
    labelList surfaceToFaceZone =
        surfaceZonesInfo::addFaceZonesToMesh
        (
            surfZones,
            namedSurfaces,
            mesh
        );

    labelList surfaceToCellZone =
        surfaceZonesInfo::addCellZonesToMesh
        (
            surfZones,
            namedSurfaces,
            mesh
        );

    // Topochange container
    polyTopoChange meshMod(mesh);

    forAll(cellToSurface, cellI)
    {
        label surfaceI = cellToSurface[cellI];

        if (surfaceI >= 0)
        {
            label zoneI = surfaceToCellZone[surfaceI];

            if (zoneI >= 0)
            {
                meshMod.setAction
                (
                    polyModifyCell
                    (
                        cellI,
                        false,          // removeFromZone
                        zoneI
                    )
                );
            }
        }
    }

    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeighbour = mesh.faceNeighbour();

    forAll(faceToSurface, faceI)
    {
        label surfaceI = faceToSurface[faceI];

        if (surfaceI < 0)
        {
            continue;
        }

        label patchID = mesh.boundaryMesh().whichPatch(faceI);

        if (mesh.isInternalFace(faceI))
        {
            label own = faceOwner[faceI];
            label nei = faceNeighbour[faceI];

            meshMod.setAction
            (
                polyModifyFace
                (
                    mesh.faces()[faceI],            // modified face
                    faceI,                          // label of face
                    own,                            // owner
                    nei,                            // neighbour
                    false,                          // face flip
                    -1,                             // patch for face
                    false,                          // remove from zone
                    surfaceToFaceZone[surfaceI],    // zone for face
                    flipMap[faceI]                  // face flip in zone
                )
            );
        }
        else if (patchID != -1 && mesh.boundaryMesh()[patchID].coupled())
        {
            label own = faceOwner[faceI];

            meshMod.setAction
            (
                polyModifyFace
                (
                    mesh.faces()[faceI],            // modified face
                    faceI,                          // label of face
                    own,                            // owner
                    -1,                             // neighbour
                    false,                          // face flip
                    patchID,                        // patch for face
                    false,                          // remove from zone
                    surfaceToFaceZone[surfaceI],    // zone for face
                    flipMap[faceI]                  // face flip in zone
                )
            );
        }
    }

    // Change the mesh (no inflation, parallel sync)
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false, true);
}


// ************************************************************************* //
