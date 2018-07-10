/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
        FatalErrorInFunction
            << "nBoundaries:" << nBoundaryFaces
            << " neiCc:" << neiCc.size()
            << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        const labelUList& faceCells = pp.faceCells();

        label bFacei = pp.start() - mesh.nInternalFaces();

        if (pp.coupled())
        {
            forAll(faceCells, i)
            {
                neiCc[bFacei] = cellCentres[faceCells[i]];
                bFacei++;
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

    forAll(patches, patchi)
    {
        // Check all coupled. Avoid using .coupled() so we also pick up AMI.
        if (isA<coupledPolyPatch>(patches[patchi]))
        {
            const coupledPolyPatch& cpp = refCast<const coupledPolyPatch>
            (
                patches[patchi]
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

    forAll(faceToSurface, facei)
    {
        if (faceToSurface[facei] == -1)
        {
            blockedFace[facei] = false;
        }
        else
        {
            blockedFace[facei] = true;
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

        label celli = mesh.findCell(insidePoint);

        if (celli != -1)
        {
            keepRegionI = cellRegion[celli];
        }
        reduce(keepRegionI, maxOp<label>());

        Info<< "    For surface " << surfName
            << " found point " << insidePoint << " in cell " << celli
            << " in global region " << keepRegionI
            << " out of " << cellRegion.nRegions() << " regions." << endl;

        if (keepRegionI == -1)
        {
            FatalErrorInFunction
                << "Point " << insidePoint
                << " is not inside the mesh." << nl
                << "Bounding box of the mesh:" << mesh.bounds()
                << exit(FatalError);
        }

        // Set all cells with this region
        forAll(cellRegion, celli)
        {
            if (cellRegion[celli] == keepRegionI)
            {
                if (cellToSurface[celli] == -2)
                {
                    cellToSurface[celli] = surfI;
                }
                else if (cellToSurface[celli] != surfI)
                {
                    WarningInFunction
                        << "Cell " << celli
                        << " at " << mesh.cellCentres()[celli]
                        << " is inside surface " << surfName
                        << " but already marked as being in zone "
                        << cellToSurface[celli] << endl
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
    labelList cellToSurface(cellCentres.size(), label(-1));

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
            FatalErrorInFunction
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

            forAll(volType, pointi)
            {
                if (cellToSurface[pointi] == -1)
                {
                    if
                    (
                        (
                            volType[pointi] == volumeType::INSIDE
                         && selectInside
                        )
                     || (
                            volType[pointi] == volumeType::OUTSIDE
                         && !selectInside
                        )
                    )
                    {
                        cellToSurface[pointi] = surfI;
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

    labelList neiFaceOwner(mesh.nFaces() - mesh.nInternalFaces(), label(-1));

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        const labelUList& faceCells = pp.faceCells();

        label bFacei = pp.start() - mesh.nInternalFaces();

        if (pp.coupled())
        {
            forAll(faceCells, i)
            {
                neiFaceOwner[bFacei] = cellToSurface[faceCells[i]];
                bFacei++;
            }
        }
    }

    syncTools::swapBoundaryFaceList(mesh, neiFaceOwner);

    forAll(faces, facei)
    {
        const label ownerSurfacei = cellToSurface[faceOwner[facei]];

        if (faceToSurface[facei] >= 0)
        {
            continue;
        }

        if (mesh.isInternalFace(facei))
        {
            const label neiSurfacei = cellToSurface[faceNeighbour[facei]];

            if
            (
                (ownerSurfacei >= 0 || neiSurfacei >= 0)
             && ownerSurfacei != neiSurfacei
            )
            {
                flipMap[facei] =
                    (
                        ownerSurfacei == max(ownerSurfacei, neiSurfacei)
                      ? false
                      : true
                    );

                faceToSurface[facei] = max(ownerSurfacei, neiSurfacei);
            }
        }
        else
        {
            label patchID = mesh.boundaryMesh().whichPatch(facei);

            if (mesh.boundaryMesh()[patchID].coupled())
            {
                const label neiSurfacei =
                    neiFaceOwner[facei - mesh.nInternalFaces()];

                if
                (
                    (ownerSurfacei >= 0 || neiSurfacei >= 0)
                 && ownerSurfacei != neiSurfacei
                )
                {
                    flipMap[facei] =
                        (
                            ownerSurfacei == max(ownerSurfacei, neiSurfacei)
                          ? false
                          : true
                        );

                    faceToSurface[facei] = max(ownerSurfacei, neiSurfacei);
                }
            }
            else
            {
                if (ownerSurfacei >= 0)
                {
                    faceToSurface[facei] = ownerSurfacei;
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
    forAll(faces, facei)
    {
        if (faceToSurface[facei] >= 0)
        {
            continue;
        }

        label patchID = mesh.boundaryMesh().whichPatch(facei);

        const label own = faceOwner[facei];

        List<pointIndexHit> surfHit;
        labelList hitSurface;

        if (mesh.isInternalFace(facei))
        {
            const label nei = faceNeighbour[facei];

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
                neiCc[facei - mesh.nInternalFaces()],
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

                vector fN = faces[facei].area(mesh.points());
                fN /= mag(fN) + small;

                if ((norm[0] & fN) < 0)
                {
                    flipMap[facei] = true;
                }
                else
                {
                    flipMap[facei] = false;
                }

                faceToSurface[facei] = hitSurface[0];
            }
        }
    }


//    labelList neiCellSurface(mesh.nFaces()-mesh.nInternalFaces());
//
//    forAll(patches, patchi)
//    {
//        const polyPatch& pp = patches[patchi];
//
//        if (pp.coupled())
//        {
//            forAll(pp, i)
//            {
//                label facei = pp.start()+i;
//                label ownSurface = cellToSurface[faceOwner[facei]];
//                neiCellSurface[facei - mesh.nInternalFaces()] = ownSurface;
//            }
//        }
//    }
//    syncTools::swapBoundaryFaceList(mesh, neiCellSurface);
//
//    forAll(patches, patchi)
//    {
//        const polyPatch& pp = patches[patchi];
//
//        if (pp.coupled())
//        {
//            forAll(pp, i)
//            {
//                label facei = pp.start()+i;
//                label ownSurface = cellToSurface[faceOwner[facei]];
//                label neiSurface =
//                    neiCellSurface[facei-mesh.nInternalFaces()];
//
//                if (faceToSurface[facei] == -1 && (ownSurface != neiSurface))
//                {
//                    // Give face the max cell zone
//                    faceToSurface[facei] =  max(ownSurface, neiSurface);
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

    forAll(cellToSurface, celli)
    {
        label surfacei = cellToSurface[celli];

        if (surfacei >= 0)
        {
            label zoneI = surfaceToCellZone[surfacei];

            if (zoneI >= 0)
            {
                meshMod.setAction
                (
                    polyModifyCell
                    (
                        celli,
                        false,          // removeFromZone
                        zoneI
                    )
                );
            }
        }
    }

    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeighbour = mesh.faceNeighbour();

    forAll(faceToSurface, facei)
    {
        label surfacei = faceToSurface[facei];

        if (surfacei < 0)
        {
            continue;
        }

        label patchID = mesh.boundaryMesh().whichPatch(facei);

        if (mesh.isInternalFace(facei))
        {
            label own = faceOwner[facei];
            label nei = faceNeighbour[facei];

            meshMod.setAction
            (
                polyModifyFace
                (
                    mesh.faces()[facei],            // modified face
                    facei,                          // label of face
                    own,                            // owner
                    nei,                            // neighbour
                    false,                          // face flip
                    -1,                             // patch for face
                    false,                          // remove from zone
                    surfaceToFaceZone[surfacei],    // zone for face
                    flipMap[facei]                  // face flip in zone
                )
            );
        }
        else if (patchID != -1 && mesh.boundaryMesh()[patchID].coupled())
        {
            label own = faceOwner[facei];

            meshMod.setAction
            (
                polyModifyFace
                (
                    mesh.faces()[facei],            // modified face
                    facei,                          // label of face
                    own,                            // owner
                    -1,                             // neighbour
                    false,                          // face flip
                    patchID,                        // patch for face
                    false,                          // remove from zone
                    surfaceToFaceZone[surfacei],    // zone for face
                    flipMap[facei]                  // face flip in zone
                )
            );
        }
    }

    // Change the mesh (no inflation, parallel sync)
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false, true);
}


// ************************************************************************* //
