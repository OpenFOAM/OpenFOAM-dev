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
#include "polyTopoChange.H"
#include "meshTools.H"
#include "polyModifyFace.H"
#include "polyModifyCell.H"
#include "polyAddFace.H"
#include "polyRemoveFace.H"
#include "polyAddPoint.H"
#include "localPointRegion.H"
#include "duplicatePoints.H"
#include "OFstream.H"
#include "regionSplit.H"
#include "removeCells.H"
#include "unitConversion.H"
#include "OBJstream.H"
#include "patchFaceOrientation.H"
#include "PatchEdgeFaceWave.H"
#include "patchEdgeFaceRegion.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::meshRefinement::createBaffle
(
    const label facei,
    const label ownPatch,
    const label neiPatch,
    polyTopoChange& meshMod
) const
{
    const face& f = mesh_.faces()[facei];
    label zoneID = mesh_.faceZones().whichZone(facei);
    bool zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh_.faceZones()[zoneID];
        zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
    }

    meshMod.setAction
    (
        polyModifyFace
        (
            f,                          // modified face
            facei,                      // label of face
            mesh_.faceOwner()[facei],   // owner
            -1,                         // neighbour
            false,                      // face flip
            ownPatch,                   // patch for face
            false,                      // remove from zone
            zoneID,                     // zone for face
            zoneFlip                    // face flip in zone
        )
    );


    label dupFacei = -1;

    if (mesh_.isInternalFace(facei))
    {
        if (neiPatch == -1)
        {
            FatalErrorInFunction
                << "No neighbour patch for internal face " << facei
                << " fc:" << mesh_.faceCentres()[facei]
                << " ownPatch:" << ownPatch << abort(FatalError);
        }

        bool reverseFlip = false;
        if (zoneID >= 0)
        {
            reverseFlip = !zoneFlip;
        }

        dupFacei = meshMod.setAction
        (
            polyAddFace
            (
                f.reverseFace(),            // modified face
                mesh_.faceNeighbour()[facei],// owner
                -1,                         // neighbour
                -1,                         // masterPointID
                -1,                         // masterEdgeID
                facei,                      // masterFaceID,
                true,                       // face flip
                neiPatch,                   // patch for face
                zoneID,                     // zone for face
                reverseFlip                 // face flip in zone
            )
        );
    }
    return dupFacei;
}


void Foam::meshRefinement::getBafflePatches
(
    const labelList& globalToMasterPatch,
    const labelList& neiLevel,
    const pointField& neiCc,
    labelList& ownPatch,
    labelList& neiPatch
) const
{
    autoPtr<OFstream> str;
    label vertI = 0;
    if (debug&OBJINTERSECTIONS)
    {
        mkDir(mesh_.time().path()/timeName());
        str.reset
        (
            new OFstream
            (
                mesh_.time().path()/timeName()/"intersections.obj"
            )
        );

        Pout<< "getBafflePatches : Writing surface intersections to file "
            << str().name() << nl << endl;
    }

    const pointField& cellCentres = mesh_.cellCentres();

    // Surfaces that need to be baffled
    const labelList surfacesToBaffle
    (
        surfaceZonesInfo::getUnnamedSurfaces(surfaces_.surfZones())
    );

    ownPatch.setSize(mesh_.nFaces());
    ownPatch = -1;
    neiPatch.setSize(mesh_.nFaces());
    neiPatch = -1;


    // Collect candidate faces
    // ~~~~~~~~~~~~~~~~~~~~~~~

    labelList testFaces(intersectedFaces());

    // Collect segments
    // ~~~~~~~~~~~~~~~~

    pointField start(testFaces.size());
    pointField end(testFaces.size());

    forAll(testFaces, i)
    {
        label facei = testFaces[i];

        label own = mesh_.faceOwner()[facei];

        if (mesh_.isInternalFace(facei))
        {
            start[i] = cellCentres[own];
            end[i] = cellCentres[mesh_.faceNeighbour()[facei]];
        }
        else
        {
            start[i] = cellCentres[own];
            end[i] = neiCc[facei-mesh_.nInternalFaces()];
        }
    }

    // Extend segments a bit
    {
        const vectorField smallVec(rootSmall*(end-start));
        start -= smallVec;
        end += smallVec;
    }


    // Do test for intersections
    // ~~~~~~~~~~~~~~~~~~~~~~~~~
    labelList surface1;
    List<pointIndexHit> hit1;
    labelList region1;
    labelList surface2;
    List<pointIndexHit> hit2;
    labelList region2;
    surfaces_.findNearestIntersection
    (
        surfacesToBaffle,
        start,
        end,

        surface1,
        hit1,
        region1,

        surface2,
        hit2,
        region2
    );

    forAll(testFaces, i)
    {
        label facei = testFaces[i];

        if (hit1[i].hit() && hit2[i].hit())
        {
            if (str.valid())
            {
                meshTools::writeOBJ(str(), start[i]);
                vertI++;
                meshTools::writeOBJ(str(), hit1[i].rawPoint());
                vertI++;
                meshTools::writeOBJ(str(), hit2[i].rawPoint());
                vertI++;
                meshTools::writeOBJ(str(), end[i]);
                vertI++;
                str()<< "l " << vertI-3 << ' ' << vertI-2 << nl;
                str()<< "l " << vertI-2 << ' ' << vertI-1 << nl;
                str()<< "l " << vertI-1 << ' ' << vertI << nl;
            }

            // Pick up the patches
            ownPatch[facei] = globalToMasterPatch
            [
                surfaces_.globalRegion(surface1[i], region1[i])
            ];
            neiPatch[facei] = globalToMasterPatch
            [
                surfaces_.globalRegion(surface2[i], region2[i])
            ];

            if (ownPatch[facei] == -1 || neiPatch[facei] == -1)
            {
                FatalErrorInFunction
                    << "problem." << abort(FatalError);
            }
        }
    }

    // No need to parallel sync since intersection data (surfaceIndex_ etc.)
    // already guaranteed to be synced...
    // However:
    // - owncc and neicc are reversed on different procs so might pick
    //   up different regions reversed? No problem. Neighbour on one processor
    //   might not be owner on the other processor but the neighbour is
    //   not used when creating baffles from proc faces.
    // - tolerances issues occasionally crop up.
    syncTools::syncFaceList(mesh_, ownPatch, maxEqOp<label>());
    syncTools::syncFaceList(mesh_, neiPatch, maxEqOp<label>());
}


Foam::Map<Foam::labelPair> Foam::meshRefinement::getZoneBafflePatches
(
    const bool allowBoundary,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch
) const
{
    Map<labelPair> bafflePatch(mesh_.nFaces()/1000);

    const PtrList<surfaceZonesInfo>& surfZones = surfaces_.surfZones();
    const faceZoneMesh& fZones = mesh_.faceZones();

    forAll(surfZones, surfI)
    {
        const word& faceZoneName = surfZones[surfI].faceZoneName();

        if (faceZoneName.size())
        {
            // Get zone
            label zoneI = fZones.findZoneID(faceZoneName);

            const faceZone& fZone = fZones[zoneI];

            // Get patch allocated for zone
            label globalRegionI = surfaces_.globalRegion(surfI, 0);
            labelPair zPatches
            (
                globalToMasterPatch[globalRegionI],
                globalToSlavePatch[globalRegionI]
            );

            Info<< "For zone " << fZone.name() << " found patches "
                << mesh_.boundaryMesh()[zPatches[0]].name() << " and "
                << mesh_.boundaryMesh()[zPatches[1]].name()
                << endl;

            forAll(fZone, i)
            {
                label facei = fZone[i];

                if (allowBoundary || mesh_.isInternalFace(facei))
                {
                    labelPair patches = zPatches;
                    if (fZone.flipMap()[i])
                    {
                       patches = reverse(patches);
                    }

                    if (!bafflePatch.insert(facei, patches))
                    {
                        FatalErrorInFunction
                            << "Face " << facei
                            << " fc:" << mesh_.faceCentres()[facei]
                            << " in zone " << fZone.name()
                            << " is in multiple zones!"
                            << abort(FatalError);
                    }
                }
            }
        }
    }
    return bafflePatch;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::createBaffles
(
    const labelList& ownPatch,
    const labelList& neiPatch
)
{
    if
    (
        ownPatch.size() != mesh_.nFaces()
     || neiPatch.size() != mesh_.nFaces()
    )
    {
        FatalErrorInFunction
            << "Illegal size :"
            << " ownPatch:" << ownPatch.size()
            << " neiPatch:" << neiPatch.size()
            << ". Should be number of faces:" << mesh_.nFaces()
            << abort(FatalError);
    }

    if (debug)
    {
        labelList syncedOwnPatch(ownPatch);
        syncTools::syncFaceList(mesh_, syncedOwnPatch, maxEqOp<label>());
        labelList syncedNeiPatch(neiPatch);
        syncTools::syncFaceList(mesh_, syncedNeiPatch, maxEqOp<label>());

        forAll(syncedOwnPatch, facei)
        {
            if
            (
                (ownPatch[facei] == -1 && syncedOwnPatch[facei] != -1)
             || (neiPatch[facei] == -1 && syncedNeiPatch[facei] != -1)
            )
            {
                FatalErrorInFunction
                    << "Non synchronised at face:" << facei
                    << " on patch:" << mesh_.boundaryMesh().whichPatch(facei)
                    << " fc:" << mesh_.faceCentres()[facei] << endl
                    << "ownPatch:" << ownPatch[facei]
                    << " syncedOwnPatch:" << syncedOwnPatch[facei]
                    << " neiPatch:" << neiPatch[facei]
                    << " syncedNeiPatch:" << syncedNeiPatch[facei]
                    << abort(FatalError);
            }
        }
    }

    // Topochange container
    polyTopoChange meshMod(mesh_);

    label nBaffles = 0;

    forAll(ownPatch, facei)
    {
        if (ownPatch[facei] != -1)
        {
            // Create baffle or repatch face. Return label of inserted baffle
            // face.
            createBaffle
            (
                facei,
                ownPatch[facei],   // owner side patch
                neiPatch[facei],   // neighbour side patch
                meshMod
            );
            nBaffles++;
        }
    }

    // Change the mesh (no inflation, parallel sync)
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

    // Update fields
    mesh_.updateMesh(map);

    // Move mesh if in inflation mode
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes.
        mesh_.clearOut();
    }


    // Reset the instance for if in overwrite mode
    mesh_.setInstance(timeName());

    //- Redo the intersections on the newly create baffle faces. Note that
    //  this changes also the cell centre positions.
    faceSet baffledFacesSet(mesh_, "baffledFacesSet", 2*nBaffles);

    const labelList& reverseFaceMap = map().reverseFaceMap();
    const labelList& faceMap = map().faceMap();

    // Pick up owner side of baffle
    forAll(ownPatch, oldFacei)
    {
        label facei = reverseFaceMap[oldFacei];

        if (ownPatch[oldFacei] != -1 && facei >= 0)
        {
            const cell& ownFaces = mesh_.cells()[mesh_.faceOwner()[facei]];

            forAll(ownFaces, i)
            {
                baffledFacesSet.insert(ownFaces[i]);
            }
        }
    }
    // Pick up neighbour side of baffle (added faces)
    forAll(faceMap, facei)
    {
        label oldFacei = faceMap[facei];

        if (oldFacei >= 0 && reverseFaceMap[oldFacei] != facei)
        {
            const cell& ownFaces = mesh_.cells()[mesh_.faceOwner()[facei]];

            forAll(ownFaces, i)
            {
                baffledFacesSet.insert(ownFaces[i]);
            }
        }
    }
    baffledFacesSet.sync(mesh_);

    updateMesh(map, baffledFacesSet.toc());

    return map;
}


void Foam::meshRefinement::checkZoneFaces() const
{
    const faceZoneMesh& fZones = mesh_.faceZones();

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            forAll(pp, i)
            {
                label facei = pp.start()+i;
                label zoneI = fZones.whichZone(facei);

                if (zoneI != -1)
                {
                    FatalErrorInFunction
                        << "face:" << facei << " on patch " << pp.name()
                        << " is in zone " << fZones[zoneI].name()
                        << exit(FatalError);
                }
            }
        }
    }
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::createZoneBaffles
(
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,
    List<labelPair>& baffles
)
{
    const labelList zonedSurfaces
    (
        surfaceZonesInfo::getNamedSurfaces(surfaces_.surfZones())
    );

    autoPtr<mapPolyMesh> map;

    // No need to sync; all processors will have all same zonedSurfaces.
    if (zonedSurfaces.size())
    {
        // Split internal faces on interface surfaces
        Info<< "Converting zoned faces into baffles ..." << endl;

        // Get faces (internal only) to be baffled. Map from face to patch
        // label.
        Map<labelPair> faceToPatch
        (
            getZoneBafflePatches
            (
                false,
                globalToMasterPatch,
                globalToSlavePatch
            )
        );

        label nZoneFaces = returnReduce(faceToPatch.size(), sumOp<label>());
        if (nZoneFaces > 0)
        {
            // Convert into labelLists
            labelList ownPatch(mesh_.nFaces(), -1);
            labelList neiPatch(mesh_.nFaces(), -1);
            forAllConstIter(Map<labelPair>, faceToPatch, iter)
            {
                ownPatch[iter.key()] = iter().first();
                neiPatch[iter.key()] = iter().second();
            }

            // Create baffles. both sides same patch.
            map = createBaffles(ownPatch, neiPatch);

            // Get pairs of faces created.
            // Just loop over faceMap and store baffle if we encounter a slave
            // face.

            baffles.setSize(faceToPatch.size());
            label baffleI = 0;

            const labelList& faceMap = map().faceMap();
            const labelList& reverseFaceMap = map().reverseFaceMap();

            forAll(faceMap, facei)
            {
                label oldFacei = faceMap[facei];

                // Does face originate from face-to-patch
                Map<labelPair>::const_iterator iter = faceToPatch.find
                (
                    oldFacei
                );

                if (iter != faceToPatch.end())
                {
                    label masterFacei = reverseFaceMap[oldFacei];
                    if (facei != masterFacei)
                    {
                        baffles[baffleI++] = labelPair(masterFacei, facei);
                    }
                }
            }

            if (baffleI != faceToPatch.size())
            {
                FatalErrorInFunction
                    << "Had " << faceToPatch.size() << " patches to create "
                    << " but encountered " << baffleI
                    << " slave faces originating from patcheable faces."
                    << abort(FatalError);
            }

            if (debug&MESH)
            {
                const_cast<Time&>(mesh_.time())++;
                Pout<< "Writing zone-baffled mesh to time " << timeName()
                    << endl;
                write
                (
                    debugType(debug),
                    writeType(writeLevel() | WRITEMESH),
                    mesh_.time().path()/"baffles"
                );
            }
        }
        Info<< "Created " << nZoneFaces << " baffles in = "
            << mesh_.time().cpuTimeIncrement() << " s\n" << nl << endl;
    }
    return map;
}


Foam::List<Foam::labelPair> Foam::meshRefinement::freeStandingBaffles
(
    const List<labelPair>& couples,
    const scalar planarAngle
) const
{
    // Done by counting the number of baffles faces per mesh edge. If edge
    // has 2 boundary faces and both are baffle faces it is the edge of a baffle
    // region.

    // All duplicate faces on edge of the patch are to be merged.
    // So we count for all edges of duplicate faces how many duplicate
    // faces use them.
    labelList nBafflesPerEdge(mesh_.nEdges(), 0);


    // This algorithm is quite tricky. We don't want to use edgeFaces and
    // also want it to run in parallel so it is now an algorithm over
    // all (boundary) faces instead.
    // We want to pick up any edges that are only used by the baffle
    // or internal faces but not by any other boundary faces. So
    // - increment count on an edge by 1 if it is used by any (uncoupled)
    //   boundary face.
    // - increment count on an edge by 1000000 if it is used by a baffle face
    // - sum in parallel
    //
    // So now any edge that is used by baffle faces only will have the
    // value 2*1000000+2*1.


    const label baffleValue = 1000000;


    // Count number of boundary faces per edge
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        // Count number of boundary faces. Discard coupled boundary faces.
        if (!pp.coupled())
        {
            label facei = pp.start();

            forAll(pp, i)
            {
                const labelList& fEdges = mesh_.faceEdges(facei);

                forAll(fEdges, fEdgeI)
                {
                    nBafflesPerEdge[fEdges[fEdgeI]]++;
                }
                facei++;
            }
        }
    }


    DynamicList<label> fe0;
    DynamicList<label> fe1;


    // Count number of duplicate boundary faces per edge
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(couples, i)
    {
        {
            label f0 = couples[i].first();
            const labelList& fEdges0 = mesh_.faceEdges(f0, fe0);
            forAll(fEdges0, fEdgeI)
            {
                nBafflesPerEdge[fEdges0[fEdgeI]] += baffleValue;
            }
        }

        {
            label f1 = couples[i].second();
            const labelList& fEdges1 = mesh_.faceEdges(f1, fe1);
            forAll(fEdges1, fEdgeI)
            {
                nBafflesPerEdge[fEdges1[fEdgeI]] += baffleValue;
            }
        }
    }

    // Add nBaffles on shared edges
    syncTools::syncEdgeList
    (
        mesh_,
        nBafflesPerEdge,
        plusEqOp<label>(),  // in-place add
        label(0)            // initial value
    );


    // Baffles which are not next to other boundaries and baffles will have
    // nBafflesPerEdge value 2*baffleValue+2*1 (from 2 boundary faces which
    // are both baffle faces)

    List<labelPair> filteredCouples(couples.size());
    label filterI = 0;

    forAll(couples, i)
    {
        const labelPair& couple = couples[i];

        if
        (
            patches.whichPatch(couple.first())
         == patches.whichPatch(couple.second())
        )
        {
            const labelList& fEdges = mesh_.faceEdges(couple.first());

            forAll(fEdges, fEdgeI)
            {
                label edgeI = fEdges[fEdgeI];

                if (nBafflesPerEdge[edgeI] == 2*baffleValue+2*1)
                {
                    filteredCouples[filterI++] = couple;
                    break;
                }
            }
        }
    }
    filteredCouples.setSize(filterI);


    label nFiltered = returnReduce(filteredCouples.size(), sumOp<label>());

    Info<< "freeStandingBaffles : detected "
        << nFiltered
        << " free-standing baffles out of "
        << returnReduce(couples.size(), sumOp<label>())
        << nl << endl;


    if (nFiltered > 0)
    {
        // Collect segments
        // ~~~~~~~~~~~~~~~~

        pointField start(filteredCouples.size());
        pointField end(filteredCouples.size());

        const pointField& cellCentres = mesh_.cellCentres();

        forAll(filteredCouples, i)
        {
            const labelPair& couple = filteredCouples[i];
            start[i] = cellCentres[mesh_.faceOwner()[couple.first()]];
            end[i] = cellCentres[mesh_.faceOwner()[couple.second()]];
        }

        // Extend segments a bit
        {
            const vectorField smallVec(rootSmall*(end-start));
            start -= smallVec;
            end += smallVec;
        }


        // Do test for intersections
        // ~~~~~~~~~~~~~~~~~~~~~~~~~
        labelList surface1;
        List<pointIndexHit> hit1;
        labelList region1;
        vectorField normal1;

        labelList surface2;
        List<pointIndexHit> hit2;
        labelList region2;
        vectorField normal2;

        surfaces_.findNearestIntersection
        (
            identity(surfaces_.surfaces().size()),
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

        //mkDir(mesh_.time().path()/timeName());
        //OBJstream str
        //(
        //    mesh_.time().path()/timeName()/"flatBaffles.obj"
        //);

        const scalar planarAngleCos = Foam::cos(degToRad(planarAngle));

        label filterI = 0;
        forAll(filteredCouples, i)
        {
            const labelPair& couple = filteredCouples[i];

            if
            (
                hit1[i].hit()
             && hit2[i].hit()
             && (
                    surface1[i] != surface2[i]
                 || hit1[i].index() != hit2[i].index()
                )
            )
            {
                // Two different hits. Check angle.
                //str.write
                //(
                //    linePointRef(hit1[i].hitPoint(), hit2[i].hitPoint()),
                //    normal1[i],
                //    normal2[i]
                //);

                if ((normal1[i]&normal2[i]) > planarAngleCos)
                {
                    // Both normals aligned
                    vector n = end[i]-start[i];
                    scalar magN = mag(n);
                    if (magN > vSmall)
                    {
                        filteredCouples[filterI++] = couple;
                    }
                }
            }
            else if (hit1[i].hit() || hit2[i].hit())
            {
                // Single hit. Do not include in freestanding baffles.
            }
        }

        filteredCouples.setSize(filterI);

        Info<< "freeStandingBaffles : detected "
            << returnReduce(filterI, sumOp<label>())
            << " planar (within " << planarAngle
            << " degrees) free-standing baffles out of "
            << nFiltered
            << nl << endl;
    }

    return filteredCouples;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::mergeBaffles
(
    const List<labelPair>& couples
)
{
    // Mesh change engine
    polyTopoChange meshMod(mesh_);

    const faceList& faces = mesh_.faces();
    const labelList& faceOwner = mesh_.faceOwner();
    const faceZoneMesh& faceZones = mesh_.faceZones();

    forAll(couples, i)
    {
        label face0 = couples[i].first();
        label face1 = couples[i].second();

        // face1 < 0 signals a coupled face that has been converted to baffle.

        label own0 = faceOwner[face0];
        label own1 = faceOwner[face1];

        if (face1 < 0 || own0 < own1)
        {
            // Use face0 as the new internal face.
            label zoneID = faceZones.whichZone(face0);
            bool zoneFlip = false;

            if (zoneID >= 0)
            {
                const faceZone& fZone = faceZones[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(face0)];
            }

            label nei = (face1 < 0 ? -1 : own1);

            meshMod.setAction(polyRemoveFace(face1));
            meshMod.setAction
            (
                polyModifyFace
                (
                    faces[face0],           // modified face
                    face0,                  // label of face being modified
                    own0,                   // owner
                    nei,                    // neighbour
                    false,                  // face flip
                    -1,                     // patch for face
                    false,                  // remove from zone
                    zoneID,                 // zone for face
                    zoneFlip                // face flip in zone
                )
            );
        }
        else
        {
            // Use face1 as the new internal face.
            label zoneID = faceZones.whichZone(face1);
            bool zoneFlip = false;

            if (zoneID >= 0)
            {
                const faceZone& fZone = faceZones[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(face1)];
            }

            meshMod.setAction(polyRemoveFace(face0));
            meshMod.setAction
            (
                polyModifyFace
                (
                    faces[face1],           // modified face
                    face1,                  // label of face being modified
                    own1,                   // owner
                    own0,                   // neighbour
                    false,                  // face flip
                    -1,                     // patch for face
                    false,                  // remove from zone
                    zoneID,                 // zone for face
                    zoneFlip                // face flip in zone
                )
            );
        }
    }

    // Change the mesh (no inflation)
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

    // Update fields
    mesh_.updateMesh(map);

    // Move mesh (since morphing does not do this)
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes.
        mesh_.clearOut();
    }

    // Reset the instance for if in overwrite mode
    mesh_.setInstance(timeName());

    // Update intersections. Recalculate intersections on merged faces since
    // this seems to give problems? Note: should not be necessary since
    // baffles preserve intersections from when they were created.
    labelList newExposedFaces(2*couples.size());
    label newI = 0;

    forAll(couples, i)
    {
        label newFace0 = map().reverseFaceMap()[couples[i].first()];
        if (newFace0 != -1)
        {
            newExposedFaces[newI++] = newFace0;
        }

        label newFace1 = map().reverseFaceMap()[couples[i].second()];
        if (newFace1 != -1)
        {
            newExposedFaces[newI++] = newFace1;
        }
    }
    newExposedFaces.setSize(newI);
    updateMesh(map, newExposedFaces);

    return map;
}


// Finds region per cell for cells inside closed named surfaces
void Foam::meshRefinement::findCellZoneGeometric
(
    const pointField& neiCc,
    const labelList& closedNamedSurfaces,   // indices of closed surfaces
    labelList& namedSurfaceIndex,           // per face index of named surface
    const labelList& surfaceToCellZone,     // cell zone index per surface

    labelList& cellToZone
) const
{
    const pointField& cellCentres = mesh_.cellCentres();
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    // Check if cell centre is inside
    labelList insideSurfaces;
    surfaces_.findInside
    (
        closedNamedSurfaces,
        cellCentres,
        insideSurfaces
    );

    forAll(insideSurfaces, celli)
    {
        if (cellToZone[celli] == -2)
        {
            label surfI = insideSurfaces[celli];

            if (surfI != -1)
            {
                cellToZone[celli] = surfaceToCellZone[surfI];
            }
        }
    }


    // Some cells with cell centres close to surface might have
    // had been put into wrong surface. Recheck with perturbed cell centre.


    // 1. Collect points

    // Count points to test.
    label nCandidates = 0;
    forAll(namedSurfaceIndex, facei)
    {
        label surfI = namedSurfaceIndex[facei];

        if (surfI != -1)
        {
            if (mesh_.isInternalFace(facei))
            {
                nCandidates += 2;
            }
            else
            {
                nCandidates += 1;
            }
        }
    }

    // Collect points.
    pointField candidatePoints(nCandidates);
    nCandidates = 0;
    forAll(namedSurfaceIndex, facei)
    {
        label surfI = namedSurfaceIndex[facei];

        if (surfI != -1)
        {
            label own = faceOwner[facei];
            const point& ownCc = cellCentres[own];

            if (mesh_.isInternalFace(facei))
            {
                label nei = faceNeighbour[facei];
                const point& neiCc = cellCentres[nei];
                // Perturbed cc
                const vector d = 1e-4*(neiCc - ownCc);
                candidatePoints[nCandidates++] = ownCc-d;
                candidatePoints[nCandidates++] = neiCc+d;
            }
            else
            {
                //const point& neiFc = mesh_.faceCentres()[facei];
                const point& neiFc = neiCc[facei-mesh_.nInternalFaces()];

                // Perturbed cc
                const vector d = 1e-4*(neiFc - ownCc);
                candidatePoints[nCandidates++] = ownCc-d;
            }
        }
    }


    // 2. Test points for inside

    surfaces_.findInside
    (
        closedNamedSurfaces,
        candidatePoints,
        insideSurfaces
    );


    // 3. Update zone information

    nCandidates = 0;
    forAll(namedSurfaceIndex, facei)
    {
        label surfI = namedSurfaceIndex[facei];

        if (surfI != -1)
        {
            label own = faceOwner[facei];

            if (mesh_.isInternalFace(facei))
            {
                label ownSurfI = insideSurfaces[nCandidates++];
                if (ownSurfI != -1)
                {
                    cellToZone[own] = surfaceToCellZone[ownSurfI];
                }

                label neiSurfI = insideSurfaces[nCandidates++];
                if (neiSurfI != -1)
                {
                    label nei = faceNeighbour[facei];

                    cellToZone[nei] = surfaceToCellZone[neiSurfI];
                }
            }
            else
            {
                label ownSurfI = insideSurfaces[nCandidates++];
                if (ownSurfI != -1)
                {
                    cellToZone[own] = surfaceToCellZone[ownSurfI];
                }
            }
        }
    }


    // Adapt the namedSurfaceIndex
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // for if any cells were not completely covered.

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        label ownZone = cellToZone[mesh_.faceOwner()[facei]];
        label neiZone = cellToZone[mesh_.faceNeighbour()[facei]];

        if (namedSurfaceIndex[facei] == -1 && (ownZone != neiZone))
        {
            // Give face the zone of max cell zone
            namedSurfaceIndex[facei] = findIndex
            (
                surfaceToCellZone,
                max(ownZone, neiZone)
            );
        }
    }

    labelList neiCellZone(mesh_.nFaces()-mesh_.nInternalFaces());
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label facei = pp.start()+i;
                label ownZone = cellToZone[mesh_.faceOwner()[facei]];
                neiCellZone[facei-mesh_.nInternalFaces()] = ownZone;
            }
        }
    }
    syncTools::swapBoundaryFaceList(mesh_, neiCellZone);

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label facei = pp.start()+i;
                label ownZone = cellToZone[mesh_.faceOwner()[facei]];
                label neiZone = neiCellZone[facei-mesh_.nInternalFaces()];

                if (namedSurfaceIndex[facei] == -1 && (ownZone != neiZone))
                {
                    // Give face the max cell zone
                    namedSurfaceIndex[facei] = findIndex
                    (
                        surfaceToCellZone,
                        max(ownZone, neiZone)
                    );
                }
            }
        }
    }

    // Sync
    syncTools::syncFaceList(mesh_, namedSurfaceIndex, maxEqOp<label>());
}


void Foam::meshRefinement::findCellZoneInsideWalk
(
    const labelList& locationSurfaces,  // indices of surfaces with inside point
    const labelList& namedSurfaceIndex, // per face index of named surface
    const labelList& surfaceToCellZone, // cell zone index per surface

    labelList& cellToZone
) const
{
    // Analyse regions. Reuse regionsplit
    boolList blockedFace(mesh_.nFaces());
    //selectSeparatedCoupledFaces(blockedFace);

    forAll(namedSurfaceIndex, facei)
    {
        if (namedSurfaceIndex[facei] == -1)
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
    regionSplit cellRegion(mesh_, blockedFace);
    blockedFace.clear();


    // Force calculation of face decomposition (used in findCell)
    (void)mesh_.tetBasePtIs();

    const PtrList<surfaceZonesInfo>& surfZones = surfaces_.surfZones();

    // For all locationSurface find the cell
    forAll(locationSurfaces, i)
    {
        label surfI = locationSurfaces[i];

        const point& insidePoint = surfZones[surfI].zoneInsidePoint();

        Info<< "For surface " << surfaces_.names()[surfI]
            << " finding inside point " << insidePoint
            << endl;

        // Find the region containing the insidePoint
        label keepRegionI = findRegion
        (
            mesh_,
            cellRegion,
            mergeDistance_*vector(1,1,1),
            insidePoint
        );

        Info<< "For surface " << surfaces_.names()[surfI]
            << " found point " << insidePoint
            << " in global region " << keepRegionI
            << " out of " << cellRegion.nRegions() << " regions." << endl;

        if (keepRegionI == -1)
        {
            FatalErrorInFunction
                << "Point " << insidePoint
                << " is not inside the mesh." << nl
                << "Bounding box of the mesh:" << mesh_.bounds()
                << exit(FatalError);
        }

        // Set all cells with this region
        forAll(cellRegion, celli)
        {
            if (cellRegion[celli] == keepRegionI)
            {
                if (cellToZone[celli] == -2)
                {
                    cellToZone[celli] = surfaceToCellZone[surfI];
                }
                else if (cellToZone[celli] != surfaceToCellZone[surfI])
                {
                    WarningInFunction
                        << "Cell " << celli
                        << " at " << mesh_.cellCentres()[celli]
                        << " is inside surface " << surfaces_.names()[surfI]
                        << " but already marked as being in zone "
                        << cellToZone[celli] << endl
                        << "This can happen if your surfaces are not"
                        << " (sufficiently) closed."
                        << endl;
                }
            }
        }
    }
}


bool Foam::meshRefinement::calcRegionToZone
(
    const label surfZoneI,
    const label ownRegion,
    const label neiRegion,

    labelList& regionToCellZone
) const
{
    bool changed = false;

    // Check whether in between different regions
    if (ownRegion != neiRegion)
    {
        // Jump. Change one of the sides to my type.

        // 1. Interface between my type and unset region.
        // Set region to keepRegion

        if (regionToCellZone[ownRegion] == -2)
        {
            if (regionToCellZone[neiRegion] == surfZoneI)
            {
                // Face between unset and my region. Put unset
                // region into keepRegion
                regionToCellZone[ownRegion] = -1;
                changed = true;
            }
            else if (regionToCellZone[neiRegion] != -2)
            {
                // Face between unset and other region.
                // Put unset region into my region
                regionToCellZone[ownRegion] = surfZoneI;
                changed = true;
            }
        }
        else if (regionToCellZone[neiRegion] == -2)
        {
            if (regionToCellZone[ownRegion] == surfZoneI)
            {
                // Face between unset and my region. Put unset
                // region into keepRegion
                regionToCellZone[neiRegion] = -1;
                changed = true;
            }
            else if (regionToCellZone[ownRegion] != -2)
            {
                // Face between unset and other region.
                // Put unset region into my region
                regionToCellZone[neiRegion] = surfZoneI;
                changed = true;
            }
        }
    }
    return changed;
}


void Foam::meshRefinement::findCellZoneTopo
(
    const point& keepPoint,
    const labelList& namedSurfaceIndex,
    const labelList& surfaceToCellZone,
    labelList& cellToZone
) const
{
    // Assumes:
    // - region containing keepPoint does not go into a cellZone
    // - all other regions can be found by crossing faces marked in
    //   namedSurfaceIndex.

    // Analyse regions. Reuse regionsplit
    boolList blockedFace(mesh_.nFaces());

    forAll(namedSurfaceIndex, facei)
    {
        if (namedSurfaceIndex[facei] == -1)
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
    regionSplit cellRegion(mesh_, blockedFace);
    blockedFace.clear();

    // Per mesh region the zone the cell should be put in.
    // -2   : not analysed yet
    // -1   : keepPoint region. Not put into any cellzone.
    // >= 0 : index of cellZone
    labelList regionToCellZone(cellRegion.nRegions(), -2);

    // See which cells already are set in the cellToZone (from geometric
    // searching) and use these to take over their zones.
    // Note: could be improved to count number of cells per region.
    forAll(cellToZone, celli)
    {
        if (cellToZone[celli] != -2)
        {
            regionToCellZone[cellRegion[celli]] = cellToZone[celli];
        }
    }


    // Find the region containing the keepPoint
    label keepRegionI = findRegion
    (
        mesh_,
        cellRegion,
        mergeDistance_*vector(1,1,1),
        keepPoint
    );

    Info<< "Found point " << keepPoint
        << " in global region " << keepRegionI
        << " out of " << cellRegion.nRegions() << " regions." << endl;

    if (keepRegionI == -1)
    {
        FatalErrorInFunction
            << "Point " << keepPoint
            << " is not inside the mesh." << nl
            << "Bounding box of the mesh:" << mesh_.bounds()
            << exit(FatalError);
    }

    // Mark default region with zone -1.
    if (regionToCellZone[keepRegionI] == -2)
    {
        regionToCellZone[keepRegionI] = -1;
    }


    // Find correspondence between cell zone and surface
    // by changing cell zone every time we cross a surface.
    while (true)
    {
        // Synchronise regionToCellZone.
        // Note:
        // - region numbers are identical on all processors
        // - keepRegion is identical ,,
        // - cellZones are identical ,,
        // This done at top of loop to account for geometric matching
        // not being synchronised.
        Pstream::listCombineGather(regionToCellZone, maxEqOp<label>());
        Pstream::listCombineScatter(regionToCellZone);


        bool changed = false;

        // Internal faces

        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            label surfI = namedSurfaceIndex[facei];

            // Connected even if no cellZone defined for surface
            if (surfI != -1)
            {
                // Calculate region to zone from cellRegions on either side
                // of internal face.
                bool changedCell = calcRegionToZone
                (
                    surfaceToCellZone[surfI],
                    cellRegion[mesh_.faceOwner()[facei]],
                    cellRegion[mesh_.faceNeighbour()[facei]],
                    regionToCellZone
                );

                changed = changed | changedCell;
            }
        }

        // Boundary faces

        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        // Get coupled neighbour cellRegion
        labelList neiCellRegion(mesh_.nFaces()-mesh_.nInternalFaces());
        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];

            if (pp.coupled())
            {
                forAll(pp, i)
                {
                    label facei = pp.start()+i;
                    neiCellRegion[facei-mesh_.nInternalFaces()] =
                        cellRegion[mesh_.faceOwner()[facei]];
                }
            }
        }
        syncTools::swapBoundaryFaceList(mesh_, neiCellRegion);

        // Calculate region to zone from cellRegions on either side of coupled
        // face.
        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];

            if (pp.coupled())
            {
                forAll(pp, i)
                {
                    label facei = pp.start()+i;

                    label surfI = namedSurfaceIndex[facei];

                    // Connected even if no cellZone defined for surface
                    if (surfI != -1)
                    {
                        bool changedCell = calcRegionToZone
                        (
                            surfaceToCellZone[surfI],
                            cellRegion[mesh_.faceOwner()[facei]],
                            neiCellRegion[facei-mesh_.nInternalFaces()],
                            regionToCellZone
                        );

                        changed = changed | changedCell;
                    }
                }
            }
        }

        if (!returnReduce(changed, orOp<bool>()))
        {
            break;
        }
    }


    forAll(regionToCellZone, regionI)
    {
        label zoneI = regionToCellZone[regionI];

        if (zoneI ==  -2)
        {
            FatalErrorInFunction
                << "For region " << regionI << " haven't set cell zone."
                << exit(FatalError);
        }
    }

    if (debug)
    {
        forAll(regionToCellZone, regionI)
        {
            Pout<< "Region " << regionI
                << " becomes cellZone:" << regionToCellZone[regionI]
                << endl;
        }
    }

    // Rework into cellToZone
    forAll(cellToZone, celli)
    {
        cellToZone[celli] = regionToCellZone[cellRegion[celli]];
    }
}


void Foam::meshRefinement::makeConsistentFaceIndex
(
    const labelList& cellToZone,
    labelList& namedSurfaceIndex
) const
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        label ownZone = cellToZone[faceOwner[facei]];
        label neiZone = cellToZone[faceNeighbour[facei]];

        if (ownZone == neiZone && namedSurfaceIndex[facei] != -1)
        {
            namedSurfaceIndex[facei] = -1;
        }
        else if (ownZone != neiZone && namedSurfaceIndex[facei] == -1)
        {
            FatalErrorInFunction
                << "Different cell zones on either side of face " << facei
                << " at " << mesh_.faceCentres()[facei]
                << " but face not marked with a surface."
                << abort(FatalError);
        }
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Get coupled neighbour cellZone
    labelList neiCellZone(mesh_.nFaces()-mesh_.nInternalFaces());
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label facei = pp.start()+i;
                neiCellZone[facei-mesh_.nInternalFaces()] =
                    cellToZone[mesh_.faceOwner()[facei]];
            }
        }
    }
    syncTools::swapBoundaryFaceList(mesh_, neiCellZone);

    // Use coupled cellZone to do check
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label facei = pp.start()+i;

                label ownZone = cellToZone[faceOwner[facei]];
                label neiZone = neiCellZone[facei-mesh_.nInternalFaces()];

                if (ownZone == neiZone && namedSurfaceIndex[facei] != -1)
                {
                    namedSurfaceIndex[facei] = -1;
                }
                else if (ownZone != neiZone && namedSurfaceIndex[facei] == -1)
                {
                    FatalErrorInFunction
                        << "Different cell zones on either side of face "
                        << facei << " at " << mesh_.faceCentres()[facei]
                        << " but face not marked with a surface."
                        << abort(FatalError);
                }
            }
        }
        else
        {
            // Unzonify boundary faces
            forAll(pp, i)
            {
                label facei = pp.start()+i;
                namedSurfaceIndex[facei] = -1;
            }
        }
    }
}


void Foam::meshRefinement::handleSnapProblems
(
    const snapParameters& snapParams,
    const bool useTopologicalSnapDetection,
    const bool removeEdgeConnectedCells,
    const scalarField& perpendicularAngle,
    const dictionary& motionDict,
    Time& runTime,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch
)
{
    Info<< nl
        << "Introducing baffles to block off problem cells" << nl
        << "----------------------------------------------" << nl
        << endl;

    labelList facePatch;
    if (useTopologicalSnapDetection)
    {
        facePatch = markFacesOnProblemCells
        (
            motionDict,
            removeEdgeConnectedCells,
            perpendicularAngle,
            globalToMasterPatch
        );
    }
    else
    {
        facePatch = markFacesOnProblemCellsGeometric(snapParams, motionDict);
    }
    Info<< "Analyzed problem cells in = "
        << runTime.cpuTimeIncrement() << " s\n" << nl << endl;

    if (debug&MESH)
    {
        faceSet problemFaces(mesh_, "problemFaces", mesh_.nFaces()/100);

        forAll(facePatch, facei)
        {
            if (facePatch[facei] != -1)
            {
                problemFaces.insert(facei);
            }
        }
        problemFaces.instance() = timeName();
        Pout<< "Dumping " << problemFaces.size()
            << " problem faces to " << problemFaces.objectPath() << endl;
        problemFaces.write();
    }

    Info<< "Introducing baffles to delete problem cells." << nl << endl;

    if (debug)
    {
        runTime++;
    }

    // Create baffles with same owner and neighbour for now.
    createBaffles(facePatch, facePatch);

    if (debug)
    {
        // Debug:test all is still synced across proc patches
        checkData();
    }
    Info<< "Created baffles in = "
        << runTime.cpuTimeIncrement() << " s\n" << nl << endl;

    printMeshInfo(debug, "After introducing baffles");

    if (debug&MESH)
    {
        Pout<< "Writing extra baffled mesh to time "
            << timeName() << endl;
        write
        (
            debugType(debug),
            writeType(writeLevel() | WRITEMESH),
            runTime.path()/"extraBaffles"
        );
        Pout<< "Dumped debug data in = "
            << runTime.cpuTimeIncrement() << " s\n" << nl << endl;
    }
}


Foam::labelList Foam::meshRefinement::freeStandingBaffleFaces
(
    const labelList& faceToZone,
    const labelList& cellToZone,
    const labelList& neiCellZone
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();


    // We want to pick up the faces to orient. These faces come in
    // two variants:
    // - faces originating from stand-alone faceZones
    //   (these will most likely have no cellZone on either side so
    //    ownZone and neiZone both -1)
    // - sticky-up faces originating from a 'bulge' in a outside of
    //   a cellZone. These will have the same cellZone on either side.
    //   How to orient these is not really clearly defined so do them
    //   same as stand-alone faceZone faces for now. (Normally these will
    //   already have been removed by the 'allowFreeStandingZoneFaces=false'
    //   default setting)

    // Note that argument neiCellZone will have -1 on uncoupled boundaries.

    DynamicList<label> faceLabels(mesh_.nFaces()/100);

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        if (faceToZone[facei] != -1)
        {
            // Free standing baffle?
            label ownZone = cellToZone[faceOwner[facei]];
            label neiZone = cellToZone[faceNeighbour[facei]];
            if (ownZone == neiZone)
            {
                faceLabels.append(facei);
            }
        }
    }
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        forAll(pp, i)
        {
            label facei = pp.start()+i;
            if (faceToZone[facei] != -1)
            {
                // Free standing baffle?
                label ownZone = cellToZone[faceOwner[facei]];
                label neiZone = neiCellZone[facei-mesh_.nInternalFaces()];
                if (ownZone == neiZone)
                {
                    faceLabels.append(facei);
                }
            }
        }
    }
    return faceLabels.shrink();
}


void Foam::meshRefinement::calcPatchNumMasterFaces
(
    const PackedBoolList& isMasterFace,
    const indirectPrimitivePatch& patch,
    labelList& nMasterFacesPerEdge
) const
{
    // Number of (master)faces per edge
    nMasterFacesPerEdge.setSize(patch.nEdges());
    nMasterFacesPerEdge = 0;

    forAll(patch.addressing(), facei)
    {
        const label meshFacei = patch.addressing()[facei];

        if (isMasterFace[meshFacei])
        {
            const labelList& fEdges = patch.faceEdges()[facei];
            forAll(fEdges, fEdgeI)
            {
                nMasterFacesPerEdge[fEdges[fEdgeI]]++;
            }
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        patch.meshEdges(mesh_.edges(), mesh_.pointEdges()),
        nMasterFacesPerEdge,
        plusEqOp<label>(),
        label(0)
    );
}


Foam::label Foam::meshRefinement::markPatchZones
(
    const indirectPrimitivePatch& patch,
    const labelList& nMasterFacesPerEdge,
    labelList& faceToZone
) const
{
    List<patchEdgeFaceRegion> allEdgeInfo(patch.nEdges());
    List<patchEdgeFaceRegion> allFaceInfo(patch.size());


    // Protect all non-manifold edges
    {
        label nProtected = 0;

        forAll(nMasterFacesPerEdge, edgeI)
        {
            if (nMasterFacesPerEdge[edgeI] > 2)
            {
                allEdgeInfo[edgeI] = -2;
                nProtected++;
            }
        }
        //Info<< "Protected from visiting "
        //    << returnReduce(nProtected, sumOp<label>())
        //    << " non-manifold edges" << nl << endl;
    }


    // Hand out zones

    DynamicList<label> changedEdges;
    DynamicList<patchEdgeFaceRegion> changedInfo;

    const scalar tol = PatchEdgeFaceWave
    <
        indirectPrimitivePatch,
        patchEdgeFaceRegion
    >::propagationTol();

    int dummyTrackData;

    const globalIndex globalFaces(patch.size());

    label facei = 0;

    label currentZoneI = 0;

    while (true)
    {
        // Pick an unset face
        label globalSeed = labelMax;
        for (; facei < allFaceInfo.size(); facei++)
        {
            if (!allFaceInfo[facei].valid(dummyTrackData))
            {
                globalSeed = globalFaces.toGlobal(facei);
                break;
            }
        }

        reduce(globalSeed, minOp<label>());

        if (globalSeed == labelMax)
        {
            break;
        }

        label proci = globalFaces.whichProcID(globalSeed);
        label seedFacei = globalFaces.toLocal(proci, globalSeed);

        //Info<< "Seeding zone " << currentZoneI
        //    << " from processor " << proci << " face " << seedFacei
        //    << endl;

        if (proci == Pstream::myProcNo())
        {
            patchEdgeFaceRegion& faceInfo = allFaceInfo[seedFacei];


            // Set face
            faceInfo = currentZoneI;

            // .. and seed its edges
            const labelList& fEdges = patch.faceEdges()[seedFacei];
            forAll(fEdges, fEdgeI)
            {
                label edgeI = fEdges[fEdgeI];

                patchEdgeFaceRegion& edgeInfo = allEdgeInfo[edgeI];

                if
                (
                    edgeInfo.updateEdge<int>
                    (
                        mesh_,
                        patch,
                        edgeI,
                        seedFacei,
                        faceInfo,
                        tol,
                        dummyTrackData
                    )
                )
                {
                    changedEdges.append(edgeI);
                    changedInfo.append(edgeInfo);
                }
            }
        }


        if (returnReduce(changedEdges.size(), sumOp<label>()) == 0)
        {
            break;
        }


        // Walk
        PatchEdgeFaceWave
        <
            indirectPrimitivePatch,
            patchEdgeFaceRegion
        > calc
        (
            mesh_,
            patch,
            changedEdges,
            changedInfo,
            allEdgeInfo,
            allFaceInfo,
            returnReduce(patch.nEdges(), sumOp<label>())
        );

        currentZoneI++;
    }


    faceToZone.setSize(patch.size());
    forAll(allFaceInfo, facei)
    {
        if (!allFaceInfo[facei].valid(dummyTrackData))
        {
            FatalErrorInFunction
                << "Problem: unvisited face " << facei
                << " at " << patch.faceCentres()[facei]
                << exit(FatalError);
        }
        faceToZone[facei] = allFaceInfo[facei].region();
    }

    return currentZoneI;
}


void Foam::meshRefinement::consistentOrientation
(
    const PackedBoolList& isMasterFace,
    const indirectPrimitivePatch& patch,
    const labelList& nMasterFacesPerEdge,
    const labelList& faceToZone,
    const Map<label>& zoneToOrientation,
    boolList& meshFlipMap
) const
{
    const polyBoundaryMesh& bm = mesh_.boundaryMesh();

    // Data on all edges and faces
    List<patchFaceOrientation> allEdgeInfo(patch.nEdges());
    List<patchFaceOrientation> allFaceInfo(patch.size());

    // Make sure we don't walk through
    // - slaves of coupled faces
    // - non-manifold edges
    {
        label nProtected = 0;

        forAll(patch.addressing(), facei)
        {
            const label meshFacei = patch.addressing()[facei];
            const label patchi = bm.whichPatch(meshFacei);

            if
            (
                patchi != -1
             && bm[patchi].coupled()
             && !isMasterFace[meshFacei]
            )
            {
                // Slave side. Mark so doesn't get visited.
                allFaceInfo[facei] = orientedSurface::NOFLIP;
                nProtected++;
            }
        }
        //Info<< "Protected from visiting "
        //    << returnReduce(nProtected, sumOp<label>())
        //    << " slaves of coupled faces" << nl << endl;
    }
    {
        label nProtected = 0;

        forAll(nMasterFacesPerEdge, edgeI)
        {
            if (nMasterFacesPerEdge[edgeI] > 2)
            {
                allEdgeInfo[edgeI] = orientedSurface::NOFLIP;
                nProtected++;
            }
        }

        Info<< "Protected from visiting "
            << returnReduce(nProtected, sumOp<label>())
            << " non-manifold edges" << nl << endl;
    }



    DynamicList<label> changedEdges;
    DynamicList<patchFaceOrientation> changedInfo;

    const scalar tol = PatchEdgeFaceWave
    <
        indirectPrimitivePatch,
        patchFaceOrientation
    >::propagationTol();

    int dummyTrackData;

    globalIndex globalFaces(patch.size());

    while (true)
    {
        // Pick an unset face
        label globalSeed = labelMax;
        forAll(allFaceInfo, facei)
        {
            if (allFaceInfo[facei] == orientedSurface::UNVISITED)
            {
                globalSeed = globalFaces.toGlobal(facei);
                break;
            }
        }

        reduce(globalSeed, minOp<label>());

        if (globalSeed == labelMax)
        {
            break;
        }

        label proci = globalFaces.whichProcID(globalSeed);
        label seedFacei = globalFaces.toLocal(proci, globalSeed);

        //Info<< "Seeding from processor " << proci << " face " << seedFacei
        //    << endl;

        if (proci == Pstream::myProcNo())
        {
            // Determine orientation of seedFace

            patchFaceOrientation& faceInfo = allFaceInfo[seedFacei];

            // Start off with correct orientation
            faceInfo = orientedSurface::NOFLIP;

            if (zoneToOrientation[faceToZone[seedFacei]] < 0)
            {
                faceInfo.flip();
            }


            const labelList& fEdges = patch.faceEdges()[seedFacei];
            forAll(fEdges, fEdgeI)
            {
                label edgeI = fEdges[fEdgeI];

                patchFaceOrientation& edgeInfo = allEdgeInfo[edgeI];

                if
                (
                    edgeInfo.updateEdge<int>
                    (
                        mesh_,
                        patch,
                        edgeI,
                        seedFacei,
                        faceInfo,
                        tol,
                        dummyTrackData
                    )
                )
                {
                    changedEdges.append(edgeI);
                    changedInfo.append(edgeInfo);
                }
            }
        }


        if (returnReduce(changedEdges.size(), sumOp<label>()) == 0)
        {
            break;
        }



        // Walk
        PatchEdgeFaceWave
        <
            indirectPrimitivePatch,
            patchFaceOrientation
        > calc
        (
            mesh_,
            patch,
            changedEdges,
            changedInfo,
            allEdgeInfo,
            allFaceInfo,
            returnReduce(patch.nEdges(), sumOp<label>())
        );
    }


    // Push master zone info over to slave (since slave faces never visited)
    {
        labelList neiStatus
        (
            mesh_.nFaces()-mesh_.nInternalFaces(),
            orientedSurface::UNVISITED
        );

        forAll(patch.addressing(), i)
        {
            const label meshFacei = patch.addressing()[i];
            if (!mesh_.isInternalFace(meshFacei))
            {
                neiStatus[meshFacei-mesh_.nInternalFaces()] =
                    allFaceInfo[i].flipStatus();
            }
        }
        syncTools::swapBoundaryFaceList(mesh_, neiStatus);

        forAll(patch.addressing(), i)
        {
            const label meshFacei = patch.addressing()[i];
            const label patchi = bm.whichPatch(meshFacei);

            if
            (
                patchi != -1
             && bm[patchi].coupled()
             && !isMasterFace[meshFacei]
            )
            {
                // Slave side. Take flipped from neighbour
                label bFacei = meshFacei-mesh_.nInternalFaces();

                if (neiStatus[bFacei] == orientedSurface::NOFLIP)
                {
                    allFaceInfo[i] = orientedSurface::FLIP;
                }
                else if (neiStatus[bFacei] == orientedSurface::FLIP)
                {
                    allFaceInfo[i] = orientedSurface::NOFLIP;
                }
                else
                {
                    FatalErrorInFunction
                        << "Incorrect status for face " << meshFacei
                        << abort(FatalError);
                }
            }
        }
    }


    // Convert to meshFlipMap and adapt faceZones

    meshFlipMap.setSize(mesh_.nFaces());
    meshFlipMap = false;

    forAll(allFaceInfo, facei)
    {
        label meshFacei = patch.addressing()[facei];

        if (allFaceInfo[facei] == orientedSurface::NOFLIP)
        {
            meshFlipMap[meshFacei] = false;
        }
        else if (allFaceInfo[facei] == orientedSurface::FLIP)
        {
            meshFlipMap[meshFacei] = true;
        }
        else
        {
            FatalErrorInFunction
                << "Problem : unvisited face " << facei
                << " centre:" << mesh_.faceCentres()[meshFacei]
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meshRefinement::baffleAndSplitMesh
(
    const bool doHandleSnapProblems,
    const snapParameters& snapParams,
    const bool useTopologicalSnapDetection,
    const bool removeEdgeConnectedCells,
    const scalarField& perpendicularAngle,
    const bool mergeFreeStanding,
    const scalar planarAngle,
    const dictionary& motionDict,
    Time& runTime,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,
    const point& keepPoint
)
{
    // Introduce baffles
    // ~~~~~~~~~~~~~~~~~

    // Split the mesh along internal faces wherever there is a pierce between
    // two cell centres.

    Info<< "Introducing baffles for "
        << returnReduce(countHits(), sumOp<label>())
        << " faces that are intersected by the surface." << nl << endl;

    // Swap neighbouring cell centres and cell level
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
    calcNeighbourData(neiLevel, neiCc);

    labelList ownPatch, neiPatch;
    getBafflePatches
    (
        globalToMasterPatch,
        neiLevel,
        neiCc,

        ownPatch,
        neiPatch
    );

    createBaffles(ownPatch, neiPatch);

    if (debug)
    {
        // Debug:test all is still synced across proc patches
        checkData();
    }

    Info<< "Created baffles in = "
        << runTime.cpuTimeIncrement() << " s\n" << nl << endl;

    printMeshInfo(debug, "After introducing baffles");

    if (debug&MESH)
    {
        Pout<< "Writing baffled mesh to time " << timeName()
            << endl;
        write
        (
            debugType(debug),
            writeType(writeLevel() | WRITEMESH),
            runTime.path()/"baffles"
        );
        Pout<< "Dumped debug data in = "
            << runTime.cpuTimeIncrement() << " s\n" << nl << endl;
    }


    // Introduce baffles to delete problem cells
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Create some additional baffles where we want surface cells removed.

    if (doHandleSnapProblems)
    {
        handleSnapProblems
        (
            snapParams,
            useTopologicalSnapDetection,
            removeEdgeConnectedCells,
            perpendicularAngle,
            motionDict,
            runTime,
            globalToMasterPatch,
            globalToSlavePatch
        );
    }


    // Select part of mesh
    // ~~~~~~~~~~~~~~~~~~~

    Info<< nl
        << "Remove unreachable sections of mesh" << nl
        << "-----------------------------------" << nl
        << endl;

    if (debug)
    {
        runTime++;
    }

    splitMeshRegions(globalToMasterPatch, globalToSlavePatch, keepPoint);

    if (debug)
    {
        // Debug:test all is still synced across proc patches
        checkData();
    }
    Info<< "Split mesh in = "
        << runTime.cpuTimeIncrement() << " s\n" << nl << endl;

    printMeshInfo(debug, "After subsetting");

    if (debug&MESH)
    {
        Pout<< "Writing subsetted mesh to time " << timeName()
            << endl;
        write
        (
            debugType(debug),
            writeType(writeLevel() | WRITEMESH),
            runTime.path()/timeName()
        );
        Pout<< "Dumped debug data in = "
            << runTime.cpuTimeIncrement() << " s\n" << nl << endl;
    }


    // Merge baffles
    // ~~~~~~~~~~~~~

    if (mergeFreeStanding)
    {
        Info<< nl
            << "Merge free-standing baffles" << nl
            << "---------------------------" << nl
            << endl;


        // List of pairs of freestanding baffle faces.
        List<labelPair> couples
        (
            freeStandingBaffles    // filter out freestanding baffles
            (
                localPointRegion::findDuplicateFacePairs(mesh_),
                planarAngle
            )
        );

        label nCouples = couples.size();
        reduce(nCouples, sumOp<label>());

        Info<< "Detected free-standing baffles : " << nCouples << endl;

        if (nCouples > 0)
        {
            // Actually merge baffles. Note: not exactly parallellized. Should
            // convert baffle faces into processor faces if they resulted
            // from them.
            mergeBaffles(couples);

            // Detect any problem cells resulting from merging of baffles
            // and delete them
            handleSnapProblems
            (
                snapParams,
                useTopologicalSnapDetection,
                removeEdgeConnectedCells,
                perpendicularAngle,
                motionDict,
                runTime,
                globalToMasterPatch,
                globalToSlavePatch
            );

            if (debug)
            {
                // Debug:test all is still synced across proc patches
                checkData();
            }
        }
        Info<< "Merged free-standing baffles in = "
            << runTime.cpuTimeIncrement() << " s\n" << nl << endl;
    }
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::splitMesh
(
    const label nBufferLayers,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,
    const point& keepPoint
)
{
    // Determine patches to put intersections into
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Swap neighbouring cell centres and cell level
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
    calcNeighbourData(neiLevel, neiCc);

    labelList ownPatch, neiPatch;
    getBafflePatches
    (
        globalToMasterPatch,
        neiLevel,
        neiCc,

        ownPatch,
        neiPatch
    );

    // Analyse regions. Reuse regionsplit
    boolList blockedFace(mesh_.nFaces(), false);

    forAll(ownPatch, facei)
    {
        if (ownPatch[facei] != -1 || neiPatch[facei] != -1)
        {
            blockedFace[facei] = true;
        }
    }
    syncTools::syncFaceList(mesh_, blockedFace, orEqOp<bool>());

    // Set region per cell based on walking
    regionSplit cellRegion(mesh_, blockedFace);
    blockedFace.clear();

    // Find the region containing the keepPoint
    const label keepRegionI = findRegion
    (
        mesh_,
        cellRegion,
        mergeDistance_*vector(1,1,1),
        keepPoint
    );

    Info<< "Found point " << keepPoint
        << " in global region " << keepRegionI
        << " out of " << cellRegion.nRegions() << " regions." << endl;

    if (keepRegionI == -1)
    {
        FatalErrorInFunction
            << "Point " << keepPoint
            << " is not inside the mesh." << nl
            << "Bounding box of the mesh:" << mesh_.bounds()
            << exit(FatalError);
    }


    // Walk out nBufferlayers from region split
    // (modifies cellRegion, ownPatch)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Takes over face patch onto points and then back to faces and cells
    // (so cell-face-point walk)

    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    // Patch for exposed faces for lack of anything sensible.
    label defaultPatch = 0;
    if (globalToMasterPatch.size())
    {
        defaultPatch = globalToMasterPatch[0];
    }

    for (label i = 0; i < nBufferLayers; i++)
    {
        // 1. From cells (via faces) to points

        labelList pointBaffle(mesh_.nPoints(), -1);

        forAll(faceNeighbour, facei)
        {
            const face& f = mesh_.faces()[facei];

            label ownRegion = cellRegion[faceOwner[facei]];
            label neiRegion = cellRegion[faceNeighbour[facei]];

            if (ownRegion == keepRegionI && neiRegion != keepRegionI)
            {
                // Note max(..) since possibly regionSplit might have split
                // off extra unreachable parts of mesh. Note: or can this only
                // happen for boundary faces?
                forAll(f, fp)
                {
                    pointBaffle[f[fp]] = max(defaultPatch, ownPatch[facei]);
                }
            }
            else if (ownRegion != keepRegionI && neiRegion == keepRegionI)
            {
                label newPatchi = neiPatch[facei];
                if (newPatchi == -1)
                {
                    newPatchi = max(defaultPatch, ownPatch[facei]);
                }
                forAll(f, fp)
                {
                    pointBaffle[f[fp]] = newPatchi;
                }
            }
        }
        for
        (
            label facei = mesh_.nInternalFaces();
            facei < mesh_.nFaces();
            facei++
        )
        {
            const face& f = mesh_.faces()[facei];

            label ownRegion = cellRegion[faceOwner[facei]];

            if (ownRegion == keepRegionI)
            {
                forAll(f, fp)
                {
                    pointBaffle[f[fp]] = max(defaultPatch, ownPatch[facei]);
                }
            }
        }

        // Sync
        syncTools::syncPointList
        (
            mesh_,
            pointBaffle,
            maxEqOp<label>(),
            label(-1)           // null value
        );


        // 2. From points back to faces

        const labelListList& pointFaces = mesh_.pointFaces();

        forAll(pointFaces, pointi)
        {
            if (pointBaffle[pointi] != -1)
            {
                const labelList& pFaces = pointFaces[pointi];

                forAll(pFaces, pFacei)
                {
                    label facei = pFaces[pFacei];

                    if (ownPatch[facei] == -1)
                    {
                        ownPatch[facei] = pointBaffle[pointi];
                    }
                }
            }
        }
        syncTools::syncFaceList(mesh_, ownPatch, maxEqOp<label>());


        // 3. From faces to cells (cellRegion) and back to faces (ownPatch)

        labelList newOwnPatch(ownPatch);

        forAll(ownPatch, facei)
        {
            if (ownPatch[facei] != -1)
            {
                label own = faceOwner[facei];

                if (cellRegion[own] != keepRegionI)
                {
                    cellRegion[own] = keepRegionI;

                    const cell& ownFaces = mesh_.cells()[own];
                    forAll(ownFaces, j)
                    {
                        if (ownPatch[ownFaces[j]] == -1)
                        {
                            newOwnPatch[ownFaces[j]] = ownPatch[facei];
                        }
                    }
                }
                if (mesh_.isInternalFace(facei))
                {
                    label nei = faceNeighbour[facei];

                    if (cellRegion[nei] != keepRegionI)
                    {
                        cellRegion[nei] = keepRegionI;

                        const cell& neiFaces = mesh_.cells()[nei];
                        forAll(neiFaces, j)
                        {
                            if (ownPatch[neiFaces[j]] == -1)
                            {
                                newOwnPatch[neiFaces[j]] = ownPatch[facei];
                            }
                        }
                    }
                }
            }
        }

        ownPatch.transfer(newOwnPatch);

        syncTools::syncFaceList(mesh_, ownPatch, maxEqOp<label>());
    }



    // Subset
    // ~~~~~~

    // Get cells to remove
    DynamicList<label> cellsToRemove(mesh_.nCells());
    forAll(cellRegion, celli)
    {
        if (cellRegion[celli] != keepRegionI)
        {
            cellsToRemove.append(celli);
        }
    }
    cellsToRemove.shrink();

    label nCellsToKeep = mesh_.nCells() - cellsToRemove.size();
    reduce(nCellsToKeep, sumOp<label>());

    Info<< "Keeping all cells in region " << keepRegionI
        << " containing point " << keepPoint << endl
        << "Selected for keeping : " << nCellsToKeep
        << " cells." << endl;


    // Remove cells
    removeCells cellRemover(mesh_);

    // Pick up patches for exposed faces
    labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
    labelList exposedPatches(exposedFaces.size());

    forAll(exposedFaces, i)
    {
        label facei = exposedFaces[i];

        if (ownPatch[facei] != -1)
        {
            exposedPatches[i] = ownPatch[facei];
        }
        else
        {
            WarningInFunction
                << "For exposed face " << facei
                << " fc:" << mesh_.faceCentres()[facei]
                << " found no patch." << endl
                << "    Taking patch " << defaultPatch
                << " instead." << endl;
            exposedPatches[i] = defaultPatch;
        }
    }

    return doRemoveCells
    (
        cellsToRemove,
        exposedFaces,
        exposedPatches,
        cellRemover
    );
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::dupNonManifoldPoints
(
    const localPointRegion& regionSide
)
{
    // Topochange container
    polyTopoChange meshMod(mesh_);

    label nNonManifPoints = returnReduce
    (
        regionSide.meshPointMap().size(),
        sumOp<label>()
    );

    Info<< "dupNonManifoldPoints : Found : " << nNonManifPoints
        << " non-manifold points (out of "
        << mesh_.globalData().nTotalPoints()
        << ')' << endl;

    // Topo change engine
    duplicatePoints pointDuplicator(mesh_);

    // Insert changes into meshMod
    pointDuplicator.setRefinement(regionSide, meshMod);

    // Change the mesh (no inflation, parallel sync)
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

    // Update fields
    mesh_.updateMesh(map);

    // Move mesh if in inflation mode
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes.
        mesh_.clearOut();
    }

    // Reset the instance for if in overwrite mode
    mesh_.setInstance(timeName());

    // Update intersections. Is mapping only (no faces created, positions stay
    // same) so no need to recalculate intersections.
    updateMesh(map, labelList(0));

    return map;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::dupNonManifoldPoints()
{
    // Analyse which points need to be duplicated
    localPointRegion regionSide(mesh_);

    return dupNonManifoldPoints(regionSide);
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::zonify
(
    const point& keepPoint,
    const bool allowFreeStandingZoneFaces
)
{
    const PtrList<surfaceZonesInfo>& surfZones = surfaces_.surfZones();

    labelList namedSurfaces(surfaceZonesInfo::getNamedSurfaces(surfZones));

    forAll(namedSurfaces, i)
    {
        label surfI = namedSurfaces[i];

        Info<< "Surface : " << surfaces_.names()[surfI] << nl
            << "    faceZone : " << surfZones[surfI].faceZoneName() << nl
            << "    cellZone : " << surfZones[surfI].cellZoneName() << endl;
    }


    // Add zones to mesh
    labelList surfaceToFaceZone =
        surfaceZonesInfo::addFaceZonesToMesh
        (
            surfZones,
            namedSurfaces,
            mesh_
        );

    labelList surfaceToCellZone =
        surfaceZonesInfo::addCellZonesToMesh
        (
            surfZones,
            namedSurfaces,
            mesh_
        );


    const pointField& cellCentres = mesh_.cellCentres();
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();


    // Swap neighbouring cell centres and cell level
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
    calcNeighbourData(neiLevel, neiCc);


    // Mark faces intersecting zoned surfaces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    // Like surfaceIndex_ but only for named surfaces.
    labelList namedSurfaceIndex(mesh_.nFaces(), -1);
    PackedBoolList posOrientation(mesh_.nFaces());

    {
        // Statistics: number of faces per faceZone
        labelList nSurfFaces(surfZones.size(), 0);

        // Note: for all internal faces? internal + coupled?
        // Since zonify is run after baffling the surfaceIndex_ on baffles is
        // not synchronised across both baffle faces. Fortunately we don't
        // do zonify baffle faces anyway (they are normal boundary faces).

        // Collect candidate faces
        // ~~~~~~~~~~~~~~~~~~~~~~~

        labelList testFaces(intersectedFaces());

        // Collect segments
        // ~~~~~~~~~~~~~~~~

        pointField start(testFaces.size());
        pointField end(testFaces.size());

        forAll(testFaces, i)
        {
            label facei = testFaces[i];

            if (mesh_.isInternalFace(facei))
            {
                start[i] = cellCentres[faceOwner[facei]];
                end[i] = cellCentres[faceNeighbour[facei]];
            }
            else
            {
                start[i] = cellCentres[faceOwner[facei]];
                end[i] = neiCc[facei-mesh_.nInternalFaces()];
            }
        }

        // Extend segments a bit
        {
            const vectorField smallVec(rootSmall*(end-start));
            start -= smallVec;
            end += smallVec;
        }


        // Do test for intersections
        // ~~~~~~~~~~~~~~~~~~~~~~~~~
        // Note that we intersect all intersected faces again. Could reuse
        // the information already in surfaceIndex_.

        labelList surface1;
        List<pointIndexHit> hit1;
        vectorField normal1;
        labelList surface2;
        List<pointIndexHit> hit2;
        vectorField normal2;
        {
            labelList region1;
            labelList region2;
            surfaces_.findNearestIntersection
            (
                namedSurfaces,
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
        }

        forAll(testFaces, i)
        {
            label facei = testFaces[i];
            const vector& area = mesh_.faceAreas()[facei];

            if (surface1[i] != -1)
            {
                // If both hit should probably choose 'nearest'
                if
                (
                    surface2[i] != -1
                 && (
                        magSqr(hit2[i].hitPoint())
                      < magSqr(hit1[i].hitPoint())
                    )
                )
                {
                    namedSurfaceIndex[facei] = surface2[i];
                    posOrientation[facei] = ((area&normal2[i]) > 0);
                    nSurfFaces[surface2[i]]++;
                }
                else
                {
                    namedSurfaceIndex[facei] = surface1[i];
                    posOrientation[facei] = ((area&normal1[i]) > 0);
                    nSurfFaces[surface1[i]]++;
                }
            }
            else if (surface2[i] != -1)
            {
                namedSurfaceIndex[facei] = surface2[i];
                posOrientation[facei] = ((area&normal2[i]) > 0);
                nSurfFaces[surface2[i]]++;
            }
        }


        // surfaceIndex might have different surfaces on both sides if
        // there happen to be a (obviously thin) surface with different
        // regions between the cell centres. If one is on a named surface
        // and the other is not this might give problems so sync.
        syncTools::syncFaceList
        (
            mesh_,
            namedSurfaceIndex,
            maxEqOp<label>()
        );

        // Print a bit
        if (debug)
        {
            forAll(nSurfFaces, surfI)
            {
                Pout<< "Surface:"
                    << surfaces_.names()[surfI]
                    << "  nZoneFaces:" << nSurfFaces[surfI] << nl;
            }
            Pout<< endl;
        }
    }


    // Put the cells into the correct zone
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Zone per cell:
    // -2 : unset
    // -1 : not in any zone
    // >=0: zoneID
    labelList cellToZone(mesh_.nCells(), -2);


    // Set using geometric test
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    // Closed surfaces with cellZone specified.
    labelList closedNamedSurfaces
    (
        surfaceZonesInfo::getClosedNamedSurfaces
        (
            surfZones,
            surfaces_.geometry(),
            surfaces_.surfaces()
        )
    );

    if (closedNamedSurfaces.size())
    {
        Info<< "Found " << closedNamedSurfaces.size()
            << " closed, named surfaces. Assigning cells in/outside"
            << " these surfaces to the corresponding cellZone."
            << nl << endl;

        findCellZoneGeometric
        (
            neiCc,
            closedNamedSurfaces,    // indices of closed surfaces
            namedSurfaceIndex,      // per face index of named surface
            surfaceToCellZone,      // cell zone index per surface

            cellToZone
        );
    }


    // Set using provided locations
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList locationSurfaces
    (
        surfaceZonesInfo::getInsidePointNamedSurfaces(surfZones)
    );

    if (locationSurfaces.size())
    {
        Info<< "Found " << locationSurfaces.size()
            << " named surfaces with a provided inside point."
            << " Assigning cells inside these surfaces"
            << " to the corresponding cellZone."
            << nl << endl;

        findCellZoneInsideWalk
        (
            locationSurfaces,       // indices of closed surfaces
            namedSurfaceIndex,      // per face index of named surface
            surfaceToCellZone,      // cell zone index per surface

            cellToZone
        );
    }


    // Set using walking
    // ~~~~~~~~~~~~~~~~~

    {
        Info<< "Walking from location-in-mesh " << keepPoint
            << " to assign cellZones "
            << "- crossing a faceZone face changes cellZone" << nl << endl;

        // Topological walk
        findCellZoneTopo
        (
            keepPoint,
            namedSurfaceIndex,
            surfaceToCellZone,

            cellToZone
        );
    }


    // Make sure namedSurfaceIndex is unset in between same cell cell zones.
    if (!allowFreeStandingZoneFaces)
    {
        Info<< "Only keeping zone faces in between different cellZones."
            << nl << endl;

        makeConsistentFaceIndex(cellToZone, namedSurfaceIndex);
    }


    //- Per face index of faceZone or -1
    labelList faceToZone(mesh_.nFaces(), -1);

    // Convert namedSurfaceIndex (index of named surfaces) to
    // actual faceZone index

    forAll(namedSurfaceIndex, facei)
    {
        label surfI = namedSurfaceIndex[facei];
        if (surfI != -1)
        {
            faceToZone[facei] = surfaceToFaceZone[surfI];
        }
    }


    // Topochange container
    polyTopoChange meshMod(mesh_);



    // Get coupled neighbour cellZone. Set to -1 on non-coupled patches.
    labelList neiCellZone;
    syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (!pp.coupled())
        {
            label bFacei = pp.start()-mesh_.nInternalFaces();
            forAll(pp, i)
            {
                neiCellZone[bFacei++] = -1;
            }
        }
    }



    // Get per face whether is it master (of a coupled set of faces)
    const PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_));



    // faceZones
    // ~~~~~~~~~
    // Faces on faceZones come in two variants:
    // - faces on the outside of a cellZone. They will be oriented to
    //   point out of the maximum cellZone.
    // - free-standing faces. These will be oriented according to the
    //   local surface normal. We do this in a two step algorithm:
    //      - do a consistent orientation
    //      - check number of faces with consistent orientation
    //      - if <0 flip the whole patch
    boolList meshFlipMap(mesh_.nFaces(), false);
    {
        // Collect all data on zone faces without cellZones on either side.
        const indirectPrimitivePatch patch
        (
            IndirectList<face>
            (
                mesh_.faces(),
                freeStandingBaffleFaces
                (
                    faceToZone,
                    cellToZone,
                    neiCellZone
                )
            ),
            mesh_.points()
        );

        label nFreeStanding = returnReduce(patch.size(), sumOp<label>());
        if (nFreeStanding > 0)
        {
            Info<< "Detected " << nFreeStanding << " free-standing zone faces"
                << endl;

            if (debug)
            {
                OBJstream str(mesh_.time().path()/"freeStanding.obj");
                str.write(patch.localFaces(), patch.localPoints(), false);
            }


            // Detect non-manifold edges
            labelList nMasterFacesPerEdge;
            calcPatchNumMasterFaces(isMasterFace, patch, nMasterFacesPerEdge);

            // Mark zones. Even a single original surface might create multiple
            // disconnected/non-manifold-connected zones
            labelList faceToConnectedZone;
            const label nZones = markPatchZones
            (
                patch,
                nMasterFacesPerEdge,
                faceToConnectedZone
            );

            Map<label> nPosOrientation(2*nZones);
            for (label zoneI = 0; zoneI < nZones; zoneI++)
            {
                nPosOrientation.insert(zoneI, 0);
            }

            // Make orientations consistent in a topological way. This just
            // checks  the first face per zone for whether nPosOrientation
            // is negative (which it never is at this point)
            consistentOrientation
            (
                isMasterFace,
                patch,
                nMasterFacesPerEdge,
                faceToConnectedZone,
                nPosOrientation,

                meshFlipMap
            );

            // Count per region the number of orientations (taking the new
            // flipMap into account)
            forAll(patch.addressing(), facei)
            {
                label meshFacei = patch.addressing()[facei];

                if (isMasterFace[meshFacei])
                {
                    label n = 1;
                    if
                    (
                        bool(posOrientation[meshFacei])
                     == meshFlipMap[meshFacei]
                    )
                    {
                        n = -1;
                    }

                    nPosOrientation.find(faceToConnectedZone[facei])() += n;
                }
            }
            Pstream::mapCombineGather(nPosOrientation, plusEqOp<label>());
            Pstream::mapCombineScatter(nPosOrientation);


            Info<< "Split " << nFreeStanding << " free-standing zone faces"
                << " into " << nZones << " disconnected regions with size"
                << " (negative denotes wrong orientation) :"
                << endl;

            for (label zoneI = 0; zoneI < nZones; zoneI++)
            {
                Info<< "    " << zoneI << "\t" << nPosOrientation[zoneI]
                    << endl;
            }
            Info<< endl;


            // Reapply with new counts (in nPosOrientation). This will cause
            // zones with a negative count to be flipped.
            consistentOrientation
            (
                isMasterFace,
                patch,
                nMasterFacesPerEdge,
                faceToConnectedZone,
                nPosOrientation,

                meshFlipMap
            );
        }
    }


    // Put the faces into the correct zone
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        label faceZoneI = faceToZone[facei];

        if (faceZoneI != -1)
        {
            // Orient face zone to have slave cells in max cell zone.
            // Note: logic to use flipMap should be consistent with logic
            //       to pick up the freeStandingBaffleFaces!

            label ownZone = cellToZone[faceOwner[facei]];
            label neiZone = cellToZone[faceNeighbour[facei]];

            bool flip;

            if (ownZone == neiZone)
            {
                // free-standing face. Use geometrically derived orientation
                flip = meshFlipMap[facei];
            }
            else
            {
                flip =
                (
                    ownZone == -1
                 || (neiZone != -1 && ownZone > neiZone)
                );
            }

            meshMod.setAction
            (
                polyModifyFace
                (
                    mesh_.faces()[facei],           // modified face
                    facei,                          // label of face
                    faceOwner[facei],               // owner
                    faceNeighbour[facei],           // neighbour
                    false,                          // face flip
                    -1,                             // patch for face
                    false,                          // remove from zone
                    faceZoneI,                      // zone for face
                    flip                            // face flip in zone
                )
            );
        }
    }


    // Set owner as no-flip
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        label facei = pp.start();

        forAll(pp, i)
        {
            label faceZoneI = faceToZone[facei];

            if (faceZoneI != -1)
            {
                label ownZone = cellToZone[faceOwner[facei]];
                label neiZone = neiCellZone[facei-mesh_.nInternalFaces()];

                bool flip;

                if (ownZone == neiZone)
                {
                    // free-standing face. Use geometrically derived orientation
                    flip = meshFlipMap[facei];
                }
                else
                {
                    flip =
                    (
                        ownZone == -1
                     || (neiZone != -1 && ownZone > neiZone)
                    );
                }

                meshMod.setAction
                (
                    polyModifyFace
                    (
                        mesh_.faces()[facei],           // modified face
                        facei,                          // label of face
                        faceOwner[facei],               // owner
                        -1,                             // neighbour
                        false,                          // face flip
                        patchi,                         // patch for face
                        false,                          // remove from zone
                        faceZoneI,                      // zone for face
                        flip                            // face flip in zone
                    )
                );
            }
            facei++;
        }
    }


    // Put the cells into the correct zone
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(cellToZone, celli)
    {
        label zoneI = cellToZone[celli];

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

    // Change the mesh (no inflation, parallel sync)
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

    // Update fields
    mesh_.updateMesh(map);

    // Move mesh if in inflation mode
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes.
        mesh_.clearOut();
    }

    // Reset the instance for if in overwrite mode
    mesh_.setInstance(timeName());

    // Print some stats (note: zones are synchronised)
    if (mesh_.cellZones().size() > 0)
    {
        Info<< "CellZones:" << endl;
        forAll(mesh_.cellZones(), zoneI)
        {
            const cellZone& cz = mesh_.cellZones()[zoneI];
            Info<< "    " << cz.name()
                << "\tsize:" << returnReduce(cz.size(), sumOp<label>())
                << endl;
        }
        Info<< endl;
    }
    if (mesh_.faceZones().size() > 0)
    {
        Info<< "FaceZones:" << endl;
        forAll(mesh_.faceZones(), zoneI)
        {
            const faceZone& fz = mesh_.faceZones()[zoneI];
            Info<< "    " << fz.name()
                << "\tsize:" << returnReduce(fz.size(), sumOp<label>())
                << endl;
        }
        Info<< endl;
    }

    // None of the faces has changed, only the zones. Still...
    updateMesh(map, labelList());

    return map;
}


// ************************************************************************* //
