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

// Repatches external face or creates baffle for internal face
// with user specified patches (might be different for both sides).
// Returns label of added face.
Foam::label Foam::meshRefinement::createBaffle
(
    const label faceI,
    const label ownPatch,
    const label neiPatch,
    polyTopoChange& meshMod
) const
{
    const face& f = mesh_.faces()[faceI];
    label zoneID = mesh_.faceZones().whichZone(faceI);
    bool zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh_.faceZones()[zoneID];
        zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
    }

    meshMod.setAction
    (
        polyModifyFace
        (
            f,                          // modified face
            faceI,                      // label of face
            mesh_.faceOwner()[faceI],   // owner
            -1,                         // neighbour
            false,                      // face flip
            ownPatch,                   // patch for face
            false,                      // remove from zone
            zoneID,                     // zone for face
            zoneFlip                    // face flip in zone
        )
    );


    label dupFaceI = -1;

    if (mesh_.isInternalFace(faceI))
    {
        if (neiPatch == -1)
        {
            FatalErrorIn
            (
                "meshRefinement::createBaffle"
                "(const label, const label, const label, polyTopoChange&)"
            )   << "No neighbour patch for internal face " << faceI
                << " fc:" << mesh_.faceCentres()[faceI]
                << " ownPatch:" << ownPatch << abort(FatalError);
        }

        bool reverseFlip = false;
        if (zoneID >= 0)
        {
            reverseFlip = !zoneFlip;
        }

        dupFaceI = meshMod.setAction
        (
            polyAddFace
            (
                f.reverseFace(),            // modified face
                mesh_.faceNeighbour()[faceI],// owner
                -1,                         // neighbour
                -1,                         // masterPointID
                -1,                         // masterEdgeID
                faceI,                      // masterFaceID,
                true,                       // face flip
                neiPatch,                   // patch for face
                zoneID,                     // zone for face
                reverseFlip                 // face flip in zone
            )
        );
    }
    return dupFaceI;
}


//// Check if we are a boundary face and normal of surface does
//// not align with test vector. In this case there'd probably be
//// a freestanding 'baffle' so we might as well not create it.
//// Note that since it is not a proper baffle we cannot detect it
//// afterwards so this code cannot be merged with the
//// filterDuplicateFaces code.
//bool Foam::meshRefinement::validBaffleTopology
//(
//    const label faceI,
//    const vector& n1,
//    const vector& n2,
//    const vector& testDir
//) const
//{
//
//    label patchI = mesh_.boundaryMesh().whichPatch(faceI);
//    if (patchI == -1 || mesh_.boundaryMesh()[patchI].coupled())
//    {
//        return true;
//    }
//    else if (mag(n1&n2) > cos(degToRad(30)))
//    {
//        // Both normals aligned. Check that test vector perpendicularish to
//        // surface normal
//        scalar magTestDir = mag(testDir);
//        if (magTestDir > VSMALL)
//        {
//            if (mag(n1&(testDir/magTestDir)) < cos(degToRad(45)))
//            {
//                //Pout<< "** disabling baffling face "
//                //    << mesh_.faceCentres()[faceI] << endl;
//                return false;
//            }
//        }
//    }
//    return true;
//}


// Determine patches for baffles on all intersected unnamed faces
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
        label faceI = testFaces[i];

        label own = mesh_.faceOwner()[faceI];

        if (mesh_.isInternalFace(faceI))
        {
            start[i] = cellCentres[own];
            end[i] = cellCentres[mesh_.faceNeighbour()[faceI]];
        }
        else
        {
            start[i] = cellCentres[own];
            end[i] = neiCc[faceI-mesh_.nInternalFaces()];
        }
    }

    // Extend segments a bit
    {
        const vectorField smallVec(Foam::sqrt(SMALL)*(end-start));
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
        label faceI = testFaces[i];

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
            ownPatch[faceI] = globalToMasterPatch
            [
                surfaces_.globalRegion(surface1[i], region1[i])
            ];
            neiPatch[faceI] = globalToMasterPatch
            [
                surfaces_.globalRegion(surface2[i], region2[i])
            ];

            if (ownPatch[faceI] == -1 || neiPatch[faceI] == -1)
            {
                FatalErrorIn("getBafflePatches(..)")
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
                label faceI = fZone[i];

                if (allowBoundary || mesh_.isInternalFace(faceI))
                {
                    labelPair patches = zPatches;
                    if (fZone.flipMap()[i])
                    {
                       patches = reverse(patches);
                    }

                    if (!bafflePatch.insert(faceI, patches))
                    {
                        FatalErrorIn("getZoneBafflePatches(..)")
                            << "Face " << faceI
                            << " fc:" << mesh_.faceCentres()[faceI]
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


// Create baffle for every face where ownPatch != -1
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
        FatalErrorIn
        (
            "meshRefinement::createBaffles"
            "(const labelList&, const labelList&)"
        )   << "Illegal size :"
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

        forAll(syncedOwnPatch, faceI)
        {
            if
            (
                (ownPatch[faceI] == -1 && syncedOwnPatch[faceI] != -1)
             || (neiPatch[faceI] == -1 && syncedNeiPatch[faceI] != -1)
            )
            {
                FatalErrorIn
                (
                    "meshRefinement::createBaffles"
                    "(const labelList&, const labelList&)"
                )   << "Non synchronised at face:" << faceI
                    << " on patch:" << mesh_.boundaryMesh().whichPatch(faceI)
                    << " fc:" << mesh_.faceCentres()[faceI] << endl
                    << "ownPatch:" << ownPatch[faceI]
                    << " syncedOwnPatch:" << syncedOwnPatch[faceI]
                    << " neiPatch:" << neiPatch[faceI]
                    << " syncedNeiPatch:" << syncedNeiPatch[faceI]
                    << abort(FatalError);
            }
        }
    }

    // Topochange container
    polyTopoChange meshMod(mesh_);

    label nBaffles = 0;

    forAll(ownPatch, faceI)
    {
        if (ownPatch[faceI] != -1)
        {
            // Create baffle or repatch face. Return label of inserted baffle
            // face.
            createBaffle
            (
                faceI,
                ownPatch[faceI],   // owner side patch
                neiPatch[faceI],   // neighbour side patch
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
    forAll(ownPatch, oldFaceI)
    {
        label faceI = reverseFaceMap[oldFaceI];

        if (ownPatch[oldFaceI] != -1 && faceI >= 0)
        {
            const cell& ownFaces = mesh_.cells()[mesh_.faceOwner()[faceI]];

            forAll(ownFaces, i)
            {
                baffledFacesSet.insert(ownFaces[i]);
            }
        }
    }
    // Pick up neighbour side of baffle (added faces)
    forAll(faceMap, faceI)
    {
        label oldFaceI = faceMap[faceI];

        if (oldFaceI >= 0 && reverseFaceMap[oldFaceI] != faceI)
        {
            const cell& ownFaces = mesh_.cells()[mesh_.faceOwner()[faceI]];

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

    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];

        if (isA<processorPolyPatch>(pp))
        {
            forAll(pp, i)
            {
                label faceI = pp.start()+i;
                label zoneI = fZones.whichZone(faceI);

                if (zoneI != -1)
                {
                    FatalErrorIn("meshRefinement::checkZoneFaces")
                        << "face:" << faceI << " on patch " << pp.name()
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

            forAll(faceMap, faceI)
            {
                label oldFaceI = faceMap[faceI];

                // Does face originate from face-to-patch
                Map<labelPair>::const_iterator iter = faceToPatch.find
                (
                    oldFaceI
                );

                if (iter != faceToPatch.end())
                {
                    label masterFaceI = reverseFaceMap[oldFaceI];
                    if (faceI != masterFaceI)
                    {
                        baffles[baffleI++] = labelPair(masterFaceI, faceI);
                    }
                }
            }

            if (baffleI != faceToPatch.size())
            {
                FatalErrorIn("meshRefinement::createZoneBaffles(..)")
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


// Extract those baffles (duplicate) faces that are on the edge of a baffle
// region. These are candidates for merging.
// Done by counting the number of baffles faces per mesh edge. If edge
// has 2 boundary faces and both are baffle faces it is the edge of a baffle
// region.
Foam::List<Foam::labelPair> Foam::meshRefinement::freeStandingBaffles
(
    const List<labelPair>& couples,
    const scalar planarAngle
) const
{
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

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        // Count number of boundary faces. Discard coupled boundary faces.
        if (!pp.coupled())
        {
            label faceI = pp.start();

            forAll(pp, i)
            {
                const labelList& fEdges = mesh_.faceEdges(faceI);

                forAll(fEdges, fEdgeI)
                {
                    nBafflesPerEdge[fEdges[fEdgeI]]++;
                }
                faceI++;
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
            const vectorField smallVec(Foam::sqrt(SMALL)*(end-start));
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
                    if (magN > VSMALL)
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


// Merge baffles. Gets pairs of faces.
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

    forAll(insideSurfaces, cellI)
    {
        if (cellToZone[cellI] == -2)
        {
            label surfI = insideSurfaces[cellI];

            if (surfI != -1)
            {
                cellToZone[cellI] = surfaceToCellZone[surfI];
            }
        }
    }


    // Some cells with cell centres close to surface might have
    // had been put into wrong surface. Recheck with perturbed cell centre.


    // 1. Collect points

    // Count points to test.
    label nCandidates = 0;
    forAll(namedSurfaceIndex, faceI)
    {
        label surfI = namedSurfaceIndex[faceI];

        if (surfI != -1)
        {
            if (mesh_.isInternalFace(faceI))
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
    forAll(namedSurfaceIndex, faceI)
    {
        label surfI = namedSurfaceIndex[faceI];

        if (surfI != -1)
        {
            label own = faceOwner[faceI];
            const point& ownCc = cellCentres[own];

            if (mesh_.isInternalFace(faceI))
            {
                label nei = faceNeighbour[faceI];
                const point& neiCc = cellCentres[nei];
                // Perturbed cc
                const vector d = 1e-4*(neiCc - ownCc);
                candidatePoints[nCandidates++] = ownCc-d;
                candidatePoints[nCandidates++] = neiCc+d;
            }
            else
            {
                //const point& neiFc = mesh_.faceCentres()[faceI];
                const point& neiFc = neiCc[faceI-mesh_.nInternalFaces()];

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
    forAll(namedSurfaceIndex, faceI)
    {
        label surfI = namedSurfaceIndex[faceI];

        if (surfI != -1)
        {
            label own = faceOwner[faceI];

            if (mesh_.isInternalFace(faceI))
            {
                label ownSurfI = insideSurfaces[nCandidates++];
                if (ownSurfI != -1)
                {
                    cellToZone[own] = surfaceToCellZone[ownSurfI];
                }

                label neiSurfI = insideSurfaces[nCandidates++];
                if (neiSurfI != -1)
                {
                    label nei = faceNeighbour[faceI];

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

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label ownZone = cellToZone[mesh_.faceOwner()[faceI]];
        label neiZone = cellToZone[mesh_.faceNeighbour()[faceI]];

        if (namedSurfaceIndex[faceI] == -1 && (ownZone != neiZone))
        {
            // Give face the zone of max cell zone
            namedSurfaceIndex[faceI] = findIndex
            (
                surfaceToCellZone,
                max(ownZone, neiZone)
            );
        }
    }

    labelList neiCellZone(mesh_.nFaces()-mesh_.nInternalFaces());
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label faceI = pp.start()+i;
                label ownZone = cellToZone[mesh_.faceOwner()[faceI]];
                neiCellZone[faceI-mesh_.nInternalFaces()] = ownZone;
            }
        }
    }
    syncTools::swapBoundaryFaceList(mesh_, neiCellZone);

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label faceI = pp.start()+i;
                label ownZone = cellToZone[mesh_.faceOwner()[faceI]];
                label neiZone = neiCellZone[faceI-mesh_.nInternalFaces()];

                if (namedSurfaceIndex[faceI] == -1 && (ownZone != neiZone))
                {
                    // Give face the max cell zone
                    namedSurfaceIndex[faceI] = findIndex
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

    forAll(namedSurfaceIndex, faceI)
    {
        if (namedSurfaceIndex[faceI] == -1)
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
            FatalErrorIn
            (
                "meshRefinement::findCellZoneInsideWalk"
                "(const labelList&, const labelList&"
                ", const labelList&, const labelList&)"
            )   << "Point " << insidePoint
                << " is not inside the mesh." << nl
                << "Bounding box of the mesh:" << mesh_.bounds()
                << exit(FatalError);
        }

        // Set all cells with this region
        forAll(cellRegion, cellI)
        {
            if (cellRegion[cellI] == keepRegionI)
            {
                if (cellToZone[cellI] == -2)
                {
                    cellToZone[cellI] = surfaceToCellZone[surfI];
                }
                else if (cellToZone[cellI] != surfaceToCellZone[surfI])
                {
                    WarningIn
                    (
                        "meshRefinement::findCellZoneInsideWalk"
                        "(const labelList&, const labelList&"
                        ", const labelList&, const labelList&)"
                    )   << "Cell " << cellI
                        << " at " << mesh_.cellCentres()[cellI]
                        << " is inside surface " << surfaces_.names()[surfI]
                        << " but already marked as being in zone "
                        << cellToZone[cellI] << endl
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

    // Check whether inbetween different regions
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


// Finds region per cell. Assumes:
// - region containing keepPoint does not go into a cellZone
// - all other regions can be found by crossing faces marked in
//   namedSurfaceIndex.
void Foam::meshRefinement::findCellZoneTopo
(
    const point& keepPoint,
    const labelList& namedSurfaceIndex,
    const labelList& surfaceToCellZone,
    labelList& cellToZone
) const
{
    // Analyse regions. Reuse regionsplit
    boolList blockedFace(mesh_.nFaces());

    forAll(namedSurfaceIndex, faceI)
    {
        if (namedSurfaceIndex[faceI] == -1)
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
    forAll(cellToZone, cellI)
    {
        if (cellToZone[cellI] != -2)
        {
            regionToCellZone[cellRegion[cellI]] = cellToZone[cellI];
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
        FatalErrorIn
        (
            "meshRefinement::findCellZoneTopo"
            "(const point&, const labelList&, const labelList&, labelList&)"
        )   << "Point " << keepPoint
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

        for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
        {
            label surfI = namedSurfaceIndex[faceI];

            // Connected even if no cellZone defined for surface
            if (surfI != -1)
            {
                // Calculate region to zone from cellRegions on either side
                // of internal face.
                bool changedCell = calcRegionToZone
                (
                    surfaceToCellZone[surfI],
                    cellRegion[mesh_.faceOwner()[faceI]],
                    cellRegion[mesh_.faceNeighbour()[faceI]],
                    regionToCellZone
                );

                changed = changed | changedCell;
            }
        }

        // Boundary faces

        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        // Get coupled neighbour cellRegion
        labelList neiCellRegion(mesh_.nFaces()-mesh_.nInternalFaces());
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (pp.coupled())
            {
                forAll(pp, i)
                {
                    label faceI = pp.start()+i;
                    neiCellRegion[faceI-mesh_.nInternalFaces()] =
                        cellRegion[mesh_.faceOwner()[faceI]];
                }
            }
        }
        syncTools::swapBoundaryFaceList(mesh_, neiCellRegion);

        // Calculate region to zone from cellRegions on either side of coupled
        // face.
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (pp.coupled())
            {
                forAll(pp, i)
                {
                    label faceI = pp.start()+i;

                    label surfI = namedSurfaceIndex[faceI];

                    // Connected even if no cellZone defined for surface
                    if (surfI != -1)
                    {
                        bool changedCell = calcRegionToZone
                        (
                            surfaceToCellZone[surfI],
                            cellRegion[mesh_.faceOwner()[faceI]],
                            neiCellRegion[faceI-mesh_.nInternalFaces()],
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
            FatalErrorIn
            (
                "meshRefinement::findCellZoneTopo"
                "(const point&, const labelList&, const labelList&, labelList&)"
            )   << "For region " << regionI << " haven't set cell zone."
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
    forAll(cellToZone, cellI)
    {
        cellToZone[cellI] = regionToCellZone[cellRegion[cellI]];
    }
}


// Make namedSurfaceIndex consistent with cellToZone - clear out any blocked
// faces inbetween same cell zone.
void Foam::meshRefinement::makeConsistentFaceIndex
(
    const labelList& cellToZone,
    labelList& namedSurfaceIndex
) const
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label ownZone = cellToZone[faceOwner[faceI]];
        label neiZone = cellToZone[faceNeighbour[faceI]];

        if (ownZone == neiZone && namedSurfaceIndex[faceI] != -1)
        {
            namedSurfaceIndex[faceI] = -1;
        }
        else if (ownZone != neiZone && namedSurfaceIndex[faceI] == -1)
        {
            FatalErrorIn("meshRefinement::zonify()")
                << "Different cell zones on either side of face " << faceI
                << " at " << mesh_.faceCentres()[faceI]
                << " but face not marked with a surface."
                << abort(FatalError);
        }
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Get coupled neighbour cellZone
    labelList neiCellZone(mesh_.nFaces()-mesh_.nInternalFaces());
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label faceI = pp.start()+i;
                neiCellZone[faceI-mesh_.nInternalFaces()] =
                    cellToZone[mesh_.faceOwner()[faceI]];
            }
        }
    }
    syncTools::swapBoundaryFaceList(mesh_, neiCellZone);

    // Use coupled cellZone to do check
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label faceI = pp.start()+i;

                label ownZone = cellToZone[faceOwner[faceI]];
                label neiZone = neiCellZone[faceI-mesh_.nInternalFaces()];

                if (ownZone == neiZone && namedSurfaceIndex[faceI] != -1)
                {
                    namedSurfaceIndex[faceI] = -1;
                }
                else if (ownZone != neiZone && namedSurfaceIndex[faceI] == -1)
                {
                    FatalErrorIn("meshRefinement::zonify()")
                        << "Different cell zones on either side of face "
                        << faceI << " at " << mesh_.faceCentres()[faceI]
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
                label faceI = pp.start()+i;
                namedSurfaceIndex[faceI] = -1;
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

        forAll(facePatch, faceI)
        {
            if (facePatch[faceI] != -1)
            {
                problemFaces.insert(faceI);
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
    const labelList& namedSurfaceIndex,
    const labelList& cellToZone,
    const labelList& neiCellZone
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();


    DynamicList<label> faceLabels(mesh_.nFaces()/20);

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label surfI = namedSurfaceIndex[faceI];
        if (surfI != -1)
        {
            // Free standing baffle?
            label ownZone = cellToZone[faceOwner[faceI]];
            label neiZone = cellToZone[faceNeighbour[faceI]];
            if (max(ownZone, neiZone) == -1)
            {
                faceLabels.append(faceI);
            }
        }
    }
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        forAll(pp, i)
        {
            label faceI = pp.start()+i;
            label surfI = namedSurfaceIndex[faceI];
            if (surfI != -1)
            {
                // Free standing baffle?
                label ownZone = cellToZone[faceOwner[faceI]];
                label neiZone = neiCellZone[faceI-mesh_.nInternalFaces()];
                if (max(ownZone, neiZone) == -1)
                {
                    faceLabels.append(faceI);
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

    forAll(patch.addressing(), faceI)
    {
        const label meshFaceI = patch.addressing()[faceI];

        if (isMasterFace[meshFaceI])
        {
            const labelList& fEdges = patch.faceEdges()[faceI];
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
        0
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

    label faceI = 0;

    label currentZoneI = 0;

    while (true)
    {
        // Pick an unset face
        label globalSeed = labelMax;
        for (; faceI < allFaceInfo.size(); faceI++)
        {
            if (!allFaceInfo[faceI].valid(dummyTrackData))
            {
                globalSeed = globalFaces.toGlobal(faceI);
                break;
            }
        }

        reduce(globalSeed, minOp<label>());

        if (globalSeed == labelMax)
        {
            break;
        }

        label procI = globalFaces.whichProcID(globalSeed);
        label seedFaceI = globalFaces.toLocal(procI, globalSeed);

        //Info<< "Seeding zone " << currentZoneI
        //    << " from processor " << procI << " face " << seedFaceI
        //    << endl;

        if (procI == Pstream::myProcNo())
        {
            patchEdgeFaceRegion& faceInfo = allFaceInfo[seedFaceI];


            // Set face
            faceInfo = currentZoneI;

            // .. and seed its edges
            const labelList& fEdges = patch.faceEdges()[seedFaceI];
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
                        seedFaceI,
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
    forAll(allFaceInfo, faceI)
    {
        if (!allFaceInfo[faceI].valid(dummyTrackData))
        {
            FatalErrorIn("meshRefinement::markPatchZones(..)")
                << "Problem: unvisited face " << faceI
                << " at " << patch.faceCentres()[faceI]
                << exit(FatalError);
        }
        faceToZone[faceI] = allFaceInfo[faceI].region();
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

        forAll(patch.addressing(), faceI)
        {
            const label meshFaceI = patch.addressing()[faceI];
            const label patchI = bm.whichPatch(meshFaceI);

            if
            (
                patchI != -1
             && bm[patchI].coupled()
             && !isMasterFace[meshFaceI]
            )
            {
                // Slave side. Mark so doesn't get visited.
                allFaceInfo[faceI] = orientedSurface::NOFLIP;
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
        forAll(allFaceInfo, faceI)
        {
            if (allFaceInfo[faceI] == orientedSurface::UNVISITED)
            {
                globalSeed = globalFaces.toGlobal(faceI);
                break;
            }
        }

        reduce(globalSeed, minOp<label>());

        if (globalSeed == labelMax)
        {
            break;
        }

        label procI = globalFaces.whichProcID(globalSeed);
        label seedFaceI = globalFaces.toLocal(procI, globalSeed);

        //Info<< "Seeding from processor " << procI << " face " << seedFaceI
        //    << endl;

        if (procI == Pstream::myProcNo())
        {
            // Determine orientation of seedFace

            patchFaceOrientation& faceInfo = allFaceInfo[seedFaceI];

            // Start off with correct orientation
            faceInfo = orientedSurface::NOFLIP;

            if (zoneToOrientation[faceToZone[seedFaceI]] < 0)
            {
                faceInfo.flip();
            }


            const labelList& fEdges = patch.faceEdges()[seedFaceI];
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
                        seedFaceI,
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
            const label meshFaceI = patch.addressing()[i];
            if (!mesh_.isInternalFace(meshFaceI))
            {
                neiStatus[meshFaceI-mesh_.nInternalFaces()] =
                    allFaceInfo[i].flipStatus();
            }
        }
        syncTools::swapBoundaryFaceList(mesh_, neiStatus);

        forAll(patch.addressing(), i)
        {
            const label meshFaceI = patch.addressing()[i];
            const label patchI = bm.whichPatch(meshFaceI);

            if
            (
                patchI != -1
             && bm[patchI].coupled()
             && !isMasterFace[meshFaceI]
            )
            {
                // Slave side. Take flipped from neighbour
                label bFaceI = meshFaceI-mesh_.nInternalFaces();

                if (neiStatus[bFaceI] == orientedSurface::NOFLIP)
                {
                    allFaceInfo[i] = orientedSurface::FLIP;
                }
                else if (neiStatus[bFaceI] == orientedSurface::FLIP)
                {
                    allFaceInfo[i] = orientedSurface::NOFLIP;
                }
                else
                {
                    FatalErrorIn("meshRefinement::consistentOrientation(..)")
                        << "Incorrect status for face " << meshFaceI
                        << abort(FatalError);
                }
            }
        }
    }


    // Convert to meshFlipMap and adapt faceZones

    meshFlipMap.setSize(mesh_.nFaces());
    meshFlipMap = false;

    forAll(allFaceInfo, faceI)
    {
        label meshFaceI = patch.addressing()[faceI];

        if (allFaceInfo[faceI] == orientedSurface::NOFLIP)
        {
            meshFlipMap[meshFaceI] = false;
        }
        else if (allFaceInfo[faceI] == orientedSurface::FLIP)
        {
            meshFlipMap[meshFaceI] = true;
        }
        else
        {
            FatalErrorIn("meshRefinement::consistentOrientation(..)")
                << "Problem : unvisited face " << faceI
                << " centre:" << mesh_.faceCentres()[meshFaceI]
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Split off unreachable areas of mesh.
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


// Split off (with optional buffer layers) unreachable areas of mesh.
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

    forAll(ownPatch, faceI)
    {
        if (ownPatch[faceI] != -1 || neiPatch[faceI] != -1)
        {
            blockedFace[faceI] = true;
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
        FatalErrorIn
        (
            "meshRefinement::splitMesh"
            "(const label, const labelList&, const point&)"
        )   << "Point " << keepPoint
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

        forAll(faceNeighbour, faceI)
        {
            const face& f = mesh_.faces()[faceI];

            label ownRegion = cellRegion[faceOwner[faceI]];
            label neiRegion = cellRegion[faceNeighbour[faceI]];

            if (ownRegion == keepRegionI && neiRegion != keepRegionI)
            {
                // Note max(..) since possibly regionSplit might have split
                // off extra unreachable parts of mesh. Note: or can this only
                // happen for boundary faces?
                forAll(f, fp)
                {
                    pointBaffle[f[fp]] = max(defaultPatch, ownPatch[faceI]);
                }
            }
            else if (ownRegion != keepRegionI && neiRegion == keepRegionI)
            {
                label newPatchI = neiPatch[faceI];
                if (newPatchI == -1)
                {
                    newPatchI = max(defaultPatch, ownPatch[faceI]);
                }
                forAll(f, fp)
                {
                    pointBaffle[f[fp]] = newPatchI;
                }
            }
        }
        for
        (
            label faceI = mesh_.nInternalFaces();
            faceI < mesh_.nFaces();
            faceI++
        )
        {
            const face& f = mesh_.faces()[faceI];

            label ownRegion = cellRegion[faceOwner[faceI]];

            if (ownRegion == keepRegionI)
            {
                forAll(f, fp)
                {
                    pointBaffle[f[fp]] = max(defaultPatch, ownPatch[faceI]);
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

        forAll(pointFaces, pointI)
        {
            if (pointBaffle[pointI] != -1)
            {
                const labelList& pFaces = pointFaces[pointI];

                forAll(pFaces, pFaceI)
                {
                    label faceI = pFaces[pFaceI];

                    if (ownPatch[faceI] == -1)
                    {
                        ownPatch[faceI] = pointBaffle[pointI];
                    }
                }
            }
        }
        syncTools::syncFaceList(mesh_, ownPatch, maxEqOp<label>());


        // 3. From faces to cells (cellRegion) and back to faces (ownPatch)

        labelList newOwnPatch(ownPatch);

        forAll(ownPatch, faceI)
        {
            if (ownPatch[faceI] != -1)
            {
                label own = faceOwner[faceI];

                if (cellRegion[own] != keepRegionI)
                {
                    cellRegion[own] = keepRegionI;

                    const cell& ownFaces = mesh_.cells()[own];
                    forAll(ownFaces, j)
                    {
                        if (ownPatch[ownFaces[j]] == -1)
                        {
                            newOwnPatch[ownFaces[j]] = ownPatch[faceI];
                        }
                    }
                }
                if (mesh_.isInternalFace(faceI))
                {
                    label nei = faceNeighbour[faceI];

                    if (cellRegion[nei] != keepRegionI)
                    {
                        cellRegion[nei] = keepRegionI;

                        const cell& neiFaces = mesh_.cells()[nei];
                        forAll(neiFaces, j)
                        {
                            if (ownPatch[neiFaces[j]] == -1)
                            {
                                newOwnPatch[neiFaces[j]] = ownPatch[faceI];
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
    forAll(cellRegion, cellI)
    {
        if (cellRegion[cellI] != keepRegionI)
        {
            cellsToRemove.append(cellI);
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
        label faceI = exposedFaces[i];

        if (ownPatch[faceI] != -1)
        {
            exposedPatches[i] = ownPatch[faceI];
        }
        else
        {
            WarningIn("meshRefinement::splitMesh(..)")
                << "For exposed face " << faceI
                << " fc:" << mesh_.faceCentres()[faceI]
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


// Find boundary points that connect to more than one cell region and
// split them.
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


// Find boundary points that connect to more than one cell region and
// split them.
Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::dupNonManifoldPoints()
{
    // Analyse which points need to be duplicated
    localPointRegion regionSide(mesh_);

    return dupNonManifoldPoints(regionSide);
}


// Zoning
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
            label faceI = testFaces[i];

            if (mesh_.isInternalFace(faceI))
            {
                start[i] = cellCentres[faceOwner[faceI]];
                end[i] = cellCentres[faceNeighbour[faceI]];
            }
            else
            {
                start[i] = cellCentres[faceOwner[faceI]];
                end[i] = neiCc[faceI-mesh_.nInternalFaces()];
            }
        }

        // Extend segments a bit
        {
            const vectorField smallVec(Foam::sqrt(SMALL)*(end-start));
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
            label faceI = testFaces[i];
            const vector& area = mesh_.faceAreas()[faceI];

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
                    namedSurfaceIndex[faceI] = surface2[i];
                    posOrientation[faceI] = ((area&normal2[i]) > 0);
                    nSurfFaces[surface2[i]]++;
                }
                else
                {
                    namedSurfaceIndex[faceI] = surface1[i];
                    posOrientation[faceI] = ((area&normal1[i]) > 0);
                    nSurfFaces[surface1[i]]++;
                }
            }
            else if (surface2[i] != -1)
            {
                namedSurfaceIndex[faceI] = surface2[i];
                posOrientation[faceI] = ((area&normal2[i]) > 0);
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


    // Make sure namedSurfaceIndex is unset inbetween same cell cell zones.
    if (!allowFreeStandingZoneFaces)
    {
        Info<< "Only keeping zone faces inbetween different cellZones."
            << nl << endl;

        makeConsistentFaceIndex(cellToZone, namedSurfaceIndex);
    }

    // Topochange container
    polyTopoChange meshMod(mesh_);



    // Get coupled neighbour cellZone. Set to -1 on non-coupled patches.
    labelList neiCellZone(mesh_.nFaces()-mesh_.nInternalFaces(), -1);
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label faceI = pp.start()+i;
                neiCellZone[faceI-mesh_.nInternalFaces()] =
                    cellToZone[mesh_.faceOwner()[faceI]];
            }
        }
    }
    syncTools::swapBoundaryFaceList(mesh_, neiCellZone);

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
                    namedSurfaceIndex,
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

            // Detect non-manifold edges
            labelList nMasterFacesPerEdge;
            calcPatchNumMasterFaces(isMasterFace, patch, nMasterFacesPerEdge);

            // Mark zones. Even a single original surface might create multiple
            // disconnected/non-manifold-connected zones
            labelList faceToZone;
            const label nZones = markPatchZones
            (
                patch,
                nMasterFacesPerEdge,
                faceToZone
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
                faceToZone,
                nPosOrientation,

                meshFlipMap
            );

            // Count per region the number of orientations (taking the new
            // flipMap into account)
            forAll(patch.addressing(), faceI)
            {
                label meshFaceI = patch.addressing()[faceI];

                if (isMasterFace[meshFaceI])
                {
                    label n = 1;
                    if
                    (
                        bool(posOrientation[meshFaceI])
                     == meshFlipMap[meshFaceI]
                    )
                    {
                        n = -1;
                    }

                    nPosOrientation.find(faceToZone[faceI])() += n;
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
                faceToZone,
                nPosOrientation,

                meshFlipMap
            );
        }
    }


    // Put the faces into the correct zone
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label surfI = namedSurfaceIndex[faceI];

        if (surfI != -1)
        {
            // Orient face zone to have slave cells in max cell zone.
            label ownZone = cellToZone[faceOwner[faceI]];
            label neiZone = cellToZone[faceNeighbour[faceI]];

            bool flip;

            label maxZone = max(ownZone, neiZone);

            if (maxZone == -1)
            {
                // free-standing face. Use geometrically derived orientation
                flip = meshFlipMap[faceI];
            }
            else if (ownZone == maxZone)
            {
                flip = false;
            }
            else
            {
                flip = true;
            }

            meshMod.setAction
            (
                polyModifyFace
                (
                    mesh_.faces()[faceI],           // modified face
                    faceI,                          // label of face
                    faceOwner[faceI],               // owner
                    faceNeighbour[faceI],           // neighbour
                    false,                          // face flip
                    -1,                             // patch for face
                    false,                          // remove from zone
                    surfaceToFaceZone[surfI],       // zone for face
                    flip                            // face flip in zone
                )
            );
        }
    }


    // Set owner as no-flip
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        label faceI = pp.start();

        forAll(pp, i)
        {
            label surfI = namedSurfaceIndex[faceI];

            if (surfI != -1)
            {
                label ownZone = cellToZone[faceOwner[faceI]];
                label neiZone = neiCellZone[faceI-mesh_.nInternalFaces()];

                bool flip;

                label maxZone = max(ownZone, neiZone);

                if (maxZone == -1)
                {
                    // free-standing face. Use geometrically derived orientation
                    flip = meshFlipMap[faceI];
                }
                else if (ownZone == neiZone)
                {
                    // Free-standing zone face or coupled boundary. Keep master
                    // face unflipped.
                    flip = !isMasterFace[faceI];
                }
                else if (neiZone == maxZone)
                {
                    flip = true;
                }
                else
                {
                    flip = false;
                }

                meshMod.setAction
                (
                    polyModifyFace
                    (
                        mesh_.faces()[faceI],           // modified face
                        faceI,                          // label of face
                        faceOwner[faceI],               // owner
                        -1,                             // neighbour
                        false,                          // face flip
                        patchI,                         // patch for face
                        false,                          // remove from zone
                        surfaceToFaceZone[surfI],       // zone for face
                        flip                            // face flip in zone
                    )
                );
            }
            faceI++;
        }
    }


    // Put the cells into the correct zone
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(cellToZone, cellI)
    {
        label zoneI = cellToZone[cellI];

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
