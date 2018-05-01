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

#include "removeFaces.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "meshTools.H"
#include "polyModifyFace.H"
#include "polyRemoveFace.H"
#include "polyRemoveCell.H"
#include "polyRemovePoint.H"
#include "syncTools.H"
#include "OFstream.H"
#include "indirectPrimitivePatch.H"
#include "Time.H"
#include "faceSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(removeFaces, 0);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Changes region of connected set of cells. Can be recursive since hopefully
// only small area of faces removed in one go.
void Foam::removeFaces::changeCellRegion
(
    const label celli,
    const label oldRegion,
    const label newRegion,
    labelList& cellRegion
) const
{
    if (cellRegion[celli] == oldRegion)
    {
        cellRegion[celli] = newRegion;

        // Step to neighbouring cells

        const labelList& cCells = mesh_.cellCells()[celli];

        forAll(cCells, i)
        {
            changeCellRegion(cCells[i], oldRegion, newRegion, cellRegion);
        }
    }
}


// Changes region of connected set of faces. Returns number of changed faces.
Foam::label Foam::removeFaces::changeFaceRegion
(
    const labelList& cellRegion,
    const boolList& removedFace,
    const labelList& nFacesPerEdge,
    const label facei,
    const label newRegion,
    const labelList& fEdges,
    labelList& faceRegion
) const
{
    label nChanged = 0;

    if (faceRegion[facei] == -1 && !removedFace[facei])
    {
        faceRegion[facei] = newRegion;

        nChanged = 1;

        // Storage for on-the-fly addressing
        DynamicList<label> fe;
        DynamicList<label> ef;

        // Step to neighbouring faces across edges that will get removed
        forAll(fEdges, i)
        {
            label edgeI = fEdges[i];

            if (nFacesPerEdge[edgeI] >= 0 && nFacesPerEdge[edgeI] <= 2)
            {
                const labelList& eFaces = mesh_.edgeFaces(edgeI, ef);

                forAll(eFaces, j)
                {
                    label nbrFacei = eFaces[j];

                    const labelList& fEdges1 = mesh_.faceEdges(nbrFacei, fe);

                    nChanged += changeFaceRegion
                    (
                        cellRegion,
                        removedFace,
                        nFacesPerEdge,
                        nbrFacei,
                        newRegion,
                        fEdges1,
                        faceRegion
                    );
                }
            }
        }
    }
    return nChanged;
}


// Mark all faces affected in any way by
// - removal of cells
// - removal of faces
// - removal of edges
// - removal of points
Foam::boolList Foam::removeFaces::getFacesAffected
(
    const labelList& cellRegion,
    const labelList& cellRegionMaster,
    const labelList& facesToRemove,
    const labelHashSet& edgesToRemove,
    const labelHashSet& pointsToRemove
) const
{
    boolList affectedFace(mesh_.nFaces(), false);

    // Mark faces affected by removal of cells
    forAll(cellRegion, celli)
    {
        label region = cellRegion[celli];

        if (region != -1 && (celli != cellRegionMaster[region]))
        {
            const labelList& cFaces = mesh_.cells()[celli];

            forAll(cFaces, cFacei)
            {
                affectedFace[cFaces[cFacei]] = true;
            }
        }
    }

    // Mark faces affected by removal of face.
    forAll(facesToRemove, i)
    {
         affectedFace[facesToRemove[i]] = true;
    }

    //  Mark faces affected by removal of edges
    forAllConstIter(labelHashSet, edgesToRemove, iter)
    {
        const labelList& eFaces = mesh_.edgeFaces(iter.key());

        forAll(eFaces, eFacei)
        {
            affectedFace[eFaces[eFacei]] = true;
        }
    }

    // Mark faces affected by removal of points
    forAllConstIter(labelHashSet, pointsToRemove, iter)
    {
        label pointi = iter.key();

        const labelList& pFaces = mesh_.pointFaces()[pointi];

        forAll(pFaces, pFacei)
        {
            affectedFace[pFaces[pFacei]] = true;
        }
    }
    return affectedFace;
}


void Foam::removeFaces::writeOBJ
(
    const indirectPrimitivePatch& fp,
    const fileName& fName
)
{
    OFstream str(fName);
    Pout<< "removeFaces::writeOBJ : Writing faces to file "
        << str.name() << endl;

    const pointField& localPoints = fp.localPoints();

    forAll(localPoints, i)
    {
        meshTools::writeOBJ(str, localPoints[i]);
    }

    const faceList& localFaces = fp.localFaces();

    forAll(localFaces, i)
    {
        const face& f = localFaces[i];

        str<< 'f';

        forAll(f, fp)
        {
            str<< ' ' << f[fp]+1;
        }
        str<< nl;
    }
}


// Inserts commands to merge faceLabels into one face.
void Foam::removeFaces::mergeFaces
(
    const labelList& cellRegion,
    const labelList& cellRegionMaster,
    const labelHashSet& pointsToRemove,
    const labelList& faceLabels,
    polyTopoChange& meshMod
) const
{
    // Construct addressing engine from faceLabels (in order of faceLabels as
    // well)
    indirectPrimitivePatch fp
    (
        IndirectList<face>
        (
            mesh_.faces(),
            faceLabels
        ),
        mesh_.points()
    );

    // Get outside vertices (in local vertex numbering)

    if (fp.edgeLoops().size() != 1)
    {
        writeOBJ(fp, mesh_.time().path()/"facesToBeMerged.obj");
        FatalErrorInFunction
            << "Cannot merge faces " << faceLabels
            << " into single face since outside vertices " << fp.edgeLoops()
            << " do not form single loop but form " << fp.edgeLoops().size()
            << " loops instead." << abort(FatalError);
    }

    const labelList& edgeLoop = fp.edgeLoops()[0];

    // Get outside vertices in order of one of the faces in faceLabels.
    // (this becomes the master face)
    // Find the first face that uses edgeLoop[0] and edgeLoop[1] as consecutive
    // vertices.

    label masterIndex = -1;
    bool reverseLoop = false;

    const labelList& pFaces = fp.pointFaces()[edgeLoop[0]];

    // Find face among pFaces which uses edgeLoop[1]
    forAll(pFaces, i)
    {
        label facei = pFaces[i];

        const face& f = fp.localFaces()[facei];

        label index1 = findIndex(f, edgeLoop[1]);

        if (index1 != -1)
        {
            // Check whether consecutive to edgeLoop[0]
            label index0 = findIndex(f, edgeLoop[0]);

            if (index0 != -1)
            {
                if (index1 == f.fcIndex(index0))
                {
                    masterIndex = facei;
                    reverseLoop = false;
                    break;
                }
                else if (index1 == f.rcIndex(index0))
                {
                    masterIndex = facei;
                    reverseLoop = true;
                    break;
                }
            }
        }
    }

    if (masterIndex == -1)
    {
        writeOBJ(fp, mesh_.time().path()/"facesToBeMerged.obj");
        FatalErrorInFunction
            << "Problem" << abort(FatalError);
    }


    // Modify the master face.
    // ~~~~~~~~~~~~~~~~~~~~~~~

    // Modify first face.
    label facei = faceLabels[masterIndex];

    label own = mesh_.faceOwner()[facei];

    if (cellRegion[own] != -1)
    {
        own = cellRegionMaster[cellRegion[own]];
    }

    label patchID, zoneID, zoneFlip;

    getFaceInfo(facei, patchID, zoneID, zoneFlip);

    label nei = -1;

    if (mesh_.isInternalFace(facei))
    {
        nei = mesh_.faceNeighbour()[facei];

        if (cellRegion[nei] != -1)
        {
            nei = cellRegionMaster[cellRegion[nei]];
        }
    }


    DynamicList<label> faceVerts(edgeLoop.size());

    forAll(edgeLoop, i)
    {
        label pointi = fp.meshPoints()[edgeLoop[i]];

        if (pointsToRemove.found(pointi))
        {
            // Pout<< "**Removing point " << pointi << " from "
            //    << edgeLoop << endl;
        }
        else
        {
            faceVerts.append(pointi);
        }
    }

    face mergedFace;
    mergedFace.transfer(faceVerts);

    if (reverseLoop)
    {
        reverse(mergedFace);
    }

    //{
    //    Pout<< "Modifying masterface " << facei
    //        << " from faces:" << faceLabels
    //        << " old verts:" << UIndirectList<face>(mesh_.faces(), faceLabels)
    //        << " for new verts:"
    //        << mergedFace
    //        << " possibly new owner " << own
    //        << " or new nei " << nei
    //        << endl;
    //}

    modFace
    (
        mergedFace,         // modified face
        facei,              // label of face being modified
        own,                // owner
        nei,                // neighbour
        false,              // face flip
        patchID,            // patch for face
        false,              // remove from zone
        zoneID,             // zone for face
        zoneFlip,           // face flip in zone

        meshMod
    );


    // Remove all but master face.
    forAll(faceLabels, patchFacei)
    {
        if (patchFacei != masterIndex)
        {
            // Pout<< "Removing face " << faceLabels[patchFacei] << endl;

            meshMod.setAction(polyRemoveFace(faceLabels[patchFacei], facei));
        }
    }
}


// Get patch, zone info for facei
void Foam::removeFaces::getFaceInfo
(
    const label facei,

    label& patchID,
    label& zoneID,
    label& zoneFlip
) const
{
    patchID = -1;

    if (!mesh_.isInternalFace(facei))
    {
        patchID = mesh_.boundaryMesh().whichPatch(facei);
    }

    zoneID = mesh_.faceZones().whichZone(facei);

    zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh_.faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
    }
}


// Return face with all pointsToRemove removed.
Foam::face Foam::removeFaces::filterFace
(
    const labelHashSet& pointsToRemove,
    const label facei
) const
{
    const face& f = mesh_.faces()[facei];

    labelList newFace(f.size(), -1);

    label newFp = 0;

    forAll(f, fp)
    {
        label vertI = f[fp];

        if (!pointsToRemove.found(vertI))
        {
            newFace[newFp++] = vertI;
        }
    }

    newFace.setSize(newFp);

    return face(newFace);
}


// Wrapper for meshMod.modifyFace. Reverses face if own>nei.
void Foam::removeFaces::modFace
(
    const face& f,
    const label masterFaceID,
    const label own,
    const label nei,
    const bool flipFaceFlux,
    const label newPatchID,
    const bool removeFromZone,
    const label zoneID,
    const bool zoneFlip,

    polyTopoChange& meshMod
) const
{
    if ((nei == -1) || (own < nei))
    {
//        if (debug)
//        {
//            Pout<< "ModifyFace (unreversed) :"
//                << "  facei:" << masterFaceID
//                << "  f:" << f
//                << "  own:" << own
//                << "  nei:" << nei
//                << "  flipFaceFlux:" << flipFaceFlux
//                << "  newPatchID:" << newPatchID
//                << "  removeFromZone:" << removeFromZone
//                << "  zoneID:" << zoneID
//                << "  zoneFlip:" << zoneFlip
//                << endl;
//        }

        meshMod.setAction
        (
            polyModifyFace
            (
                f,              // modified face
                masterFaceID,   // label of face being modified
                own,            // owner
                nei,            // neighbour
                flipFaceFlux,   // face flip
                newPatchID,     // patch for face
                removeFromZone, // remove from zone
                zoneID,         // zone for face
                zoneFlip        // face flip in zone
            )
        );
    }
    else
    {
//        if (debug)
//        {
//            Pout<< "ModifyFace (!reversed) :"
//                << "  facei:" << masterFaceID
//                << "  f:" << f.reverseFace()
//                << "  own:" << nei
//                << "  nei:" << own
//                << "  flipFaceFlux:" << flipFaceFlux
//                << "  newPatchID:" << newPatchID
//                << "  removeFromZone:" << removeFromZone
//                << "  zoneID:" << zoneID
//                << "  zoneFlip:" << zoneFlip
//                << endl;
//        }

        meshMod.setAction
        (
            polyModifyFace
            (
                f.reverseFace(),// modified face
                masterFaceID,   // label of face being modified
                nei,            // owner
                own,            // neighbour
                flipFaceFlux,   // face flip
                newPatchID,     // patch for face
                removeFromZone, // remove from zone
                zoneID,         // zone for face
                zoneFlip        // face flip in zone
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::removeFaces::removeFaces
(
    const polyMesh& mesh,
    const scalar minCos
)
:
    mesh_(mesh),
    minCos_(minCos)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Removing face connects cells. This function works out a consistent set of
// cell regions.
// - returns faces to remove. Can be extended with additional faces
//   (if owner would become neighbour)
// - sets cellRegion to -1 or to region number
// - regionMaster contains for every region number a master cell.
Foam::label Foam::removeFaces::compatibleRemoves
(
    const labelList& facesToRemove,
    labelList& cellRegion,
    labelList& regionMaster,
    labelList& newFacesToRemove
) const
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    cellRegion.setSize(mesh_.nCells());
    cellRegion = -1;

    regionMaster.setSize(mesh_.nCells());
    regionMaster = -1;

    label nRegions = 0;

    forAll(facesToRemove, i)
    {
        label facei = facesToRemove[i];

        if (!mesh_.isInternalFace(facei))
        {
            FatalErrorInFunction
                << "Not internal face:" << facei << abort(FatalError);
        }


        label own = faceOwner[facei];
        label nei = faceNeighbour[facei];

        label region0 = cellRegion[own];
        label region1 = cellRegion[nei];

        if (region0 == -1)
        {
            if (region1 == -1)
            {
                // Create new region
                cellRegion[own] = nRegions;
                cellRegion[nei] = nRegions;

                // Make owner (lowest numbered!) the master of the region
                regionMaster[nRegions] = own;
                nRegions++;
            }
            else
            {
                // Add owner to neighbour region
                cellRegion[own] = region1;
                // See if owner becomes the master of the region
                regionMaster[region1] = min(own, regionMaster[region1]);
            }
        }
        else
        {
            if (region1 == -1)
            {
                // Add neighbour to owner region
                cellRegion[nei] = region0;
                // nei is higher numbered than own so guaranteed not lower
                // than master of region0.
            }
            else if (region0 != region1)
            {
                // Both have regions. Keep lowest numbered region and master.
                label freedRegion = -1;
                label keptRegion = -1;

                if (region0 < region1)
                {
                    changeCellRegion
                    (
                        nei,
                        region1,    // old region
                        region0,    // new region
                        cellRegion
                    );

                    keptRegion = region0;
                    freedRegion = region1;
                }
                else if (region1 < region0)
                {
                    changeCellRegion
                    (
                        own,
                        region0,    // old region
                        region1,    // new region
                        cellRegion
                    );

                    keptRegion = region1;
                    freedRegion = region0;
                }

                label master0 = regionMaster[region0];
                label master1 = regionMaster[region1];

                regionMaster[freedRegion] = -1;
                regionMaster[keptRegion] = min(master0, master1);
            }
        }
    }

    regionMaster.setSize(nRegions);


    // Various checks
    // - master is lowest numbered in any region
    // - regions have more than 1 cell
    {
        labelList nCells(regionMaster.size(), 0);

        forAll(cellRegion, celli)
        {
            label r = cellRegion[celli];

            if (r != -1)
            {
                nCells[r]++;

                if (celli < regionMaster[r])
                {
                    FatalErrorInFunction
                        << "Not lowest numbered : cell:" << celli
                        << " region:" << r
                        << " regionmaster:" << regionMaster[r]
                        << abort(FatalError);
                }
            }
        }

        forAll(nCells, region)
        {
            if (nCells[region] == 1)
            {
                FatalErrorInFunction
                    << "Region " << region
                    << " has only " << nCells[region] << " cells in it"
                    << abort(FatalError);
            }
        }
    }


    // Count number of used regions
    label nUsedRegions = 0;

    forAll(regionMaster, i)
    {
        if (regionMaster[i] != -1)
        {
            nUsedRegions++;
        }
    }

    // Recreate facesToRemove to be consistent with the cellRegions.
    DynamicList<label> allFacesToRemove(facesToRemove.size());

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        label own = faceOwner[facei];
        label nei = faceNeighbour[facei];

        if (cellRegion[own] != -1 && cellRegion[own] == cellRegion[nei])
        {
            // Both will become the same cell so add face to list of faces
            // to be removed.
            allFacesToRemove.append(facei);
        }
    }

    newFacesToRemove.transfer(allFacesToRemove);

    return nUsedRegions;
}


void Foam::removeFaces::setRefinement
(
    const labelList& faceLabels,
    const labelList& cellRegion,
    const labelList& cellRegionMaster,
    polyTopoChange& meshMod
) const
{
    if (debug)
    {
        faceSet facesToRemove(mesh_, "facesToRemove", faceLabels);
        Pout<< "Writing faces to remove to faceSet " << facesToRemove.name()
            << endl;
        facesToRemove.write();
    }

    // Make map of all faces to be removed
    boolList removedFace(mesh_.nFaces(), false);

    forAll(faceLabels, i)
    {
        label facei = faceLabels[i];

        if (!mesh_.isInternalFace(facei))
        {
            FatalErrorInFunction
                << "Face to remove is not internal face:" << facei
                << abort(FatalError);
        }

        removedFace[facei] = true;
    }


    // Edges to be removed
    // ~~~~~~~~~~~~~~~~~~~


    // Edges to remove
    labelHashSet edgesToRemove(faceLabels.size());

    // Per face the region it is in. -1 for removed faces, -2 for regions
    // consisting of single face only.
    labelList faceRegion(mesh_.nFaces(), -1);

    // Number of connected face regions
    label nRegions = 0;

    // Storage for on-the-fly addressing
    DynamicList<label> fe;
    DynamicList<label> ef;


    {
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        // Usage of edges by non-removed faces.
        // See below about initialization.
        labelList nFacesPerEdge(mesh_.nEdges(), -1);

        // Count usage of edges by non-removed faces.
        forAll(faceLabels, i)
        {
            label facei = faceLabels[i];

            const labelList& fEdges = mesh_.faceEdges(facei, fe);

            forAll(fEdges, i)
            {
                label edgeI = fEdges[i];

                if (nFacesPerEdge[edgeI] == -1)
                {
                    nFacesPerEdge[edgeI] = mesh_.edgeFaces(edgeI, ef).size()-1;
                }
                else
                {
                    nFacesPerEdge[edgeI]--;
                }
            }
        }

        // Count usage for edges not on faces-to-be-removed.
        // Note that this only needs to be done for possibly coupled edges
        // so we could choose to loop only over boundary faces and use faceEdges
        // of those.

        forAll(mesh_.edges(), edgeI)
        {
            if (nFacesPerEdge[edgeI] == -1)
            {
                // Edge not yet handled in loop above so is not used by any
                // face to be removed.

                const labelList& eFaces = mesh_.edgeFaces(edgeI, ef);

                if (eFaces.size() > 2)
                {
                    nFacesPerEdge[edgeI] = eFaces.size();
                }
                else if (eFaces.size() == 2)
                {
                    // nFacesPerEdge already -1 so do nothing.
                }
                else
                {
                    const edge& e = mesh_.edges()[edgeI];

                    FatalErrorInFunction
                        << "Problem : edge has too few face neighbours:"
                        << eFaces << endl
                        << "edge:" << edgeI
                        << " vertices:" << e
                        << " coords:" << mesh_.points()[e[0]]
                        << mesh_.points()[e[1]]
                        << abort(FatalError);
                }
            }
        }



        if (debug)
        {
            OFstream str(mesh_.time().path()/"edgesWithTwoFaces.obj");
            Pout<< "Dumping edgesWithTwoFaces to " << str.name() << endl;
            label vertI = 0;

            forAll(nFacesPerEdge, edgeI)
            {
                if (nFacesPerEdge[edgeI] == 2)
                {
                    // Edge will get removed.
                    const edge& e = mesh_.edges()[edgeI];

                    meshTools::writeOBJ(str, mesh_.points()[e[0]]);
                    vertI++;
                    meshTools::writeOBJ(str, mesh_.points()[e[1]]);
                    vertI++;
                    str<< "l " << vertI-1 << ' ' << vertI << nl;
                }
            }
        }


        // Now all unaffected edges will have labelMax, all affected edges the
        // number of unremoved faces.

        // Filter for edges in between two remaining boundary faces that
        // make too big an angle.
        forAll(nFacesPerEdge, edgeI)
        {
            if (nFacesPerEdge[edgeI] == 2)
            {
                // See if they are two boundary faces
                label f0 = -1;
                label f1 = -1;

                const labelList& eFaces = mesh_.edgeFaces(edgeI, ef);

                forAll(eFaces, i)
                {
                    label facei = eFaces[i];

                    if (!removedFace[facei] && !mesh_.isInternalFace(facei))
                    {
                        if (f0 == -1)
                        {
                            f0 = facei;
                        }
                        else
                        {
                            f1 = facei;
                            break;
                        }
                    }
                }

                if (f0 != -1 && f1 != -1)
                {
                    // Edge has two boundary faces remaining.
                    // See if should be merged.

                    label patch0 = patches.whichPatch(f0);
                    label patch1 = patches.whichPatch(f1);

                    if (patch0 != patch1)
                    {
                        // Different patches. Do not merge edge.
                        WarningInFunction
                            << "not merging faces " << f0 << " and "
                            << f1 << " across patch boundary edge " << edgeI
                            << endl;

                        // Mark so it gets preserved
                        nFacesPerEdge[edgeI] = 3;
                    }
                    else if (minCos_ < 1 && minCos_ > -1)
                    {
                        const polyPatch& pp0 = patches[patch0];
                        const vectorField& n0 = pp0.faceNormals();

                        if
                        (
                            mag
                            (
                                n0[f0 - pp0.start()]
                              & n0[f1 - pp0.start()]
                            )
                            < minCos_
                        )
                        {
                            WarningInFunction
                                << "not merging faces " << f0 << " and "
                                << f1 << " across edge " << edgeI
                                << endl;

                            // Angle between two remaining faces too large.
                            // Mark so it gets preserved
                            nFacesPerEdge[edgeI] = 3;
                        }
                    }
                }
                else if (f0 != -1 || f1 != -1)
                {
                    const edge& e = mesh_.edges()[edgeI];

                    // Only found one boundary face. Problem.
                    FatalErrorInFunction
                        << "Problem : edge would have one boundary face"
                        << " and one internal face using it." << endl
                        << "Your remove pattern is probably incorrect." << endl
                        << "edge:" << edgeI
                        << " nFaces:" << nFacesPerEdge[edgeI]
                        << " vertices:" << e
                        << " coords:" << mesh_.points()[e[0]]
                        << mesh_.points()[e[1]]
                        << " face0:" << f0
                        << " face1:" << f1
                        << abort(FatalError);
                }
            }
        }



        // Check locally (before synchronizing) for strangeness
        forAll(nFacesPerEdge, edgeI)
        {
            if (nFacesPerEdge[edgeI] == 1)
            {
                const edge& e = mesh_.edges()[edgeI];

                FatalErrorInFunction
                    << "Problem : edge would get 1 face using it only"
                    << " edge:" << edgeI
                    << " nFaces:" << nFacesPerEdge[edgeI]
                    << " vertices:" << e
                    << " coords:" << mesh_.points()[e[0]]
                    << ' ' << mesh_.points()[e[1]]
                    << abort(FatalError);
            }
            // Could check here for boundary edge with <=1 faces remaining.
        }


        // Synchronize edge usage. This is to make sure that both sides remove
        // (or not remove) an edge on the boundary at the same time.
        //
        // Coupled edges (edge0, edge1 are opposite each other)
        // a. edge not on face to be removed, edge has >= 3 faces
        // b.  ,,                             edge has 2 faces
        // c. edge has >= 3 remaining faces
        // d. edge has 2 remaining faces (assume angle>minCos already handled)
        //
        // - a + a: do not remove edge
        // - a + b: do not remove edge
        // - a + c: do not remove edge
        // - a + d: do not remove edge
        //
        // - b + b: do not remove edge
        // - b + c: do not remove edge
        // - b + d: remove edge
        //
        // - c + c: do not remove edge
        // - c + d: do not remove edge
        // - d + d: remove edge
        //
        //
        // So code situation a. with >= 3
        //                   b. with -1
        //                   c. with >=3
        //                   d. with 2
        // then do max and check result.
        //
        // a+a : max(3,3) = 3. do not remove
        // a+b : max(3,-1) = 3. do not remove
        // a+c : max(3,3) = 3. do not remove
        // a+d : max(3,2) = 3. do not remove
        //
        // b+b : max(-1,-1) = -1. do not remove
        // b+c : max(-1,3) = 3. do not remove
        // b+d : max(-1,2) = 2. remove
        //
        // c+c : max(3,3) = 3. do not remove
        // c+d : max(3,2) = 3. do not remove
        //
        // d+d : max(2,2) = 2. remove

        syncTools::syncEdgeList
        (
            mesh_,
            nFacesPerEdge,
            maxEqOp<label>(),
            labelMin                // guaranteed to be overridden by maxEqOp
        );

        // Convert to labelHashSet
        forAll(nFacesPerEdge, edgeI)
        {
            if (nFacesPerEdge[edgeI] == 0)
            {
                // 0: edge not used anymore.
                edgesToRemove.insert(edgeI);
            }
            else if (nFacesPerEdge[edgeI] == 1)
            {
                // 1: illegal. Tested above.
            }
            else if (nFacesPerEdge[edgeI] == 2)
            {
                // 2: merge faces.
                edgesToRemove.insert(edgeI);
            }
        }

        if (debug)
        {
            OFstream str(mesh_.time().path()/"edgesToRemove.obj");
            Pout<< "Dumping edgesToRemove to " << str.name() << endl;
            label vertI = 0;

            forAllConstIter(labelHashSet, edgesToRemove, iter)
            {
                // Edge will get removed.
                const edge& e = mesh_.edges()[iter.key()];

                meshTools::writeOBJ(str, mesh_.points()[e[0]]);
                vertI++;
                meshTools::writeOBJ(str, mesh_.points()[e[1]]);
                vertI++;
                str<< "l " << vertI-1 << ' ' << vertI << nl;
            }
        }


        // Walk to fill faceRegion with faces that will be connected across
        // edges that will be removed.

        label startFacei = 0;

        while (true)
        {
            // Find unset region.
            for (; startFacei < mesh_.nFaces(); startFacei++)
            {
                if (faceRegion[startFacei] == -1 && !removedFace[startFacei])
                {
                    break;
                }
            }

            if (startFacei == mesh_.nFaces())
            {
                break;
            }

            // Start walking face-edge-face, crossing edges that will get
            // removed. Every thus connected region will get single region
            // number.
            label nRegion = changeFaceRegion
            (
                cellRegion,
                removedFace,
                nFacesPerEdge,
                startFacei,
                nRegions,
                mesh_.faceEdges(startFacei, fe),
                faceRegion
            );

            if (nRegion < 1)
            {
                FatalErrorInFunction << "Problem" << abort(FatalError);
            }
            else if (nRegion == 1)
            {
                // Reset face to be single region.
                faceRegion[startFacei] = -2;
            }
            else
            {
                nRegions++;
            }
        }


        // Check we're deciding the same on both sides. Since the regioning
        // is done based on nFacesPerEdge (which is synced) this should
        // indeed be the case.

        labelList nbrFaceRegion(faceRegion);

        syncTools::swapFaceList
        (
            mesh_,
            nbrFaceRegion
        );

        labelList toNbrRegion(nRegions, -1);

        for
        (
            label facei = mesh_.nInternalFaces();
            facei < mesh_.nFaces();
            facei++
        )
        {
            // Get the neighbouring region.
            label nbrRegion = nbrFaceRegion[facei];
            label myRegion = faceRegion[facei];

            if (myRegion <= -1 || nbrRegion <= -1)
            {
                if (nbrRegion != myRegion)
                {
                    FatalErrorInFunction
                        << "Inconsistent face region across coupled patches."
                        << endl
                        << "This side has for facei:" << facei
                        << " region:" << myRegion << endl
                        << "The other side has region:" << nbrRegion
                        << endl
                        << "(region -1 means face is to be deleted)"
                        << abort(FatalError);
                }
            }
            else if (toNbrRegion[myRegion] == -1)
            {
                // First visit of region. Store correspondence.
                toNbrRegion[myRegion] = nbrRegion;
            }
            else
            {
                // Second visit of this region.
                if (toNbrRegion[myRegion] != nbrRegion)
                {
                    FatalErrorInFunction
                        << "Inconsistent face region across coupled patches."
                        << endl
                        << "This side has for facei:" << facei
                        << " region:" << myRegion
                        << " with coupled neighbouring regions:"
                        << toNbrRegion[myRegion] << " and "
                        << nbrRegion
                        << abort(FatalError);
                }
            }
        }
    }

    // if (debug)
    //{
    //    labelListList regionToFaces(invertOneToMany(nRegions, faceRegion));
    //
    //    forAll(regionToFaces, regionI)
    //    {
    //        Pout<< "    " << regionI << " faces:" << regionToFaces[regionI]
    //            << endl;
    //    }
    //}


    // Points to be removed
    // ~~~~~~~~~~~~~~~~~~~~

    labelHashSet pointsToRemove(4*faceLabels.size());


    // Per point count the number of unremoved edges. Store the ones that
    // are only used by 2 unremoved edges.
    {
        // Usage of points by non-removed edges.
        labelList nEdgesPerPoint(mesh_.nPoints());

        const labelListList& pointEdges = mesh_.pointEdges();

        forAll(pointEdges, pointi)
        {
            nEdgesPerPoint[pointi] = pointEdges[pointi].size();
        }

        forAllConstIter(labelHashSet, edgesToRemove, iter)
        {
            // Edge will get removed.
            const edge& e = mesh_.edges()[iter.key()];

            forAll(e, i)
            {
                nEdgesPerPoint[e[i]]--;
            }
        }

        // Check locally (before synchronizing) for strangeness
        forAll(nEdgesPerPoint, pointi)
        {
            if (nEdgesPerPoint[pointi] == 1)
            {
                FatalErrorInFunction
                    << "Problem : point would get 1 edge using it only."
                    << " pointi:" << pointi
                    << " coord:" << mesh_.points()[pointi]
                    << abort(FatalError);
            }
        }

        // Synchronize point usage. This is to make sure that both sides remove
        // (or not remove) a point on the boundary at the same time.
        syncTools::syncPointList
        (
            mesh_,
            nEdgesPerPoint,
            maxEqOp<label>(),
            labelMin
        );

        forAll(nEdgesPerPoint, pointi)
        {
            if (nEdgesPerPoint[pointi] == 0)
            {
                pointsToRemove.insert(pointi);
            }
            else if (nEdgesPerPoint[pointi] == 1)
            {
                // Already checked before
            }
            else if (nEdgesPerPoint[pointi] == 2)
            {
                // Remove point and merge edges.
                pointsToRemove.insert(pointi);
            }
        }
    }


    if (debug)
    {
        OFstream str(mesh_.time().path()/"pointsToRemove.obj");
        Pout<< "Dumping pointsToRemove to " << str.name() << endl;

        forAllConstIter(labelHashSet, pointsToRemove, iter)
        {
            meshTools::writeOBJ(str, mesh_.points()[iter.key()]);
        }
    }


    // All faces affected in any way
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Get all faces affected in any way by removal of points/edges/faces/cells
    boolList affectedFace
    (
        getFacesAffected
        (
            cellRegion,
            cellRegionMaster,
            faceLabels,
            edgesToRemove,
            pointsToRemove
        )
    );

    //
    // Now we know
    // - faceLabels         : faces to remove (sync since no boundary faces)
    // - cellRegion/Master  : cells to remove (sync since cells)
    // - pointsToRemove     : points to remove (sync)
    // - faceRegion         : connected face region of faces to be merged (sync)
    // - affectedFace       : faces with points removed and/or owner/neighbour
    //                        changed (non sync)


    // Start modifying mesh and keep track of faces changed.


    // Do all removals
    // ~~~~~~~~~~~~~~~

    // Remove split faces.
    forAll(faceLabels, labelI)
    {
        label facei = faceLabels[labelI];

        // Remove face if not yet uptodate (which is never; but want to be
        // consistent with rest of face removals/modifications)
        if (affectedFace[facei])
        {
            affectedFace[facei] = false;

            meshMod.setAction(polyRemoveFace(facei, -1));
        }
    }


    // Remove points.
    forAllConstIter(labelHashSet, pointsToRemove, iter)
    {
        label pointi = iter.key();

        meshMod.setAction(polyRemovePoint(pointi, -1));
    }


    // Remove cells.
    forAll(cellRegion, celli)
    {
        label region = cellRegion[celli];

        if (region != -1 && (celli != cellRegionMaster[region]))
        {
            meshMod.setAction(polyRemoveCell(celli, cellRegionMaster[region]));
        }
    }



    // Merge faces across edges to be merged
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Invert faceRegion so we get region to faces.
    {
        labelListList regionToFaces(invertOneToMany(nRegions, faceRegion));

        forAll(regionToFaces, regionI)
        {
            const labelList& rFaces = regionToFaces[regionI];

            if (rFaces.size() <= 1)
            {
                FatalErrorInFunction
                    << "Region:" << regionI
                    << " contains only faces " << rFaces
                    << abort(FatalError);
            }

            // rFaces[0] is master, rest gets removed.
            mergeFaces
            (
                cellRegion,
                cellRegionMaster,
                pointsToRemove,
                rFaces,
                meshMod
            );

            forAll(rFaces, i)
            {
                affectedFace[rFaces[i]] = false;
            }
        }
    }


    // Remaining affected faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    // Check if any remaining faces have not been updated for new slave/master
    // or points removed.
    forAll(affectedFace, facei)
    {
        if (affectedFace[facei])
        {
            affectedFace[facei] = false;

            face f(filterFace(pointsToRemove, facei));

            label own = mesh_.faceOwner()[facei];

            if (cellRegion[own] != -1)
            {
                own = cellRegionMaster[cellRegion[own]];
            }

            label patchID, zoneID, zoneFlip;

            getFaceInfo(facei, patchID, zoneID, zoneFlip);

            label nei = -1;

            if (mesh_.isInternalFace(facei))
            {
                nei = mesh_.faceNeighbour()[facei];

                if (cellRegion[nei] != -1)
                {
                    nei = cellRegionMaster[cellRegion[nei]];
                }
            }

//            if (debug)
//            {
//                Pout<< "Modifying " << facei
//                    << " old verts:" << mesh_.faces()[facei]
//                    << " for new verts:" << f
//                    << " or for new owner " << own << " or for new nei "
//                    << nei
//                    << endl;
//            }

            modFace
            (
                f,                  // modified face
                facei,              // label of face being modified
                own,                // owner
                nei,                // neighbour
                false,              // face flip
                patchID,            // patch for face
                false,              // remove from zone
                zoneID,             // zone for face
                zoneFlip,           // face flip in zone

                meshMod
            );
        }
    }
}


// ************************************************************************* //
