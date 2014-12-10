/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    const label cellI,
    const label oldRegion,
    const label newRegion,
    labelList& cellRegion
) const
{
    if (cellRegion[cellI] == oldRegion)
    {
        cellRegion[cellI] = newRegion;

        // Step to neighbouring cells

        const labelList& cCells = mesh_.cellCells()[cellI];

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
    const label faceI,
    const label newRegion,
    const labelList& fEdges,
    labelList& faceRegion
) const
{
    label nChanged = 0;

    if (faceRegion[faceI] == -1 && !removedFace[faceI])
    {
        faceRegion[faceI] = newRegion;

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
                    label nbrFaceI = eFaces[j];

                    const labelList& fEdges1 = mesh_.faceEdges(nbrFaceI, fe);

                    nChanged += changeFaceRegion
                    (
                        cellRegion,
                        removedFace,
                        nFacesPerEdge,
                        nbrFaceI,
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
    forAll(cellRegion, cellI)
    {
        label region = cellRegion[cellI];

        if (region != -1 && (cellI != cellRegionMaster[region]))
        {
            const labelList& cFaces = mesh_.cells()[cellI];

            forAll(cFaces, cFaceI)
            {
                affectedFace[cFaces[cFaceI]] = true;
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

        forAll(eFaces, eFaceI)
        {
            affectedFace[eFaces[eFaceI]] = true;
        }
    }

    // Mark faces affected by removal of points
    forAllConstIter(labelHashSet, pointsToRemove, iter)
    {
        label pointI = iter.key();

        const labelList& pFaces = mesh_.pointFaces()[pointI];

        forAll(pFaces, pFaceI)
        {
            affectedFace[pFaces[pFaceI]] = true;
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
        FatalErrorIn("removeFaces::mergeFaces")
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
        label faceI = pFaces[i];

        const face& f = fp.localFaces()[faceI];

        label index1 = findIndex(f, edgeLoop[1]);

        if (index1 != -1)
        {
            // Check whether consecutive to edgeLoop[0]
            label index0 = findIndex(f, edgeLoop[0]);

            if (index0 != -1)
            {
                if (index1 == f.fcIndex(index0))
                {
                    masterIndex = faceI;
                    reverseLoop = false;
                    break;
                }
                else if (index1 == f.rcIndex(index0))
                {
                    masterIndex = faceI;
                    reverseLoop = true;
                    break;
                }
            }
        }
    }

    if (masterIndex == -1)
    {
        writeOBJ(fp, mesh_.time().path()/"facesToBeMerged.obj");
        FatalErrorIn("removeFaces::mergeFaces")
            << "Problem" << abort(FatalError);
    }


    // Modify the master face.
    // ~~~~~~~~~~~~~~~~~~~~~~~

    // Modify first face.
    label faceI = faceLabels[masterIndex];

    label own = mesh_.faceOwner()[faceI];

    if (cellRegion[own] != -1)
    {
        own = cellRegionMaster[cellRegion[own]];
    }

    label patchID, zoneID, zoneFlip;

    getFaceInfo(faceI, patchID, zoneID, zoneFlip);

    label nei = -1;

    if (mesh_.isInternalFace(faceI))
    {
        nei = mesh_.faceNeighbour()[faceI];

        if (cellRegion[nei] != -1)
        {
            nei = cellRegionMaster[cellRegion[nei]];
        }
    }


    DynamicList<label> faceVerts(edgeLoop.size());

    forAll(edgeLoop, i)
    {
        label pointI = fp.meshPoints()[edgeLoop[i]];

        if (pointsToRemove.found(pointI))
        {
            //Pout<< "**Removing point " << pointI << " from "
            //    << edgeLoop << endl;
        }
        else
        {
            faceVerts.append(pointI);
        }
    }

    face mergedFace;
    mergedFace.transfer(faceVerts);

    if (reverseLoop)
    {
        reverse(mergedFace);
    }

    //{
    //    Pout<< "Modifying masterface " << faceI
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
        faceI,              // label of face being modified
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
    forAll(faceLabels, patchFaceI)
    {
        if (patchFaceI != masterIndex)
        {
            //Pout<< "Removing face " << faceLabels[patchFaceI] << endl;

            meshMod.setAction(polyRemoveFace(faceLabels[patchFaceI], faceI));
        }
    }
}


// Get patch, zone info for faceI
void Foam::removeFaces::getFaceInfo
(
    const label faceI,

    label& patchID,
    label& zoneID,
    label& zoneFlip
) const
{
    patchID = -1;

    if (!mesh_.isInternalFace(faceI))
    {
        patchID = mesh_.boundaryMesh().whichPatch(faceI);
    }

    zoneID = mesh_.faceZones().whichZone(faceI);

    zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh_.faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
    }
}


// Return face with all pointsToRemove removed.
Foam::face Foam::removeFaces::filterFace
(
    const labelHashSet& pointsToRemove,
    const label faceI
) const
{
    const face& f = mesh_.faces()[faceI];

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
//                << "  faceI:" << masterFaceID
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
//                << "  faceI:" << masterFaceID
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
        label faceI = facesToRemove[i];

        if (!mesh_.isInternalFace(faceI))
        {
            FatalErrorIn
            (
                "removeFaces::compatibleRemoves(const labelList&"
                ", labelList&, labelList&, labelList&)"
            )   << "Not internal face:" << faceI << abort(FatalError);
        }


        label own = faceOwner[faceI];
        label nei = faceNeighbour[faceI];

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

        forAll(cellRegion, cellI)
        {
            label r = cellRegion[cellI];

            if (r != -1)
            {
                nCells[r]++;

                if (cellI < regionMaster[r])
                {
                    FatalErrorIn
                    (
                        "removeFaces::compatibleRemoves(const labelList&"
                        ", labelList&, labelList&, labelList&)"
                    )   << "Not lowest numbered : cell:" << cellI
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
                FatalErrorIn
                (
                    "removeFaces::compatibleRemoves(const labelList&"
                    ", labelList&, labelList&, labelList&)"
                )   << "Region " << region
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

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label own = faceOwner[faceI];
        label nei = faceNeighbour[faceI];

        if (cellRegion[own] != -1 && cellRegion[own] == cellRegion[nei])
        {
            // Both will become the same cell so add face to list of faces
            // to be removed.
            allFacesToRemove.append(faceI);
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
        label faceI = faceLabels[i];

        if (!mesh_.isInternalFace(faceI))
        {
            FatalErrorIn
            (
                "removeFaces::setRefinement(const labelList&"
                ", const labelList&, const labelList&, polyTopoChange&)"
            )   << "Face to remove is not internal face:" << faceI
                << abort(FatalError);
        }

        removedFace[faceI] = true;
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
            label faceI = faceLabels[i];

            const labelList& fEdges = mesh_.faceEdges(faceI, fe);

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

                    FatalErrorIn("removeFaces::setRefinement")
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

        // Filter for edges inbetween two remaining boundary faces that
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
                    label faceI = eFaces[i];

                    if (!removedFace[faceI] && !mesh_.isInternalFace(faceI))
                    {
                        if (f0 == -1)
                        {
                            f0 = faceI;
                        }
                        else
                        {
                            f1 = faceI;
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
                        WarningIn("removeFaces::setRefinement")
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
                            WarningIn("removeFaces::setRefinement")
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
                    FatalErrorIn("removeFaces::setRefinement")
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

                FatalErrorIn("removeFaces::setRefinement")
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

        label startFaceI = 0;

        while (true)
        {
            // Find unset region.
            for (; startFaceI < mesh_.nFaces(); startFaceI++)
            {
                if (faceRegion[startFaceI] == -1 && !removedFace[startFaceI])
                {
                    break;
                }
            }

            if (startFaceI == mesh_.nFaces())
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
                startFaceI,
                nRegions,
                mesh_.faceEdges(startFaceI, fe),
                faceRegion
            );

            if (nRegion < 1)
            {
                FatalErrorIn("setRefinement") << "Problem" << abort(FatalError);
            }
            else if (nRegion == 1)
            {
                // Reset face to be single region.
                faceRegion[startFaceI] = -2;
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
            label faceI = mesh_.nInternalFaces();
            faceI < mesh_.nFaces();
            faceI++
        )
        {
            // Get the neighbouring region.
            label nbrRegion = nbrFaceRegion[faceI];
            label myRegion = faceRegion[faceI];

            if (myRegion <= -1 || nbrRegion <= -1)
            {
                if (nbrRegion != myRegion)
                {
                    FatalErrorIn("removeFaces::setRefinement")
                        << "Inconsistent face region across coupled patches."
                        << endl
                        << "This side has for faceI:" << faceI
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
                    FatalErrorIn("removeFaces::setRefinement")
                        << "Inconsistent face region across coupled patches."
                        << endl
                        << "This side has for faceI:" << faceI
                        << " region:" << myRegion
                        << " with coupled neighbouring regions:"
                        << toNbrRegion[myRegion] << " and "
                        << nbrRegion
                        << abort(FatalError);
                }
            }
        }
    }

    //if (debug)
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

        forAll(pointEdges, pointI)
        {
            nEdgesPerPoint[pointI] = pointEdges[pointI].size();
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
        forAll(nEdgesPerPoint, pointI)
        {
            if (nEdgesPerPoint[pointI] == 1)
            {
                FatalErrorIn("removeFaces::setRefinement")
                    << "Problem : point would get 1 edge using it only."
                    << " pointI:" << pointI
                    << " coord:" << mesh_.points()[pointI]
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

        forAll(nEdgesPerPoint, pointI)
        {
            if (nEdgesPerPoint[pointI] == 0)
            {
                pointsToRemove.insert(pointI);
            }
            else if (nEdgesPerPoint[pointI] == 1)
            {
                // Already checked before
            }
            else if (nEdgesPerPoint[pointI] == 2)
            {
                // Remove point and merge edges.
                pointsToRemove.insert(pointI);
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
        label faceI = faceLabels[labelI];

        // Remove face if not yet uptodate (which is never; but want to be
        // consistent with rest of face removals/modifications)
        if (affectedFace[faceI])
        {
            affectedFace[faceI] = false;

            meshMod.setAction(polyRemoveFace(faceI, -1));
        }
    }


    // Remove points.
    forAllConstIter(labelHashSet, pointsToRemove, iter)
    {
        label pointI = iter.key();

        meshMod.setAction(polyRemovePoint(pointI, -1));
    }


    // Remove cells.
    forAll(cellRegion, cellI)
    {
        label region = cellRegion[cellI];

        if (region != -1 && (cellI != cellRegionMaster[region]))
        {
            meshMod.setAction(polyRemoveCell(cellI, cellRegionMaster[region]));
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
                FatalErrorIn("setRefinement")
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
    forAll(affectedFace, faceI)
    {
        if (affectedFace[faceI])
        {
            affectedFace[faceI] = false;

            face f(filterFace(pointsToRemove, faceI));

            label own = mesh_.faceOwner()[faceI];

            if (cellRegion[own] != -1)
            {
                own = cellRegionMaster[cellRegion[own]];
            }

            label patchID, zoneID, zoneFlip;

            getFaceInfo(faceI, patchID, zoneID, zoneFlip);

            label nei = -1;

            if (mesh_.isInternalFace(faceI))
            {
                nei = mesh_.faceNeighbour()[faceI];

                if (cellRegion[nei] != -1)
                {
                    nei = cellRegionMaster[cellRegion[nei]];
                }
            }

//            if (debug)
//            {
//                Pout<< "Modifying " << faceI
//                    << " old verts:" << mesh_.faces()[faceI]
//                    << " for new verts:" << f
//                    << " or for new owner " << own << " or for new nei "
//                    << nei
//                    << endl;
//            }

            modFace
            (
                f,                  // modified face
                faceI,              // label of face being modified
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
