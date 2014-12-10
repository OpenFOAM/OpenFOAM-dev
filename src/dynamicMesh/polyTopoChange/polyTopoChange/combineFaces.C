/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "combineFaces.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "polyRemoveFace.H"
#include "polyAddFace.H"
#include "polyModifyFace.H"
#include "polyRemovePoint.H"
#include "polyAddPoint.H"
#include "syncTools.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(combineFaces, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Test face for (almost) convexeness. Allows certain convexity before
// complaining.
// For every two consecutive edges calculate the normal. If it differs too
// much from the average face normal complain.
bool Foam::combineFaces::convexFace
(
    const scalar minConcaveCos,
    const pointField& points,
    const face& f
)
{
    // Get outwards pointing normal of f.
    vector n = f.normal(points);
    n /= mag(n);

    // Get edge from f[0] to f[size-1];
    vector ePrev(points[f.first()] - points[f.last()]);
    scalar magEPrev = mag(ePrev);
    ePrev /= magEPrev + VSMALL;

    forAll(f, fp0)
    {
        // Get vertex after fp
        label fp1 = f.fcIndex(fp0);

        // Normalized vector between two consecutive points
        vector e10(points[f[fp1]] - points[f[fp0]]);
        scalar magE10 = mag(e10);
        e10 /= magE10 + VSMALL;

        if (magEPrev > SMALL && magE10 > SMALL)
        {
            vector edgeNormal = ePrev ^ e10;

            if ((edgeNormal & n) < 0)
            {
                // Concave. Check angle.
                if ((ePrev & e10) < minConcaveCos)
                {
                    return false;
                }
            }
        }

        ePrev = e10;
        magEPrev = magE10;
    }

    // Not a single internal angle is concave so face is convex.
    return true;
}


// Determines if set of faces is valid to collapse into single face.
bool Foam::combineFaces::validFace
(
    const scalar minConcaveCos,
    const indirectPrimitivePatch& bigFace
)
{
    // Get outside vertices (in local vertex numbering)
    const labelListList& edgeLoops = bigFace.edgeLoops();

    if (edgeLoops.size() > 1)
    {
        return false;
    }

    bool isNonManifold = bigFace.checkPointManifold(false, NULL);
    if (isNonManifold)
    {
        return false;
    }

    // Check for convexness
    face f(getOutsideFace(bigFace));

    return convexFace(minConcaveCos, bigFace.points(), f);
}


void Foam::combineFaces::regioniseFaces
(
    const scalar minCos,
    const label cellI,
    const labelList& cEdges,
    Map<label>& faceRegion
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(cEdges, i)
    {
        label edgeI = cEdges[i];

        label f0, f1;
        meshTools::getEdgeFaces(mesh_, cellI, edgeI, f0, f1);

        label p0 = patches.whichPatch(f0);
        label p1 = patches.whichPatch(f1);

        // Face can be merged if
        // - same non-coupled patch
        // - small angle
        if (p0 != -1 && p0 == p1 && !patches[p0].coupled())
        {
            vector f0Normal = mesh_.faceAreas()[f0];
            f0Normal /= mag(f0Normal);
            vector f1Normal = mesh_.faceAreas()[f1];
            f1Normal /= mag(f1Normal);

            if ((f0Normal&f1Normal) > minCos)
            {
                Map<label>::const_iterator f0Fnd = faceRegion.find(f0);

                label region0 = -1;
                if (f0Fnd != faceRegion.end())
                {
                    region0 = f0Fnd();
                }

                Map<label>::const_iterator f1Fnd = faceRegion.find(f1);

                label region1 = -1;
                if (f1Fnd != faceRegion.end())
                {
                    region1 = f1Fnd();
                }

                if (region0 == -1)
                {
                    if (region1 == -1)
                    {
                        label useRegion = faceRegion.size();
                        faceRegion.insert(f0, useRegion);
                        faceRegion.insert(f1, useRegion);
                    }
                    else
                    {
                        faceRegion.insert(f0, region1);
                    }
                }
                else
                {
                    if (region1 == -1)
                    {
                        faceRegion.insert(f1, region0);
                    }
                    else if (region0 != region1)
                    {
                        // Merge the two regions
                        label useRegion = min(region0, region1);
                        label freeRegion = max(region0, region1);

                        forAllIter(Map<label>, faceRegion, iter)
                        {
                            if (iter() == freeRegion)
                            {
                                iter() = useRegion;
                            }
                        }
                    }
                }
            }
        }
    }
}


bool Foam::combineFaces::faceNeighboursValid
(
    const label cellI,
    const Map<label>& faceRegion
) const
{
    if (faceRegion.size() <= 1)
    {
        return true;
    }

    const cell& cFaces = mesh_.cells()[cellI];

    DynamicList<label> storage;

    // Test for face collapsing to edge since too many neighbours merged.
    forAll(cFaces, cFaceI)
    {
        label faceI = cFaces[cFaceI];

        if (!faceRegion.found(faceI))
        {
            const labelList& fEdges = mesh_.faceEdges(faceI, storage);

            // Count number of remaining faces neighbouring faceI. This has
            // to be 3 or more.

            // Unregioned neighbouring faces
            DynamicList<label> neighbourFaces(cFaces.size());
            // Regioned neighbouring faces
            labelHashSet neighbourRegions;

            forAll(fEdges, i)
            {
                label edgeI = fEdges[i];
                label nbrI = meshTools::otherFace(mesh_, cellI, faceI, edgeI);

                Map<label>::const_iterator iter = faceRegion.find(nbrI);

                if (iter == faceRegion.end())
                {
                    if (findIndex(neighbourFaces, nbrI) == -1)
                    {
                        neighbourFaces.append(nbrI);
                    }
                }
                else
                {
                    neighbourRegions.insert(iter());
                }
            }

            if ((neighbourFaces.size()+neighbourRegions.size()) < 3)
            {
                return false;
            }
        }
    }
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::combineFaces::combineFaces
(
    const polyMesh& mesh,
    const bool undoable
)
:
    mesh_(mesh),
    undoable_(undoable),
    masterFace_(0),
    faceSetsVertices_(0),
    savedPointLabels_(0),
    savedPoints_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelListList Foam::combineFaces::getMergeSets
(
    const scalar featureCos,
    const scalar minConcaveCos,
    const labelHashSet& boundaryCells
) const
{
    // Lists of faces that can be merged.
    DynamicList<labelList> allFaceSets(boundaryCells.size() / 10);
    // Storage for on-the-fly cell-edge addressing.
    DynamicList<label> storage;

    // On all cells regionise the faces
    forAllConstIter(labelHashSet, boundaryCells, iter)
    {
        label cellI = iter.key();

        const cell& cFaces = mesh_.cells()[cellI];

        const labelList& cEdges = mesh_.cellEdges(cellI, storage);

        // Region per face
        Map<label> faceRegion(cFaces.size());
        regioniseFaces(featureCos, cellI, cEdges, faceRegion);

        // Now we have in faceRegion for every face the region with planar
        // face sharing the same region. We now check whether the resulting
        // sets cause a face
        // - to become a set of edges since too many faces are merged.
        // - to become convex

        if (faceNeighboursValid(cellI, faceRegion))
        {
            // Create region-to-faces addressing
            Map<labelList> regionToFaces(faceRegion.size());

            forAllConstIter(Map<label>, faceRegion, iter)
            {
                label faceI = iter.key();
                label region = iter();

                Map<labelList>::iterator regionFnd = regionToFaces.find(region);

                if (regionFnd != regionToFaces.end())
                {
                    labelList& setFaces = regionFnd();
                    label sz = setFaces.size();
                    setFaces.setSize(sz+1);
                    setFaces[sz] = faceI;
                }
                else
                {
                    regionToFaces.insert(region, labelList(1, faceI));
                }
            }

            // For every set check if it forms a valid convex face

            forAllConstIter(Map<labelList>, regionToFaces, iter)
            {
                // Make face out of setFaces
                indirectPrimitivePatch bigFace
                (
                    IndirectList<face>
                    (
                        mesh_.faces(),
                        iter()
                    ),
                    mesh_.points()
                );

                // Only store if -only one outside loop -which forms convex face
                if (validFace(minConcaveCos, bigFace))
                {
                    allFaceSets.append(iter());
                }
            }
        }
    }

    return allFaceSets.shrink();
}


Foam::labelListList Foam::combineFaces::getMergeSets
(
    const scalar featureCos,
    const scalar minConcaveCos
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Pick up all cells on boundary
    labelHashSet boundaryCells(mesh_.nFaces()-mesh_.nInternalFaces());

    forAll(patches, patchI)
    {
        const polyPatch& patch = patches[patchI];

        if (!patch.coupled())
        {
            forAll(patch, i)
            {
                boundaryCells.insert(mesh_.faceOwner()[patch.start()+i]);
            }
        }
    }

    return getMergeSets(featureCos, minConcaveCos, boundaryCells);
}


// Gets outside edgeloop as a face
// - in same order as faces
// - in mesh vertex labels
Foam::face Foam::combineFaces::getOutsideFace
(
    const indirectPrimitivePatch& fp
)
{
    if (fp.edgeLoops().size() != 1)
    {
        FatalErrorIn
        (
            "combineFaces::getOutsideFace(const indirectPrimitivePatch&)"
        )   << "Multiple outside loops:" << fp.edgeLoops()
            << abort(FatalError);
    }

    // Get first boundary edge. Since guaranteed one edgeLoop when in here this
    // edge must be on it.
    label bEdgeI = fp.nInternalEdges();

    const edge& e = fp.edges()[bEdgeI];

    const labelList& eFaces = fp.edgeFaces()[bEdgeI];

    if (eFaces.size() != 1)
    {
        FatalErrorIn
        (
            "combineFaces::getOutsideFace(const indirectPrimitivePatch&)"
        )   << "boundary edge:" << bEdgeI
            << " points:" << fp.meshPoints()[e[0]]
            << ' ' << fp.meshPoints()[e[1]]
            << " on indirectPrimitivePatch has " << eFaces.size()
            << " faces using it" << abort(FatalError);
    }


    // Outside loop
    const labelList& outsideLoop = fp.edgeLoops()[0];


    // Get order of edge e in outside loop
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // edgeLoopConsistent : true  = edge in same order as outsideloop
    //                      false = opposite order
    bool edgeLoopConsistent = false;

    {
        label index0 = findIndex(outsideLoop, e[0]);
        label index1 = findIndex(outsideLoop, e[1]);

        if (index0 == -1 || index1 == -1)
        {
            FatalErrorIn
            (
                "combineFaces::getOutsideFace"
                "(const indirectPrimitivePatch&)"
            )   << "Cannot find boundary edge:" << e
                << " points:" << fp.meshPoints()[e[0]]
                << ' ' << fp.meshPoints()[e[1]]
                << " in edgeLoop:" << outsideLoop << abort(FatalError);
        }
        else if (index1 == outsideLoop.fcIndex(index0))
        {
            edgeLoopConsistent = true;
        }
        else if (index0 == outsideLoop.fcIndex(index1))
        {
            edgeLoopConsistent = false;
        }
        else
        {
            FatalErrorIn
            (
                "combineFaces::getOutsideFace"
                "(const indirectPrimitivePatch&)"
            )   << "Cannot find boundary edge:" << e
                << " points:" << fp.meshPoints()[e[0]]
                << ' ' << fp.meshPoints()[e[1]]
                << " on consecutive points in edgeLoop:"
                << outsideLoop << abort(FatalError);
        }
    }


    // Get face in local vertices
    const face& localF = fp.localFaces()[eFaces[0]];

    // Get order of edge in localF
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // faceEdgeConsistent : true  = edge in same order as localF
    //                      false = opposite order
    bool faceEdgeConsistent = false;

    {
        // Find edge in face.
        label index = findIndex(fp.faceEdges()[eFaces[0]], bEdgeI);

        if (index == -1)
        {
            FatalErrorIn
            (
                "combineFaces::getOutsideFace"
                "(const indirectPrimitivePatch&)"
            )   << "Cannot find boundary edge:" << e
                << " points:" << fp.meshPoints()[e[0]]
                << ' ' << fp.meshPoints()[e[1]]
                << " in face:" << eFaces[0]
                << " edges:" << fp.faceEdges()[eFaces[0]]
                << abort(FatalError);
        }
        else if (localF[index] == e[0] && localF.nextLabel(index) == e[1])
        {
            faceEdgeConsistent = true;
        }
        else if (localF[index] == e[1] && localF.nextLabel(index) == e[0])
        {
            faceEdgeConsistent = false;
        }
        else
        {
            FatalErrorIn
            (
                "combineFaces::getOutsideFace"
                "(const indirectPrimitivePatch&)"
            )   << "Cannot find boundary edge:" << e
                << " points:" << fp.meshPoints()[e[0]]
                << ' ' << fp.meshPoints()[e[1]]
                << " in face:" << eFaces[0] << " verts:" << localF
                << abort(FatalError);
        }
    }

    // Get face in mesh points.
    face meshFace(renumber(fp.meshPoints(), outsideLoop));

    if (faceEdgeConsistent != edgeLoopConsistent)
    {
        reverse(meshFace);
    }
    return meshFace;
}


void Foam::combineFaces::setRefinement
(
    const labelListList& faceSets,
    polyTopoChange& meshMod
)
{
    if (undoable_)
    {
        masterFace_.setSize(faceSets.size());
        faceSetsVertices_.setSize(faceSets.size());
        forAll(faceSets, setI)
        {
            const labelList& setFaces = faceSets[setI];

            masterFace_[setI] = setFaces[0];
            faceSetsVertices_[setI] = UIndirectList<face>
            (
                mesh_.faces(),
                setFaces
            );
        }
    }

    // Running count of number of faces using a point
    labelList nPointFaces(mesh_.nPoints(), 0);

    const labelListList& pointFaces = mesh_.pointFaces();

    forAll(pointFaces, pointI)
    {
        nPointFaces[pointI] = pointFaces[pointI].size();
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(faceSets, setI)
    {
        const labelList& setFaces = faceSets[setI];

        // Check
        if (debug)
        {
            forAll(setFaces, i)
            {
                label patchI = patches.whichPatch(setFaces[i]);

                if (patchI == -1 || patches[patchI].coupled())
                {
                    FatalErrorIn
                    (
                        "combineFaces::setRefinement"
                        "(const bool, const labelListList&"
                        ", polyTopoChange&)"
                    )   << "Can only merge non-coupled boundary faces"
                        << " but found internal or coupled face:"
                        << setFaces[i] << " in set " << setI
                        << abort(FatalError);
                }
            }
        }

        // Make face out of setFaces
        indirectPrimitivePatch bigFace
        (
            IndirectList<face>
            (
                mesh_.faces(),
                setFaces
            ),
            mesh_.points()
        );

        // Get outside vertices (in local vertex numbering)
        const labelListList& edgeLoops = bigFace.edgeLoops();

        if (edgeLoops.size() != 1)
        {
            FatalErrorIn
            (
                "combineFaces::setRefinement"
                "(const bool, const labelListList&, polyTopoChange&)"
            )   << "Faces to-be-merged " << setFaces
                << " do not form a single big face." << nl
                << abort(FatalError);
        }


        // Delete all faces except master, whose face gets modified.

        // Modify master face
        // ~~~~~~~~~~~~~~~~~~

        label masterFaceI = setFaces[0];

        // Get outside face in mesh vertex labels
        face outsideFace(getOutsideFace(bigFace));

        label zoneID = mesh_.faceZones().whichZone(masterFaceI);

        bool zoneFlip = false;

        if (zoneID >= 0)
        {
            const faceZone& fZone = mesh_.faceZones()[zoneID];

            zoneFlip = fZone.flipMap()[fZone.whichFace(masterFaceI)];
        }

        label patchI = mesh_.boundaryMesh().whichPatch(masterFaceI);

        meshMod.setAction
        (
            polyModifyFace
            (
                outsideFace,                    // modified face
                masterFaceI,                    // label of face being modified
                mesh_.faceOwner()[masterFaceI], // owner
                -1,                             // neighbour
                false,                          // face flip
                patchI,                         // patch for face
                false,                          // remove from zone
                zoneID,                         // zone for face
                zoneFlip                        // face flip in zone
            )
        );


        // Delete all non-master faces
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        for (label i = 1; i < setFaces.size(); i++)
        {
            meshMod.setAction(polyRemoveFace(setFaces[i]));
        }


        // Mark unused points for deletion
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Remove count of all faces combined
        forAll(setFaces, i)
        {
            const face& f = mesh_.faces()[setFaces[i]];

            forAll(f, fp)
            {
                nPointFaces[f[fp]]--;
            }
        }
        // But keep count on new outside face
        forAll(outsideFace, fp)
        {
            nPointFaces[outsideFace[fp]]++;
        }
    }


    // If one or more processors use the point it needs to be kept.
    syncTools::syncPointList
    (
        mesh_,
        nPointFaces,
        plusEqOp<label>(),
        label(0)            // null value
    );

    // Remove all unused points. Store position if undoable.
    if (!undoable_)
    {
        forAll(nPointFaces, pointI)
        {
            if (nPointFaces[pointI] == 0)
            {
                meshMod.setAction(polyRemovePoint(pointI));
            }
        }
    }
    else
    {
        // Count removed points
        label n = 0;
        forAll(nPointFaces, pointI)
        {
            if (nPointFaces[pointI] == 0)
            {
                n++;
            }
        }

        savedPointLabels_.setSize(n);
        savedPoints_.setSize(n);
        Map<label> meshToSaved(2*n);

        // Remove points and store position
        n = 0;
        forAll(nPointFaces, pointI)
        {
            if (nPointFaces[pointI] == 0)
            {
                meshMod.setAction(polyRemovePoint(pointI));

                savedPointLabels_[n] = pointI;
                savedPoints_[n] = mesh_.points()[pointI];

                meshToSaved.insert(pointI, n);
                n++;
            }
        }

        // Update stored vertex labels. Negative indices index into local points
        forAll(faceSetsVertices_, setI)
        {
            faceList& setFaces = faceSetsVertices_[setI];

            forAll(setFaces, i)
            {
                face& f = setFaces[i];

                forAll(f, fp)
                {
                    label pointI = f[fp];

                    if (nPointFaces[pointI] == 0)
                    {
                        f[fp] = -meshToSaved[pointI]-1;
                    }
                }
            }
        }
    }
}


void Foam::combineFaces::updateMesh(const mapPolyMesh& map)
{
    if (undoable_)
    {
        // Master face just renumbering of point labels
        inplaceRenumber(map.reverseFaceMap(), masterFace_);

        // Stored faces refer to backed-up vertices (not changed)
        // and normal vertices (need to be renumbered)
        forAll(faceSetsVertices_, setI)
        {
            faceList& faces = faceSetsVertices_[setI];

            forAll(faces, i)
            {
                // Note: faces[i] can have negative elements.
                face& f = faces[i];

                forAll(f, fp)
                {
                    label pointI = f[fp];

                    if (pointI >= 0)
                    {
                        f[fp] = map.reversePointMap()[pointI];

                        if (f[fp] < 0)
                        {
                            FatalErrorIn
                            (
                                "combineFaces::updateMesh"
                                "(const mapPolyMesh&)"
                            )   << "In set " << setI << " at position " << i
                                << " with master face "
                                << masterFace_[setI] << nl
                                << "the points of the slave face " << faces[i]
                                << " don't exist anymore."
                                << abort(FatalError);
                        }
                    }
                }
            }
        }
    }
}


// Note that faces on coupled patches are never combined (or at least
// getMergeSets never picks them up. Thus the points on them will never get
// deleted since still used by the coupled faces. So never a risk of undoing
// things (faces/points) on coupled patches.
void Foam::combineFaces::setUnrefinement
(
    const labelList& masterFaces,
    polyTopoChange& meshMod,
    Map<label>& restoredPoints,
    Map<label>& restoredFaces,
    Map<label>& restoredCells
)
{
    if (!undoable_)
    {
        FatalErrorIn
        (
            "combineFaces::setUnrefinement"
            "(const labelList&, polyTopoChange&"
            ", Map<label>&, Map<label>&, Map<label>&)"
        )   << "Can only call setUnrefinement if constructed with"
            << " unrefinement capability." << exit(FatalError);
    }


    // Restored points
    labelList addedPoints(savedPoints_.size(), -1);

    // Invert set-to-master-face
    Map<label> masterToSet(masterFace_.size());

    forAll(masterFace_, setI)
    {
        if (masterFace_[setI] >= 0)
        {
            masterToSet.insert(masterFace_[setI], setI);
        }
    }

    forAll(masterFaces, i)
    {
        label masterFaceI = masterFaces[i];

        Map<label>::const_iterator iter = masterToSet.find(masterFaceI);

        if (iter == masterToSet.end())
        {
            FatalErrorIn
            (
                "combineFaces::setUnrefinement"
                "(const labelList&, polyTopoChange&"
                ", Map<label>&, Map<label>&, Map<label>&)"
            )   << "Master face " << masterFaceI
                << " is not the master of one of the merge sets"
                << " or has already been merged"
                << abort(FatalError);
        }

        label setI = iter();


        // Update faces of the merge setI for reintroduced vertices
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        faceList& faces = faceSetsVertices_[setI];

        if (faces.empty())
        {
            FatalErrorIn
            (
                "combineFaces::setUnrefinement"
                "(const labelList&, polyTopoChange&"
                ", Map<label>&, Map<label>&, Map<label>&)"
            )   << "Set " << setI << " with master face " << masterFaceI
                << " has already been merged." << abort(FatalError);
        }

        forAll(faces, i)
        {
            face& f = faces[i];

            forAll(f, fp)
            {
                label pointI = f[fp];

                if (pointI < 0)
                {
                    label localI = -pointI-1;

                    if (addedPoints[localI] == -1)
                    {
                        // First occurrence of saved point. Reintroduce point
                        addedPoints[localI] = meshMod.setAction
                        (
                            polyAddPoint
                            (
                                savedPoints_[localI],   // point
                                -1,                     // master point
                                -1,                     // zone for point
                                true                    // supports a cell
                            )
                        );
                        restoredPoints.insert
                        (
                            addedPoints[localI],        // current point label
                            savedPointLabels_[localI]   // point label when it
                                                        // was stored
                        );
                    }
                    f[fp] = addedPoints[localI];
                }
            }
        }


        // Restore
        // ~~~~~~~

        label own = mesh_.faceOwner()[masterFaceI];
        label zoneID = mesh_.faceZones().whichZone(masterFaceI);
        bool zoneFlip = false;
        if (zoneID >= 0)
        {
            const faceZone& fZone = mesh_.faceZones()[zoneID];
            zoneFlip = fZone.flipMap()[fZone.whichFace(masterFaceI)];
        }
        label patchI = mesh_.boundaryMesh().whichPatch(masterFaceI);

        if (mesh_.boundaryMesh()[patchI].coupled())
        {
            FatalErrorIn
            (
                "combineFaces::setUnrefinement"
                "(const labelList&, polyTopoChange&"
                ", Map<label>&, Map<label>&, Map<label>&)"
            )   << "Master face " << masterFaceI << " is on coupled patch "
                << mesh_.boundaryMesh()[patchI].name()
                << abort(FatalError);
        }

        //Pout<< "Restoring new master face " << masterFaceI
        //    << " to vertices " << faces[0] << endl;

        // Modify the master face.
        meshMod.setAction
        (
            polyModifyFace
            (
                faces[0],                       // original face
                masterFaceI,                    // label of face
                own,                            // owner
                -1,                             // neighbour
                false,                          // face flip
                patchI,                         // patch for face
                false,                          // remove from zone
                zoneID,                         // zone for face
                zoneFlip                        // face flip in zone
            )
        );
        restoredFaces.insert(masterFaceI, masterFaceI);

        // Add the previously removed faces
        for (label i = 1; i < faces.size(); i++)
        {
            //Pout<< "Restoring removed face with vertices " << faces[i]
            //    << endl;

            label faceI = meshMod.setAction
            (
                polyAddFace
                (
                    faces[i],        // vertices
                    own,                    // owner,
                    -1,                     // neighbour,
                    -1,                     // masterPointID,
                    -1,                     // masterEdgeID,
                    masterFaceI,             // masterFaceID,
                    false,                  // flipFaceFlux,
                    patchI,                 // patchID,
                    zoneID,                 // zoneID,
                    zoneFlip                // zoneFlip
                )
            );
            restoredFaces.insert(faceI, masterFaceI);
        }

        // Clear out restored set
        faceSetsVertices_[setI].clear();
        masterFace_[setI] = -1;
    }
}


// ************************************************************************* //
