/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "meshCutAndRemove.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "polyAddFace.H"
#include "polyAddPoint.H"
#include "polyRemovePoint.H"
#include "polyRemoveFace.H"
#include "polyModifyFace.H"
#include "cellCuts.H"
#include "polyTopoChangeMap.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(meshCutAndRemove, 0);
}


// * * * * * * * * * * * * * Private Static Functions  * * * * * * * * * * * //

// Returns -1 or index in elems1 of first shared element.
Foam::label Foam::meshCutAndRemove::firstCommon
(
    const labelList& elems1,
    const labelList& elems2
)
{
    forAll(elems1, elemI)
    {
        label index1 = findIndex(elems2, elems1[elemI]);

        if (index1 != -1)
        {
            return index1;
        }
    }
    return -1;
}


// Check if twoCuts at two consecutive position in cuts.
bool Foam::meshCutAndRemove::isIn
(
    const edge& twoCuts,
    const labelList& cuts
)
{
    label index = findIndex(cuts, twoCuts[0]);

    if (index == -1)
    {
        return false;
    }

    return
    (
        cuts[cuts.fcIndex(index)] == twoCuts[1]
     || cuts[cuts.rcIndex(index)] == twoCuts[1]
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::meshCutAndRemove::findCutCell
(
    const cellCuts& cuts,
    const labelList& cellLabels
) const
{
    forAll(cellLabels, labelI)
    {
        label celli = cellLabels[labelI];

        if (cuts.cellLoops()[celli].size())
        {
            return celli;
        }
    }
    return -1;
}


Foam::label Foam::meshCutAndRemove::findInternalFacePoint
(
    const labelList& pointLabels
) const
{
    forAll(pointLabels, labelI)
    {
        label pointi = pointLabels[labelI];

        const labelList& pFaces = mesh().pointFaces()[pointi];

        forAll(pFaces, pFacei)
        {
            label facei = pFaces[pFacei];

            if (mesh().isInternalFace(facei))
            {
                return pointi;
            }
        }
    }

    if (pointLabels.empty())
    {
        FatalErrorInFunction
            << "Empty pointLabels" << abort(FatalError);
    }

    return -1;
}


Foam::label Foam::meshCutAndRemove::findPatchFacePoint
(
    const face& f,
    const label exposedPatchi
) const
{
    const labelListList& pointFaces = mesh().pointFaces();
    const polyBoundaryMesh& patches = mesh().boundaryMesh();

    forAll(f, fp)
    {
        label pointi = f[fp];

        if (pointi < mesh().nPoints())
        {
            const labelList& pFaces = pointFaces[pointi];

            forAll(pFaces, i)
            {
                if (patches.whichPatch(pFaces[i]) == exposedPatchi)
                {
                    return pointi;
                }
            }
        }
    }
    return -1;
}


void Foam::meshCutAndRemove::faceCells
(
    const cellCuts& cuts,
    const label exposedPatchi,
    const label facei,
    label& own,
    label& nei,
    label& patchID
) const
{
    const labelListList& anchorPts = cuts.cellAnchorPoints();
    const labelListList& cellLoops = cuts.cellLoops();

    const face& f = mesh().faces()[facei];

    own = mesh().faceOwner()[facei];

    if (cellLoops[own].size() && firstCommon(f, anchorPts[own]) == -1)
    {
        // owner has been split and this is the removed part.
        own = -1;
    }

    nei = -1;

    if (mesh().isInternalFace(facei))
    {
        nei = mesh().faceNeighbour()[facei];

        if (cellLoops[nei].size() && firstCommon(f, anchorPts[nei]) == -1)
        {
            nei = -1;
        }
    }

    patchID = mesh().boundaryMesh().whichPatch(facei);

    if (patchID == -1 && (own == -1 || nei == -1))
    {
        // Face was internal but becomes external
        patchID = exposedPatchi;
    }
}


void Foam::meshCutAndRemove::getZoneInfo
(
    const label facei,
    label& zoneID,
    bool& zoneFlip
) const
{
    zoneID = mesh().faceZones().whichZone(facei);

    zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh().faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
    }
}


void Foam::meshCutAndRemove::addFace
(
    polyTopoChange& meshMod,
    const label facei,
    const label masterPointi,
    const face& newFace,
    const label own,
    const label nei,
    const label patchID
)
{
    label zoneID;
    bool zoneFlip;

    getZoneInfo(facei, zoneID, zoneFlip);

    if ((nei == -1) || (own != -1 && own < nei))
    {
        // Ordering ok.
        if (debug & 2)
        {
            Pout<< "Adding face " << newFace
                << " with new owner:" << own
                << " with new neighbour:" << nei
                << " patchID:" << patchID
                << " anchor:" << masterPointi
                << " zoneID:" << zoneID
                << " zoneFlip:" << zoneFlip
                << endl;
        }

        meshMod.setAction
        (
            polyAddFace
            (
                newFace,                    // face
                own,                        // owner
                nei,                        // neighbour
                masterPointi,               // master point
                -1,                         // master edge
                -1,                         // master face for addition
                false,                      // flux flip
                patchID,                    // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face zone flip
            )
        );
    }
    else
    {
        // Reverse owner/neighbour
        if (debug & 2)
        {
            Pout<< "Adding (reversed) face " << newFace.reverseFace()
                << " with new owner:" << nei
                << " with new neighbour:" << own
                << " patchID:" << patchID
                << " anchor:" << masterPointi
                << " zoneID:" << zoneID
                << " zoneFlip:" << zoneFlip
                << endl;
        }

        meshMod.setAction
        (
            polyAddFace
            (
                newFace.reverseFace(),      // face
                nei,                        // owner
                own,                        // neighbour
                masterPointi,               // master point
                -1,                         // master edge
                -1,                         // master face for addition
                false,                      // flux flip
                patchID,                    // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face zone flip
            )
        );
    }
}


// Modifies existing facei for either new owner/neighbour or new face points.
void Foam::meshCutAndRemove::modFace
(
    polyTopoChange& meshMod,
    const label facei,
    const face& newFace,
    const label own,
    const label nei,
    const label patchID
)
{
    label zoneID;
    bool zoneFlip;

    getZoneInfo(facei, zoneID, zoneFlip);

    if
    (
        (own != mesh().faceOwner()[facei])
     || (
            mesh().isInternalFace(facei)
         && (nei != mesh().faceNeighbour()[facei])
        )
     || (newFace != mesh().faces()[facei])
    )
    {
        if (debug & 2)
        {
            Pout<< "Modifying face " << facei
                << " old vertices:" << mesh().faces()[facei]
                << " new vertices:" << newFace
                << " new owner:" << own
                << " new neighbour:" << nei
                << " new patch:" << patchID
                << " new zoneID:" << zoneID
                << " new zoneFlip:" << zoneFlip
                << endl;
        }

        if ((nei == -1) || (own != -1 && own < nei))
        {
            meshMod.setAction
            (
                polyModifyFace
                (
                    newFace,            // modified face
                    facei,              // label of face being modified
                    own,                // owner
                    nei,                // neighbour
                    false,              // face flip
                    patchID,            // patch for face
                    false,              // remove from zone
                    zoneID,             // zone for face
                    zoneFlip            // face flip in zone
                )
            );
        }
        else
        {
            meshMod.setAction
            (
                polyModifyFace
                (
                    newFace.reverseFace(),  // modified face
                    facei,                  // label of face being modified
                    nei,                    // owner
                    own,                    // neighbour
                    false,                  // face flip
                    patchID,                // patch for face
                    false,                  // remove from zone
                    zoneID,                 // zone for face
                    zoneFlip                // face flip in zone
                )
            );
        }
    }
}


void Foam::meshCutAndRemove::copyFace
(
    const face& f,
    const label startFp,
    const label endFp,
    face& newFace
) const
{
    label fp = startFp;

    label newFp = 0;

    while (fp != endFp)
    {
        newFace[newFp++] = f[fp];

        fp = (fp + 1) % f.size();
    }
    newFace[newFp] = f[fp];
}


// Actually split face in two along splitEdge v0, v1 (the two vertices in new
// vertex numbering). Generates faces in same ordering
// as original face. Replaces cutEdges by the points introduced on them
// (addedPoints_).
void Foam::meshCutAndRemove::splitFace
(
    const face& f,
    const label v0,
    const label v1,

    face& f0,
    face& f1
) const
{
    // Check if we find any new vertex which is part of the splitEdge.
    label startFp = findIndex(f, v0);

    if (startFp == -1)
    {
        FatalErrorInFunction
            << "Cannot find vertex (new numbering) " << v0
            << " on face " << f
            << abort(FatalError);
    }

    label endFp = findIndex(f, v1);

    if (endFp == -1)
    {
        FatalErrorInFunction
            << "Cannot find vertex (new numbering) " << v1
            << " on face " << f
            << abort(FatalError);
    }


    f0.setSize((endFp + 1 + f.size() - startFp) % f.size());
    f1.setSize(f.size() - f0.size() + 2);

    copyFace(f, startFp, endFp, f0);
    copyFace(f, endFp, startFp, f1);
}


Foam::face Foam::meshCutAndRemove::addEdgeCutsToFace(const label facei) const
{
    const face& f = mesh().faces()[facei];

    face newFace(2 * f.size());

    label newFp = 0;

    forAll(f, fp)
    {
        // Duplicate face vertex.
        newFace[newFp++] = f[fp];

        // Check if edge has been cut.
        label fp1 = f.fcIndex(fp);

        HashTable<label, edge, Hash<edge>>::const_iterator fnd =
            addedPoints_.find(edge(f[fp], f[fp1]));

        if (fnd != addedPoints_.end())
        {
            // edge has been cut. Introduce new vertex.
            newFace[newFp++] = fnd();
        }
    }

    newFace.setSize(newFp);

    return newFace;
}


// Walk loop (loop of cuts) across circumference of celli. Returns face in
// new vertices.
// Note: tricky bit is that it can use existing edges which have been split.
Foam::face Foam::meshCutAndRemove::loopToFace
(
    const label celli,
    const labelList& loop
) const
{
    face newFace(2*loop.size());

    label newFacei = 0;

    forAll(loop, fp)
    {
        label cut = loop[fp];

        if (isEdge(cut))
        {
            label edgeI = getEdge(cut);

            const edge& e = mesh().edges()[edgeI];

            label vertI = addedPoints_[e];

            newFace[newFacei++] = vertI;
        }
        else
        {
            // cut is vertex.
            label vertI = getVertex(cut);

            newFace[newFacei++] = vertI;

            label nextCut = loop[loop.fcIndex(fp)];

            if (!isEdge(nextCut))
            {
                // From vertex to vertex -> cross cut only if no existing edge.
                label nextVertI = getVertex(nextCut);

                label edgeI = meshTools::findEdge(mesh(), vertI, nextVertI);

                if (edgeI != -1)
                {
                    // Existing edge. Insert split-edge point if any.
                    HashTable<label, edge, Hash<edge>>::const_iterator fnd =
                        addedPoints_.find(mesh().edges()[edgeI]);

                    if (fnd != addedPoints_.end())
                    {
                        newFace[newFacei++] = fnd();
                    }
                }
            }
        }
    }
    newFace.setSize(newFacei);

    return newFace;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::meshCutAndRemove::meshCutAndRemove(const polyMesh& mesh)
:
    edgeVertex(mesh),
    addedFaces_(),
    addedPoints_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meshCutAndRemove::setRefinement
(
    const label exposedPatchi,
    const cellCuts& cuts,
    const labelList& cutPatch,
    polyTopoChange& meshMod
)
{
    // Clear and size maps here since mesh size will change.
    addedFaces_.clear();
    addedFaces_.resize(cuts.nLoops());

    addedPoints_.clear();
    addedPoints_.resize(cuts.nLoops());

    if (cuts.nLoops() == 0)
    {
        return;
    }

    const labelListList& anchorPts = cuts.cellAnchorPoints();
    const labelListList& cellLoops = cuts.cellLoops();
    const polyBoundaryMesh& patches = mesh().boundaryMesh();

    if (exposedPatchi < 0 || exposedPatchi >= patches.size())
    {
        FatalErrorInFunction
            << "Illegal exposed patch " << exposedPatchi
            << abort(FatalError);
    }


    //
    // Add new points along cut edges.
    //

    forAll(cuts.edgeIsCut(), edgeI)
    {
        if (cuts.edgeIsCut()[edgeI])
        {
            const edge& e = mesh().edges()[edgeI];

            // Check if there is any cell using this edge.
            if (debug && findCutCell(cuts, mesh().edgeCells()[edgeI]) == -1)
            {
                FatalErrorInFunction
                    << "Problem: cut edge but none of the cells using it is\n"
                    << "edge:" << edgeI << " verts:" << e
                    << abort(FatalError);
            }

            // One of the edge end points should be master point of nbCelli.
            label masterPointi = e.start();

            const point& v0 = mesh().points()[e.start()];
            const point& v1 = mesh().points()[e.end()];

            scalar weight = cuts.edgeWeight()[edgeI];

            point newPt = weight*v1 + (1.0-weight)*v0;

            label addedPointi =
                meshMod.setAction
                (
                    polyAddPoint
                    (
                        newPt,              // point
                        masterPointi,       // master point
                        -1,                 // zone for point
                        true                // supports a cell
                    )
                );

            // Store on (hash of) edge.
            addedPoints_.insert(e, addedPointi);

            if (debug & 2)
            {
                Pout<< "Added point " << addedPointi
                    << " to vertex "
                    << masterPointi << " of edge " << edgeI
                    << " vertices " << e << endl;
            }
        }
    }


    //
    // Remove all points that will not be used anymore
    //
    {
        boolList usedPoint(mesh().nPoints(), false);

        forAll(cellLoops, celli)
        {
            const labelList& loop = cellLoops[celli];

            if (loop.size())
            {
                // Cell is cut. Uses only anchor points and loop itself.
                forAll(loop, fp)
                {
                    label cut = loop[fp];

                    if (!isEdge(cut))
                    {
                        usedPoint[getVertex(cut)] = true;
                    }
                }

                const labelList& anchors = anchorPts[celli];

                forAll(anchors, i)
                {
                    usedPoint[anchors[i]] = true;
                }
            }
            else
            {
                // Cell is not cut so use all its points
                const labelList& cPoints = mesh().cellPoints()[celli];

                forAll(cPoints, i)
                {
                    usedPoint[cPoints[i]] = true;
                }
            }
        }


        // Check
        const Map<edge>& faceSplitCut = cuts.faceSplitCut();

        forAllConstIter(Map<edge>, faceSplitCut, iter)
        {
            const edge& fCut = iter();

            forAll(fCut, i)
            {
                label cut = fCut[i];

                if (!isEdge(cut))
                {
                    label pointi = getVertex(cut);

                    if (!usedPoint[pointi])
                    {
                        FatalErrorInFunction
                            << "Problem: faceSplitCut not used by any loop"
                            << " or cell anchor point"
                            << "face:" << iter.key() << " point:" << pointi
                            << " coord:" << mesh().points()[pointi]
                            << abort(FatalError);
                    }
                }
            }
        }

        forAll(cuts.pointIsCut(), pointi)
        {
            if (cuts.pointIsCut()[pointi])
            {
                if (!usedPoint[pointi])
                {
                    FatalErrorInFunction
                        << "Problem: point is marked as cut but"
                        << " not used by any loop"
                        << " or cell anchor point"
                        << "point:" << pointi
                        << " coord:" << mesh().points()[pointi]
                        << abort(FatalError);
                }
            }
        }


        // Remove unused points.
        forAll(usedPoint, pointi)
        {
            if (!usedPoint[pointi])
            {
                meshMod.setAction(polyRemovePoint(pointi));

                if (debug & 2)
                {
                    Pout<< "Removing unused point " << pointi << endl;
                }
            }
        }
    }


    //
    // For all cut cells add an internal or external face
    //

    forAll(cellLoops, celli)
    {
        const labelList& loop = cellLoops[celli];

        if (loop.size())
        {
            if (cutPatch[celli] < 0 || cutPatch[celli] >= patches.size())
            {
                FatalErrorInFunction
                    << "Illegal patch " << cutPatch[celli]
                    << " provided for cut cell " << celli
                    << abort(FatalError);
            }

            //
            // Convert loop (=list of cuts) into proper face.
            // cellCuts sets orientation is towards anchor side so reverse.
            //
            face newFace(loopToFace(celli, loop));

            reverse(newFace);

            // Pick any anchor point on cell
            label masterPointi = findPatchFacePoint(newFace, exposedPatchi);

            label addedFacei =
                meshMod.setAction
                (
                    polyAddFace
                    (
                        newFace,                // face
                        celli,                  // owner
                        -1,                     // neighbour
                        masterPointi,           // master point
                        -1,                     // master edge
                        -1,                     // master face for addition
                        false,                  // flux flip
                        cutPatch[celli],        // patch for face
                        -1,                     // zone for face
                        false                   // face zone flip
                    )
                );

            addedFaces_.insert(celli, addedFacei);

            if (debug & 2)
            {
                Pout<< "Added splitting face " << newFace << " index:"
                    << addedFacei << " from masterPoint:" << masterPointi
                    << " to owner " << celli << " with anchors:"
                    << anchorPts[celli]
                    << " from Loop:";

                // Gets edgeweights of loop
                scalarField weights(loop.size());
                forAll(loop, i)
                {
                    label cut = loop[i];

                    weights[i] =
                    (
                        isEdge(cut)
                      ? cuts.edgeWeight()[getEdge(cut)]
                      : -great
                    );
                }
                writeCuts(Pout, loop, weights);
                Pout<< endl;
            }
        }
    }


    //
    // Modify faces to use only anchorpoints and loop points
    // (so throw away part without anchorpoints)
    //


    // Maintain whether face has been updated (for -split edges
    // -new owner/neighbour)
    boolList faceUptodate(mesh().nFaces(), false);

    const Map<edge>& faceSplitCuts = cuts.faceSplitCut();

    forAllConstIter(Map<edge>, faceSplitCuts, iter)
    {
        label facei = iter.key();

        // Renumber face to include split edges.
        face newFace(addEdgeCutsToFace(facei));

        // Edge splitting the face. Convert edge to new vertex numbering.
        const edge& splitEdge = iter();

        label cut0 = splitEdge[0];

        label v0;
        if (isEdge(cut0))
        {
            label edgeI = getEdge(cut0);
            v0 = addedPoints_[mesh().edges()[edgeI]];
        }
        else
        {
            v0 = getVertex(cut0);
        }

        label cut1 = splitEdge[1];
        label v1;
        if (isEdge(cut1))
        {
            label edgeI = getEdge(cut1);
            v1 = addedPoints_[mesh().edges()[edgeI]];
        }
        else
        {
            v1 = getVertex(cut1);
        }

        // Split face along the elements of the splitEdge.
        face f0, f1;
        splitFace(newFace, v0, v1, f0, f1);

        label own = mesh().faceOwner()[facei];

        label nei = -1;

        if (mesh().isInternalFace(facei))
        {
            nei = mesh().faceNeighbour()[facei];
        }

        if (debug & 2)
        {
            Pout<< "Split face " << mesh().faces()[facei]
                << " own:" << own << " nei:" << nei
                << " into f0:" << f0
                << " and f1:" << f1 << endl;
        }


        // Check which cell using face uses anchorPoints (so is kept)
        // and which one doesn't (gets removed)

        // Bit tricky. We have to know whether this faceSplit splits owner/
        // neighbour or both. Even if cell is cut we have to make sure this is
        // the one that cuts it (this face cut might not be the one splitting
        // the cell)
        // The face f gets split into two parts, f0 and f1.
        // Each of these can have a different owner and or neighbour.

        const face& f = mesh().faces()[facei];

        label f0Own = -1;
        label f1Own = -1;

        if (cellLoops[own].empty())
        {
            // Owner side is not split so keep both halves.
            f0Own = own;
            f1Own = own;
        }
        else if (isIn(splitEdge, cellLoops[own]))
        {
            // Owner is cut by this splitCut. See which of f0, f1 gets
            // preserved and becomes owner, and which gets removed.
            if (firstCommon(f0, anchorPts[own]) != -1)
            {
                // f0 preserved so f1 gets deleted
                f0Own = own;
                f1Own = -1;
            }
            else
            {
                f0Own = -1;
                f1Own = own;
            }
        }
        else
        {
            // Owner not cut by this splitCut but by another.
            // Check on original face whether
            // use anchorPts.
            if (firstCommon(f, anchorPts[own]) != -1)
            {
                // both f0 and f1 owner side preserved
                f0Own = own;
                f1Own = own;
            }
            else
            {
                // both f0 and f1 owner side removed
                f0Own = -1;
                f1Own = -1;
            }
        }


        label f0Nei = -1;
        label f1Nei = -1;

        if (nei != -1)
        {
            if (cellLoops[nei].empty())
            {
                f0Nei = nei;
                f1Nei = nei;
            }
            else if (isIn(splitEdge, cellLoops[nei]))
            {
                // Neighbour is cut by this splitCut. So anchor part of it
                // gets kept, non-anchor bit gets removed. See which of f0, f1
                // connects to which part.

                if (firstCommon(f0, anchorPts[nei]) != -1)
                {
                    f0Nei = nei;
                    f1Nei = -1;
                }
                else
                {
                    f0Nei = -1;
                    f1Nei = nei;
                }
            }
            else
            {
                // neighbour not cut by this splitCut. Check on original face
                // whether use anchorPts.

                if (firstCommon(f, anchorPts[nei]) != -1)
                {
                    f0Nei = nei;
                    f1Nei = nei;
                }
                else
                {
                    // both f0 and f1 on neighbour side removed
                    f0Nei = -1;
                    f1Nei = -1;
                }
            }
        }


        if (debug & 2)
        {
            Pout<< "f0 own:" << f0Own << " nei:" << f0Nei
                << "  f1 own:" << f1Own << " nei:" << f1Nei
                << endl;
        }


        // If faces were internal but now become external set a patch.
        // If they were external already keep the patch.
        label patchID = patches.whichPatch(facei);

        if (patchID == -1)
        {
            patchID = exposedPatchi;
        }


        // Do as much as possible by modifying facei. Delay any remove
        // face. Keep track of whether facei has been used.

        bool modifiedFacei = false;

        if (f0Own == -1)
        {
            if (f0Nei != -1)
            {
                // f0 becomes external face (note:modFace will reverse face)
                modFace(meshMod, facei, f0, f0Own, f0Nei, patchID);
                modifiedFacei = true;
            }
        }
        else
        {
            if (f0Nei == -1)
            {
                // f0 becomes external face
                modFace(meshMod, facei, f0, f0Own, f0Nei, patchID);
                modifiedFacei = true;
            }
            else
            {
                // f0 stays internal face.
                modFace(meshMod, facei, f0, f0Own, f0Nei, -1);
                modifiedFacei = true;
            }
        }


        // f1 is added face (if at all)

        if (f1Own == -1)
        {
            if (f1Nei == -1)
            {
                // f1 not needed.
            }
            else
            {
                // f1 becomes external face (note:modFace will reverse face)
                if (!modifiedFacei)
                {
                    modFace(meshMod, facei, f1, f1Own, f1Nei, patchID);
                    modifiedFacei = true;
                }
                else
                {
                    label masterPointi = findPatchFacePoint(f1, patchID);

                    addFace
                    (
                        meshMod,
                        facei,          // face for zone info
                        masterPointi,   // inflation point
                        f1,             // vertices of face
                        f1Own,
                        f1Nei,
                        patchID         // patch for new face
                    );
                }
            }
        }
        else
        {
            if (f1Nei == -1)
            {
                // f1 becomes external face
                if (!modifiedFacei)
                {
                    modFace(meshMod, facei, f1, f1Own, f1Nei, patchID);
                    modifiedFacei = true;
                }
                else
                {
                    label masterPointi = findPatchFacePoint(f1, patchID);

                    addFace
                    (
                        meshMod,
                        facei,
                        masterPointi,
                        f1,
                        f1Own,
                        f1Nei,
                        patchID
                    );
                }
            }
            else
            {
                // f1 is internal face.
                if (!modifiedFacei)
                {
                    modFace(meshMod, facei, f1, f1Own, f1Nei, -1);
                    modifiedFacei = true;
                }
                else
                {
                    label masterPointi = findPatchFacePoint(f1, -1);

                    addFace(meshMod, facei, masterPointi, f1, f1Own, f1Nei, -1);
                }
            }
        }

        if (f0Own == -1 && f0Nei == -1 && !modifiedFacei)
        {
            meshMod.setAction(polyRemoveFace(facei));

            if (debug & 2)
            {
                Pout<< "Removed face " << facei << endl;
            }
        }

        faceUptodate[facei] = true;
    }


    //
    // Faces that have not been split but just appended to. Are guaranteed
    // to be reachable from an edgeCut.
    //

    const boolList& edgeIsCut = cuts.edgeIsCut();

    forAll(edgeIsCut, edgeI)
    {
        if (edgeIsCut[edgeI])
        {
            const labelList& eFaces = mesh().edgeFaces()[edgeI];

            forAll(eFaces, i)
            {
                label facei = eFaces[i];

                if (!faceUptodate[facei])
                {
                    // So the face has not been split itself (i.e. its owner
                    // or neighbour have not been split) so it only
                    // borders by edge a cell which has been split.

                    // Get (new or original) owner and neighbour of facei
                    label own, nei, patchID;
                    faceCells(cuts, exposedPatchi, facei, own, nei, patchID);


                    if (own == -1 && nei == -1)
                    {
                        meshMod.setAction(polyRemoveFace(facei));

                        if (debug & 2)
                        {
                            Pout<< "Removed face " << facei << endl;
                        }
                    }
                    else
                    {
                        // Renumber face to include split edges.
                        face newFace(addEdgeCutsToFace(facei));

                        if (debug & 2)
                        {
                            Pout<< "Added edge cuts to face " << facei
                                << " f:" << mesh().faces()[facei]
                                << " newFace:" << newFace << endl;
                        }

                        modFace
                        (
                            meshMod,
                            facei,
                            newFace,
                            own,
                            nei,
                            patchID
                        );
                    }

                    faceUptodate[facei] = true;
                }
            }
        }
    }


    //
    // Remove any faces on the non-anchor side of a split cell.
    // Note: could loop through all cut cells only and check their faces but
    //       looping over all faces is cleaner and probably faster for dense
    //       cut patterns.

    const faceList& faces = mesh().faces();

    forAll(faces, facei)
    {
        if (!faceUptodate[facei])
        {
            // Get (new or original) owner and neighbour of facei
            label own, nei, patchID;
            faceCells(cuts, exposedPatchi, facei, own, nei, patchID);

            if (own == -1 && nei == -1)
            {
                meshMod.setAction(polyRemoveFace(facei));

                if (debug & 2)
                {
                    Pout<< "Removed face " << facei << endl;
                }
            }
            else
            {
                modFace(meshMod, facei, faces[facei], own, nei, patchID);
            }

            faceUptodate[facei] = true;
        }
    }

    if (debug)
    {
        Pout<< "meshCutAndRemove:" << nl
            << "    cells split:" << cuts.nLoops() << nl
            << "    faces added:" << addedFaces_.size() << nl
            << "    points added on edges:" << addedPoints_.size() << nl
            << endl;
    }
}


void Foam::meshCutAndRemove::topoChange(const polyTopoChangeMap& map)
{
    // Update stored labels for mesh change.
    {
        Map<label> newAddedFaces(addedFaces_.size());

        forAllConstIter(Map<label>, addedFaces_, iter)
        {
            label celli = iter.key();
            label newCelli = map.reverseCellMap()[celli];

            label addedFacei = iter();

            label newAddedFacei = map.reverseFaceMap()[addedFacei];

            if ((newCelli >= 0) && (newAddedFacei >= 0))
            {
                if
                (
                    (debug & 2)
                 && (newCelli != celli || newAddedFacei != addedFacei)
                )
                {
                    Pout<< "meshCutAndRemove::topoChange :"
                        << " updating addedFace for cell " << celli
                        << " from " << addedFacei
                        << " to " << newAddedFacei
                        << endl;
                }
                newAddedFaces.insert(newCelli, newAddedFacei);
            }
        }

        // Copy
        addedFaces_.transfer(newAddedFaces);
    }

    {
        HashTable<label, edge, Hash<edge>> newAddedPoints(addedPoints_.size());

        for
        (
            HashTable<label, edge, Hash<edge>>::const_iterator iter =
                addedPoints_.begin();
            iter != addedPoints_.end();
            ++iter
        )
        {
            const edge& e = iter.key();

            label newStart = map.reversePointMap()[e.start()];

            label newEnd = map.reversePointMap()[e.end()];

            label addedPointi = iter();

            label newAddedPointi = map.reversePointMap()[addedPointi];

            if ((newStart >= 0) && (newEnd >= 0) && (newAddedPointi >= 0))
            {
                edge newE = edge(newStart, newEnd);

                if
                (
                    (debug & 2)
                 && (e != newE || newAddedPointi != addedPointi)
                )
                {
                    Pout<< "meshCutAndRemove::topoChange :"
                        << " updating addedPoints for edge " << e
                        << " from " << addedPointi
                        << " to " << newAddedPointi
                        << endl;
                }

                newAddedPoints.insert(newE, newAddedPointi);
            }
        }

        // Copy
        addedPoints_.transfer(newAddedPoints);
    }
}


// ************************************************************************* //
