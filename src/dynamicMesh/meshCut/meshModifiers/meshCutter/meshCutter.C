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

#include "meshCutter.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "cellCuts.H"
#include "mapPolyMesh.H"
#include "meshTools.H"
#include "polyModifyFace.H"
#include "polyAddPoint.H"
#include "polyAddFace.H"
#include "polyAddCell.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(meshCutter, 0);
}


// * * * * * * * * * * * * * Private Static Functions  * * * * * * * * * * * //

// Returns true if lst1 and lst2 share elements
bool Foam::meshCutter::uses(const labelList& elems1, const labelList& elems2)
{
    forAll(elems1, elemI)
    {
        if (findIndex(elems2, elems1[elemI]) != -1)
        {
            return true;
        }
    }
    return false;
}


// Check if twoCuts at two consecutive position in cuts.
bool Foam::meshCutter::isIn
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

// Returns the cell in cellLabels that is cut. Or -1.
Foam::label Foam::meshCutter::findCutCell
(
    const cellCuts& cuts,
    const labelList& cellLabels
) const
{
    forAll(cellLabels, labelI)
    {
        label cellI = cellLabels[labelI];

        if (cuts.cellLoops()[cellI].size())
        {
            return cellI;
        }
    }
    return -1;
}


//- Returns first pointI in pointLabels that uses an internal
//  face. Used to find point to inflate cell/face from (has to be
//  connected to internal face). Returns -1 (so inflate from nothing) if
//  none found.
Foam::label Foam::meshCutter::findInternalFacePoint
(
    const labelList& pointLabels
) const
{
    forAll(pointLabels, labelI)
    {
        label pointI = pointLabels[labelI];

        const labelList& pFaces = mesh().pointFaces()[pointI];

        forAll(pFaces, pFaceI)
        {
            label faceI = pFaces[pFaceI];

            if (mesh().isInternalFace(faceI))
            {
                return pointI;
            }
        }
    }

    if (pointLabels.empty())
    {
        FatalErrorIn("meshCutter::findInternalFacePoint(const labelList&)")
            << "Empty pointLabels" << abort(FatalError);
    }

    return -1;
}


// Get new owner and neighbour of face. Checks anchor points to see if
// need to get original or added cell.
void Foam::meshCutter::faceCells
(
    const cellCuts& cuts,
    const label faceI,
    label& own,
    label& nei
) const
{
    const labelListList& anchorPts = cuts.cellAnchorPoints();
    const labelListList& cellLoops = cuts.cellLoops();

    const face& f = mesh().faces()[faceI];

    own = mesh().faceOwner()[faceI];

    if (cellLoops[own].size() && uses(f, anchorPts[own]))
    {
        own = addedCells_[own];
    }

    nei = -1;

    if (mesh().isInternalFace(faceI))
    {
        nei = mesh().faceNeighbour()[faceI];

        if (cellLoops[nei].size() && uses(f, anchorPts[nei]))
        {
            nei = addedCells_[nei];
        }
    }
}


void Foam::meshCutter::getFaceInfo
(
    const label faceI,
    label& patchID,
    label& zoneID,
    label& zoneFlip
) const
{
    patchID = -1;

    if (!mesh().isInternalFace(faceI))
    {
        patchID = mesh().boundaryMesh().whichPatch(faceI);
    }

    zoneID = mesh().faceZones().whichZone(faceI);

    zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh().faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
    }
}


// Adds a face on top of existing faceI.
void Foam::meshCutter::addFace
(
    polyTopoChange& meshMod,
    const label faceI,
    const face& newFace,
    const label own,
    const label nei
)
{
    label patchID, zoneID, zoneFlip;

    getFaceInfo(faceI, patchID, zoneID, zoneFlip);

    if ((nei == -1) || (own < nei))
    {
        // Ordering ok.
        if (debug & 2)
        {
            Pout<< "Adding face " << newFace
                << " with new owner:" << own
                << " with new neighbour:" << nei
                << " patchID:" << patchID
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
                -1,                         // master point
                -1,                         // master edge
                faceI,                      // master face for addition
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
                -1,                         // master point
                -1,                         // master edge
                faceI,                      // master face for addition
                false,                      // flux flip
                patchID,                    // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face zone flip
            )
        );
    }
}


// Modifies existing faceI for either new owner/neighbour or new face points.
void Foam::meshCutter::modFace
(
    polyTopoChange& meshMod,
    const label faceI,
    const face& newFace,
    const label own,
    const label nei
)
{
    label patchID, zoneID, zoneFlip;

    getFaceInfo(faceI, patchID, zoneID, zoneFlip);

    if
    (
        (own != mesh().faceOwner()[faceI])
     || (
            mesh().isInternalFace(faceI)
         && (nei != mesh().faceNeighbour()[faceI])
        )
     || (newFace != mesh().faces()[faceI])
    )
    {
        if (debug & 2)
        {
            Pout<< "Modifying face " << faceI
                << " old vertices:" << mesh().faces()[faceI]
                << " new vertices:" << newFace
                << " new owner:" << own
                << " new neighbour:" << nei
                << " new zoneID:" << zoneID
                << " new zoneFlip:" << zoneFlip
                << endl;
        }

        if ((nei == -1) || (own < nei))
        {
            meshMod.setAction
            (
                polyModifyFace
                (
                    newFace,            // modified face
                    faceI,              // label of face being modified
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
                    faceI,                  // label of face being modified
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


// Copies face starting from startFp up to and including endFp.
void Foam::meshCutter::copyFace
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
void Foam::meshCutter::splitFace
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
        FatalErrorIn
        (
            "meshCutter::splitFace"
            ", const face&, const label, const label, face&, face&)"
        )   << "Cannot find vertex (new numbering) " << v0
            << " on face " << f
            << abort(FatalError);
    }

    label endFp = findIndex(f, v1);

    if (endFp == -1)
    {
        FatalErrorIn
        (
            "meshCutter::splitFace("
            ", const face&, const label, const label, face&, face&)"
        )   << "Cannot find vertex (new numbering) " << v1
            << " on face " << f
            << abort(FatalError);
    }


    f0.setSize((endFp + 1 + f.size() - startFp) % f.size());
    f1.setSize(f.size() - f0.size() + 2);

    copyFace(f, startFp, endFp, f0);
    copyFace(f, endFp, startFp, f1);
}


// Adds additional vertices (from edge cutting) to face. Used for faces which
// are not split but still might use edge that has been cut.
Foam::face Foam::meshCutter::addEdgeCutsToFace(const label faceI) const
{
    const face& f = mesh().faces()[faceI];

    face newFace(2 * f.size());

    label newFp = 0;

    forAll(f, fp)
    {
        // Duplicate face vertex .
        newFace[newFp++] = f[fp];

        // Check if edge has been cut.
        label fp1 = f.fcIndex(fp);

        HashTable<label, edge, Hash<edge> >::const_iterator fnd =
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


// Walk loop (loop of cuts) across circumference of cellI. Returns face in
// new vertices.
// Note: tricky bit is that it can use existing edges which have been split.
Foam::face Foam::meshCutter::loopToFace
(
    const label cellI,
    const labelList& loop
) const
{
    face newFace(2*loop.size());

    label newFaceI = 0;

    forAll(loop, fp)
    {
        label cut = loop[fp];

        if (isEdge(cut))
        {
            label edgeI = getEdge(cut);

            const edge& e = mesh().edges()[edgeI];

            label vertI = addedPoints_[e];

            newFace[newFaceI++] = vertI;
        }
        else
        {
            // cut is vertex.
            label vertI = getVertex(cut);

            newFace[newFaceI++] = vertI;

            label nextCut = loop[loop.fcIndex(fp)];

            if (!isEdge(nextCut))
            {
                // From vertex to vertex -> cross cut only if no existing edge.
                label nextVertI = getVertex(nextCut);

                label edgeI = meshTools::findEdge(mesh(), vertI, nextVertI);

                if (edgeI != -1)
                {
                    // Existing edge. Insert split-edge point if any.
                    HashTable<label, edge, Hash<edge> >::const_iterator fnd =
                        addedPoints_.find(mesh().edges()[edgeI]);

                    if (fnd != addedPoints_.end())
                    {
                        newFace[newFaceI++] = fnd();
                    }
                }
            }
        }
    }
    newFace.setSize(newFaceI);

    return newFace;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::meshCutter::meshCutter(const polyMesh& mesh)
:
    edgeVertex(mesh),
    addedCells_(),
    addedFaces_(),
    addedPoints_()

{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshCutter::~meshCutter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meshCutter::setRefinement
(
    const cellCuts& cuts,
    polyTopoChange& meshMod
)
{
    // Clear and size maps here since mesh size will change.
    addedCells_.clear();
    addedCells_.resize(cuts.nLoops());

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
                FatalErrorIn
                (
                    "meshCutter::setRefinement(const cellCuts&"
                    ", polyTopoChange&)"
                )   << "Problem: cut edge but none of the cells using it is\n"
                    << "edge:" << edgeI << " verts:" << e
                    << abort(FatalError);
            }

            // One of the edge end points should be master point of nbCellI.
            label masterPointI = e.start();

            const point& v0 = mesh().points()[e.start()];
            const point& v1 = mesh().points()[e.end()];

            scalar weight = cuts.edgeWeight()[edgeI];

            point newPt = weight*v1 + (1.0-weight)*v0;

            label addedPointI =
                meshMod.setAction
                (
                    polyAddPoint
                    (
                        newPt,              // point
                        masterPointI,       // master point
                        -1,                 // zone for point
                        true                // supports a cell
                    )
                );

            // Store on (hash of) edge.
            addedPoints_.insert(e, addedPointI);

            if (debug & 2)
            {
                Pout<< "Added point " << addedPointI
                    << " to vertex "
                    << masterPointI << " of edge " << edgeI
                    << " vertices " << e << endl;
            }
        }
    }

    //
    // Add cells (on 'anchor' side of cell)
    //

    forAll(cellLoops, cellI)
    {
        if (cellLoops[cellI].size())
        {
            // Add a cell to the existing cell
            label addedCellI =
                meshMod.setAction
                (
                    polyAddCell
                    (
                        -1,                 // master point
                        -1,                 // master edge
                        -1,                 // master face
                        cellI,              // master cell
                        -1                  // zone for cell
                    )
                );

            addedCells_.insert(cellI, addedCellI);

            if (debug & 2)
            {
                Pout<< "Added cell " << addedCells_[cellI] << " to cell "
                    << cellI << endl;
            }
        }
    }


    //
    // For all cut cells add an internal face
    //

    forAll(cellLoops, cellI)
    {
        const labelList& loop = cellLoops[cellI];

        if (loop.size())
        {
            //
            // Convert loop (=list of cuts) into proper face.
            // Orientation should already be ok. (done by cellCuts)
            //
            face newFace(loopToFace(cellI, loop));

            // Pick any anchor point on cell
            label masterPointI = findInternalFacePoint(anchorPts[cellI]);

            label addedFaceI =
                meshMod.setAction
                (
                    polyAddFace
                    (
                        newFace,                // face
                        cellI,                  // owner
                        addedCells_[cellI],     // neighbour
                        masterPointI,           // master point
                        -1,                     // master edge
                        -1,                     // master face for addition
                        false,                  // flux flip
                        -1,                     // patch for face
                        -1,                     // zone for face
                        false                   // face zone flip
                    )
                );

            addedFaces_.insert(cellI, addedFaceI);

            if (debug & 2)
            {
                // Gets edgeweights of loop
                scalarField weights(loop.size());
                forAll(loop, i)
                {
                    label cut = loop[i];

                    weights[i] =
                    (
                        isEdge(cut)
                      ? cuts.edgeWeight()[getEdge(cut)]
                      : -GREAT
                    );
                }

                Pout<< "Added splitting face " << newFace << " index:"
                    << addedFaceI
                    << " to owner " << cellI
                    << " neighbour " << addedCells_[cellI]
                    << " from Loop:";
                writeCuts(Pout, loop, weights);
                Pout<< endl;
            }
        }
    }


    //
    // Modify faces on the outside and create new ones
    // (in effect split old faces into two)
    //

    // Maintain whether face has been updated (for -split edges
    // -new owner/neighbour)
    boolList faceUptodate(mesh().nFaces(), false);

    const Map<edge>& faceSplitCuts = cuts.faceSplitCut();

    forAllConstIter(Map<edge>, faceSplitCuts, iter)
    {
        label faceI = iter.key();

        // Renumber face to include split edges.
        face newFace(addEdgeCutsToFace(faceI));

        // Edge splitting the face. Convert cuts to new vertex numbering.
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

        label own = mesh().faceOwner()[faceI];

        label nei = -1;

        if (mesh().isInternalFace(faceI))
        {
            nei = mesh().faceNeighbour()[faceI];
        }

        if (debug & 2)
        {
            Pout<< "Split face " << mesh().faces()[faceI]
                << " own:" << own << " nei:" << nei
                << " into f0:" << f0
                << " and f1:" << f1 << endl;
        }

        // Check which face uses anchorPoints (connects to addedCell)
        // and which one doesn't (connects to original cell)

        // Bit tricky. We have to know whether this faceSplit splits owner/
        // neighbour or both. Even if cell is cut we have to make sure this is
        // the one that cuts it (this face cut might not be the one splitting
        // the cell)

        const face& f = mesh().faces()[faceI];

        label f0Owner = -1;
        label f1Owner = -1;

        if (cellLoops[own].empty())
        {
            f0Owner = own;
            f1Owner = own;
        }
        else if (isIn(splitEdge, cellLoops[own]))
        {
            // Owner is cut by this splitCut. See which of f0, f1 gets
            // owner, which gets addedCells_[owner]
            if (uses(f0, anchorPts[own]))
            {
                f0Owner = addedCells_[own];
                f1Owner = own;
            }
            else
            {
                f0Owner = own;
                f1Owner = addedCells_[own];
            }
        }
        else
        {
            // Owner not cut by this splitCut. Check on original face whether
            // use anchorPts.
            if (uses(f, anchorPts[own]))
            {
                label newCellI = addedCells_[own];
                f0Owner = newCellI;
                f1Owner = newCellI;
            }
            else
            {
                f0Owner = own;
                f1Owner = own;
            }
        }


        label f0Neighbour = -1;
        label f1Neighbour = -1;

        if (nei != -1)
        {
            if (cellLoops[nei].empty())
            {
                f0Neighbour = nei;
                f1Neighbour = nei;
            }
            else if (isIn(splitEdge, cellLoops[nei]))
            {
                // Neighbour is cut by this splitCut. See which of f0, f1
                // gets which neighbour/addedCells_[neighbour]
                if (uses(f0, anchorPts[nei]))
                {
                    f0Neighbour = addedCells_[nei];
                    f1Neighbour = nei;
                }
                else
                {
                    f0Neighbour = nei;
                    f1Neighbour = addedCells_[nei];
                }
            }
            else
            {
                // neighbour not cut by this splitCut. Check on original face
                // whether use anchorPts.
                if (uses(f, anchorPts[nei]))
                {
                    label newCellI = addedCells_[nei];
                    f0Neighbour = newCellI;
                    f1Neighbour = newCellI;
                }
                else
                {
                    f0Neighbour = nei;
                    f1Neighbour = nei;
                }
            }
        }

        // f0 is the added face, f1 the modified one
        addFace(meshMod, faceI, f0, f0Owner, f0Neighbour);

        modFace(meshMod, faceI, f1, f1Owner, f1Neighbour);

        faceUptodate[faceI] = true;
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
                label faceI = eFaces[i];

                if (!faceUptodate[faceI])
                {
                    // Renumber face to include split edges.
                    face newFace(addEdgeCutsToFace(faceI));

                    if (debug & 2)
                    {
                        Pout<< "Added edge cuts to face " << faceI
                            << " f:" << mesh().faces()[faceI]
                            << " newFace:" << newFace << endl;
                    }

                    // Get (new or original) owner and neighbour of faceI
                    label own, nei;
                    faceCells(cuts, faceI, own, nei);

                    modFace(meshMod, faceI, newFace, own, nei);

                    faceUptodate[faceI] = true;
                }
            }
        }
    }


    //
    // Correct any original faces on split cell for new neighbour/owner
    //

    forAll(cellLoops, cellI)
    {
        if (cellLoops[cellI].size())
        {
            const labelList& cllFaces = mesh().cells()[cellI];

            forAll(cllFaces, cllFaceI)
            {
                label faceI = cllFaces[cllFaceI];

                if (!faceUptodate[faceI])
                {
                    // Update face with new owner/neighbour (if any)
                    const face& f = mesh().faces()[faceI];

                    if (debug && (f != addEdgeCutsToFace(faceI)))
                    {
                        FatalErrorIn
                        (
                            "meshCutter::setRefinement(const cellCuts&"
                            ", polyTopoChange&)"
                        )   << "Problem: edges added to face which does "
                            << " not use a marked cut" << endl
                            << "faceI:" << faceI << endl
                            << "face:" << f << endl
                            << "newFace:" << addEdgeCutsToFace(faceI)
                            << abort(FatalError);
                    }

                    // Get (new or original) owner and neighbour of faceI
                    label own, nei;
                    faceCells(cuts, faceI, own, nei);

                    modFace
                    (
                        meshMod,
                        faceI,
                        f,
                        own,
                        nei
                    );

                    faceUptodate[faceI] = true;
                }
            }
        }
    }

    if (debug)
    {
        Pout<< "meshCutter:" << nl
            << "    cells split:" << addedCells_.size() << nl
            << "    faces added:" << addedFaces_.size() << nl
            << "    points added on edges:" << addedPoints_.size() << nl
            << endl;
    }
}


void Foam::meshCutter::updateMesh(const mapPolyMesh& morphMap)
{
    // Update stored labels for mesh change.

    {
        // Create copy since new label might (temporarily) clash with existing
        // key.
        Map<label> newAddedCells(addedCells_.size());

        forAllConstIter(Map<label>, addedCells_, iter)
        {
            label cellI = iter.key();
            label newCellI = morphMap.reverseCellMap()[cellI];

            label addedCellI = iter();

            label newAddedCellI = morphMap.reverseCellMap()[addedCellI];

            if (newCellI >= 0 && newAddedCellI >= 0)
            {
                if
                (
                    (debug & 2)
                 && (newCellI != cellI || newAddedCellI != addedCellI)
                )
                {
                    Pout<< "meshCutter::updateMesh :"
                        << " updating addedCell for cell " << cellI
                        << " from " << addedCellI
                        << " to " << newAddedCellI << endl;
                }
                newAddedCells.insert(newCellI, newAddedCellI);
            }
        }

        // Copy
        addedCells_.transfer(newAddedCells);
    }

    {
        Map<label> newAddedFaces(addedFaces_.size());

        forAllConstIter(Map<label>, addedFaces_, iter)
        {
            label cellI = iter.key();
            label newCellI = morphMap.reverseCellMap()[cellI];

            label addedFaceI = iter();

            label newAddedFaceI = morphMap.reverseFaceMap()[addedFaceI];

            if ((newCellI >= 0) && (newAddedFaceI >= 0))
            {
                if
                (
                    (debug & 2)
                 && (newCellI != cellI || newAddedFaceI != addedFaceI)
                )
                {
                    Pout<< "meshCutter::updateMesh :"
                        << " updating addedFace for cell " << cellI
                        << " from " << addedFaceI
                        << " to " << newAddedFaceI
                        << endl;
                }
                newAddedFaces.insert(newCellI, newAddedFaceI);
            }
        }

        // Copy
        addedFaces_.transfer(newAddedFaces);
    }

    {
        HashTable<label, edge, Hash<edge> > newAddedPoints(addedPoints_.size());

        for
        (
            HashTable<label, edge, Hash<edge> >::const_iterator iter =
                addedPoints_.begin();
            iter != addedPoints_.end();
            ++iter
        )
        {
            const edge& e = iter.key();

            label newStart = morphMap.reversePointMap()[e.start()];

            label newEnd = morphMap.reversePointMap()[e.end()];

            label addedPointI = iter();

            label newAddedPointI = morphMap.reversePointMap()[addedPointI];

            if ((newStart >= 0) && (newEnd >= 0) && (newAddedPointI >= 0))
            {
                edge newE = edge(newStart, newEnd);

                if
                (
                    (debug & 2)
                 && (e != newE || newAddedPointI != addedPointI)
                )
                {
                    Pout<< "meshCutter::updateMesh :"
                        << " updating addedPoints for edge " << e
                        << " from " << addedPointI
                        << " to " << newAddedPointI
                        << endl;
                }

                newAddedPoints.insert(newE, newAddedPointI);
            }
        }

        // Copy
        addedPoints_.transfer(newAddedPoints);
    }
}


// ************************************************************************* //
