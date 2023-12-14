/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
#include "polyTopoChangeMap.H"
#include "meshTools.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(meshCutter, 0);
}


// * * * * * * * * * * * * * Private Static Functions  * * * * * * * * * * * //

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

Foam::label Foam::meshCutter::findCutCell
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


Foam::label Foam::meshCutter::findInternalFacePoint
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


void Foam::meshCutter::faceCells
(
    const cellCuts& cuts,
    const label facei,
    label& own,
    label& nei
) const
{
    const labelListList& anchorPts = cuts.cellAnchorPoints();
    const labelListList& cellLoops = cuts.cellLoops();

    const face& f = mesh().faces()[facei];

    own = mesh().faceOwner()[facei];

    if (cellLoops[own].size() && uses(f, anchorPts[own]))
    {
        own = addedCells_[own];
    }

    nei = -1;

    if (mesh().isInternalFace(facei))
    {
        nei = mesh().faceNeighbour()[facei];

        if (cellLoops[nei].size() && uses(f, anchorPts[nei]))
        {
            nei = addedCells_[nei];
        }
    }
}


void Foam::meshCutter::getFaceInfo
(
    const label facei,
    label& patchID,
    label& zoneID,
    label& zoneFlip
) const
{
    patchID = -1;

    if (!mesh().isInternalFace(facei))
    {
        patchID = mesh().boundaryMesh().whichPatch(facei);
    }

    zoneID = mesh().faceZones().whichZone(facei);

    zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh().faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
    }
}


void Foam::meshCutter::addFace
(
    polyTopoChange& meshMod,
    const label facei,
    const face& newFace,
    const label own,
    const label nei
)
{
    label patchID, zoneID, zoneFlip;

    getFaceInfo(facei, patchID, zoneID, zoneFlip);

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

        meshMod.addFace
        (
            newFace,                    // face
            own,                        // owner
            nei,                        // neighbour
            -1,                         // master point
            -1,                         // master edge
            facei,                      // master face for addition
            false,                      // flux flip
            patchID,                    // patch for face
            zoneID,                     // zone for face
            zoneFlip                    // face zone flip
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

        meshMod.addFace
        (
            newFace.reverseFace(),      // face
            nei,                        // owner
            own,                        // neighbour
            -1,                         // master point
            -1,                         // master edge
            facei,                      // master face for addition
            false,                      // flux flip
            patchID,                    // patch for face
            zoneID,                     // zone for face
            zoneFlip                    // face zone flip
        );
    }
}


void Foam::meshCutter::modFace
(
    polyTopoChange& meshMod,
    const label facei,
    const face& newFace,
    const label own,
    const label nei
)
{
    label patchID, zoneID, zoneFlip;

    getFaceInfo(facei, patchID, zoneID, zoneFlip);

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
                << " new zoneID:" << zoneID
                << " new zoneFlip:" << zoneFlip
                << endl;
        }

        if ((nei == -1) || (own < nei))
        {
            meshMod.modifyFace
            (
                newFace,            // modified face
                facei,              // label of face being modified
                own,                // owner
                nei,                // neighbour
                false,              // face flip
                patchID,            // patch for face
                zoneID,             // zone for face
                zoneFlip            // face flip in zone
            );
        }
        else
        {
            meshMod.modifyFace
            (
                newFace.reverseFace(),  // modified face
                facei,                  // label of face being modified
                nei,                    // owner
                own,                    // neighbour
                false,                  // face flip
                patchID,                // patch for face
                zoneID,                 // zone for face
                zoneFlip                // face flip in zone
            );
        }
    }
}


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


Foam::face Foam::meshCutter::addEdgeCutsToFace(const label facei) const
{
    const face& f = mesh().faces()[facei];

    face newFace(2 * f.size());

    label newFp = 0;

    forAll(f, fp)
    {
        // Duplicate face vertex .
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


Foam::face Foam::meshCutter::loopToFace
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

    if (returnReduce(cuts.nLoops(), sumOp<label>()) == 0)
    {
        return;
    }

    const labelListList& anchorPts = cuts.cellAnchorPoints();
    const labelListList& cellLoops = cuts.cellLoops();

    if (debug)
    {
        // Check that any edge is cut only if any cell using it is cut
        boolList edgeOnCutCell(mesh().nEdges(), false);
        forAll(cuts.cellLoops(), celli)
        {
            if (cuts.cellLoops()[celli].size())
            {
                const labelList& cEdges = mesh().cellEdges(celli);
                forAll(cEdges, i)
                {
                    edgeOnCutCell[cEdges[i]] = true;
                }
            }
        }
        syncTools::syncEdgeList(mesh(), edgeOnCutCell, orEqOp<bool>(), false);

        forAll(cuts.edgeIsCut(), edgeI)
        {
            if (cuts.edgeIsCut()[edgeI] && !edgeOnCutCell[edgeI])
            {
                const edge& e = mesh().edges()[edgeI];

                WarningInFunction
                    << "Problem: cut edge but none of the cells using"
                    << " it is cut\n"
                    << "edge:" << edgeI << " verts:" << e
                    << " at:" << e.line(mesh().points())
                    << endl;    // abort(FatalError);
            }
        }
    }


    //
    // Add new points along cut edges.
    //

    forAll(cuts.edgeIsCut(), edgeI)
    {
        if (cuts.edgeIsCut()[edgeI])
        {
            const edge& e = mesh().edges()[edgeI];

            // One of the edge end points should be master point of nbCelli.
            label masterPointi = e.start();

            const point& v0 = mesh().points()[e.start()];
            const point& v1 = mesh().points()[e.end()];

            scalar weight = cuts.edgeWeight()[edgeI];

            point newPt = weight*v1 + (1.0-weight)*v0;

            label addedPointi = meshMod.addPoint
            (
                newPt,              // point
                masterPointi,       // master point
                -1,                 // zone for point
                true                // supports a cell
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
    // Add cells (on 'anchor' side of cell)
    //

    forAll(cellLoops, celli)
    {
        if (cellLoops[celli].size())
        {
            // Add a cell to the existing cell
            label addedCelli = meshMod.addCell
            (
                -1,                 // master point
                -1,                 // master edge
                -1,                 // master face
                celli,              // master cell
                mesh().cellZones().whichZone(celli) // zone for cell
            );

            addedCells_.insert(celli, addedCelli);

            if (debug & 2)
            {
                Pout<< "Added cell " << addedCells_[celli] << " to cell "
                    << celli << endl;
            }
        }
    }


    //
    // For all cut cells add an internal face
    //

    forAll(cellLoops, celli)
    {
        const labelList& loop = cellLoops[celli];

        if (loop.size())
        {
            // Convert loop (=list of cuts) into proper face.
            // Orientation should already be ok. (done by cellCuts)
            //
            face newFace(loopToFace(celli, loop));

            // Pick any anchor point on cell
            label masterPointi = findInternalFacePoint(anchorPts[celli]);

            label addedFacei =
                meshMod.addFace
            (
                newFace,                // face
                celli,                  // owner
                addedCells_[celli],     // neighbour
                masterPointi,           // master point
                -1,                     // master edge
                -1,                     // master face for addition
                false,                  // flux flip
                -1,                     // patch for face
                -1,                     // zone for face
                false                   // face zone flip
            );

            addedFaces_.insert(celli, addedFacei);

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
                      : -great
                    );
                }

                Pout<< "Added splitting face " << newFace << " index:"
                    << addedFacei
                    << " to owner " << celli
                    << " neighbour " << addedCells_[celli]
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
        label facei = iter.key();

        // Renumber face to include split edges.
        face newFace(addEdgeCutsToFace(facei));

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

        // Check which face uses anchorPoints (connects to addedCell)
        // and which one doesn't (connects to original cell)

        // Bit tricky. We have to know whether this faceSplit splits owner/
        // neighbour or both. Even if cell is cut we have to make sure this is
        // the one that cuts it (this face cut might not be the one splitting
        // the cell)

        const face& f = mesh().faces()[facei];

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
                label newCelli = addedCells_[own];
                f0Owner = newCelli;
                f1Owner = newCelli;
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
                    label newCelli = addedCells_[nei];
                    f0Neighbour = newCelli;
                    f1Neighbour = newCelli;
                }
                else
                {
                    f0Neighbour = nei;
                    f1Neighbour = nei;
                }
            }
        }

        // f0 is the added face, f1 the modified one
        addFace(meshMod, facei, f0, f0Owner, f0Neighbour);

        modFace(meshMod, facei, f1, f1Owner, f1Neighbour);

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
                    // Renumber face to include split edges.
                    face newFace(addEdgeCutsToFace(facei));

                    if (debug & 2)
                    {
                        Pout<< "Added edge cuts to face " << facei
                            << " f:" << mesh().faces()[facei]
                            << " newFace:" << newFace << endl;
                    }

                    // Get (new or original) owner and neighbour of facei
                    label own, nei;
                    faceCells(cuts, facei, own, nei);

                    modFace(meshMod, facei, newFace, own, nei);

                    faceUptodate[facei] = true;
                }
            }
        }
    }


    //
    // Correct any original faces on split cell for new neighbour/owner
    //

    forAll(cellLoops, celli)
    {
        if (cellLoops[celli].size())
        {
            const labelList& cllFaces = mesh().cells()[celli];

            forAll(cllFaces, cllFacei)
            {
                label facei = cllFaces[cllFacei];

                if (!faceUptodate[facei])
                {
                    // Update face with new owner/neighbour (if any)
                    const face& f = mesh().faces()[facei];

                    if (debug && (f != addEdgeCutsToFace(facei)))
                    {
                        FatalErrorInFunction
                            << "Problem: edges added to face which does "
                            << " not use a marked cut" << endl
                            << "facei:" << facei << endl
                            << "face:" << f << endl
                            << "newFace:" << addEdgeCutsToFace(facei)
                            << abort(FatalError);
                    }

                    // Get (new or original) owner and neighbour of facei
                    label own, nei;
                    faceCells(cuts, facei, own, nei);

                    modFace
                    (
                        meshMod,
                        facei,
                        f,
                        own,
                        nei
                    );

                    faceUptodate[facei] = true;
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


void Foam::meshCutter::topoChange(const polyTopoChangeMap& map)
{
    // Update stored labels for mesh change.

    {
        // Create copy since new label might (temporarily) clash with existing
        // key.
        Map<label> newAddedCells(addedCells_.size());

        forAllConstIter(Map<label>, addedCells_, iter)
        {
            label celli = iter.key();
            label newCelli = map.reverseCellMap()[celli];

            label addedCelli = iter();

            label newAddedCelli = map.reverseCellMap()[addedCelli];

            if (newCelli >= 0 && newAddedCelli >= 0)
            {
                if
                (
                    (debug & 2)
                 && (newCelli != celli || newAddedCelli != addedCelli)
                )
                {
                    Pout<< "meshCutter::topoChange :"
                        << " updating addedCell for cell " << celli
                        << " from " << addedCelli
                        << " to " << newAddedCelli << endl;
                }
                newAddedCells.insert(newCelli, newAddedCelli);
            }
        }

        // Copy
        addedCells_.transfer(newAddedCells);
    }

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
                    Pout<< "meshCutter::topoChange :"
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
                    Pout<< "meshCutter::topoChange :"
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
