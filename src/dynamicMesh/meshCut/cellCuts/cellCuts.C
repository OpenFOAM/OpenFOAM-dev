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

#include "cellCuts.H"
#include "polyMesh.H"
#include "Time.H"
#include "ListOps.H"
#include "cellLooper.H"
#include "refineCell.H"
#include "meshTools.H"
#include "geomCellLooper.H"
#include "OFstream.H"
#include "plane.H"
#include "syncTools.H"
#include "dummyTransform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellCuts, 0);

    //- Template specialization for pTraits<edge> so we can use syncTools
    //  functionality
    template<>
    class pTraits<edge>
    {
    public:

        //- Component type
        typedef edge cmptType;
    };
}


// * * * * * * * * * * * * * Private Static Functions  * * * * * * * * * * * //

Foam::label Foam::cellCuts::findPartIndex
(
    const labelList& elems,
    const label nElems,
    const label val
)
{
    for (label i = 0; i < nElems; i++)
    {
        if (elems[i] == val)
        {
            return i;
        }
    }
    return -1;
}


Foam::boolList Foam::cellCuts::expand
(
    const label size,
    const labelList& labels
)
{
    boolList result(size, false);

    forAll(labels, labelI)
    {
        result[labels[labelI]] = true;
    }
    return result;
}


Foam::scalarField Foam::cellCuts::expand
(
    const label size,
    const labelList& labels,
    const scalarField& weights
)
{
    scalarField result(size, -great);

    forAll(labels, labelI)
    {
        result[labels[labelI]] = weights[labelI];
    }
    return result;
}


Foam::label Foam::cellCuts::firstUnique
(
    const labelList& lst,
    const Map<label>& map
)
{
    forAll(lst, i)
    {
        if (!map.found(lst[i]))
        {
            return i;
        }
    }
    return -1;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cellCuts::syncProc()
{
    if (!Pstream::parRun())
    {
        return;
    }

    syncTools::syncPointList(mesh(), pointIsCut_, orEqOp<bool>(), false);
    syncTools::syncEdgeList(mesh(), edgeIsCut_, orEqOp<bool>(), false);
    syncTools::syncEdgeList(mesh(), edgeWeight_, maxEqOp<scalar>(), -great);

    {
        const label nBnd = mesh().nFaces()-mesh().nInternalFaces();

        // Convert faceSplitCut into face-local data: vertex and edge w.r.t.
        // vertex 0: (since this is same on both sides)
        //
        //      Sending-side vertex  Receiving-side vertex
        //      0                   0
        //      1                   3
        //      2                   2
        //      3                   1
        //
        //      Sending-side edge    Receiving side edge
        //      0-1                  3-0
        //      1-2                  2-3
        //      2-3                  1-2
        //      3-0                  0-1
        //
        // Encoding is as index:
        //      0  : not set
        //      >0 : vertex, vertex is index-1
        //      <0 : edge, edge is -index-1

        edgeList relCuts(nBnd, edge(0, 0));

        const polyBoundaryMesh& pbm = mesh().boundaryMesh();

        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];

            if (isA<processorPolyPatch>(pp) || isA<cyclicPolyPatch>(pp))
            {
                forAll(pp, i)
                {
                    label facei = pp.start()+i;
                    label bFacei = facei-mesh().nInternalFaces();

                    const Map<edge>::const_iterator iter =
                        faceSplitCut_.find(facei);
                    if (iter != faceSplitCut_.end())
                    {
                        const face& f = mesh().faces()[facei];
                        const labelList& fEdges = mesh().faceEdges()[facei];
                        const edge& cuts = iter();

                        forAll(cuts, i)
                        {
                            if (isEdge(cuts[i]))
                            {
                                label edgei = getEdge(cuts[i]);
                                label index = findIndex(fEdges, edgei);
                                relCuts[bFacei][i] = -index-1;
                            }
                            else
                            {
                                label index = findIndex(f, getVertex(cuts[i]));
                                relCuts[bFacei][i] = index+1;
                            }
                        }
                    }
                }
            }
        }

        // Exchange
        syncTools::syncBoundaryFaceList
        (
            mesh(),
            relCuts,
            eqOp<edge>(),
            dummyTransform()
        );

        // Convert relCuts back into mesh based data
        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];

            if (isA<processorPolyPatch>(pp) || isA<cyclicPolyPatch>(pp))
            {
                forAll(pp, i)
                {
                    label facei = pp.start()+i;
                    label bFacei = facei-mesh().nInternalFaces();

                    const edge& relCut = relCuts[bFacei];
                    if (relCut != edge(0, 0))
                    {
                        const face& f = mesh().faces()[facei];
                        const labelList& fEdges = mesh().faceEdges()[facei];

                        edge absoluteCut(0, 0);
                        forAll(relCut, i)
                        {
                            if (relCut[i] < 0)
                            {
                                label oppFp = -relCut[i]-1;
                                label fp = f.size()-1-oppFp;
                                absoluteCut[i] = edgeToEVert(fEdges[fp]);
                            }
                            else
                            {
                                label oppFp = relCut[i]-1;
                                label fp = f.size()-1-oppFp;
                                absoluteCut[i] = vertToEVert(f[fp]);
                            }
                        }

                        if
                        (
                           !faceSplitCut_.insert(facei, absoluteCut)
                         && faceSplitCut_[facei] != absoluteCut
                        )
                        {
                            FatalErrorInFunction
                                << "Cut " << faceSplitCut_[facei]
                                << " on face " << mesh().faceCentres()[facei]
                                << " of coupled patch " << pp.name()
                                << " is not consistent with coupled cut "
                                << absoluteCut
                                << exit(FatalError);
                        }
                    }
                }
            }
        }
    }
}


void Foam::cellCuts::writeUncutOBJ
(
    const fileName& dir,
    const label celli
) const
{
    // Cell edges
    OFstream cutsStream(dir / "cell_" + name(celli) + ".obj");

    Pout<< "Writing cell for time " <<  mesh().time().timeName()
        << " to " << cutsStream.name() << nl;

    meshTools::writeOBJ
    (
        cutsStream,
        mesh().cells(),
        mesh().faces(),
        mesh().points(),
        labelList(1, celli)
    );

    // Loop cutting cell in two
    OFstream cutStream(dir / "cellCuts_" + name(celli) + ".obj");

    Pout<< "Writing raw cuts on cell for time " <<  mesh().time().timeName()
        << " to " << cutStream.name() << nl;

    const labelList& cPoints = mesh().cellPoints()[celli];

    forAll(cPoints, i)
    {
        label pointi = cPoints[i];
        if (pointIsCut_[pointi])
        {
            meshTools::writeOBJ(cutStream, mesh().points()[pointi]);
        }
    }

    const pointField& pts = mesh().points();

    const labelList& cEdges = mesh().cellEdges()[celli];

    forAll(cEdges, i)
    {
        label edgeI = cEdges[i];

        if (edgeIsCut_[edgeI])
        {
            const edge& e = mesh().edges()[edgeI];

            const scalar w = edgeWeight_[edgeI];

            meshTools::writeOBJ(cutStream, w*pts[e[1]] + (1-w)*pts[e[0]]);
        }
    }
}


void Foam::cellCuts::writeOBJ
(
    const fileName& dir,
    const label celli,
    const pointField& loopPoints,
    const labelList& anchors
) const
{
    // Cell edges
    OFstream cutsStream(dir / "cell_" + name(celli) + ".obj");

    Pout<< "Writing cell for time " <<  mesh().time().timeName()
        << " to " << cutsStream.name() << nl;

    meshTools::writeOBJ
    (
        cutsStream,
        mesh().cells(),
        mesh().faces(),
        mesh().points(),
        labelList(1, celli)
    );


    // Loop cutting cell in two
    OFstream loopStream(dir / "cellLoop_" + name(celli) + ".obj");

    Pout<< "Writing loop for time " <<  mesh().time().timeName()
        << " to " << loopStream.name() << nl;

    label vertI = 0;

    writeOBJ(loopStream, loopPoints, vertI);


    // Anchors for cell
    OFstream anchorStream(dir / "anchors_" + name(celli) + ".obj");

    Pout<< "Writing anchors for time " <<  mesh().time().timeName()
        << " to " << anchorStream.name() << endl;

    forAll(anchors, i)
    {
        meshTools::writeOBJ(anchorStream, mesh().points()[anchors[i]]);
    }
}


Foam::label Foam::cellCuts::edgeEdgeToFace
(
    const label celli,
    const label edgeA,
    const label edgeB
) const
{
    const labelList& cFaces = mesh().cells()[celli];

    forAll(cFaces, cFacei)
    {
        label facei = cFaces[cFacei];

        const labelList& fEdges = mesh().faceEdges()[facei];

        if
        (
            findIndex(fEdges, edgeA) != -1
         && findIndex(fEdges, edgeB) != -1
        )
        {
           return facei;
        }
    }

    // Coming here means the loop is illegal since the two edges
    // are not shared by a face. We just mark loop as invalid and continue.

    WarningInFunction
        << "cellCuts : Cannot find face on cell "
        << celli << " that has both edges " << edgeA << ' ' << edgeB << endl
        << "faces : " << cFaces << endl
        << "edgeA : " << mesh().edges()[edgeA] << endl
        << "edgeB : " << mesh().edges()[edgeB] << endl
        << "Marking the loop across this cell as invalid" << endl;

    return -1;
}


Foam::label Foam::cellCuts::edgeVertexToFace
(
    const label celli,
    const label edgeI,
    const label vertI
) const
{
    const labelList& cFaces = mesh().cells()[celli];

    forAll(cFaces, cFacei)
    {
        label facei = cFaces[cFacei];

        const face& f = mesh().faces()[facei];

        const labelList& fEdges = mesh().faceEdges()[facei];

        if
        (
            findIndex(fEdges, edgeI) != -1
         && findIndex(f, vertI) != -1
        )
        {
           return facei;
        }
    }

    WarningInFunction
        << "cellCuts : Cannot find face on cell "
        << celli << " that has both edge " << edgeI << " and vertex "
        << vertI << endl
        << "faces : " << cFaces << endl
        << "edge : " << mesh().edges()[edgeI] << endl
        << "Marking the loop across this cell as invalid" << endl;

    return -1;
}


Foam::label Foam::cellCuts::vertexVertexToFace
(
    const label celli,
    const label vertA,
    const label vertB
) const
{
    const labelList& cFaces = mesh().cells()[celli];

    forAll(cFaces, cFacei)
    {
        label facei = cFaces[cFacei];

        const face& f = mesh().faces()[facei];

        if
        (
            findIndex(f, vertA) != -1
         && findIndex(f, vertB) != -1
        )
        {
           return facei;
        }
    }

    WarningInFunction
        << "cellCuts : Cannot find face on cell "
        << celli << " that has vertex " << vertA << " and vertex "
        << vertB << endl
        << "faces : " << cFaces << endl
        << "Marking the loop across this cell as invalid" << endl;

    return -1;
}


void Foam::cellCuts::calcFaceCuts() const
{
    if (faceCutsPtr_.valid())
    {
        FatalErrorInFunction
            << "faceCuts already calculated" << abort(FatalError);
    }

    const faceList& faces = mesh().faces();

    faceCutsPtr_.reset(new labelListList(mesh().nFaces()));
    labelListList& faceCuts = faceCutsPtr_();

    for (label facei = 0; facei < mesh().nFaces(); facei++)
    {
        const face& f = faces[facei];

        // Big enough storage (possibly all points and all edges cut). Shrink
        // later on.
        labelList& cuts = faceCuts[facei];

        cuts.setSize(2*f.size());

        label cutI = 0;

        // Do point-edge-point walk over face and collect all cuts.
        // Problem is that we want to start from one of the endpoints of a
        // string of connected cuts; we don't want to start somewhere in the
        // middle.

        // Pass1: find first point cut not preceded by a cut.
        label startFp = -1;

        forAll(f, fp)
        {
            label v0 = f[fp];

            if (pointIsCut_[v0])
            {
                label vMin1 = f[f.rcIndex(fp)];

                if
                (
                    !pointIsCut_[vMin1]
                 && !edgeIsCut_[findEdge(facei, v0, vMin1)]
                )
                {
                    cuts[cutI++] = vertToEVert(v0);
                    startFp = f.fcIndex(fp);
                    break;
                }
            }
        }

        // Pass2: first edge cut not preceded by point cut
        if (startFp == -1)
        {
            forAll(f, fp)
            {
                label fp1 = f.fcIndex(fp);

                label v0 = f[fp];
                label v1 = f[fp1];

                label edgeI = findEdge(facei, v0, v1);

                if (edgeIsCut_[edgeI] && !pointIsCut_[v0])
                {
                    cuts[cutI++] = edgeToEVert(edgeI);
                    startFp = fp1;
                    break;
                }
            }
        }

        // Pass3: nothing found so far. Either face is not cut at all or all
        // vertices are cut. Start from 0.
        if (startFp == -1)
        {
            startFp = 0;
        }

        // Store all cuts starting from startFp;
        label fp = startFp;

        bool allVerticesCut = true;

        forAll(f, i)
        {
            label fp1 = f.fcIndex(fp);

            // Get the three items: current vertex, next vertex and edge
            // in between
            label v0 = f[fp];
            label v1 = f[fp1];
            label edgeI = findEdge(facei, v0, v1);

            if (pointIsCut_[v0])
            {
                cuts[cutI++] = vertToEVert(v0);
            }
            else
            {
                // For check if all vertices have been cut (= illegal)
                allVerticesCut = false;
            }

            if (edgeIsCut_[edgeI])
            {
                cuts[cutI++] = edgeToEVert(edgeI);
            }

            fp = fp1;
        }


        if (allVerticesCut)
        {
            WarningInFunction
                << "Face " << facei << " vertices " << f
                << " has all its vertices cut. Not cutting face." << endl;

            cutI = 0;
        }

        // Remove duplicate starting point
        if (cutI > 1 && cuts[cutI-1] == cuts[0])
        {
            cutI--;
        }
        cuts.setSize(cutI);
    }
}


Foam::label Foam::cellCuts::findEdge
(
    const label facei,
    const label v0,
    const label v1
) const
{
    const edgeList& edges = mesh().edges();

    const labelList& fEdges = mesh().faceEdges()[facei];

    forAll(fEdges, i)
    {
        const edge& e = edges[fEdges[i]];

        if
        (
            (e[0] == v0 && e[1] == v1)
         || (e[0] == v1 && e[1] == v0)
        )
        {
            return fEdges[i];
        }
    }

    return -1;
}


Foam::label Foam::cellCuts::loopFace
(
    const label celli,
    const labelList& loop
) const
{
    const cell& cFaces = mesh().cells()[celli];

    forAll(cFaces, cFacei)
    {
        label facei = cFaces[cFacei];

        const labelList& fEdges = mesh().faceEdges()[facei];
        const face& f = mesh().faces()[facei];

        bool allOnFace = true;

        forAll(loop, i)
        {
            label cut = loop[i];

            if (isEdge(cut))
            {
                if (findIndex(fEdges, getEdge(cut)) == -1)
                {
                    // Edge not on face. Skip face.
                    allOnFace = false;
                    break;
                }
            }
            else
            {
                if (findIndex(f, getVertex(cut)) == -1)
                {
                    // Vertex not on face. Skip face.
                    allOnFace = false;
                    break;
                }
            }
        }

        if (allOnFace)
        {
            // Found face where all elements of loop are on the face.
            return facei;
        }
    }
    return -1;
}


bool Foam::cellCuts::walkPoint
(
    const label celli,
    const label startCut,

    const label exclude0,
    const label exclude1,

    const label otherCut,

    label& nVisited,
    labelList& visited
) const
{
    label vertI = getVertex(otherCut);

    const labelList& pFaces = mesh().pointFaces()[vertI];

    forAll(pFaces, pFacei)
    {
        label otherFacei = pFaces[pFacei];

        if
        (
            otherFacei != exclude0
         && otherFacei != exclude1
         && meshTools::faceOnCell(mesh(), celli, otherFacei)
        )
        {
            label oldNVisited = nVisited;

            bool foundLoop =
                walkCell
                (
                    celli,
                    startCut,
                    otherFacei,
                    otherCut,
                    nVisited,
                    visited
                );

            if (foundLoop)
            {
                return true;
            }

            // No success. Restore state and continue
            nVisited = oldNVisited;
        }
    }
    return false;
}


bool Foam::cellCuts::crossEdge
(
    const label celli,
    const label startCut,
    const label facei,
    const label otherCut,

    label& nVisited,
    labelList& visited
) const
{
    // Cross edge to other face
    label edgeI = getEdge(otherCut);

    label otherFacei = meshTools::otherFace(mesh(), celli, facei, edgeI);

    // Store old state
    label oldNVisited = nVisited;

    // Recurse to otherCut
    bool foundLoop =
        walkCell
        (
            celli,
            startCut,
            otherFacei,
            otherCut,
            nVisited,
            visited
        );

    if (foundLoop)
    {
        return true;
    }
    else
    {
        // No success. Restore state (i.e. backtrack)
        nVisited = oldNVisited;

        return false;
    }
}


bool Foam::cellCuts::addCut
(
    const label celli,
    const label cut,
    label& nVisited,
    labelList& visited
) const
{
    // Check for duplicate cuts.
    if (findPartIndex(visited, nVisited, cut) != -1)
    {
        // Truncate (copy of) visited for ease of printing.
        labelList truncVisited(visited);
        truncVisited.setSize(nVisited);

        Pout<< "For cell " << celli << " : trying to add duplicate cut " << cut;
        labelList cuts(1, cut);
        writeCuts(Pout, cuts, loopWeights(cuts));

        Pout<< " to path:";
        writeCuts(Pout, truncVisited, loopWeights(truncVisited));
        Pout<< endl;

        return false;
    }
    else
    {
        visited[nVisited++] = cut;

        return true;
    }
}


bool Foam::cellCuts::walkFace
(
    const label celli,
    const label startCut,
    const label facei,
    const label cut,

    label& lastCut,
    label& beforeLastCut,
    label& nVisited,
    labelList& visited
) const
{
    const labelList& fCuts = faceCuts()[facei];

    if (fCuts.size() < 2)
    {
        return false;
    }

    // Easy case : two cuts.
    if (fCuts.size() == 2)
    {
        if (fCuts[0] == cut)
        {
            if (!addCut(celli, cut, nVisited, visited))
            {
                return false;
            }

            beforeLastCut = cut;
            lastCut = fCuts[1];

            return true;
        }
        else
        {
            if (!addCut(celli, cut, nVisited, visited))
            {
                return false;
            }

            beforeLastCut = cut;
            lastCut = fCuts[0];

            return true;
        }
    }

    // Harder case: more than 2 cuts on face.
    // Should be start or end of string of cuts. Store all but last cut.
    if (fCuts[0] == cut)
    {
        // Walk forward
        for (label i = 0; i < fCuts.size()-1; i++)
        {
            if (!addCut(celli, fCuts[i], nVisited, visited))
            {
                return false;
            }
        }
        beforeLastCut = fCuts[fCuts.size()-2];
        lastCut = fCuts[fCuts.size()-1];
    }
    else if (fCuts[fCuts.size()-1] == cut)
    {
        for (label i = fCuts.size()-1; i >= 1; --i)
        {
            if (!addCut(celli, fCuts[i], nVisited, visited))
            {
                return false;
            }
        }
        beforeLastCut = fCuts[1];
        lastCut = fCuts[0];
    }
    else
    {
        WarningInFunction
            << "In middle of cut. cell:" << celli << " face:" << facei
            << " cuts:" << fCuts << " current cut:" << cut << endl;

        return false;
    }

    return true;
}



bool Foam::cellCuts::walkCell
(
    const label celli,
    const label startCut,
    const label facei,
    const label cut,

    label& nVisited,
    labelList& visited
) const
{
    // Walk across face, storing cuts. Return last two cuts visited.
    label lastCut = -1;
    label beforeLastCut = -1;


    if (debug & 2)
    {
        Pout<< "For cell:" << celli << " walked across face " << facei
            << " from cut ";
        labelList cuts(1, cut);
        writeCuts(Pout, cuts, loopWeights(cuts));
        Pout<< endl;
    }

    bool validWalk = walkFace
    (
        celli,
        startCut,
        facei,
        cut,

        lastCut,
        beforeLastCut,
        nVisited,
        visited
    );

    if (!validWalk)
    {
        return false;
    }

    if (debug & 2)
    {
        Pout<< "    to last cut ";
        labelList cuts(1, lastCut);
        writeCuts(Pout, cuts, loopWeights(cuts));
        Pout<< endl;
    }

    // Check if starting point reached.
    if (lastCut == startCut)
    {
        if (nVisited >= 3)
        {
            if (debug & 2)
            {
                // Truncate (copy of) visited for ease of printing.
                labelList truncVisited(visited);
                truncVisited.setSize(nVisited);

                Pout<< "For cell " << celli << " : found closed path:";
                writeCuts(Pout, truncVisited, loopWeights(truncVisited));
                Pout<< " closed by " << lastCut << endl;
            }

            return true;
        }
        else
        {
            // Loop but too short. Probably folded back on itself.
            return false;
        }
    }


    // Check last cut and one before that to determine type.

    // From beforeLastCut to lastCut:
    //  - from edge to edge
    //      (- always not along existing edge)
    //  - from edge to vertex
    //      - not along existing edge
    //      - along existing edge: not allowed?
    //  - from vertex to edge
    //      - not along existing edge
    //      - along existing edge. Not allowed. See above.
    //  - from vertex to vertex
    //      - not along existing edge
    //      - along existing edge

    if (isEdge(beforeLastCut))
    {
        if (isEdge(lastCut))
        {
            // beforeLastCut=edge, lastCut=edge.

            // Cross lastCut (=edge) into face not facei
            return crossEdge
            (
                celli,
                startCut,
                facei,
                lastCut,
                nVisited,
                visited
            );
        }
        else
        {
            // beforeLastCut=edge, lastCut=vertex.

            // Go from lastCut to all connected faces (except facei)
            return walkPoint
            (
                celli,
                startCut,
                facei,          // exclude0: don't cross back on current face
                -1,             // exclude1
                lastCut,
                nVisited,
                visited
            );
        }
    }
    else
    {
        if (isEdge(lastCut))
        {
            // beforeLastCut=vertex, lastCut=edge.
            return crossEdge
            (
                celli,
                startCut,
                facei,
                lastCut,
                nVisited,
                visited
            );
        }
        else
        {
            // beforeLastCut=vertex, lastCut=vertex. Check if along existing
            // edge.
            label edgeI =
                findEdge
                (
                    facei,
                    getVertex(beforeLastCut),
                    getVertex(lastCut)
                );

            if (edgeI != -1)
            {
                // Cut along existing edge. So is in fact on two faces.
                // Get faces on both sides of the edge to make
                // sure we don't fold back on to those.

                label f0, f1;
                meshTools::getEdgeFaces(mesh(), celli, edgeI, f0, f1);

                return walkPoint
                (
                    celli,
                    startCut,
                    f0,
                    f1,
                    lastCut,
                    nVisited,
                    visited
                );

            }
            else
            {
                // Cross cut across face.
                return walkPoint
                (
                    celli,
                    startCut,
                    facei,      // face to exclude
                    -1,         // face to exclude
                    lastCut,
                    nVisited,
                    visited
                );
            }
        }
    }
}


void Foam::cellCuts::calcCellLoops(const labelList& cutCells)
{
    // Determine for every cut cell the loop (= face) it is cut by. Done by
    // starting from a cut edge or cut vertex and walking across faces, from
    // cut to cut, until starting cut hit.
    // If multiple loops are possible across a cell circumference takes the
    // first one found.

    // Calculate cuts per face.
    const labelListList& allFaceCuts = faceCuts();

    // Per cell the number of faces with valid cuts. Is used as quick
    // rejection to see if cell can be cut.
    labelList nCutFaces(mesh().nCells(), 0);

    forAll(allFaceCuts, facei)
    {
        const labelList& fCuts = allFaceCuts[facei];

        if (fCuts.size() == mesh().faces()[facei].size())
        {
            // Too many cuts on face. WalkCell would get very upset so disable.
            nCutFaces[mesh().faceOwner()[facei]] = labelMin;

            if (mesh().isInternalFace(facei))
            {
                nCutFaces[mesh().faceNeighbour()[facei]] = labelMin;
            }
        }
        else if (fCuts.size() >= 2)
        {
            // Could be valid cut. Update count for owner and neighbour.
            nCutFaces[mesh().faceOwner()[facei]]++;

            if (mesh().isInternalFace(facei))
            {
                nCutFaces[mesh().faceNeighbour()[facei]]++;
            }
        }
    }


    // Stack of visited cuts (nVisited used as stack pointer)
    // Size big enough.
    labelList visited(mesh().nPoints());

    forAll(cutCells, i)
    {
        label celli = cutCells[i];

        bool validLoop = false;

        // Quick rejection: has enough faces that are cut?
        if (nCutFaces[celli] >= 3)
        {
            const labelList& cFaces = mesh().cells()[celli];

            if (debug & 2)
            {
                Pout<< "cell:" << celli << " cut faces:" << endl;
                forAll(cFaces, i)
                {
                    label facei = cFaces[i];
                    const labelList& fCuts = allFaceCuts[facei];

                    Pout<< "    face:" << facei << " cuts:";
                    writeCuts(Pout, fCuts, loopWeights(fCuts));
                    Pout<< endl;
                }
            }

            label nVisited = 0;

            // Determine the first cut face to start walking from.
            forAll(cFaces, cFacei)
            {
                label facei = cFaces[cFacei];

                const labelList& fCuts = allFaceCuts[facei];

                // Take first or last cut of multiple on face.
                // Note that in calcFaceCuts
                // we have already made sure this is the start or end of a
                // string of cuts.
                if (fCuts.size() >= 2)
                {
                    // Try walking from start of fCuts.
                    nVisited = 0;

                    if (debug & 2)
                    {
                        Pout<< "cell:" << celli
                            << " start walk at face:" << facei
                            << " cut:";
                        labelList cuts(1, fCuts[0]);
                        writeCuts(Pout, cuts, loopWeights(cuts));
                        Pout<< endl;
                    }

                    validLoop =
                        walkCell
                        (
                            celli,
                            fCuts[0],
                            facei,
                            fCuts[0],

                            nVisited,
                            visited
                        );

                    if (validLoop)
                    {
                        break;
                    }

                    // No need to try and walk from end of cuts since
                    // all paths already tried by walkCell.
                }
            }

            if (validLoop)
            {
                // Copy nVisited elements out of visited (since visited is
                // never truncated for efficiency reasons)

                labelList& loop = cellLoops_[celli];

                loop.setSize(nVisited);

                for (label i = 0; i < nVisited; i++)
                {
                    loop[i] = visited[i];
                }
            }
            else
            {
                // Invalid loop. Leave cellLoops_[celli] zero size which
                // flags this.
                Pout<< "calcCellLoops(const labelList&) : did not find valid"
                    << " loop for cell " << celli << endl;
                // Dump cell and cuts on cell.
                writeUncutOBJ(".", celli);

                cellLoops_[celli].setSize(0);
            }
        }
        else
        {
            //Pout<< "calcCellLoops(const labelList&) : did not find valid"
            //    << " loop for cell " << celli << " since not enough cut faces"
            //    << endl;
            cellLoops_[celli].setSize(0);
        }
    }
}


void Foam::cellCuts::walkEdges
(
    const label celli,
    const label pointi,
    const label status,

    Map<label>& edgeStatus,
    Map<label>& pointStatus
) const
{
    if (pointStatus.insert(pointi, status))
    {
        // First visit to pointi

        const labelList& pEdges = mesh().pointEdges()[pointi];

        forAll(pEdges, pEdgeI)
        {
            label edgeI = pEdges[pEdgeI];

            if
            (
                meshTools::edgeOnCell(mesh(), celli, edgeI)
             && edgeStatus.insert(edgeI, status)
            )
            {
                // First visit to edgeI so recurse.

                label v2 = mesh().edges()[edgeI].otherVertex(pointi);

                walkEdges(celli, v2, status, edgeStatus, pointStatus);
            }
        }
    }
}


Foam::labelList Foam::cellCuts::nonAnchorPoints
(
    const labelList& cellPoints,
    const labelList& anchorPoints,
    const labelList& loop
) const
{
    labelList newElems(cellPoints.size());
    label newElemI = 0;

    forAll(cellPoints, i)
    {
        label pointi = cellPoints[i];

        if
        (
            findIndex(anchorPoints, pointi) == -1
         && findIndex(loop, vertToEVert(pointi)) == -1
        )
        {
            newElems[newElemI++] = pointi;
        }
    }

    newElems.setSize(newElemI);

    return newElems;
}


bool Foam::cellCuts::loopAnchorConsistent
(
    const label celli,
    const pointField& loopPts,
    const labelList& anchorPoints
) const
{
    // Create identity face for ease of calculation of normal etc.
    face f(identity(loopPts.size()));

    vector normal = f.normal(loopPts);
    point ctr = f.centre(loopPts);


    // Get average position of anchor points.
    vector avg(Zero);

    forAll(anchorPoints, ptI)
    {
        avg += mesh().points()[anchorPoints[ptI]];
    }
    avg /= anchorPoints.size();


    if (((avg - ctr) & normal) > 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::cellCuts::calcAnchors
(
    const label celli,
    const labelList& loop,
    const pointField& loopPts,

    labelList& anchorPoints
) const
{
    const edgeList& edges = mesh().edges();

    const labelList& cPoints = mesh().cellPoints()[celli];
    const labelList& cEdges = mesh().cellEdges()[celli];
    const cell& cFaces = mesh().cells()[celli];

    // Points on loop

    // Status of point:
    // 0 - on loop
    // 1 - point set 1
    // 2 - point set 2
    Map<label> pointStatus(2*cPoints.size());
    Map<label> edgeStatus(2*cEdges.size());

    // Mark loop vertices
    forAll(loop, i)
    {
        label cut = loop[i];

        if (isEdge(cut))
        {
            edgeStatus.insert(getEdge(cut), 0);
        }
        else
        {
            pointStatus.insert(getVertex(cut), 0);
        }
    }
    // Since edges between two cut vertices have not been marked, mark them
    // explicitly
    forAll(cEdges, i)
    {
        label edgeI = cEdges[i];
        const edge& e = edges[edgeI];

        if (pointStatus.found(e[0]) && pointStatus.found(e[1]))
        {
            edgeStatus.insert(edgeI, 0);
        }
    }


    // Find uncut starting vertex
    label uncutIndex = firstUnique(cPoints, pointStatus);

    if (uncutIndex == -1)
    {
        WarningInFunction
            << "Invalid loop " << loop << " for cell " << celli << endl
            << "Can not find point on cell which is not cut by loop."
            << endl;

        writeOBJ(".", celli, loopPts, labelList(0));

        return false;
    }

    // Walk unset vertices and edges and mark with 1 in pointStatus, edgeStatus
    walkEdges(celli, cPoints[uncutIndex], 1, edgeStatus, pointStatus);

    // Find new uncut starting vertex
    uncutIndex = firstUnique(cPoints, pointStatus);

    if (uncutIndex == -1)
    {
        // All vertices either in loop or in anchor. So split is along single
        // face.
        WarningInFunction
            << "Invalid loop " << loop << " for cell " << celli << endl
            << "All vertices of cell are either in loop or in anchor set"
            << endl;

        writeOBJ(".", celli, loopPts, labelList(0));

        return false;
    }

    // Walk unset vertices and edges and mark with 2. These are the
    // pointset 2.
    walkEdges(celli, cPoints[uncutIndex], 2, edgeStatus, pointStatus);

    // Collect both sets in lists.
    DynamicList<label> connectedPoints(cPoints.size());
    DynamicList<label> otherPoints(cPoints.size());

    forAllConstIter(Map<label>, pointStatus, iter)
    {
        if (iter() == 1)
        {
            connectedPoints.append(iter.key());
        }
        else if (iter() == 2)
        {
            otherPoints.append(iter.key());
        }
    }
    connectedPoints.shrink();
    otherPoints.shrink();

    // Check that all points have been used.
    uncutIndex = firstUnique(cPoints, pointStatus);

    if (uncutIndex != -1)
    {
        WarningInFunction
            << "Invalid loop " << loop << " for cell " << celli
            << " since it splits the cell into more than two cells" << endl;

        writeOBJ(".", celli, loopPts, connectedPoints);

        return false;
    }


    // Check that both parts (connectedPoints, otherPoints) have enough faces.
    labelHashSet connectedFaces(2*cFaces.size());
    labelHashSet otherFaces(2*cFaces.size());

    forAllConstIter(Map<label>, pointStatus, iter)
    {
        label pointi = iter.key();

        const labelList& pFaces = mesh().pointFaces()[pointi];

        if (iter() == 1)
        {
            forAll(pFaces, pFacei)
            {
                if (meshTools::faceOnCell(mesh(), celli, pFaces[pFacei]))
                {
                    connectedFaces.insert(pFaces[pFacei]);
                }
            }
        }
        else if (iter() == 2)
        {
            forAll(pFaces, pFacei)
            {
                if (meshTools::faceOnCell(mesh(), celli, pFaces[pFacei]))
                {
                    otherFaces.insert(pFaces[pFacei]);
                }
            }
        }
    }

    if (connectedFaces.size() < 3)
    {
        WarningInFunction
            << "Invalid loop " << loop << " for cell " << celli
            << " since would have too few faces on one side." << nl
            << "All faces:" << cFaces << endl;

        writeOBJ(".", celli, loopPts, connectedPoints);

        return false;
    }

    if (otherFaces.size() < 3)
    {
        WarningInFunction
            << "Invalid loop " << loop << " for cell " << celli
            << " since would have too few faces on one side." << nl
            << "All faces:" << cFaces << endl;

        writeOBJ(".", celli, loopPts, otherPoints);

        return false;
    }


    // Check that faces are split into two regions and not more.
    // When walking across the face the only transition of pointStatus is
    // from set1 to loop to set2 or back. Not allowed is from set1 to loop to
    // set1.
    {
        forAll(cFaces, i)
        {
            label facei = cFaces[i];

            const face& f = mesh().faces()[facei];

            bool hasSet1 = false;
            bool hasSet2 = false;

            label prevStat = pointStatus[f[0]];

            forAll(f, fp)
            {
                label v0 = f[fp];
                label pStat = pointStatus[v0];

                if (pStat == prevStat)
                {
                }
                else if (pStat == 0)
                {
                    // Loop.
                }
                else if (pStat == 1)
                {
                    if (hasSet1)
                    {
                        // Second occurrence of set1.
                        WarningInFunction
                            << "Invalid loop " << loop << " for cell " << celli
                            << " since face " << f << " would be split into"
                            << " more than two faces" << endl;

                        writeOBJ(".", celli, loopPts, otherPoints);

                        return false;
                    }

                    hasSet1 = true;
                }
                else if (pStat == 2)
                {
                    if (hasSet2)
                    {
                        // Second occurrence of set1.
                        WarningInFunction
                            << "Invalid loop " << loop << " for cell " << celli
                            << " since face " << f << " would be split into"
                            << " more than two faces" << endl;

                        writeOBJ(".", celli, loopPts, otherPoints);

                        return false;
                    }

                    hasSet2 = true;
                }
                else
                {
                    FatalErrorInFunction
                        << abort(FatalError);
                }

                prevStat = pStat;


                label v1 = f.nextLabel(fp);
                label edgeI = findEdge(facei, v0, v1);

                label eStat = edgeStatus[edgeI];

                if (eStat == prevStat)
                {
                }
                else if (eStat == 0)
                {
                    // Loop.
                }
                else if (eStat == 1)
                {
                    if (hasSet1)
                    {
                        // Second occurrence of set1.
                        WarningInFunction
                            << "Invalid loop " << loop << " for cell " << celli
                            << " since face " << f << " would be split into"
                            << " more than two faces" << endl;

                        writeOBJ(".", celli, loopPts, otherPoints);

                        return false;
                    }

                    hasSet1 = true;
                }
                else if (eStat == 2)
                {
                    if (hasSet2)
                    {
                        // Second occurrence of set1.
                        WarningInFunction
                            << "Invalid loop " << loop << " for cell " << celli
                            << " since face " << f << " would be split into"
                            << " more than two faces" << endl;

                        writeOBJ(".", celli, loopPts, otherPoints);

                        return false;
                    }

                    hasSet2 = true;
                }
                prevStat = eStat;
            }
        }
    }




    // Check which one of point sets to use.
    bool loopOk = loopAnchorConsistent(celli, loopPts, connectedPoints);

    //if (debug)
    {
        // Additional check: are non-anchor points on other side?
        bool otherLoopOk = loopAnchorConsistent(celli, loopPts, otherPoints);

        if (loopOk == otherLoopOk)
        {
            // Both sets of points are supposedly on the same side as the
            // loop normal. Oops.

            WarningInFunction
                << " For cell:" << celli
                << " achorpoints and nonanchorpoints are geometrically"
                << " on same side!" << endl
                << "cellPoints:" << cPoints << endl
                << "loop:" << loop << endl
                << "anchors:" << connectedPoints << endl
                << "otherPoints:" << otherPoints << endl;

            writeOBJ(".", celli, loopPts, connectedPoints);
        }
    }

    if (loopOk)
    {
        // connectedPoints on 'outside' of loop so these are anchor points
        anchorPoints.transfer(connectedPoints);
        connectedPoints.clear();
    }
    else
    {
        anchorPoints.transfer(otherPoints);
        otherPoints.clear();
    }
    return true;
}


Foam::pointField Foam::cellCuts::loopPoints
(
    const labelList& loop,
    const scalarField& loopWeights
) const
{
    pointField loopPts(loop.size());

    forAll(loop, fp)
    {
        loopPts[fp] = coord(loop[fp], loopWeights[fp]);
    }
    return loopPts;
}


Foam::scalarField Foam::cellCuts::loopWeights(const labelList& loop) const
{
    scalarField weights(loop.size());

    forAll(loop, fp)
    {
        label cut = loop[fp];

        if (isEdge(cut))
        {
            label edgeI = getEdge(cut);

            weights[fp] = edgeWeight_[edgeI];
        }
        else
        {
            weights[fp] = -great;
        }
    }
    return weights;
}


bool Foam::cellCuts::validEdgeLoop
(
    const labelList& loop,
    const scalarField& loopWeights
) const
{
    forAll(loop, fp)
    {
        label cut = loop[fp];

        if (isEdge(cut))
        {
            label edgeI = getEdge(cut);

            // Check: cut compatible only if can be snapped to existing one.
            if (edgeIsCut_[edgeI])
            {
                scalar edgeLen = mesh().edges()[edgeI].mag(mesh().points());

                if
                (
                    mag(loopWeights[fp] - edgeWeight_[edgeI])
                  > geomCellLooper::snapTol()*edgeLen
                )
                {
                    // Positions differ too much->would create two cuts.
                    return false;
                }
            }
        }
    }
    return true;
}


Foam::label Foam::cellCuts::countFaceCuts
(
    const label facei,
    const labelList& loop
) const
{
    // Includes cuts through vertices and through edges.
    // Assumes that if edge is cut both in edgeIsCut and in loop that the
    // position of the cut is the same.

    label nCuts = 0;

    // Count cut vertices
    const face& f = mesh().faces()[facei];

    forAll(f, fp)
    {
        label vertI = f[fp];

        // Vertex already cut or mentioned in current loop.
        if
        (
            pointIsCut_[vertI]
         || (findIndex(loop, vertToEVert(vertI)) != -1)
        )
        {
            nCuts++;
        }
    }

    // Count cut edges.
    const labelList& fEdges = mesh().faceEdges()[facei];

    forAll(fEdges, fEdgeI)
    {
        label edgeI = fEdges[fEdgeI];

        // Edge already cut or mentioned in current loop.
        if
        (
            edgeIsCut_[edgeI]
         || (findIndex(loop, edgeToEVert(edgeI)) != -1)
        )
        {
            nCuts++;
        }
    }

    return nCuts;
}


bool Foam::cellCuts::conservativeValidLoop
(
    const label celli,
    const labelList& loop
) const
{

    if (loop.size() < 2)
    {
        return false;
    }

    forAll(loop, cutI)
    {
        if (isEdge(loop[cutI]))
        {
            label edgeI = getEdge(loop[cutI]);

            if (edgeIsCut_[edgeI])
            {
                // edge compatibility already checked.
            }
            else
            {
                // Quick rejection: vertices of edge itself cannot be cut.
                const edge& e = mesh().edges()[edgeI];

                if (pointIsCut_[e.start()] || pointIsCut_[e.end()])
                {
                    return false;
                }


                // Check faces using this edge
                const labelList& eFaces = mesh().edgeFaces()[edgeI];

                forAll(eFaces, eFacei)
                {
                    label nCuts = countFaceCuts(eFaces[eFacei], loop);

                    if (nCuts > 2)
                    {
                        return false;
                    }
                }
            }
        }
        else
        {
            // Vertex cut

            label vertI = getVertex(loop[cutI]);

            if (!pointIsCut_[vertI])
            {
                // New cut through vertex.

                // Check edges using vertex.
                const labelList& pEdges = mesh().pointEdges()[vertI];

                forAll(pEdges, pEdgeI)
                {
                    label edgeI = pEdges[pEdgeI];

                    if (edgeIsCut_[edgeI])
                    {
                        return false;
                    }
                }

                // Check faces using vertex.
                const labelList& pFaces = mesh().pointFaces()[vertI];

                forAll(pFaces, pFacei)
                {
                    label nCuts = countFaceCuts(pFaces[pFacei], loop);

                    if (nCuts > 2)
                    {
                        return false;
                    }
                }
            }
        }
    }

    // All ok.
    return true;
}


bool Foam::cellCuts::validLoop
(
    const label celli,
    const labelList& loop,
    const scalarField& loopWeights,

    Map<edge>& newFaceSplitCut,
    labelList& anchorPoints
) const
{
    // Determine compatibility of loop with existing cut pattern. Does not use
    // derived cut-addressing (faceCuts), only pointIsCut, edgeIsCut.
    // Adds any cross-cuts found to newFaceSplitCut and sets cell points on
    // one side of the loop in anchorPoints.

    if (loop.size() < 2)
    {
        return false;
    }

    if (debug & 4)
    {
        // Allow as fallback the 'old' loop checking where only a single
        // cut per face is allowed.
        if (!conservativeValidLoop(celli, loop))
        {
            Info << "Invalid conservative loop: " << loop << endl;
            return  false;
        }
    }

    forAll(loop, fp)
    {
        label cut = loop[fp];
        label nextCut = loop[(fp+1) % loop.size()];

        // Label (if any) of face cut (so cut not along existing edge)
        label meshFacei = -1;

        if (isEdge(cut))
        {
            label edgeI = getEdge(cut);

            // Look one cut ahead to find if it is along existing edge.

            if (isEdge(nextCut))
            {
                // From edge to edge -> cross cut
                label nextEdgeI = getEdge(nextCut);

                // Find face and mark as to be split.
                meshFacei = edgeEdgeToFace(celli, edgeI, nextEdgeI);

                if (meshFacei == -1)
                {
                    // Can't find face using both cut edges.
                    return false;
                }
            }
            else
            {
                // From edge to vertex -> cross cut only if no existing edge.

                label nextVertI = getVertex(nextCut);
                const edge& e = mesh().edges()[edgeI];

                if (e.start() != nextVertI && e.end() != nextVertI)
                {
                    // New edge. Find face and mark as to be split.
                    meshFacei = edgeVertexToFace(celli, edgeI, nextVertI);

                    if (meshFacei == -1)
                    {
                        // Can't find face. Ilegal.
                        return false;
                    }
                }
            }
        }
        else
        {
            // Cut is vertex.
            label vertI = getVertex(cut);

            if (isEdge(nextCut))
            {
                // From vertex to edge -> cross cut only if no existing edge.
                label nextEdgeI = getEdge(nextCut);

                const edge& nextE = mesh().edges()[nextEdgeI];

                if (nextE.start() != vertI && nextE.end() != vertI)
                {
                    // Should be cross cut. Find face.
                    meshFacei = edgeVertexToFace(celli, nextEdgeI, vertI);

                    if (meshFacei == -1)
                    {
                        return false;
                    }
                }
            }
            else
            {
                // From vertex to vertex -> cross cut only if no existing edge.
                label nextVertI = getVertex(nextCut);

                if (meshTools::findEdge(mesh(), vertI, nextVertI) == -1)
                {
                    // New edge. Find face.
                    meshFacei = vertexVertexToFace(celli, vertI, nextVertI);

                    if (meshFacei == -1)
                    {
                        return false;
                    }
                }
            }
        }

        if (meshFacei != -1)
        {
            // meshFacei is cut across along cut-nextCut (so not along existing
            // edge). Check if this is compatible with existing pattern.
            edge cutEdge(cut, nextCut);

            Map<edge>::const_iterator iter = faceSplitCut_.find(meshFacei);

            if (iter == faceSplitCut_.end())
            {
                // Face not yet cut so insert.
                newFaceSplitCut.insert(meshFacei, cutEdge);
            }
            else
            {
                // Face already cut. Ok if same edge.
                if (iter() != cutEdge)
                {
                    return false;
                }
            }
        }
    }

    // Is there a face on which all cuts are?
    label faceContainingLoop = loopFace(celli, loop);

    if (faceContainingLoop != -1)
    {
        WarningInFunction
            << "Found loop on cell " << celli << " with all points"
            << " on face " << faceContainingLoop << endl;

        //writeOBJ(".", celli, loopPoints(loop, loopWeights), labelList(0));

        return false;
    }

    // Calculate anchor points
    // Final success is determined by whether anchor points can be determined.
    return calcAnchors
    (
        celli,
        loop,
        loopPoints(loop, loopWeights),
        anchorPoints
    );
}


void Foam::cellCuts::setFromCellLoops()
{
    // 'Uncut' edges/vertices that are not used in loops.
    pointIsCut_ = false;

    edgeIsCut_ = false;

    faceSplitCut_.clear();

    forAll(cellLoops_, celli)
    {
        const labelList& loop = cellLoops_[celli];

        if (loop.size())
        {
            // Storage for cross-face cuts
            Map<edge> faceSplitCuts(loop.size());
            // Storage for points on one side of cell.
            labelList anchorPoints;

            if
            (
               !validLoop
                (
                    celli,
                    loop,
                    loopWeights(loop),
                    faceSplitCuts,
                    anchorPoints
                )
            )
            {
                //writeOBJ(".", celli, loopPoints(celli), anchorPoints);

                WarningInFunction
                    << "Illegal loop " << loop
                    << " when recreating cut-addressing"
                    << " from existing cellLoops for cell " << celli
                    << endl;

                cellLoops_[celli].setSize(0);
                cellAnchorPoints_[celli].setSize(0);
            }
            else
            {
                // Copy anchor points.
                cellAnchorPoints_[celli].transfer(anchorPoints);

                // Copy faceSplitCuts into overall faceSplit info.
                forAllConstIter(Map<edge>, faceSplitCuts, iter)
                {
                    faceSplitCut_.insert(iter.key(), iter());
                }

                // Update edgeIsCut, pointIsCut information
                forAll(loop, cutI)
                {
                    label cut = loop[cutI];

                    if (isEdge(cut))
                    {
                        edgeIsCut_[getEdge(cut)] = true;
                    }
                    else
                    {
                        pointIsCut_[getVertex(cut)] = true;
                    }
                }
            }
        }
    }

    // Reset edge weights
    forAll(edgeIsCut_, edgeI)
    {
        if (!edgeIsCut_[edgeI])
        {
            // Weight not used. Set to illegal value to make any use fall over.
            edgeWeight_[edgeI] = -great;
        }
    }
}


bool Foam::cellCuts::setFromCellLoop
(
    const label celli,
    const labelList& loop,
    const scalarField& loopWeights
)
{
    // Update basic cut information from single cellLoop. Returns true if loop
    // was valid.

    // Dump loop for debugging.
    if (debug)
    {
        OFstream str("last_cell.obj");

        str<< "# edges of cell " << celli << nl;

        meshTools::writeOBJ
        (
            str,
            mesh().cells(),
            mesh().faces(),
            mesh().points(),
            labelList(1, celli)
        );


        OFstream loopStr("last_loop.obj");

        loopStr<< "# looppoints for cell " << celli << nl;

        pointField pointsOfLoop = loopPoints(loop, loopWeights);

        forAll(pointsOfLoop, i)
        {
            meshTools::writeOBJ(loopStr, pointsOfLoop[i]);
        }

        str << 'l';

        forAll(pointsOfLoop, i)
        {
            str << ' ' << i + 1;
        }
        str << ' ' << 1 << nl;
    }

    bool okLoop = false;

    if (validEdgeLoop(loop, loopWeights))
    {
        // Storage for cuts across face
        Map<edge> faceSplitCuts(loop.size());
        // Storage for points on one side of cell
        labelList anchorPoints;

        okLoop =
            validLoop(celli, loop, loopWeights, faceSplitCuts, anchorPoints);

        if (okLoop)
        {
            // Valid loop. Copy cellLoops and anchorPoints
            cellLoops_[celli] = loop;
            cellAnchorPoints_[celli].transfer(anchorPoints);

            // Copy split cuts
            forAllConstIter(Map<edge>, faceSplitCuts, iter)
            {
                faceSplitCut_.insert(iter.key(), iter());
            }


            // Update edgeIsCut, pointIsCut information
            forAll(loop, cutI)
            {
                label cut = loop[cutI];

                if (isEdge(cut))
                {
                    label edgeI = getEdge(cut);

                    edgeIsCut_[edgeI] = true;

                    edgeWeight_[edgeI] = loopWeights[cutI];
                }
                else
                {
                    label vertI = getVertex(cut);

                    pointIsCut_[vertI] = true;
                }
            }
        }
    }

    return okLoop;
}


void Foam::cellCuts::setFromCellLoops
(
    const labelList& cellLabels,
    const labelListList& cellLoops,
    const List<scalarField>& cellLoopWeights
)
{
    // 'Uncut' edges/vertices that are not used in loops.
    pointIsCut_ = false;

    edgeIsCut_ = false;

    forAll(cellLabels, cellLabelI)
    {
        label celli = cellLabels[cellLabelI];

        const labelList& loop = cellLoops[cellLabelI];

        if (loop.size())
        {
            const scalarField& loopWeights = cellLoopWeights[cellLabelI];

            if (setFromCellLoop(celli, loop, loopWeights))
            {
                // Valid loop. Call above will have upated all already.
            }
            else
            {
                // Clear cellLoops
                cellLoops_[celli].setSize(0);
            }
        }
    }
}


void Foam::cellCuts::setFromCellCutter
(
    const cellLooper& cellCutter,
    const List<refineCell>& refCells
)
{
    // 'Uncut' edges/vertices that are not used in loops.
    pointIsCut_ = false;

    edgeIsCut_ = false;

    // storage for loop of cuts (cut vertices and/or cut edges)
    labelList cellLoop;
    scalarField cellLoopWeights;

    // For debugging purposes
    DynamicList<label> invalidCutCells(2);
    DynamicList<labelList> invalidCutLoops(2);
    DynamicList<scalarField> invalidCutLoopWeights(2);


    forAll(refCells, refCelli)
    {
        const refineCell& refCell = refCells[refCelli];

        label celli = refCell.cellNo();

        const vector& refDir = refCell.direction();


        // Cut cell. Determines cellLoop and cellLoopWeights
        bool goodCut =
            cellCutter.cut
            (
                refDir,
                celli,

                pointIsCut_,
                edgeIsCut_,
                edgeWeight_,

                cellLoop,
                cellLoopWeights
            );

        // Check whether edge refinement is on a per face basis compatible with
        // current pattern.
        if (goodCut)
        {
            if (setFromCellLoop(celli, cellLoop, cellLoopWeights))
            {
                // Valid loop. Will have updated all info already.
            }
            else
            {
                cellLoops_[celli].setSize(0);

                WarningInFunction
                    << "Found loop on cell " << celli
                    << " that resulted in an unexpected bad cut." << nl
                    << "    Suggestions:" << nl
                    << "      - Turn on the debug switch for 'cellCuts' to get"
                    << " geometry files that identify this cell." << nl
                    << "      - Also keep in mind to check the defined"
                    << " reference directions, as these are most likely the"
                    << " origin of the problem."
                    << nl << endl;

                // Discarded by validLoop
                if (debug)
                {
                    invalidCutCells.append(celli);
                    invalidCutLoops.append(cellLoop);
                    invalidCutLoopWeights.append(cellLoopWeights);
                }
            }
        }
        else
        {
            // Clear cellLoops
            cellLoops_[celli].setSize(0);
        }
    }

    if (debug && invalidCutCells.size())
    {
        invalidCutCells.shrink();
        invalidCutLoops.shrink();
        invalidCutLoopWeights.shrink();

        fileName cutsFile("invalidLoopCells.obj");

        Pout<< "cellCuts : writing inValidLoops cells to " << cutsFile << endl;

        OFstream cutsStream(cutsFile);

        meshTools::writeOBJ
        (
            cutsStream,
            mesh().cells(),
            mesh().faces(),
            mesh().points(),
            invalidCutCells
        );

        fileName loopsFile("invalidLoops.obj");

        Pout<< "cellCuts : writing inValidLoops loops to " << loopsFile << endl;

        OFstream loopsStream(loopsFile);

        label vertI = 0;

        forAll(invalidCutLoops, i)
        {
            writeOBJ
            (
                loopsStream,
                loopPoints(invalidCutLoops[i], invalidCutLoopWeights[i]),
                vertI
            );
        }
    }
}


void Foam::cellCuts::setFromCellCutter
(
    const cellLooper& cellCutter,
    const labelList& cellLabels,
    const List<plane>& cellCutPlanes
)
{
    // 'Uncut' edges/vertices that are not used in loops.
    pointIsCut_ = false;

    edgeIsCut_ = false;

    // storage for loop of cuts (cut vertices and/or cut edges)
    labelList cellLoop;
    scalarField cellLoopWeights;

    // For debugging purposes
    DynamicList<label> invalidCutCells(2);
    DynamicList<labelList> invalidCutLoops(2);
    DynamicList<scalarField> invalidCutLoopWeights(2);


    forAll(cellLabels, i)
    {
        label celli = cellLabels[i];

        // Cut cell. Determines cellLoop and cellLoopWeights
        bool goodCut =
            cellCutter.cut
            (
                cellCutPlanes[i],
                celli,

                pointIsCut_,
                edgeIsCut_,
                edgeWeight_,

                cellLoop,
                cellLoopWeights
            );

        // Check whether edge refinement is on a per face basis compatible with
        // current pattern.
        if (goodCut)
        {
            if (setFromCellLoop(celli, cellLoop, cellLoopWeights))
            {
                // Valid loop. Will have updated all info already.
            }
            else
            {
                cellLoops_[celli].setSize(0);

                // Discarded by validLoop
                if (debug)
                {
                    invalidCutCells.append(celli);
                    invalidCutLoops.append(cellLoop);
                    invalidCutLoopWeights.append(cellLoopWeights);
                }
            }
        }
        else
        {
            // Clear cellLoops
            cellLoops_[celli].setSize(0);
        }
    }

    if (debug && invalidCutCells.size())
    {
        invalidCutCells.shrink();
        invalidCutLoops.shrink();
        invalidCutLoopWeights.shrink();

        fileName cutsFile("invalidLoopCells.obj");

        Pout<< "cellCuts : writing inValidLoops cells to " << cutsFile << endl;

        OFstream cutsStream(cutsFile);

        meshTools::writeOBJ
        (
            cutsStream,
            mesh().cells(),
            mesh().faces(),
            mesh().points(),
            invalidCutCells
        );

        fileName loopsFile("invalidLoops.obj");

        Pout<< "cellCuts : writing inValidLoops loops to " << loopsFile << endl;

        OFstream loopsStream(loopsFile);

        label vertI = 0;

        forAll(invalidCutLoops, i)
        {
            writeOBJ
            (
                loopsStream,
                loopPoints(invalidCutLoops[i], invalidCutLoopWeights[i]),
                vertI
            );
        }
    }
}


void Foam::cellCuts::orientPlanesAndLoops()
{
    // Determine anchorPoints if not yet done by validLoop.
    forAll(cellLoops_, celli)
    {
        const labelList& loop = cellLoops_[celli];

        if (loop.size() && cellAnchorPoints_[celli].empty())
        {
            // Leave anchor points empty if illegal loop.
            calcAnchors
            (
                celli,
                loop,
                loopPoints(celli),
                cellAnchorPoints_[celli]
            );
        }
    }

    if (debug & 2)
    {
        Pout<< "cellAnchorPoints:" << endl;
    }
    forAll(cellAnchorPoints_, celli)
    {
        if (cellLoops_[celli].size())
        {
            if (cellAnchorPoints_[celli].empty())
            {
                FatalErrorInFunction
                    << "No anchor points for cut cell " << celli << endl
                    << "cellLoop:" << cellLoops_[celli] << abort(FatalError);
            }

            if (debug & 2)
            {
                Pout<< "    cell:" << celli << " anchored at "
                    << cellAnchorPoints_[celli] << endl;
            }
        }
    }

    // Calculate number of valid cellLoops
    nLoops_ = 0;

    forAll(cellLoops_, celli)
    {
        if (cellLoops_[celli].size())
        {
            nLoops_++;
        }
    }
}


void Foam::cellCuts::calcLoopsAndAddressing(const labelList& cutCells)
{
    // Sanity check on weights
    forAll(edgeIsCut_, edgeI)
    {
        if (edgeIsCut_[edgeI])
        {
            scalar weight = edgeWeight_[edgeI];

            if (weight < 0 || weight > 1)
            {
                FatalErrorInFunction
                    << "Weight out of range [0,1]. Edge " << edgeI
                    << " verts:" << mesh().edges()[edgeI]
                    << " weight:" << weight << abort(FatalError);
            }
        }
        else
        {
            // Weight not used. Set to illegal value to make any use fall over.
            edgeWeight_[edgeI] = -great;
        }
    }


    // Calculate faces that split cells in two
    calcCellLoops(cutCells);

    if (debug & 2)
    {
        Pout<< "-- cellLoops --" << endl;
        forAll(cellLoops_, celli)
        {
            const labelList& loop = cellLoops_[celli];

            if (loop.size())
            {
                Pout<< "cell:" << celli << "  ";
                writeCuts(Pout, loop, loopWeights(loop));
                Pout<< endl;
            }
        }
    }

    // Redo basic cut information (pointIsCut, edgeIsCut, faceSplitCut)
    // using cellLoop only.
    setFromCellLoops();
}


void Foam::cellCuts::check() const
{
    // Check weights for unsnapped values
    forAll(edgeIsCut_, edgeI)
    {
        if (edgeIsCut_[edgeI])
        {
            if
            (
                edgeWeight_[edgeI] < geomCellLooper::snapTol()
             || edgeWeight_[edgeI] > (1 - geomCellLooper::snapTol())
            )
            {
                // Should have been snapped.
                WarningInFunction
                    << "edge:" << edgeI << " vertices:"
                    << mesh().edges()[edgeI]
                    << " weight:" << edgeWeight_[edgeI] << " should have been"
                    << " snapped to one of its endpoints"
                    << endl;
            }
        }
        else
        {
            if (edgeWeight_[edgeI] > - 1)
            {
                FatalErrorInFunction
                    << "edge:" << edgeI << " vertices:"
                    << mesh().edges()[edgeI]
                    << " weight:" << edgeWeight_[edgeI] << " is not cut but"
                    << " its weight is not marked invalid"
                    << abort(FatalError);
            }
        }
    }

    // Check that all elements of cellloop are registered
    forAll(cellLoops_, celli)
    {
        const labelList& loop = cellLoops_[celli];

        forAll(loop, i)
        {
            label cut = loop[i];

            if
            (
                (isEdge(cut) && !edgeIsCut_[getEdge(cut)])
             || (!isEdge(cut) && !pointIsCut_[getVertex(cut)])
            )
            {
                labelList cuts(1, cut);
                writeCuts(Pout, cuts, loopWeights(cuts));

                FatalErrorInFunction
                    << "cell:" << celli << " loop:"
                    << loop
                    << " cut:" << cut << " is not marked as cut"
                    << abort(FatalError);
            }
        }
    }

    // Check that no elements of cell loop are anchor point.
    forAll(cellLoops_, celli)
    {
        const labelList& anchors = cellAnchorPoints_[celli];

        const labelList& loop = cellLoops_[celli];

        if (loop.size() && anchors.empty())
        {
            FatalErrorInFunction
                << "cell:" << celli << " loop:" << loop
                << " has no anchor points"
                << abort(FatalError);
        }


        forAll(loop, i)
        {
            label cut = loop[i];

            if
            (
                !isEdge(cut)
             && findIndex(anchors, getVertex(cut)) != -1
            )
            {
                FatalErrorInFunction
                    << "cell:" << celli << " loop:" << loop
                    << " anchor points:" << anchors
                    << " anchor:" << getVertex(cut) << " is part of loop"
                    << abort(FatalError);
            }
        }
    }


    // Check that cut faces have a neighbour that is cut.
    boolList nbrCellIsCut;
    {
        boolList cellIsCut(mesh().nCells(), false);
        forAll(cellLoops_, celli)
        {
            cellIsCut[celli] = cellLoops_[celli].size();
        }
        syncTools::swapBoundaryCellList(mesh(), cellIsCut, nbrCellIsCut);
    }

    forAllConstIter(Map<edge>, faceSplitCut_, iter)
    {
        label facei = iter.key();

        if (mesh().isInternalFace(facei))
        {
            label own = mesh().faceOwner()[facei];
            label nei = mesh().faceNeighbour()[facei];

            if (cellLoops_[own].empty() && cellLoops_[nei].empty())
            {
                FatalErrorInFunction
                    << "Internal face:" << facei << " cut by " << iter()
                    << " has owner:" << own
                    << " and neighbour:" << nei
                    << " that are both uncut"
                    << abort(FatalError);
            }
        }
        else
        {
            label bFacei = facei - mesh().nInternalFaces();

            label own = mesh().faceOwner()[facei];

            if (cellLoops_[own].empty() && !nbrCellIsCut[bFacei])
            {
                FatalErrorInFunction
                    << "Boundary face:" << facei << " cut by " << iter()
                    << " has owner:" << own
                    << " that is uncut"
                    << abort(FatalError);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellCuts::cellCuts
(
    const polyMesh& mesh,
    const labelList& cutCells,
    const labelList& meshVerts,
    const labelList& meshEdges,
    const scalarField& meshEdgeWeights
)
:
    edgeVertex(mesh),
    pointIsCut_(expand(mesh.nPoints(), meshVerts)),
    edgeIsCut_(expand(mesh.nEdges(), meshEdges)),
    edgeWeight_(expand(mesh.nEdges(), meshEdges, meshEdgeWeights)),
    faceSplitCut_(cutCells.size()),
    cellLoops_(mesh.nCells()),
    nLoops_(-1),
    cellAnchorPoints_(mesh.nCells())
{
    if (debug)
    {
        Pout<< "cellCuts : constructor from cut verts and edges" << endl;
    }

    calcLoopsAndAddressing(cutCells);

    // Calculate planes and flip cellLoops if necessary
    orientPlanesAndLoops();

    if (debug)
    {
        check();
    }

    clearOut();

    if (debug)
    {
        Pout<< "cellCuts : leaving constructor from cut verts and edges"
            << endl;
    }
}


Foam::cellCuts::cellCuts
(
    const polyMesh& mesh,
    const labelList& meshVerts,
    const labelList& meshEdges,
    const scalarField& meshEdgeWeights
)
:
    edgeVertex(mesh),
    pointIsCut_(expand(mesh.nPoints(), meshVerts)),
    edgeIsCut_(expand(mesh.nEdges(), meshEdges)),
    edgeWeight_(expand(mesh.nEdges(), meshEdges, meshEdgeWeights)),
    faceSplitCut_(mesh.nFaces()/10 + 1),
    cellLoops_(mesh.nCells()),
    nLoops_(-1),
    cellAnchorPoints_(mesh.nCells())
{
    // Construct from pattern of cuts. Finds out itself which cells are cut.
    // (can go wrong if e.g. all neighbours of cell are refined)

    if (debug)
    {
        Pout<< "cellCuts : constructor from cellLoops" << endl;
    }

    calcLoopsAndAddressing(identity(mesh.nCells()));

    // Adds cuts on other side of coupled boundaries
    syncProc();

    // Calculate planes and flip cellLoops if necessary
    orientPlanesAndLoops();

    if (debug)
    {
        check();
    }

    clearOut();

    if (debug)
    {
        Pout<< "cellCuts : leaving constructor from cellLoops" << endl;
    }
}


Foam::cellCuts::cellCuts
(
    const polyMesh& mesh,
    const labelList& cellLabels,
    const labelListList& cellLoops,
    const List<scalarField>& cellEdgeWeights
)
:
    edgeVertex(mesh),
    pointIsCut_(mesh.nPoints(), false),
    edgeIsCut_(mesh.nEdges(), false),
    edgeWeight_(mesh.nEdges(), -great),
    faceSplitCut_(cellLabels.size()),
    cellLoops_(mesh.nCells()),
    nLoops_(-1),
    cellAnchorPoints_(mesh.nCells())
{
    if (debug)
    {
        Pout<< "cellCuts : constructor from cellLoops" << endl;
    }

    // Update pointIsCut, edgeIsCut, faceSplitCut from cell loops.
    // Makes sure cuts are consistent
    setFromCellLoops(cellLabels, cellLoops, cellEdgeWeights);

    // Adds cuts on other side of coupled boundaries
    syncProc();

    // Calculate planes and flip cellLoops if necessary
    orientPlanesAndLoops();

    if (debug)
    {
        check();
    }

    clearOut();

    if (debug)
    {
        Pout<< "cellCuts : leaving constructor from cellLoops" << endl;
    }
}


Foam::cellCuts::cellCuts
(
    const polyMesh& mesh,
    const cellLooper& cellCutter,
    const List<refineCell>& refCells
)
:
    edgeVertex(mesh),
    pointIsCut_(mesh.nPoints(), false),
    edgeIsCut_(mesh.nEdges(), false),
    edgeWeight_(mesh.nEdges(), -great),
    faceSplitCut_(refCells.size()),
    cellLoops_(mesh.nCells()),
    nLoops_(-1),
    cellAnchorPoints_(mesh.nCells())
{
    if (debug)
    {
        Pout<< "cellCuts : constructor from cellCutter" << endl;
    }

    // Update pointIsCut, edgeIsCut, faceSplitCut from cell loops.
    // Makes sure cuts are consistent
    setFromCellCutter(cellCutter, refCells);

    // Adds cuts on other side of coupled boundaries
    syncProc();

    // Calculate planes and flip cellLoops if necessary
    orientPlanesAndLoops();

    if (debug)
    {
        check();
    }

    clearOut();

    if (debug)
    {
        Pout<< "cellCuts : leaving constructor from cellCutter" << endl;
    }
}


Foam::cellCuts::cellCuts
(
    const polyMesh& mesh,
    const cellLooper& cellCutter,
    const labelList& cellLabels,
    const List<plane>& cutPlanes
)
:
    edgeVertex(mesh),
    pointIsCut_(mesh.nPoints(), false),
    edgeIsCut_(mesh.nEdges(), false),
    edgeWeight_(mesh.nEdges(), -great),
    faceSplitCut_(cellLabels.size()),
    cellLoops_(mesh.nCells()),
    nLoops_(-1),
    cellAnchorPoints_(mesh.nCells())
{
    if (debug)
    {
        Pout<< "cellCuts : constructor from cellCutter with prescribed plane"
            << endl;
    }

    // Update pointIsCut, edgeIsCut, faceSplitCut from cell loops.
    // Makes sure cuts are consistent
    setFromCellCutter(cellCutter, cellLabels, cutPlanes);

    // Adds cuts on other side of coupled boundaries
    syncProc();

    // Calculate planes and flip cellLoops if necessary
    orientPlanesAndLoops();

    if (debug)
    {
        check();
    }

    clearOut();

    if (debug)
    {
        Pout<< "cellCuts : leaving constructor from cellCutter with prescribed"
            << " plane" << endl;
    }
}


Foam::cellCuts::cellCuts
(
    const polyMesh& mesh,
    const boolList& pointIsCut,
    const boolList& edgeIsCut,
    const scalarField& edgeWeight,
    const Map<edge>& faceSplitCut,
    const labelListList& cellLoops,
    const label nLoops,
    const labelListList& cellAnchorPoints
)
:
    edgeVertex(mesh),
    pointIsCut_(pointIsCut),
    edgeIsCut_(edgeIsCut),
    edgeWeight_(edgeWeight),
    faceSplitCut_(faceSplitCut),
    cellLoops_(cellLoops),
    nLoops_(nLoops),
    cellAnchorPoints_(cellAnchorPoints)
{
    if (debug)
    {
        Pout<< "cellCuts : constructor from components" << endl;
        Pout<< "cellCuts : leaving constructor from components" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellCuts::~cellCuts()
{
    clearOut();
}


void Foam::cellCuts::clearOut()
{
    faceCutsPtr_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::pointField Foam::cellCuts::loopPoints(const label celli) const
{
    const labelList& loop = cellLoops_[celli];

    pointField loopPts(loop.size());

    forAll(loop, fp)
    {
        label cut = loop[fp];

        if (isEdge(cut))
        {
            label edgeI = getEdge(cut);

            loopPts[fp] = coord(cut, edgeWeight_[edgeI]);
        }
        else
        {
            loopPts[fp] = coord(cut, -great);
        }
    }
    return loopPts;
}


void Foam::cellCuts::flip(const label celli)
{
    labelList& loop = cellLoops_[celli];

    reverse(loop);

    // Reverse anchor point set.
    cellAnchorPoints_[celli] =
        nonAnchorPoints
        (
            mesh().cellPoints()[celli],
            cellAnchorPoints_[celli],
            loop
        );
}


void Foam::cellCuts::flipLoopOnly(const label celli)
{
    labelList& loop = cellLoops_[celli];

    reverse(loop);
}


void Foam::cellCuts::writeOBJ
(
    Ostream& os,
    const pointField& loopPts,
    label& vertI
) const
{
    label startVertI = vertI;

    forAll(loopPts, fp)
    {
        const point& pt = loopPts[fp];

        os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;

        vertI++;
    }

    os  << 'f';
    forAll(loopPts, fp)
    {
        os  << ' ' << startVertI + fp + 1;
    }
    os  << endl;
}


void Foam::cellCuts::writeOBJ(Ostream& os) const
{
    label vertI = 0;

    forAll(cellLoops_, celli)
    {
        writeOBJ(os, loopPoints(celli), vertI);
    }
}


void Foam::cellCuts::writeCellOBJ(const fileName& dir, const label celli) const
{
    const labelList& anchors = cellAnchorPoints_[celli];

    writeOBJ(dir, celli, loopPoints(celli), anchors);
}


// ************************************************************************* //
