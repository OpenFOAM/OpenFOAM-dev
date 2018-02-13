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

#include "topoCellLooper.H"
#include "cellFeatures.H"
#include "polyMesh.H"
#include "unitConversion.H"
#include "DynamicList.H"
#include "ListOps.H"
#include "meshTools.H"
#include "hexMatcher.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
   defineTypeNameAndDebug(topoCellLooper, 0);
   addToRunTimeSelectionTable(cellLooper, topoCellLooper, word);
}

// Angle for polys to be considered splitHexes.
const Foam::scalar Foam::topoCellLooper::featureCos = Foam::cos(degToRad(10.0));


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// In-memory truncate a list
template<class T>
void Foam::topoCellLooper::subsetList
(
    const label startI,
    const label freeI,
    DynamicList<T>& lst
)
{
    if (startI == 0)
    {
        // Truncate (setSize decides itself not to do anything if nothing
        // changed)
        if (freeI < 0)
        {
            FatalErrorInFunction
                << "  lst:" << lst << abort(FatalError);
        }
        lst.setCapacity(freeI);
    }
    else
    {
        // Shift elements down
        label newI = 0;
        for (label elemI = startI; elemI < freeI; elemI++)
        {
            lst[newI++] = lst[elemI];
        }

        if ((freeI - startI) < 0)
        {
            FatalErrorInFunction
                << "  lst:" << lst << abort(FatalError);
        }

        lst.setCapacity(freeI - startI);
    }
}


// Returns label of edge nFeaturePts away from startEdge (in the direction of
// startVertI) and not counting non-featurePoints.
//
// When stepping to this face it can happen in 3 different ways:
//
//  --|------
//    |
//  1 |   0
//    |A
//    |
//    |
//  --|------
//
//  A: jumping from face0 to face1 across edge A.
//      startEdge != -1
//      startVert == -1
//
//  --|------
//    |
//  1 |   0
//    +B
//    |
//    |
//  --|------
//
//  B: jumping from face0 to face1 across (non-feature) point B
//      startEdge == -1
//      startVert != -1
//
//  --|------
//  0 |   1
//    |C
//  --+
//    |
//    |
//  --|------
//
//  C: jumping from face0 to face1 across (feature) point C on edge.
//      startEdge != -1
//      startVert != -1
//
void Foam::topoCellLooper::walkFace
(
    const cellFeatures& features,
    const label facei,
    const label startEdgeI,
    const label startVertI,
    const label nFeaturePts,

    label& edgeI,
    label& vertI
) const
{
    const labelList& fEdges = mesh().faceEdges()[facei];

    edgeI = startEdgeI;

    vertI = startVertI;

    // Number of feature points crossed so far
    label nVisited = 0;

    if (vertI == -1)
    {
        // Started on edge. Go to one of its endpoints.
        vertI = mesh().edges()[edgeI].start();

        if (features.isFeatureVertex(facei, vertI))
        {
            nVisited++;
        }
    }

    if ((edgeI == -1) || !meshTools::edgeOnFace(mesh(), facei, edgeI))
    {
        // Either edge is not set or not on current face.  Just take one of
        // the edges on this face as starting edge.
        edgeI = getFirstVertEdge(facei, vertI);
    }

    // Now we should have starting edge on face and a vertex on that edge.

    do
    {

        edgeI = meshTools::otherEdge(mesh(), fEdges, edgeI, vertI);

        if (nVisited == nFeaturePts)
        {
            break;
        }

        vertI = mesh().edges()[edgeI].otherVertex(vertI);

        if (features.isFeatureVertex(facei, vertI))
        {
            nVisited++;
        }
    }
    while (true);
}


// Returns list of vertices on 'superEdge' i.e. list of edges connected by
// non-feature points. First and last are feature points, ones in between are
// not.
Foam::labelList Foam::topoCellLooper::getSuperEdge
(
    const cellFeatures& features,
    const label facei,
    const label startEdgeI,
    const label startVertI
) const
{
    const labelList& fEdges = mesh().faceEdges()[facei];

    labelList superVerts(fEdges.size());
    label superVertI = 0;


    label edgeI = startEdgeI;

    label vertI = startVertI;

    superVerts[superVertI++] = vertI;

    label prevEdgeI = -1;

    do
    {
        vertI = mesh().edges()[edgeI].otherVertex(vertI);

        superVerts[superVertI++] = vertI;

        prevEdgeI = edgeI;

        edgeI = meshTools::otherEdge(mesh(), fEdges, edgeI, vertI);
    }
    while (!features.isFeaturePoint(prevEdgeI, edgeI));

    superVerts.setSize(superVertI);

    return superVerts;
}


// Return non-feature edge from cells' edges that is most perpendicular
// to refinement direction.
Foam::label Foam::topoCellLooper::getAlignedNonFeatureEdge
(
    const vector& refDir,
    const label celli,
    const cellFeatures& features
) const
{
    const labelList& cEdges = mesh().cellEdges()[celli];

    const point& ctr = mesh().cellCentres()[celli];

    label cutEdgeI = -1;
    scalar maxCos = -great;

    forAll(cEdges, cEdgeI)
    {
        label edgeI = cEdges[cEdgeI];

        if (!features.isFeatureEdge(edgeI))
        {
            const edge& e = mesh().edges()[edgeI];

            // Get plane spanned by e.start, e.end and cell centre.
            vector e0 = mesh().points()[e.start()] - ctr;
            vector e1 = mesh().points()[e.end()] - ctr;

            vector n = e0 ^ e1;
            n /= mag(n);

            scalar cosAngle = mag(refDir & n);

            if (cosAngle > maxCos)
            {
                maxCos = cosAngle;

                cutEdgeI = edgeI;
            }
        }
    }

    return cutEdgeI;
}


// Starts from edge and vertex on edge on face (or neighbouring face)
// and steps either to existing vertex (vertI != -1) or to edge (vertI == -1)
// by walking point-edge and crossing nFeats featurePoints.
void Foam::topoCellLooper::walkAcrossFace
(
    const cellFeatures& features,
    const label facei,
    const label startEdgeI,
    const label startVertI,
    const label nFeats,

    label& edgeI,
    label& vertI
) const
{
    label oppositeVertI = -1;
    label oppositeEdgeI = -1;

    // Go to oppositeEdge and a vertex on that.
    walkFace
    (
        features,
        facei,
        startEdgeI,
        startVertI,
        nFeats,

        oppositeEdgeI,
        oppositeVertI
    );

    // Loop over super edge to find internal points if there are any.

    labelList superEdge =
        getSuperEdge
        (
            features,
            facei,
            oppositeEdgeI,
            oppositeVertI
        );

    label sz = superEdge.size();

    if (sz == 2)
    {
        // No non-feature point in between feature points.
        // Mark edge.

        vertI = -1;
        edgeI = oppositeEdgeI;
    }
    else if (sz == 3)
    {
        vertI = superEdge[1];
        edgeI = -1;
    }
    else
    {
        //Should choose acc. to geometry!
        label index = sz/2;

        if (debug)
        {
            Pout<< "    Don't know what to do. Stepped to non-feature point "
                << "at index " << index << " in superEdge:" << superEdge
                << endl;
        }

        vertI = superEdge[index];
        edgeI = -1;
    }
}


// Walks cell circumference. Updates face, edge, vertex.
//
// Position on face is given by:
//
//  vertI == -1, facei != -1, edgeI != -1
//      on edge of face. Cross edge to neighbouring face.
//
//  vertI != -1, edgeI != -1, facei == -1
//      coming from edge onto vertex vertI. Need to step to one
//      of the faces not using edgeI.
//
//  vertI != -1, edgeI == -1, facei != -1
//      coming from vertex on side of face. Step to one of the faces
//      using vertI but not facei
void Foam::topoCellLooper::walkSplitHex
(
    const label celli,
    const cellFeatures& features,
    const label fromFacei,
    const label fromEdgeI,
    const label fromVertI,

    DynamicList<label>& loop,
    DynamicList<scalar>& loopWeights
) const
{
    // Work vars giving position on cell
    label facei = fromFacei;
    label edgeI = fromEdgeI;
    label vertI = fromVertI;

    do
    {
        if (debug)
        {
            Pout<< "Entering walk with : cell:" << celli << " face:" << facei;
            if (facei != -1)
            {
                Pout<< " verts:" << mesh().faces()[facei];
            }
            Pout<< " edge:" << edgeI;
            if (edgeI != -1)
            {
                Pout<< " verts:" << mesh().edges()[edgeI];
            }
            Pout<< " vert:" << vertI << endl;
        }

        label startLoop = -1;

        if
        (
            (vertI != -1)
         && (
                (startLoop =
                    findIndex
                    (
                        loop,
                        vertToEVert(vertI)
                    )
                )
            != -1
            )
        )
        {
            // Breaking walk since vertI already cut
            label firstFree = loop.size();

            subsetList(startLoop, firstFree, loop);
            subsetList(startLoop, firstFree, loopWeights);

            break;
        }
        if
        (
            (edgeI != -1)
         && (
                (startLoop =
                    findIndex
                    (
                        loop,
                        edgeToEVert(edgeI)
                    )
                )
             != -1
            )
        )
        {
            // Breaking walk since edgeI already cut
            label firstFree = loop.size();

            subsetList(startLoop, firstFree, loop);
            subsetList(startLoop, firstFree, loopWeights);

            break;
        }


        if (vertI == -1)
        {
            // On edge
            if (edgeI == -1)
            {
                FatalErrorInFunction
                    << abort(FatalError);
            }

            loop.append(edgeToEVert(edgeI));
            loopWeights.append(0.5);

            // Cross edge to next face
            facei = meshTools::otherFace(mesh(), celli, facei, edgeI);

            if (debug)
            {
                Pout<< "    stepped across edge " << mesh().edges()[edgeI]
                    << " to face " << facei << " verts:"
                    << mesh().faces()[facei] << endl;
            }

            label nextEdgeI = -1;
            label nextVertI = -1;

            // Walk face along its edges
            walkAcrossFace
            (
                features,
                facei,
                edgeI,
                vertI,
                2,

                nextEdgeI,
                nextVertI
            );

            edgeI = nextEdgeI;
            vertI = nextVertI;
        }
        else
        {
            // On vertex.

            loop.append(vertToEVert(vertI));
            loopWeights.append(-great);

            if (edgeI == -1)
            {
                // Normal vertex on edge of face. Get edges connected to it
                // which are not on facei.
                labelList nextEdges = getVertEdgesNonFace
                (
                    celli,
                    facei,
                    vertI
                );

                if (nextEdges.empty())
                {
                    // Cross to other face (there is only one since no edges)
                    const labelList& pFaces = mesh().pointFaces()[vertI];

                    forAll(pFaces, pFacei)
                    {
                        label thisFacei = pFaces[pFacei];

                        if
                        (
                            (thisFacei != facei)
                         && meshTools::faceOnCell(mesh(), celli, thisFacei)
                        )
                        {
                            facei = thisFacei;
                            break;
                        }
                    }

                    if (debug)
                    {
                        Pout<< "    stepped from non-edge vertex " << vertI
                            << " to face " << facei << " verts:"
                            << mesh().faces()[facei]
                            << " since candidate edges:" << nextEdges << endl;
                    }

                    label nextEdgeI = -1;
                    label nextVertI = -1;

                    walkAcrossFace
                    (
                        features,
                        facei,
                        edgeI,
                        vertI,
                        2,          // 2 vertices to cross

                        nextEdgeI,
                        nextVertI
                    );

                    edgeI = nextEdgeI;
                    vertI = nextVertI;
                }
                else if (nextEdges.size() == 1)
                {
                    // One edge. Go along this one.
                    edgeI = nextEdges[0];

                    if (debug)
                    {
                        Pout<< "    stepped from non-edge vertex " << vertI
                            << " along edge " << edgeI << " verts:"
                            << mesh().edges()[edgeI]
                            << " out of candidate edges:"
                            << nextEdges << endl;
                    }

                    vertI = mesh().edges()[edgeI].otherVertex(vertI);

                    facei = -1;
                }
                else
                {
                    // Multiple edges to choose from. Pick any one.
                    // (ideally should be geometric)

                    label index = nextEdges.size()/2;

                    edgeI = nextEdges[index];

                    if (debug)
                    {
                        Pout<< "    stepped from non-edge vertex " << vertI
                            << " along edge " << edgeI << " verts:"
                            << mesh().edges()[edgeI]
                            << " out of candidate edges:" << nextEdges << endl;
                    }

                    vertI = mesh().edges()[edgeI].otherVertex(vertI);

                    facei = -1;
                }
            }
            else
            {
                // Get faces connected to startVertI but not startEdgeI
                labelList nextFaces =
                    getVertFacesNonEdge
                    (
                        celli,
                        edgeI,
                        vertI
                    );

                if (nextFaces.size() == 1)
                {
                    // Only one face to cross.
                    facei = nextFaces[0];

                    label nextEdgeI = -1;
                    label nextVertI = -1;

                    walkAcrossFace
                    (
                        features,
                        facei,
                        edgeI,
                        vertI,
                        2,          // 2 vertices to cross

                        nextEdgeI,
                        nextVertI
                    );

                    edgeI = nextEdgeI;
                    vertI = nextVertI;
                }
                else if (nextFaces.size() == 2)
                {
                    // Split face. Get edge in between.
                    facei = -1;

                    edgeI =
                        meshTools::getSharedEdge
                        (
                            mesh(),
                            nextFaces[0],
                            nextFaces[1]
                        );

                    vertI = mesh().edges()[edgeI].otherVertex(vertI);
                }
                else
                {
                    FatalErrorInFunction
                        << "Choosing from more than "
                        << "two candidates:" << nextFaces
                        << " when coming from vertex " << vertI << " on cell "
                        << celli << abort(FatalError);
                }
            }
        }

        if (debug)
        {
            Pout<< "Walked to : face:" << facei;
            if (facei != -1)
            {
                Pout<< " verts:" << mesh().faces()[facei];
            }
            Pout<< " edge:" << edgeI;
            if (edgeI != -1)
            {
                Pout<< " verts:" << mesh().edges()[edgeI];
            }
            Pout<< " vert:" << vertI << endl;
        }
    }
    while (true);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::topoCellLooper::topoCellLooper(const polyMesh& mesh)
:
    hexCellLooper(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::topoCellLooper::~topoCellLooper()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::topoCellLooper::cut
(
    const vector& refDir,
    const label celli,
    const boolList& vertIsCut,
    const boolList& edgeIsCut,
    const scalarField& edgeWeight,

    labelList& loop,
    scalarField& loopWeights
) const
{
    if (mesh().cellShapes()[celli].model() == hex_)
    {
        // Let parent handle hex case.
        return
            hexCellLooper::cut
            (
                refDir,
                celli,
                vertIsCut,
                edgeIsCut,
                edgeWeight,
                loop,
                loopWeights
            );
    }
    else
    {
        cellFeatures superCell(mesh(), featureCos, celli);

        if (hexMatcher().isA(superCell.faces()))
        {
            label edgeI =
                getAlignedNonFeatureEdge
                (
                    refDir,
                    celli,
                    superCell
                );

            label vertI = -1;

            label facei = -1;

            if (edgeI != -1)
            {
                // Found non-feature edge. Start walking from vertex on edge.
                vertI = mesh().edges()[edgeI].start();
            }
            else
            {
                // No 'matching' non-feature edge found on cell. Get starting
                // normal i.e. feature edge.
                edgeI = getMisAlignedEdge(refDir, celli);

                // Get any face using edge
                label face0;
                label face1;
                meshTools::getEdgeFaces(mesh(), celli, edgeI, face0, face1);

                facei = face0;
            }

            label nEstCuts = 2*mesh().cells()[celli].size();

            DynamicList<label> localLoop(nEstCuts);
            DynamicList<scalar> localLoopWeights(nEstCuts);

            walkSplitHex
            (
                celli,
                superCell,
                facei,
                edgeI,
                vertI,

                localLoop,
                localLoopWeights
            );

            if (localLoop.size() <=2)
            {
                return false;
            }
            else
            {
                loop.transfer(localLoop);
                loopWeights.transfer(localLoopWeights);

                return true;
            }
        }
        else
        {
            // Let parent handle poly case.
            return hexCellLooper::cut
            (
                refDir,
                celli,
                vertIsCut,
                edgeIsCut,
                edgeWeight,
                loop,
                loopWeights
            );
        }
    }
}


bool Foam::topoCellLooper::cut
(
    const plane& cutPlane,
    const label celli,
    const boolList& vertIsCut,
    const boolList& edgeIsCut,
    const scalarField& edgeWeight,

    labelList& loop,
    scalarField& loopWeights
) const
{
    // Let parent handle cut with plane.
    return
        hexCellLooper::cut
        (
            cutPlane,
            celli,
            vertIsCut,
            edgeIsCut,
            edgeWeight,

            loop,
            loopWeights
        );
}


// ************************************************************************* //
