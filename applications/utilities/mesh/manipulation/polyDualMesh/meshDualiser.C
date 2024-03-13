/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "meshDualiser.H"
#include "meshTools.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "polyTopoChangeMap.H"
#include "edgeFaceCirculator.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(meshDualiser, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Dump state so far.
void Foam::meshDualiser::dumpPolyTopoChange
(
    const polyTopoChange& meshMod,
    const fileName& prefix
)
{
    OFstream str1(prefix + "Faces.obj");
    OFstream str2(prefix + "Edges.obj");

    Info<< "Dumping current polyTopoChange. Faces to " << str1.name()
        << " , points and edges to " << str2.name() << endl;

    const DynamicList<point>& points = meshMod.points();

    forAll(points, pointi)
    {
        meshTools::writeOBJ(str1, points[pointi]);
        meshTools::writeOBJ(str2, points[pointi]);
    }

    const DynamicList<face>& faces = meshMod.faces();

    forAll(faces, facei)
    {
        const face& f = faces[facei];

        str1<< 'f';
        forAll(f, fp)
        {
            str1<< ' ' << f[fp]+1;
        }
        str1<< nl;

        str2<< 'l';
        forAll(f, fp)
        {
            str2<< ' ' << f[fp]+1;
        }
        str2<< ' ' << f[0]+1 << nl;
    }
}


Foam::label Foam::meshDualiser::findDualCell
(
    const label celli,
    const label pointi
) const
{
    const labelList& dualCells = pointToDualCells_[pointi];

    if (dualCells.size() == 1)
    {
        return dualCells[0];
    }
    else
    {
        label index = findIndex(mesh_.pointCells()[pointi], celli);

        return dualCells[index];
    }
}


void Foam::meshDualiser::generateDualBoundaryEdges
(
    const PackedBoolList& isBoundaryEdge,
    const label pointi,
    polyTopoChange& meshMod
)
{
    const labelList& pEdges = mesh_.pointEdges()[pointi];

    forAll(pEdges, pEdgeI)
    {
        label edgeI = pEdges[pEdgeI];

        if (edgeToDualPoint_[edgeI] == -1 && isBoundaryEdge.get(edgeI) == 1)
        {
            const edge& e = mesh_.edges()[edgeI];

            edgeToDualPoint_[edgeI] = meshMod.addPoint
            (
                e.centre(mesh_.points()),
                pointi, // masterPoint
                true    // inCell
            );
        }
    }
}


// Return true if point on face has same dual cells on both owner and neighbour
// sides.
bool Foam::meshDualiser::sameDualCell
(
    const label facei,
    const label pointi
) const
{
    if (!mesh_.isInternalFace(facei))
    {
        FatalErrorInFunction
            << "face:" << facei << " is not internal face."
            << abort(FatalError);
    }

    label own = mesh_.faceOwner()[facei];
    label nei = mesh_.faceNeighbour()[facei];

    return findDualCell(own, pointi) == findDualCell(nei, pointi);
}


Foam::label Foam::meshDualiser::addInternalFace
(
    const label masterPointi,
    const label masterEdgeI,
    const label masterFacei,

    const bool edgeOrder,
    const label dualCell0,
    const label dualCell1,
    const DynamicList<label>& verts,
    polyTopoChange& meshMod
) const
{
    face newFace(verts);

    if (edgeOrder != (dualCell0 < dualCell1))
    {
        reverse(newFace);
    }

    label dualFacei;

    if (dualCell0 < dualCell1)
    {
        dualFacei = meshMod.addFace
        (
            newFace,
            dualCell0,      // own
            dualCell1,      // nei
            masterFacei,    // masterFaceID
            false,          // flipFaceFlux
            -1,             // patchID
            -1,             // zoneID
            false           // zoneFlip
        );

        // pointField dualPoints(meshMod.points());
        // vector n(newFace.normal(dualPoints));
        // n /= mag(n);
        // Pout<< "Generated internal dualFace:" << dualFacei
        //    << " verts:" << newFace
        //    << " points:" << UIndirectList<point>(meshMod.points(), newFace)()
        //    << " n:" << n
        //    << " between dualowner:" << dualCell0
        //    << " dualneighbour:" << dualCell1
        //    << endl;
    }
    else
    {
        dualFacei = meshMod.addFace
        (
            newFace,
            dualCell1,      // own
            dualCell0,      // nei
            masterFacei,    // masterFaceID
            false,          // flipFaceFlux
            -1,             // patchID
            -1,             // zoneID
            false           // zoneFlip
        );

        // pointField dualPoints(meshMod.points());
        // vector n(newFace.normal(dualPoints));
        // n /= mag(n);
        // Pout<< "Generated internal dualFace:" << dualFacei
        //    << " verts:" << newFace
        //    << " points:" << UIndirectList<point>(meshMod.points(), newFace)()
        //    << " n:" << n
        //    << " between dualowner:" << dualCell1
        //    << " dualneighbour:" << dualCell0
        //    << endl;
    }
    return dualFacei;
}


Foam::label Foam::meshDualiser::addBoundaryFace
(
    const label masterPointi,
    const label masterEdgeI,
    const label masterFacei,

    const label dualCelli,
    const label patchi,
    const DynamicList<label>& verts,
    polyTopoChange& meshMod
) const
{
    face newFace(verts);

    label dualFacei = meshMod.addFace
    (
        newFace,
        dualCelli,      // own
        -1,             // nei
        masterFacei,    // masterFaceID
        false,          // flipFaceFlux
        patchi,         // patchID
        -1,             // zoneID
        false           // zoneFlip
    );

    // pointField dualPoints(meshMod.points());
    // vector n(newFace.normal(dualPoints));
    // n /= mag(n);
    // Pout<< "Generated boundary dualFace:" << dualFacei
    //    << " verts:" << newFace
    //    << " points:" << UIndirectList<point>(meshMod.points(), newFace)()
    //    << " n:" << n
    //    << " on dualowner:" << dualCelli
    //    << endl;
    return dualFacei;
}


// Walks around edgeI.
// splitFace=true : creates multiple faces
// splitFace=false: creates single face if same dual cells on both sides,
//                  multiple faces otherwise.
void Foam::meshDualiser::createFacesAroundEdge
(
    const bool splitFace,
    const PackedBoolList& isBoundaryEdge,
    const label edgeI,
    const label startFacei,
    polyTopoChange& meshMod,
    boolList& doneEFaces
) const
{
    const edge& e = mesh_.edges()[edgeI];
    const labelList& eFaces = mesh_.edgeFaces()[edgeI];

    label fp = edgeFaceCirculator::getMinIndex
    (
        mesh_.faces()[startFacei],
        e[0],
        e[1]
    );

    edgeFaceCirculator ie
    (
        mesh_,
        startFacei, // face
        true,       // ownerSide
        fp,         // fp
        isBoundaryEdge.get(edgeI) == 1  // isBoundaryEdge
    );
    ie.setCanonical();

    bool edgeOrder = ie.sameOrder(e[0], e[1]);
    label startFaceLabel = ie.faceLabel();

    // Pout<< "At edge:" << edgeI << " verts:" << e
    //    << " points:" << mesh_.points()[e[0]] << mesh_.points()[e[1]]
    //    << " started walking at face:" << ie.faceLabel()
    //    << " verts:" << mesh_.faces()[ie.faceLabel()]
    //    << " edgeOrder:" << edgeOrder
    //    << " in direction of cell:" << ie.cellLabel()
    //    << endl;

    // Walk and collect face.
    DynamicList<label> verts(100);

    if (edgeToDualPoint_[edgeI] != -1)
    {
        verts.append(edgeToDualPoint_[edgeI]);
    }
    if (faceToDualPoint_[ie.faceLabel()] != -1)
    {
        doneEFaces[findIndex(eFaces, ie.faceLabel())] = true;
        verts.append(faceToDualPoint_[ie.faceLabel()]);
    }
    if (cellToDualPoint_[ie.cellLabel()] != -1)
    {
        verts.append(cellToDualPoint_[ie.cellLabel()]);
    }

    label currentDualCell0 = findDualCell(ie.cellLabel(), e[0]);
    label currentDualCell1 = findDualCell(ie.cellLabel(), e[1]);

    ++ie;

    while (true)
    {
        label facei = ie.faceLabel();

        // Mark face as visited.
        doneEFaces[findIndex(eFaces, facei)] = true;

        if (faceToDualPoint_[facei] != -1)
        {
            verts.append(faceToDualPoint_[facei]);
        }

        label celli = ie.cellLabel();

        if (celli == -1)
        {
            // At ending boundary face. We've stored the face point above
            // so this is the whole face.
            break;
        }


        label dualCell0 = findDualCell(celli, e[0]);
        label dualCell1 = findDualCell(celli, e[1]);

        // Generate face. (always if splitFace=true; only if needed to
        // separate cells otherwise)
        if
        (
            splitFace
         || (
                dualCell0 != currentDualCell0
             || dualCell1 != currentDualCell1
            )
        )
        {
            // Close current face.
            addInternalFace
            (
                -1,         // masterPointi
                edgeI,      // masterEdgeI
                -1,         // masterFacei
                edgeOrder,
                currentDualCell0,
                currentDualCell1,
                verts.shrink(),
                meshMod
            );

            // Restart
            currentDualCell0 = dualCell0;
            currentDualCell1 = dualCell1;

            verts.clear();
            if (edgeToDualPoint_[edgeI] != -1)
            {
                verts.append(edgeToDualPoint_[edgeI]);
            }
            if (faceToDualPoint_[facei] != -1)
            {
                verts.append(faceToDualPoint_[facei]);
            }
        }

        if (cellToDualPoint_[celli] != -1)
        {
            verts.append(cellToDualPoint_[celli]);
        }

        ++ie;

        if (ie == ie.end())
        {
            // Back at start face (for internal edge only). See if this needs
            // adding.
            if (isBoundaryEdge.get(edgeI) == 0)
            {
                label startDual = faceToDualPoint_[startFaceLabel];

                if (startDual != -1 && findIndex(verts, startDual) == -1)
                {
                    verts.append(startDual);
                }
            }
            break;
        }
    }

    verts.shrink();
    addInternalFace
    (
        -1,         // masterPointi
        edgeI,      // masterEdgeI
        -1,         // masterFacei
        edgeOrder,
        currentDualCell0,
        currentDualCell1,
        verts,
        meshMod
    );
}


// Walks around circumference of facei. Creates single face. Gets given
// starting (feature) edge to start from. Returns ending edge. (all edges
// in form of index in faceEdges)
void Foam::meshDualiser::createFaceFromInternalFace
(
    const label facei,
    label& fp,
    polyTopoChange& meshMod
) const
{
    const face& f = mesh_.faces()[facei];
    const labelList& fEdges = mesh_.faceEdges()[facei];
    label own = mesh_.faceOwner()[facei];
    label nei = mesh_.faceNeighbour()[facei];

    // Pout<< "createFaceFromInternalFace : At face:" << facei
    //    << " verts:" << f
    //    << " points:" << UIndirectList<point>(mesh_.points(), f)()
    //    << " started walking at edge:" << fEdges[fp]
    //    << " verts:" << mesh_.edges()[fEdges[fp]]
    //    << endl;


    // Walk and collect face.
    DynamicList<label> verts(100);

    verts.append(faceToDualPoint_[facei]);
    verts.append(edgeToDualPoint_[fEdges[fp]]);

    // Step to vertex after edge mid
    fp = f.fcIndex(fp);

    // Get cells on either side of face at that point
    label currentDualCell0 = findDualCell(own, f[fp]);
    label currentDualCell1 = findDualCell(nei, f[fp]);

    forAll(f, i)
    {
        // Check vertex
        if (pointToDualPoint_[f[fp]] != -1)
        {
            verts.append(pointToDualPoint_[f[fp]]);
        }

        // Edge between fp and fp+1
        label edgeI = fEdges[fp];

        if (edgeToDualPoint_[edgeI] != -1)
        {
            verts.append(edgeToDualPoint_[edgeI]);
        }

        // Next vertex on edge
        label nextFp = f.fcIndex(fp);

        // Get dual cells on nextFp to check whether face needs closing.
        label dualCell0 = findDualCell(own, f[nextFp]);
        label dualCell1 = findDualCell(nei, f[nextFp]);

        if (dualCell0 != currentDualCell0 || dualCell1 != currentDualCell1)
        {
            // Check: make sure that there is a midpoint on the edge.
            if (edgeToDualPoint_[edgeI] == -1)
            {
                FatalErrorInFunction
                    << "face:" << facei << " verts:" << f
                    << " points:" << UIndirectList<point>(mesh_.points(), f)()
                    << " no feature edge between " << f[fp]
                    << " and " << f[nextFp] << " although have different"
                    << " dual cells." << endl
                    << "point " << f[fp] << " has dual cells "
                    << currentDualCell0 << " and " << currentDualCell1
                    << " ; point "<< f[nextFp] << " has dual cells "
                    << dualCell0 << " and " << dualCell1
                    << abort(FatalError);
            }


            // Close current face.
            verts.shrink();
            addInternalFace
            (
                -1,         // masterPointi
                -1,         // masterEdgeI
                facei,      // masterFacei
                true,       // edgeOrder,
                currentDualCell0,
                currentDualCell1,
                verts,
                meshMod
            );
            break;
        }

        fp = nextFp;
    }
}


// Given a point on a face converts the faces around the point.
// (pointFaces()). Gets starting face and marks off visited faces in donePFaces.
void Foam::meshDualiser::createFacesAroundBoundaryPoint
(
    const label patchi,
    const label patchPointi,
    const label startFacei,
    polyTopoChange& meshMod,
    boolList& donePFaces            // pFaces visited
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const polyPatch& pp = patches[patchi];
    const labelList& pFaces = pp.pointFaces()[patchPointi];
    const labelList& own = mesh_.faceOwner();

    label pointi = pp.meshPoints()[patchPointi];

    if (pointToDualPoint_[pointi] == -1)
    {
        // Not a feature point. Loop over all connected
        // pointFaces.

        // Starting face
        label facei = startFacei;

        DynamicList<label> verts(4);

        while (true)
        {
            label index = findIndex(pFaces, facei-pp.start());

            // Has face been visited already?
            if (donePFaces[index])
            {
                break;
            }
            donePFaces[index] = true;

            // Insert face centre
            verts.append(faceToDualPoint_[facei]);

            label dualCelli = findDualCell(own[facei], pointi);

            // Get the edge before the patchPointi
            const face& f = mesh_.faces()[facei];
            label fp = findIndex(f, pointi);
            label prevFp = f.rcIndex(fp);
            label edgeI = mesh_.faceEdges()[facei][prevFp];

            if (edgeToDualPoint_[edgeI] != -1)
            {
                verts.append(edgeToDualPoint_[edgeI]);
            }

            // Get next boundary face (whilst staying on edge).
            edgeFaceCirculator circ
            (
                mesh_,
                facei,
                true,   // ownerSide
                prevFp, // index of edge in face
                true    // isBoundaryEdge
            );

            do
            {
                ++circ;
            }
            while (mesh_.isInternalFace(circ.faceLabel()));

            // Step to next face
            facei = circ.faceLabel();

            if (facei < pp.start() || facei >= pp.start()+pp.size())
            {
                FatalErrorInFunction
                    << "Walked from face on patch:" << patchi
                    << " to face:" << facei
                    << " fc:" << mesh_.faceCentres()[facei]
                    << " on patch:" << patches.whichPatch(facei)
                    << abort(FatalError);
            }

            // Check if different cell.
            if (dualCelli != findDualCell(own[facei], pointi))
            {
                FatalErrorInFunction
                    << "Different dual cells but no feature edge"
                    << " in between point:" << pointi
                    << " coord:" << mesh_.points()[pointi]
                    << abort(FatalError);
            }
        }

        verts.shrink();

        label dualCelli = findDualCell(own[facei], pointi);

        // Bit dodgy: create dualface from the last face (instead of from
        // the central point). This will also use the original faceZone to
        // put the new face (which might span multiple original faces) in.

        addBoundaryFace
        (
            // pointi,     // masterPointi
            -1,         // masterPointi
            -1,         // masterEdgeI
            facei,      // masterFacei
            dualCelli,
            patchi,
            verts,
            meshMod
        );
    }
    else
    {
        label facei = startFacei;

        // Storage for face
        DynamicList<label> verts(mesh_.faces()[facei].size());

        // Starting point.
        verts.append(pointToDualPoint_[pointi]);

        // Find edge between pointi and next point on face.
        const labelList& fEdges = mesh_.faceEdges()[facei];
        label nextEdgeI = fEdges[findIndex(mesh_.faces()[facei], pointi)];
        if (edgeToDualPoint_[nextEdgeI] != -1)
        {
            verts.append(edgeToDualPoint_[nextEdgeI]);
        }

        do
        {
            label index = findIndex(pFaces, facei-pp.start());

            // Has face been visited already?
            if (donePFaces[index])
            {
                break;
            }
            donePFaces[index] = true;

            // Face centre
            verts.append(faceToDualPoint_[facei]);

            // Find edge before pointi on facei
            const labelList& fEdges = mesh_.faceEdges()[facei];
            const face& f = mesh_.faces()[facei];
            label prevFp = f.rcIndex(findIndex(f, pointi));
            label edgeI = fEdges[prevFp];

            if (edgeToDualPoint_[edgeI] != -1)
            {
                // Feature edge. Close any face so far. Note: uses face to
                // create dualFace from. Could use pointi instead.
                verts.append(edgeToDualPoint_[edgeI]);
                addBoundaryFace
                (
                    -1,     // masterPointi
                    -1,     // masterEdgeI
                    facei,  // masterFacei
                    findDualCell(own[facei], pointi),
                    patchi,
                    verts.shrink(),
                    meshMod
                );
                verts.clear();

                verts.append(pointToDualPoint_[pointi]);
                verts.append(edgeToDualPoint_[edgeI]);
            }

            // Cross edgeI to next boundary face
            edgeFaceCirculator circ
            (
                mesh_,
                facei,
                true,   // ownerSide
                prevFp, // index of edge in face
                true    // isBoundaryEdge
            );

            do
            {
                ++circ;
            }
            while (mesh_.isInternalFace(circ.faceLabel()));

            // Step to next face. Quit if not on same patch.
            facei = circ.faceLabel();
        }
        while
        (
            facei != startFacei
         && facei >= pp.start()
         && facei < pp.start()+pp.size()
        );

        if (verts.size() > 2)
        {
            // Note: face created from face, not from pointi
            addBoundaryFace
            (
                -1,             // masterPointi
                -1,             // masterEdgeI
                startFacei,     // masterFacei
                findDualCell(own[facei], pointi),
                patchi,
                verts.shrink(),
                meshMod
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshDualiser::meshDualiser(const polyMesh& mesh)
:
    mesh_(mesh),
    pointToDualCells_(mesh_.nPoints()),
    pointToDualPoint_(mesh_.nPoints(), -1),
    cellToDualPoint_(mesh_.nCells()),
    faceToDualPoint_(mesh_.nFaces(), -1),
    edgeToDualPoint_(mesh_.nEdges(), -1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meshDualiser::setRefinement
(
    const bool splitFace,
    const labelList& featureFaces,
    const labelList& featureEdges,
    const labelList& singleCellFeaturePoints,
    const labelList& multiCellFeaturePoints,
    polyTopoChange& meshMod
)
{
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();
    const vectorField& cellCentres = mesh_.cellCentres();

    // Mark boundary edges and points.
    // (Note: in 1.4.2 we can use the built-in mesh point ordering
    //  facility instead)
    PackedBoolList isBoundaryEdge(mesh_.nEdges());
    for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); facei++)
    {
        const labelList& fEdges = mesh_.faceEdges()[facei];

        forAll(fEdges, i)
        {
            isBoundaryEdge.set(fEdges[i], 1);
        }
    }


    if (splitFace)
    {
        // This is a special mode where whenever we are walking around an edge
        // every area through a cell becomes a separate dualface. So two
        // dual cells will probably have more than one dualface between them!
        // This mode implies that
        // - all faces have to be feature faces since there has to be a
        //   dualpoint at the face centre.
        // - all edges have to be feature edges ,,
        boolList featureFaceSet(mesh_.nFaces(), false);
        forAll(featureFaces, i)
        {
            featureFaceSet[featureFaces[i]] = true;
        }
        label facei = findIndex(featureFaceSet, false);

        if (facei != -1)
        {
            FatalErrorInFunction
                << "In split-face-mode (splitFace=true) but not all faces"
                << " marked as feature faces." << endl
                << "First conflicting face:" << facei
                << " centre:" << mesh_.faceCentres()[facei]
                << abort(FatalError);
        }

        boolList featureEdgeSet(mesh_.nEdges(), false);
        forAll(featureEdges, i)
        {
            featureEdgeSet[featureEdges[i]] = true;
        }
        label edgeI = findIndex(featureEdgeSet, false);

        if (edgeI != -1)
        {
            const edge& e = mesh_.edges()[edgeI];
            FatalErrorInFunction
                << "In split-face-mode (splitFace=true) but not all edges"
                << " marked as feature edges." << endl
                << "First conflicting edge:" << edgeI
                << " verts:" << e
                << " coords:" << mesh_.points()[e[0]] << mesh_.points()[e[1]]
                << abort(FatalError);
        }
    }
    else
    {
        // Check that all boundary faces are feature faces.

        boolList featureFaceSet(mesh_.nFaces(), false);
        forAll(featureFaces, i)
        {
            featureFaceSet[featureFaces[i]] = true;
        }
        for
        (
            label facei = mesh_.nInternalFaces();
            facei < mesh_.nFaces();
            facei++
        )
        {
            if (!featureFaceSet[facei])
            {
                FatalErrorInFunction
                    << "Not all boundary faces marked as feature faces."
                    << endl
                    << "First conflicting face:" << facei
                    << " centre:" << mesh_.faceCentres()[facei]
                    << abort(FatalError);
            }
        }
    }




    // Start creating cells, points, and faces (in that order)


    // 1. Mark which cells to create
    // Mostly every point becomes one cell but sometimes (for feature points)
    // all cells surrounding a feature point become cells. Also a non-manifold
    // point can create two cells! So a dual cell is uniquely defined by a
    // mesh point + cell (as in pointCells index)

    // 2. Mark which face centres to create

    // 3. Internal faces can now consist of
    //      - only cell centres of walk around edge
    //      - cell centres + face centres of walk around edge
    //      - same but now other side is not a single cell

    // 4. Boundary faces (or internal faces between cell zones!) now consist of
    //      - walk around boundary point.



    autoPtr<OFstream> dualCcStr;
    if (debug)
    {
        dualCcStr.reset(new OFstream("dualCc.obj"));
        Pout<< "Dumping centres of dual cells to " << dualCcStr().name()
            << endl;
    }


    // Dual cells (from points)
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    // pointToDualCells_[pointi]
    // - single entry : all cells surrounding point all become the same
    //                  cell.
    // - multiple entries: in order of pointCells.


    // feature points that become single cell
    forAll(singleCellFeaturePoints, i)
    {
        label pointi = singleCellFeaturePoints[i];

        pointToDualPoint_[pointi] = meshMod.addPoint
        (
            mesh_.points()[pointi],
            pointi,                                 // masterPoint
            true                                    // inCell
        );

        // Generate single cell
        pointToDualCells_[pointi].setSize(1);
        pointToDualCells_[pointi][0] = meshMod.addCell
        (
            -1,     // masterCellID,
            -1      // zoneID
        );
        if (dualCcStr.valid())
        {
            meshTools::writeOBJ(dualCcStr(), mesh_.points()[pointi]);
        }
    }

    // feature points that become multiple cells
    forAll(multiCellFeaturePoints, i)
    {
        label pointi = multiCellFeaturePoints[i];

        if (pointToDualCells_[pointi].size() > 0)
        {
            FatalErrorInFunction
                << "Point " << pointi << " at:" << mesh_.points()[pointi]
                << " is both in singleCellFeaturePoints"
                << " and multiCellFeaturePoints."
                << abort(FatalError);
        }

        pointToDualPoint_[pointi] = meshMod.addPoint
        (
            mesh_.points()[pointi],
            pointi,                                 // masterPoint
            true                                    // inCell
        );

        // Create dualcell for every cell connected to dual point

        const labelList& pCells = mesh_.pointCells()[pointi];

        pointToDualCells_[pointi].setSize(pCells.size());

        forAll(pCells, pCelli)
        {
            pointToDualCells_[pointi][pCelli] = meshMod.addCell
            (
                -1,                                         // masterCellID
                -1
            );
            if (dualCcStr.valid())
            {
                meshTools::writeOBJ
                (
                    dualCcStr(),
                    0.5*(mesh_.points()[pointi]+cellCentres[pCells[pCelli]])
                );
            }
        }
    }
    // Normal points
    forAll(mesh_.points(), pointi)
    {
        if (pointToDualCells_[pointi].empty())
        {
            pointToDualCells_[pointi].setSize(1);
            pointToDualCells_[pointi][0] = meshMod.addCell
            (
                -1,     // masterCellID,
                -1      // zoneID
            );

            if (dualCcStr.valid())
            {
                meshTools::writeOBJ(dualCcStr(), mesh_.points()[pointi]);
            }
        }
    }


    // Dual points (from cell centres, feature faces, feature edges)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(cellToDualPoint_, celli)
    {
        cellToDualPoint_[celli] = meshMod.addPoint
        (
            cellCentres[celli],
            mesh_.faces()[mesh_.cells()[celli][0]][0],     // masterPoint
            true    // inCell
        );
    }

    // From face to dual point

    forAll(featureFaces, i)
    {
        label facei = featureFaces[i];

        faceToDualPoint_[facei] = meshMod.addPoint
        (
            mesh_.faceCentres()[facei],
            mesh_.faces()[facei][0],     // masterPoint
            true    // inCell
        );
    }
    // Detect whether different dual cells on either side of a face. This
    // would necessitate having a dual face built from the face and thus a
    // dual point at the face centre.
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        if (faceToDualPoint_[facei] == -1)
        {
            const face& f = mesh_.faces()[facei];

            forAll(f, fp)
            {
                label ownDualCell = findDualCell(own[facei], f[fp]);
                label neiDualCell = findDualCell(nei[facei], f[fp]);

                if (ownDualCell != neiDualCell)
                {
                    faceToDualPoint_[facei] = meshMod.addPoint
                    (
                        mesh_.faceCentres()[facei],
                        f[fp],  // masterPoint
                        true    // inCell
                    );

                    break;
                }
            }
        }
    }

    // From edge to dual point

    forAll(featureEdges, i)
    {
        label edgeI = featureEdges[i];

        const edge& e = mesh_.edges()[edgeI];

        edgeToDualPoint_[edgeI] = meshMod.addPoint
        (
            e.centre(mesh_.points()),
            e[0],   // masterPoint
            true    // inCell
        );
    }

    // Detect whether different dual cells on either side of an edge. This
    // would necessitate having a dual face built perpendicular to the edge
    // and thus a dual point at the mid of the edge.
    // Note: not really true - the face can be built without the edge centre!
    const labelListList& edgeCells = mesh_.edgeCells();

    forAll(edgeCells, edgeI)
    {
       if (edgeToDualPoint_[edgeI] == -1)
       {
            const edge& e = mesh_.edges()[edgeI];

            // We need a point on the edge if not all cells on both sides
            // are the same.

            const labelList& eCells = mesh_.edgeCells()[edgeI];

            label dualE0 = findDualCell(eCells[0], e[0]);
            label dualE1 = findDualCell(eCells[0], e[1]);

            for (label i = 1; i < eCells.size(); i++)
            {
                label newDualE0 = findDualCell(eCells[i], e[0]);

                if (dualE0 != newDualE0)
                {
                    edgeToDualPoint_[edgeI] = meshMod.addPoint
                    (
                        e.centre(mesh_.points()),
                        e[0],                               // masterPoint
                        true                                // inCell
                    );

                    break;
                }

                label newDualE1 = findDualCell(eCells[i], e[1]);

                if (dualE1 != newDualE1)
                {
                    edgeToDualPoint_[edgeI] = meshMod.addPoint
                    (
                        e.centre(mesh_.points()),
                        e[1],   // masterPoint
                        true    // inCell
                    );

                    break;
                }
            }
        }
    }

    // Make sure all boundary edges emanating from feature points are
    // feature edges as well.
    forAll(singleCellFeaturePoints, i)
    {
        generateDualBoundaryEdges
        (
            isBoundaryEdge,
            singleCellFeaturePoints[i],
            meshMod
        );
    }
    forAll(multiCellFeaturePoints, i)
    {
        generateDualBoundaryEdges
        (
            isBoundaryEdge,
            multiCellFeaturePoints[i],
            meshMod
        );
    }


    // Check for duplicate points
    if (debug)
    {
        dumpPolyTopoChange(meshMod, "generatedPoints_");
    }


    // Now we have all points and cells
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  - pointToDualCells_ : per point a single dualCell or multiple dualCells
    //  - pointToDualPoint_ : per point -1 or the dual point at the coordinate
    //  - edgeToDualPoint_  : per edge -1 or the edge centre
    //  - faceToDualPoint_  : per face -1 or the face centre
    //  - cellToDualPoint_  : per cell the cell centre
    // Now we have to walk all edges and construct faces. Either single face
    // per edge or multiple (-if nonmanifold edge -if different dualcells)

    const edgeList& edges = mesh_.edges();

    forAll(edges, edgeI)
    {
        const labelList& eFaces = mesh_.edgeFaces()[edgeI];

        boolList doneEFaces(eFaces.size(), false);

        forAll(eFaces, i)
        {
            if (!doneEFaces[i])
            {
                // We found a face that hasn't yet been visited. This might
                // happen for non-manifold edges where a single edge can
                // become multiple faces.

                label startFacei = eFaces[i];

                // Pout<< "Walking edge:" << edgeI
                //    << " points:" << mesh_.points()[e[0]]
                //    << mesh_.points()[e[1]]
                //    << " startFace:" << startFacei
                //    << " at:" << mesh_.faceCentres()[startFacei]
                //    << endl;

                createFacesAroundEdge
                (
                    splitFace,
                    isBoundaryEdge,
                    edgeI,
                    startFacei,
                    meshMod,
                    doneEFaces
                );
            }
        }
    }

    if (debug)
    {
        dumpPolyTopoChange(meshMod, "generatedFacesFromEdges_");
    }

    // Create faces from feature faces. These can be internal or external faces.
    // - feature face : centre needs to be included.
    //      - single cells on either side: triangulate
    //      - multiple cells: create single face between unique cell pair. Only
    //                        create face where cells differ on either side.
    // - non-feature face : in between cell zones.
    forAll(faceToDualPoint_, facei)
    {
        if (faceToDualPoint_[facei] != -1 && mesh_.isInternalFace(facei))
        {
            const face& f = mesh_.faces()[facei];
            const labelList& fEdges = mesh_.faceEdges()[facei];

            // Starting edge
            label fp = 0;

            do
            {
                // Find edge that is in dual mesh and where the
                // next point (fp+1) has different dual cells on either side.
                bool foundStart = false;

                do
                {
                    if
                    (
                        edgeToDualPoint_[fEdges[fp]] != -1
                    && !sameDualCell(facei, f.nextLabel(fp))
                    )
                    {
                        foundStart = true;
                        break;
                    }
                    fp = f.fcIndex(fp);
                }
                while (fp != 0);

                if (!foundStart)
                {
                    break;
                }

                // Walk from edge fp and generate a face.
                createFaceFromInternalFace
                (
                    facei,
                    fp,
                    meshMod
                );
            }
            while (fp != 0);
        }
    }

    if (debug)
    {
        dumpPolyTopoChange(meshMod, "generatedFacesFromFeatFaces_");
    }


    // Create boundary faces. Every boundary point has one or more dualcells.
    // These need to be closed.
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        const labelListList& pointFaces = pp.pointFaces();

        forAll(pointFaces, patchPointi)
        {
            const labelList& pFaces = pointFaces[patchPointi];

            boolList donePFaces(pFaces.size(), false);

            forAll(pFaces, i)
            {
                if (!donePFaces[i])
                {
                    // Starting face
                    label startFacei = pp.start()+pFaces[i];

                    // Pout<< "Walking around point:" << pointi
                    //    << " coord:" << mesh_.points()[pointi]
                    //    << " on patch:" << patchi
                    //    << " startFace:" << startFacei
                    //    << " at:" << mesh_.faceCentres()[startFacei]
                    //    << endl;

                    createFacesAroundBoundaryPoint
                    (
                        patchi,
                        patchPointi,
                        startFacei,
                        meshMod,
                        donePFaces            // pFaces visited
                    );
                }
            }
        }
    }

    if (debug)
    {
        dumpPolyTopoChange(meshMod, "generatedFacesFromBndFaces_");
    }
}


// ************************************************************************* //
