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

#include "meshDualiser.H"
#include "meshTools.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "edgeFaceCirculator.H"
#include "mergePoints.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(meshDualiser, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::meshDualiser::checkPolyTopoChange(const polyTopoChange& meshMod)
{
    // Assume no removed points
    pointField points(meshMod.points().size());
    forAll(meshMod.points(), i)
    {
        points[i] = meshMod.points()[i];
    }

    labelList oldToNew;
    label nUnique = mergePoints
    (
        points,
        1e-6,
        false,
        oldToNew
    );

    if (nUnique < points.size())
    {
        labelListList newToOld(invertOneToMany(nUnique, oldToNew));

        forAll(newToOld, newI)
        {
            if (newToOld[newI].size() != 1)
            {
                FatalErrorIn
                (
                    "meshDualiser::checkPolyTopoChange(const polyTopoChange&)"
                )   << "duplicate verts:" << newToOld[newI]
                    << " coords:"
                    << UIndirectList<point>(points, newToOld[newI])()
                    << abort(FatalError);
            }
        }
    }
}


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

    forAll(points, pointI)
    {
        meshTools::writeOBJ(str1, points[pointI]);
        meshTools::writeOBJ(str2, points[pointI]);
    }

    const DynamicList<face>& faces = meshMod.faces();

    forAll(faces, faceI)
    {
        const face& f = faces[faceI];

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


//- Given cell and point on mesh finds the corresponding dualCell. Most
//  points become only one cell but the feature points might be split.
Foam::label Foam::meshDualiser::findDualCell
(
    const label cellI,
    const label pointI
) const
{
    const labelList& dualCells = pointToDualCells_[pointI];

    if (dualCells.size() == 1)
    {
        return dualCells[0];
    }
    else
    {
        label index = findIndex(mesh_.pointCells()[pointI], cellI);

        return dualCells[index];
    }
}


// Helper function to generate dualpoints on all boundary edges emanating
// from (boundary & feature) point
void Foam::meshDualiser::generateDualBoundaryEdges
(
    const PackedBoolList& isBoundaryEdge,
    const label pointI,
    polyTopoChange& meshMod
)
{
    const labelList& pEdges = mesh_.pointEdges()[pointI];

    forAll(pEdges, pEdgeI)
    {
        label edgeI = pEdges[pEdgeI];

        if (edgeToDualPoint_[edgeI] == -1 && isBoundaryEdge.get(edgeI) == 1)
        {
            const edge& e = mesh_.edges()[edgeI];

            edgeToDualPoint_[edgeI] = meshMod.addPoint
            (
                e.centre(mesh_.points()),
                pointI, // masterPoint
                -1,     // zoneID
                true    // inCell
            );
        }
    }
}


// Return true if point on face has same dual cells on both owner and neighbour
// sides.
bool Foam::meshDualiser::sameDualCell
(
    const label faceI,
    const label pointI
) const
{
    if (!mesh_.isInternalFace(faceI))
    {
        FatalErrorIn("sameDualCell(const label, const label)")
            << "face:" << faceI << " is not internal face."
            << abort(FatalError);
    }

    label own = mesh_.faceOwner()[faceI];
    label nei = mesh_.faceNeighbour()[faceI];

    return findDualCell(own, pointI) == findDualCell(nei, pointI);
}


Foam::label Foam::meshDualiser::addInternalFace
(
    const label masterPointI,
    const label masterEdgeI,
    const label masterFaceI,

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

    if (debug)
    {
        pointField facePoints(meshMod.points(), newFace);

        labelList oldToNew;
        label nUnique = mergePoints
        (
            facePoints,
            1e-6,
            false,
            oldToNew
        );

        if (nUnique < facePoints.size())
        {
            FatalErrorIn("addInternalFace(..)")
                << "verts:" << verts << " newFace:" << newFace
                << " face points:" << facePoints
                << abort(FatalError);
        }
    }


    label zoneID = -1;
    bool zoneFlip = false;
    if (masterFaceI != -1)
    {
        zoneID = mesh_.faceZones().whichZone(masterFaceI);

        if (zoneID != -1)
        {
            const faceZone& fZone = mesh_.faceZones()[zoneID];

            zoneFlip = fZone.flipMap()[fZone.whichFace(masterFaceI)];
        }
    }

    label dualFaceI;

    if (dualCell0 < dualCell1)
    {
        dualFaceI = meshMod.addFace
        (
            newFace,
            dualCell0,      // own
            dualCell1,      // nei
            masterPointI,   // masterPointID
            masterEdgeI,    // masterEdgeID
            masterFaceI,    // masterFaceID
            false,          // flipFaceFlux
            -1,             // patchID
            zoneID,         // zoneID
            zoneFlip        // zoneFlip
        );

        //pointField dualPoints(meshMod.points());
        //vector n(newFace.normal(dualPoints));
        //n /= mag(n);
        //Pout<< "Generated internal dualFace:" << dualFaceI
        //    << " verts:" << newFace
        //    << " points:" << UIndirectList<point>(meshMod.points(), newFace)()
        //    << " n:" << n
        //    << " between dualowner:" << dualCell0
        //    << " dualneigbour:" << dualCell1
        //    << endl;
    }
    else
    {
        dualFaceI = meshMod.addFace
        (
            newFace,
            dualCell1,      // own
            dualCell0,      // nei
            masterPointI,   // masterPointID
            masterEdgeI,    // masterEdgeID
            masterFaceI,    // masterFaceID
            false,          // flipFaceFlux
            -1,             // patchID
            zoneID,         // zoneID
            zoneFlip        // zoneFlip
        );

        //pointField dualPoints(meshMod.points());
        //vector n(newFace.normal(dualPoints));
        //n /= mag(n);
        //Pout<< "Generated internal dualFace:" << dualFaceI
        //    << " verts:" << newFace
        //    << " points:" << UIndirectList<point>(meshMod.points(), newFace)()
        //    << " n:" << n
        //    << " between dualowner:" << dualCell1
        //    << " dualneigbour:" << dualCell0
        //    << endl;
    }
    return dualFaceI;
}


Foam::label Foam::meshDualiser::addBoundaryFace
(
    const label masterPointI,
    const label masterEdgeI,
    const label masterFaceI,

    const label dualCellI,
    const label patchI,
    const DynamicList<label>& verts,
    polyTopoChange& meshMod
) const
{
    face newFace(verts);

    label zoneID = -1;
    bool zoneFlip = false;
    if (masterFaceI != -1)
    {
        zoneID = mesh_.faceZones().whichZone(masterFaceI);

        if (zoneID != -1)
        {
            const faceZone& fZone = mesh_.faceZones()[zoneID];

            zoneFlip = fZone.flipMap()[fZone.whichFace(masterFaceI)];
        }
    }

    label dualFaceI = meshMod.addFace
    (
        newFace,
        dualCellI,      // own
        -1,             // nei
        masterPointI,   // masterPointID
        masterEdgeI,    // masterEdgeID
        masterFaceI,    // masterFaceID
        false,          // flipFaceFlux
        patchI,         // patchID
        zoneID,         // zoneID
        zoneFlip        // zoneFlip
    );

    //pointField dualPoints(meshMod.points());
    //vector n(newFace.normal(dualPoints));
    //n /= mag(n);
    //Pout<< "Generated boundary dualFace:" << dualFaceI
    //    << " verts:" << newFace
    //    << " points:" << UIndirectList<point>(meshMod.points(), newFace)()
    //    << " n:" << n
    //    << " on dualowner:" << dualCellI
    //    << endl;
    return dualFaceI;
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
    const label startFaceI,
    polyTopoChange& meshMod,
    boolList& doneEFaces
) const
{
    const edge& e = mesh_.edges()[edgeI];
    const labelList& eFaces = mesh_.edgeFaces()[edgeI];

    label fp = edgeFaceCirculator::getMinIndex
    (
        mesh_.faces()[startFaceI],
        e[0],
        e[1]
    );

    edgeFaceCirculator ie
    (
        mesh_,
        startFaceI, // face
        true,       // ownerSide
        fp,         // fp
        isBoundaryEdge.get(edgeI) == 1  // isBoundaryEdge
    );
    ie.setCanonical();

    bool edgeOrder = ie.sameOrder(e[0], e[1]);
    label startFaceLabel = ie.faceLabel();

    //Pout<< "At edge:" << edgeI << " verts:" << e
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
        label faceI = ie.faceLabel();

        // Mark face as visited.
        doneEFaces[findIndex(eFaces, faceI)] = true;

        if (faceToDualPoint_[faceI] != -1)
        {
            verts.append(faceToDualPoint_[faceI]);
        }

        label cellI = ie.cellLabel();

        if (cellI == -1)
        {
            // At ending boundary face. We've stored the face point above
            // so this is the whole face.
            break;
        }


        label dualCell0 = findDualCell(cellI, e[0]);
        label dualCell1 = findDualCell(cellI, e[1]);

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
                -1,         // masterPointI
                edgeI,      // masterEdgeI
                -1,         // masterFaceI
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
            if (faceToDualPoint_[faceI] != -1)
            {
                verts.append(faceToDualPoint_[faceI]);
            }
        }

        if (cellToDualPoint_[cellI] != -1)
        {
            verts.append(cellToDualPoint_[cellI]);
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
        -1,         // masterPointI
        edgeI,      // masterEdgeI
        -1,         // masterFaceI
        edgeOrder,
        currentDualCell0,
        currentDualCell1,
        verts,
        meshMod
    );
}


// Walks around circumference of faceI. Creates single face. Gets given
// starting (feature) edge to start from. Returns ending edge. (all edges
// in form of index in faceEdges)
void Foam::meshDualiser::createFaceFromInternalFace
(
    const label faceI,
    label& fp,
    polyTopoChange& meshMod
) const
{
    const face& f = mesh_.faces()[faceI];
    const labelList& fEdges = mesh_.faceEdges()[faceI];
    label own = mesh_.faceOwner()[faceI];
    label nei = mesh_.faceNeighbour()[faceI];

    //Pout<< "createFaceFromInternalFace : At face:" << faceI
    //    << " verts:" << f
    //    << " points:" << UIndirectList<point>(mesh_.points(), f)()
    //    << " started walking at edge:" << fEdges[fp]
    //    << " verts:" << mesh_.edges()[fEdges[fp]]
    //    << endl;


    // Walk and collect face.
    DynamicList<label> verts(100);

    verts.append(faceToDualPoint_[faceI]);
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
                FatalErrorIn("createFacesFromInternalFace(..)")
                    << "face:" << faceI << " verts:" << f
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
                -1,         // masterPointI
                -1,         // masterEdgeI
                faceI,      // masterFaceI
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
    const label patchI,
    const label patchPointI,
    const label startFaceI,
    polyTopoChange& meshMod,
    boolList& donePFaces            // pFaces visited
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const polyPatch& pp = patches[patchI];
    const labelList& pFaces = pp.pointFaces()[patchPointI];
    const labelList& own = mesh_.faceOwner();

    label pointI = pp.meshPoints()[patchPointI];

    if (pointToDualPoint_[pointI] == -1)
    {
        // Not a feature point. Loop over all connected
        // pointFaces.

        // Starting face
        label faceI = startFaceI;

        DynamicList<label> verts(4);

        while (true)
        {
            label index = findIndex(pFaces, faceI-pp.start());

            // Has face been visited already?
            if (donePFaces[index])
            {
                break;
            }
            donePFaces[index] = true;

            // Insert face centre
            verts.append(faceToDualPoint_[faceI]);

            label dualCellI = findDualCell(own[faceI], pointI);

            // Get the edge before the patchPointI
            const face& f = mesh_.faces()[faceI];
            label fp = findIndex(f, pointI);
            label prevFp = f.rcIndex(fp);
            label edgeI = mesh_.faceEdges()[faceI][prevFp];

            if (edgeToDualPoint_[edgeI] != -1)
            {
                verts.append(edgeToDualPoint_[edgeI]);
            }

            // Get next boundary face (whilst staying on edge).
            edgeFaceCirculator circ
            (
                mesh_,
                faceI,
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
            faceI = circ.faceLabel();

            if (faceI < pp.start() || faceI >= pp.start()+pp.size())
            {
                FatalErrorIn
                (
                    "meshDualiser::createFacesAroundBoundaryPoint(..)"
                )   << "Walked from face on patch:" << patchI
                    << " to face:" << faceI
                    << " fc:" << mesh_.faceCentres()[faceI]
                    << " on patch:" << patches.whichPatch(faceI)
                    << abort(FatalError);
            }

            // Check if different cell.
            if (dualCellI != findDualCell(own[faceI], pointI))
            {
                FatalErrorIn
                (
                    "meshDualiser::createFacesAroundBoundaryPoint(..)"
                )   << "Different dual cells but no feature edge"
                    << " inbetween point:" << pointI
                    << " coord:" << mesh_.points()[pointI]
                    << abort(FatalError);
            }
        }

        verts.shrink();

        label dualCellI = findDualCell(own[faceI], pointI);

        //Bit dodgy: create dualface from the last face (instead of from
        // the central point). This will also use the original faceZone to
        // put the new face (which might span multiple original faces) in.

        addBoundaryFace
        (
            //pointI,     // masterPointI
            -1,         // masterPointI
            -1,         // masterEdgeI
            faceI,      // masterFaceI
            dualCellI,
            patchI,
            verts,
            meshMod
        );
    }
    else
    {
        label faceI = startFaceI;

        // Storage for face
        DynamicList<label> verts(mesh_.faces()[faceI].size());

        // Starting point.
        verts.append(pointToDualPoint_[pointI]);

        // Find edge between pointI and next point on face.
        const labelList& fEdges = mesh_.faceEdges()[faceI];
        label nextEdgeI = fEdges[findIndex(mesh_.faces()[faceI], pointI)];
        if (edgeToDualPoint_[nextEdgeI] != -1)
        {
            verts.append(edgeToDualPoint_[nextEdgeI]);
        }

        do
        {
            label index = findIndex(pFaces, faceI-pp.start());

            // Has face been visited already?
            if (donePFaces[index])
            {
                break;
            }
            donePFaces[index] = true;

            // Face centre
            verts.append(faceToDualPoint_[faceI]);

            // Find edge before pointI on faceI
            const labelList& fEdges = mesh_.faceEdges()[faceI];
            const face& f = mesh_.faces()[faceI];
            label prevFp = f.rcIndex(findIndex(f, pointI));
            label edgeI = fEdges[prevFp];

            if (edgeToDualPoint_[edgeI] != -1)
            {
                // Feature edge. Close any face so far. Note: uses face to
                // create dualFace from. Could use pointI instead.
                verts.append(edgeToDualPoint_[edgeI]);
                addBoundaryFace
                (
                    -1,     // masterPointI
                    -1,     // masterEdgeI
                    faceI,  // masterFaceI
                    findDualCell(own[faceI], pointI),
                    patchI,
                    verts.shrink(),
                    meshMod
                );
                verts.clear();

                verts.append(pointToDualPoint_[pointI]);
                verts.append(edgeToDualPoint_[edgeI]);
            }

            // Cross edgeI to next boundary face
            edgeFaceCirculator circ
            (
                mesh_,
                faceI,
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
            faceI = circ.faceLabel();
        }
        while
        (
            faceI != startFaceI
         && faceI >= pp.start()
         && faceI < pp.start()+pp.size()
        );

        if (verts.size() > 2)
        {
            // Note: face created from face, not from pointI
            addBoundaryFace
            (
                -1,             // masterPointI
                -1,             // masterEdgeI
                startFaceI,     // masterFaceI
                findDualCell(own[faceI], pointI),
                patchI,
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
    for (label faceI = mesh_.nInternalFaces(); faceI < mesh_.nFaces(); faceI++)
    {
        const labelList& fEdges = mesh_.faceEdges()[faceI];

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
        label faceI = findIndex(featureFaceSet, false);

        if (faceI != -1)
        {
            FatalErrorIn
            (
                "meshDualiser::setRefinement"
                "(const labelList&, const labelList&, const labelList&"
                ", const labelList, polyTopoChange&)"
            )   << "In split-face-mode (splitFace=true) but not all faces"
                << " marked as feature faces." << endl
                << "First conflicting face:" << faceI
                << " centre:" << mesh_.faceCentres()[faceI]
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
            FatalErrorIn
            (
                "meshDualiser::setRefinement"
                "(const labelList&, const labelList&, const labelList&"
                ", const labelList, polyTopoChange&)"
            )   << "In split-face-mode (splitFace=true) but not all edges"
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
            label faceI = mesh_.nInternalFaces();
            faceI < mesh_.nFaces();
            faceI++
        )
        {
            if (!featureFaceSet[faceI])
            {
                FatalErrorIn
                (
                    "meshDualiser::setRefinement"
                    "(const labelList&, const labelList&, const labelList&"
                    ", const labelList, polyTopoChange&)"
                )   << "Not all boundary faces marked as feature faces."
                    << endl
                    << "First conflicting face:" << faceI
                    << " centre:" << mesh_.faceCentres()[faceI]
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

    // pointToDualCells_[pointI]
    // - single entry : all cells surrounding point all become the same
    //                  cell.
    // - multiple entries: in order of pointCells.


    // feature points that become single cell
    forAll(singleCellFeaturePoints, i)
    {
        label pointI = singleCellFeaturePoints[i];

        pointToDualPoint_[pointI] = meshMod.addPoint
        (
            mesh_.points()[pointI],
            pointI,                                 // masterPoint
            mesh_.pointZones().whichZone(pointI),   // zoneID
            true                                    // inCell
        );

        // Generate single cell
        pointToDualCells_[pointI].setSize(1);
        pointToDualCells_[pointI][0] = meshMod.addCell
        (
            pointI, //masterPointID,
            -1,     //masterEdgeID,
            -1,     //masterFaceID,
            -1,     //masterCellID,
            -1      //zoneID
        );
        if (dualCcStr.valid())
        {
            meshTools::writeOBJ(dualCcStr(), mesh_.points()[pointI]);
        }
    }

    // feature points that become multiple cells
    forAll(multiCellFeaturePoints, i)
    {
        label pointI = multiCellFeaturePoints[i];

        if (pointToDualCells_[pointI].size() > 0)
        {
            FatalErrorIn
            (
                "meshDualiser::setRefinement"
                "(const labelList&, const labelList&, const labelList&"
                ", const labelList, polyTopoChange&)"
            )   << "Point " << pointI << " at:" << mesh_.points()[pointI]
                << " is both in singleCellFeaturePoints"
                << " and multiCellFeaturePoints."
                << abort(FatalError);
        }

        pointToDualPoint_[pointI] = meshMod.addPoint
        (
            mesh_.points()[pointI],
            pointI,                                 // masterPoint
            mesh_.pointZones().whichZone(pointI),   // zoneID
            true                                    // inCell
        );

        // Create dualcell for every cell connected to dual point

        const labelList& pCells = mesh_.pointCells()[pointI];

        pointToDualCells_[pointI].setSize(pCells.size());

        forAll(pCells, pCellI)
        {
            pointToDualCells_[pointI][pCellI] = meshMod.addCell
            (
                pointI,                                     //masterPointID
                -1,                                         //masterEdgeID
                -1,                                         //masterFaceID
                -1,                                         //masterCellID
                mesh_.cellZones().whichZone(pCells[pCellI]) //zoneID
            );
            if (dualCcStr.valid())
            {
                meshTools::writeOBJ
                (
                    dualCcStr(),
                    0.5*(mesh_.points()[pointI]+cellCentres[pCells[pCellI]])
                );
            }
        }
    }
    // Normal points
    forAll(mesh_.points(), pointI)
    {
        if (pointToDualCells_[pointI].empty())
        {
            pointToDualCells_[pointI].setSize(1);
            pointToDualCells_[pointI][0] = meshMod.addCell
            (
                pointI, //masterPointID,
                -1,     //masterEdgeID,
                -1,     //masterFaceID,
                -1,     //masterCellID,
                -1      //zoneID
            );

            if (dualCcStr.valid())
            {
                meshTools::writeOBJ(dualCcStr(), mesh_.points()[pointI]);
            }
        }
    }


    // Dual points (from cell centres, feature faces, feature edges)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(cellToDualPoint_, cellI)
    {
        cellToDualPoint_[cellI] = meshMod.addPoint
        (
            cellCentres[cellI],
            mesh_.faces()[mesh_.cells()[cellI][0]][0],     // masterPoint
            -1,     // zoneID
            true    // inCell
        );
    }

    // From face to dual point

    forAll(featureFaces, i)
    {
        label faceI = featureFaces[i];

        faceToDualPoint_[faceI] = meshMod.addPoint
        (
            mesh_.faceCentres()[faceI],
            mesh_.faces()[faceI][0],     // masterPoint
            -1,     // zoneID
            true    // inCell
        );
    }
    // Detect whether different dual cells on either side of a face. This
    // would neccesitate having a dual face built from the face and thus a
    // dual point at the face centre.
    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        if (faceToDualPoint_[faceI] == -1)
        {
            const face& f = mesh_.faces()[faceI];

            forAll(f, fp)
            {
                label ownDualCell = findDualCell(own[faceI], f[fp]);
                label neiDualCell = findDualCell(nei[faceI], f[fp]);

                if (ownDualCell != neiDualCell)
                {
                    faceToDualPoint_[faceI] = meshMod.addPoint
                    (
                        mesh_.faceCentres()[faceI],
                        f[fp],  // masterPoint
                        -1,     // zoneID
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
            -1,     // zoneID
            true    // inCell
        );
    }

    // Detect whether different dual cells on either side of an edge. This
    // would neccesitate having a dual face built perpendicular to the edge
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
                        mesh_.pointZones().whichZone(e[0]), // zoneID
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
                        mesh_.pointZones().whichZone(e[1]), // zoneID
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
        checkPolyTopoChange(meshMod);
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

                label startFaceI = eFaces[i];

                //Pout<< "Walking edge:" << edgeI
                //    << " points:" << mesh_.points()[e[0]]
                //    << mesh_.points()[e[1]]
                //    << " startFace:" << startFaceI
                //    << " at:" << mesh_.faceCentres()[startFaceI]
                //    << endl;

                createFacesAroundEdge
                (
                    splitFace,
                    isBoundaryEdge,
                    edgeI,
                    startFaceI,
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
    // - non-feature face : inbetween cell zones.
    forAll(faceToDualPoint_, faceI)
    {
        if (faceToDualPoint_[faceI] != -1 && mesh_.isInternalFace(faceI))
        {
            const face& f = mesh_.faces()[faceI];
            const labelList& fEdges = mesh_.faceEdges()[faceI];

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
                    && !sameDualCell(faceI, f.nextLabel(fp))
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
                    faceI,
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

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        const labelListList& pointFaces = pp.pointFaces();

        forAll(pointFaces, patchPointI)
        {
            const labelList& pFaces = pointFaces[patchPointI];

            boolList donePFaces(pFaces.size(), false);

            forAll(pFaces, i)
            {
                if (!donePFaces[i])
                {
                    // Starting face
                    label startFaceI = pp.start()+pFaces[i];

                    //Pout<< "Walking around point:" << pointI
                    //    << " coord:" << mesh_.points()[pointI]
                    //    << " on patch:" << patchI
                    //    << " startFace:" << startFaceI
                    //    << " at:" << mesh_.faceCentres()[startFaceI]
                    //    << endl;

                    createFacesAroundBoundaryPoint
                    (
                        patchI,
                        patchPointI,
                        startFaceI,
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
