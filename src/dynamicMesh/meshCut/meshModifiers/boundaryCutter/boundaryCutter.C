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

#include "boundaryCutter.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "polyAddCell.H"
#include "polyAddFace.H"
#include "polyAddPoint.H"
#include "polyModifyFace.H"
#include "polyModifyPoint.H"
#include "mapPolyMesh.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(boundaryCutter, 0);
}


// * * * * * * * * * * * * * Private Static Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::boundaryCutter::getFaceInfo
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


// Adds additional vertices (from edge cutting) to face. Used for faces which
// are not split but still might use edge that has been cut.
Foam::face Foam::boundaryCutter::addEdgeCutsToFace
(
    const label faceI,
    const Map<labelList>& edgeToAddedPoints
) const
{
    const edgeList& edges = mesh_.edges();
    const face& f = mesh_.faces()[faceI];
    const labelList& fEdges = mesh_.faceEdges()[faceI];

    // Storage for face
    DynamicList<label> newFace(2 * f.size());

    forAll(f, fp)
    {
        // Duplicate face vertex .
        newFace.append(f[fp]);

        // Check if edge has been cut.
        label v1 = f.nextLabel(fp);

        label edgeI = meshTools::findEdge(edges, fEdges, f[fp], v1);

        Map<labelList>::const_iterator fnd = edgeToAddedPoints.find(edgeI);

        if (fnd != edgeToAddedPoints.end())
        {
            // edge has been cut. Introduce new vertices. Check order.
            const labelList& addedPoints = fnd();

            if (edges[edgeI].start() == f[fp])
            {
                // Introduce in same order.
                forAll(addedPoints, i)
                {
                    newFace.append(addedPoints[i]);
                }
            }
            else
            {
                // Introduce in opposite order.
                forAllReverse(addedPoints, i)
                {
                    newFace.append(addedPoints[i]);
                }
            }
        }
    }

    face returnFace;
    returnFace.transfer(newFace);

    if (debug)
    {
        Pout<< "addEdgeCutsToFace:" << nl
            << "    from : " << f << nl
            << "    to   : " << returnFace << endl;
    }

    return returnFace;
}


void Foam::boundaryCutter::addFace
(
    const label faceI,
    const face& newFace,

    bool& modifiedFace,     // have we already 'used' faceI
    polyTopoChange& meshMod
) const
{
    // Information about old face
    label patchID, zoneID, zoneFlip;
    getFaceInfo(faceI, patchID, zoneID, zoneFlip);
    label own = mesh_.faceOwner()[faceI];
    label masterPoint = mesh_.faces()[faceI][0];

    if (!modifiedFace)
    {
        meshMod.setAction
        (
            polyModifyFace
            (
                newFace,       // face
                faceI,
                own,           // owner
                -1,            // neighbour
                false,         // flux flip
                patchID,       // patch for face
                false,         // remove from zone
                zoneID,        // zone for face
                zoneFlip       // face zone flip
            )
        );

        modifiedFace = true;
    }
    else
    {
        meshMod.setAction
        (
            polyAddFace
            (
                newFace,       // face
                own,           // owner
                -1,            // neighbour
                masterPoint,   // master point
                -1,            // master edge
                -1,            // master face for addition
                false,         // flux flip
                patchID,       // patch for face
                zoneID,        // zone for face
                zoneFlip       // face zone flip
            )
        );
    }
}



// Splits a face using the cut edges and modified points
bool Foam::boundaryCutter::splitFace
(
    const label faceI,
    const Map<point>& pointToPos,
    const Map<labelList>& edgeToAddedPoints,
    polyTopoChange& meshMod
) const
{
    const edgeList& edges = mesh_.edges();
    const face& f = mesh_.faces()[faceI];
    const labelList& fEdges = mesh_.faceEdges()[faceI];

    // Count number of split edges and total number of splits.
    label nSplitEdges = 0;
    label nModPoints = 0;
    label nTotalSplits = 0;

    forAll(f, fp)
    {
        if (pointToPos.found(f[fp]))
        {
            nModPoints++;
            nTotalSplits++;
        }

        // Check if edge has been cut.
        label nextV = f.nextLabel(fp);

        label edgeI = meshTools::findEdge(edges, fEdges, f[fp], nextV);

        Map<labelList>::const_iterator fnd = edgeToAddedPoints.find(edgeI);

        if (fnd != edgeToAddedPoints.end())
        {
            nSplitEdges++;
            nTotalSplits += fnd().size();
        }
    }

    if (debug)
    {
        Pout<< "Face:" << faceI
            << " nModPoints:" << nModPoints
            << " nSplitEdges:" << nSplitEdges
            << " nTotalSplits:" << nTotalSplits << endl;
    }

    if (nSplitEdges == 0 && nModPoints == 0)
    {
        FatalErrorIn("boundaryCutter::splitFace") << "Problem : face:" << faceI
            << " nSplitEdges:" << nSplitEdges
            << " nTotalSplits:" << nTotalSplits
            << abort(FatalError);
        return false;
    }
    else if (nSplitEdges + nModPoints == 1)
    {
        // single or multiple cuts on a single edge or single modified point
        // Dont cut and let caller handle this.
        Warning << "Face " << faceI << " has only one edge cut " << endl;
        return false;
    }
    else
    {
        // So guaranteed to have two edges cut or points modified. Split face:
        // - find starting cut
        // - walk to next cut. Make face
        // - loop until face done.

        // Information about old face
        label patchID, zoneID, zoneFlip;
        getFaceInfo(faceI, patchID, zoneID, zoneFlip);

        // Get face with new points on cut edges for ease of looping
        face extendedFace(addEdgeCutsToFace(faceI, edgeToAddedPoints));

        // Find first added point. This is the starting vertex for splitting.
        label startFp = -1;

        forAll(extendedFace, fp)
        {
            if (extendedFace[fp] >= mesh_.nPoints())
            {
                startFp = fp;
                break;
            }
        }

        if (startFp == -1)
        {
            // No added point. Maybe there is a modified point?
            forAll(extendedFace, fp)
            {
                if (pointToPos.found(extendedFace[fp]))
                {
                    startFp = fp;
                    break;
                }
            }
        }

        if (startFp == -1)
        {
            FatalErrorIn("boundaryCutter::splitFace")
                << "Problem" << abort(FatalError);
        }

        // Have we already modified existing face (first face gets done
        // as modification; all following ones as polyAddFace)
        bool modifiedFace = false;

        // Example face:
        //    +--+
        //   /   |
        //  /    |
        // +     +
        //  \    |
        //   \   |
        //    +--+
        //
        // Needs to get split into:
        // - three 'side' faces a,b,c
        // - one middle face d
        //    +--+
        //   /|\A|
        //  / | \|
        // + C|D +
        //  \ | /|
        //   \|/B|
        //    +--+


        // Storage for new face
        DynamicList<label> newFace(extendedFace.size());

        label fp = startFp;

        forAll(extendedFace, i)
        {
            label pointI = extendedFace[fp];

            newFace.append(pointI);

            if
            (
                newFace.size() > 2
             && (
                    pointI >= mesh_.nPoints()
                 || pointToPos.found(pointI)
                )
            )
            {
                // Enough vertices to create a face from.
                face tmpFace;
                tmpFace.transfer(newFace);

                // Add face tmpFace
                addFace(faceI, tmpFace, modifiedFace, meshMod);

                // Starting point is also the starting point for the new face
                newFace.append(extendedFace[startFp]);
                newFace.append(extendedFace[fp]);
            }

            fp = (fp+1) % extendedFace.size();
        }

        // Check final face.
        if (newFace.size() > 2)
        {
            // Enough vertices to create a face from.
            face tmpFace;
            tmpFace.transfer(newFace);

            // Add face tmpFace
            addFace(faceI, tmpFace, modifiedFace, meshMod);
        }

        // Split something
        return true;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::boundaryCutter::boundaryCutter(const polyMesh& mesh)
:
    mesh_(mesh),
    edgeAddedPoints_(),
    faceAddedPoint_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::boundaryCutter::~boundaryCutter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::boundaryCutter::setRefinement
(
    const Map<point>& pointToPos,
    const Map<List<point> >& edgeToCuts,
    const Map<labelPair>& faceToSplit,
    const Map<point>& faceToFeaturePoint,
    polyTopoChange& meshMod
)
{
    // Clear and size maps here since mesh size will change.
    edgeAddedPoints_.clear();

    faceAddedPoint_.clear();
    faceAddedPoint_.resize(faceToFeaturePoint.size());


    //
    // Points that just need to be moved
    // Note: could just as well be handled outside of setRefinement.
    //

    forAllConstIter(Map<point>, pointToPos, iter)
    {
        meshMod.setAction
        (
            polyModifyPoint
            (
                iter.key(), // point
                iter(),     // position
                false,      // no zone
                -1,         // zone for point
                true        // supports a cell
            )
        );
    }


    //
    // Add new points along cut edges.
    //

    // Map from edge label to sorted list of points
    Map<labelList> edgeToAddedPoints(edgeToCuts.size());

    forAllConstIter(Map<List<point> >, edgeToCuts, iter)
    {
        label edgeI = iter.key();

        const edge& e = mesh_.edges()[edgeI];

        // Sorted (from start to end) list of cuts on edge
        const List<point>& cuts = iter();

        forAll(cuts, cutI)
        {
            // point on feature to move to
            const point& featurePoint = cuts[cutI];

            label addedPointI =
                meshMod.setAction
                (
                    polyAddPoint
                    (
                        featurePoint,               // point
                        e.start(),                  // master point
                        -1,                         // zone for point
                        true                        // supports a cell
                    )
                );

            Map<labelList>::iterator fnd = edgeToAddedPoints.find(edgeI);

            if (fnd != edgeToAddedPoints.end())
            {
                labelList& addedPoints = fnd();

                label sz = addedPoints.size();
                addedPoints.setSize(sz+1);
                addedPoints[sz] = addedPointI;
            }
            else
            {
                edgeToAddedPoints.insert(edgeI, labelList(1, addedPointI));
            }

            if (debug)
            {
                Pout<< "Added point " << addedPointI << " for edge " << edgeI
                    << " with cuts:" << edgeToAddedPoints[edgeI] << endl;
            }
        }
    }


    //
    // Introduce feature points.
    //

    forAllConstIter(Map<point>, faceToFeaturePoint, iter)
    {
        label faceI = iter.key();

        const face& f = mesh_.faces()[faceI];

        if (faceToSplit.found(faceI))
        {
            FatalErrorIn("boundaryCutter::setRefinement")
                << "Face " << faceI << " vertices " << f
                << " is both marked for face-centre decomposition and"
                << " diagonal splitting."
                << abort(FatalError);
        }

        if (mesh_.isInternalFace(faceI))
        {
            FatalErrorIn("boundaryCutter::setRefinement")
                << "Face " << faceI << " vertices " << f
                << " is not an external face. Cannot split it"
                << abort(FatalError);
        }

        label addedPointI =
            meshMod.setAction
            (
                polyAddPoint
                (
                    iter(), // point
                    f[0],   // master point
                    -1,     // zone for point
                    true    // supports a cell
                )
            );
        faceAddedPoint_.insert(faceI, addedPointI);

        if (debug)
        {
            Pout<< "Added point " << addedPointI << " for feature point "
                << iter() << " on face " << faceI << " with centre "
                << mesh_.faceCentres()[faceI] << endl;
        }
    }


    //
    // Split or retriangulate faces
    //


    // Maintain whether face has been updated (for -split edges
    // -new owner/neighbour)
    boolList faceUptodate(mesh_.nFaces(), false);


    // Triangulate faces containing feature points
    forAllConstIter(Map<label>, faceAddedPoint_, iter)
    {
        label faceI = iter.key();

        // Get face with new points on cut edges.
        face newFace(addEdgeCutsToFace(faceI, edgeToAddedPoints));

        label addedPointI = iter();

        // Information about old face
        label patchID, zoneID, zoneFlip;
        getFaceInfo(faceI, patchID, zoneID, zoneFlip);
        label own = mesh_.faceOwner()[faceI];
        label masterPoint = mesh_.faces()[faceI][0];

        // Triangulate face around mid point

        face tri(3);

        forAll(newFace, fp)
        {
            label nextV = newFace.nextLabel(fp);

            tri[0] = newFace[fp];
            tri[1] = nextV;
            tri[2] = addedPointI;

            if (fp == 0)
            {
                // Modify the existing face.
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        tri,                        // face
                        faceI,
                        own,                        // owner
                        -1,                         // neighbour
                        false,                      // flux flip
                        patchID,                    // patch for face
                        false,                      // remove from zone
                        zoneID,                     // zone for face
                        zoneFlip                    // face zone flip
                    )
                );
            }
            else
            {
                // Add additional faces
                meshMod.setAction
                (
                    polyAddFace
                    (
                        tri,                        // face
                        own,                        // owner
                        -1,                         // neighbour
                        masterPoint,                // master point
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

        faceUptodate[faceI] = true;
    }


    // Diagonally split faces
    forAllConstIter(Map<labelPair>, faceToSplit, iter)
    {
        label faceI = iter.key();

        const face& f = mesh_.faces()[faceI];

        if (faceAddedPoint_.found(faceI))
        {
            FatalErrorIn("boundaryCutter::setRefinement")
                << "Face " << faceI << " vertices " << f
                << " is both marked for face-centre decomposition and"
                << " diagonal splitting."
                << abort(FatalError);
        }


        // Get face with new points on cut edges.
        face newFace(addEdgeCutsToFace(faceI, edgeToAddedPoints));

        // Information about old face
        label patchID, zoneID, zoneFlip;
        getFaceInfo(faceI, patchID, zoneID, zoneFlip);
        label own = mesh_.faceOwner()[faceI];
        label masterPoint = mesh_.faces()[faceI][0];

        // Split face from one side of diagonal to other.
        const labelPair& diag = iter();

        label fp0 = findIndex(newFace, f[diag[0]]);
        label fp1 = findIndex(newFace, f[diag[1]]);

        if (fp0 == -1 || fp1 == -1 || fp0 == fp1)
        {
            FatalErrorIn("boundaryCutter::setRefinement")
                << "Problem : Face " << faceI << " vertices " << f
                << " newFace:" << newFace << " diagonal:" << f[diag[0]]
                << ' ' << f[diag[1]]
                << abort(FatalError);
        }

        // Replace existing face by newFace from fp0 to fp1 and add new one
        // from fp1 to fp0.

        DynamicList<label> newVerts(newFace.size());

        // Get vertices from fp0 to (and including) fp1
        label fp = fp0;

        do
        {
            newVerts.append(newFace[fp]);

            fp = (fp == newFace.size()-1 ? 0 : fp+1);

        } while (fp != fp1);

        newVerts.append(newFace[fp1]);


        // Modify the existing face.
        meshMod.setAction
        (
            polyModifyFace
            (
                face(newVerts.shrink()),    // face
                faceI,
                own,                        // owner
                -1,                         // neighbour
                false,                      // flux flip
                patchID,                    // patch for face
                false,                      // remove from zone
                zoneID,                     // zone for face
                zoneFlip                    // face zone flip
            )
        );


        newVerts.clear();

        // Get vertices from fp1 to (and including) fp0

        do
        {
            newVerts.append(newFace[fp]);

            fp = (fp == newFace.size()-1 ? 0 : fp+1);

        } while (fp != fp0);

        newVerts.append(newFace[fp0]);

        // Add additional face
        meshMod.setAction
        (
            polyAddFace
            (
                face(newVerts.shrink()),    // face
                own,                        // owner
                -1,                         // neighbour
                masterPoint,                // master point
                -1,                         // master edge
                -1,                         // master face for addition
                false,                      // flux flip
                patchID,                    // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face zone flip
            )
        );

        faceUptodate[faceI] = true;
    }


    // Split external faces without feature point but using cut edges.
    // Does right handed walk but not really.
    forAllConstIter(Map<labelList>, edgeToAddedPoints, iter)
    {
        label edgeI = iter.key();

        const labelList& eFaces = mesh_.edgeFaces()[edgeI];

        forAll(eFaces, i)
        {
            label faceI = eFaces[i];

            if (!faceUptodate[faceI] && !mesh_.isInternalFace(faceI))
            {
                // Is external face so split
                if (splitFace(faceI, pointToPos, edgeToAddedPoints, meshMod))
                {
                    // Successfull split
                    faceUptodate[faceI] = true;
                }
            }
        }
    }


    // Add cut edges (but don't split) any other faces using any cut edge.
    // These can be external faces where splitFace hasn't cut them or
    // internal faces.
    forAllConstIter(Map<labelList>, edgeToAddedPoints, iter)
    {
        label edgeI = iter.key();

        const labelList& eFaces = mesh_.edgeFaces()[edgeI];

        forAll(eFaces, i)
        {
            label faceI = eFaces[i];

            if (!faceUptodate[faceI])
            {
                // Renumber face to include split edges.
                face newFace(addEdgeCutsToFace(faceI, edgeToAddedPoints));

                label own = mesh_.faceOwner()[faceI];

                label nei = -1;

                if (mesh_.isInternalFace(faceI))
                {
                    nei = mesh_.faceNeighbour()[faceI];
                }

                label patchID, zoneID, zoneFlip;
                getFaceInfo(faceI, patchID, zoneID, zoneFlip);

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

                faceUptodate[faceI] = true;
            }
        }
    }

    // Convert edge to points storage from edge labels (not preserved)
    // to point labels
    edgeAddedPoints_.resize(edgeToCuts.size());

    forAllConstIter(Map<labelList>, edgeToAddedPoints, iter)
    {
        edgeAddedPoints_.insert(mesh_.edges()[iter.key()], iter());
    }
}


void Foam::boundaryCutter::updateMesh(const mapPolyMesh& morphMap)
{
    // Update stored labels for mesh change.

    //
    // Do faceToAddedPoint
    //

    {
        // Create copy since we're deleting entries.
        Map<label> newAddedPoints(faceAddedPoint_.size());

        forAllConstIter(Map<label>, faceAddedPoint_, iter)
        {
            label oldFaceI = iter.key();

            label newFaceI = morphMap.reverseFaceMap()[oldFaceI];

            label oldPointI = iter();

            label newPointI = morphMap.reversePointMap()[oldPointI];

            if (newFaceI >= 0 && newPointI >= 0)
            {
                newAddedPoints.insert(newFaceI, newPointI);
            }
        }

        // Copy
        faceAddedPoint_.transfer(newAddedPoints);
    }


    //
    // Do edgeToAddedPoints
    //


    {
        // Create copy since we're deleting entries
        HashTable<labelList, edge, Hash<edge> >
            newEdgeAddedPoints(edgeAddedPoints_.size());

        for
        (
            HashTable<labelList, edge, Hash<edge> >::const_iterator iter =
                edgeAddedPoints_.begin();
            iter != edgeAddedPoints_.end();
            ++iter
        )
        {
            const edge& e = iter.key();

            label newStart = morphMap.reversePointMap()[e.start()];

            label newEnd = morphMap.reversePointMap()[e.end()];

            if (newStart >= 0 && newEnd >= 0)
            {
                const labelList& addedPoints = iter();

                labelList newAddedPoints(addedPoints.size());
                label newI = 0;

                forAll(addedPoints, i)
                {
                    label newAddedPointI =
                        morphMap.reversePointMap()[addedPoints[i]];

                    if (newAddedPointI >= 0)
                    {
                        newAddedPoints[newI++] = newAddedPointI;
                    }
                }
                if (newI > 0)
                {
                    newAddedPoints.setSize(newI);

                    edge newE = edge(newStart, newEnd);

                    newEdgeAddedPoints.insert(newE, newAddedPoints);
                }
            }
        }

        // Copy
        edgeAddedPoints_.transfer(newEdgeAddedPoints);
    }
}


// ************************************************************************* //
