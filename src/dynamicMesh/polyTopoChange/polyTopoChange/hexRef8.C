/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "hexRef8.H"

#include "polyMesh.H"
#include "polyTopoChange.H"
#include "meshTools.H"
#include "polyAddFace.H"
#include "polyAddPoint.H"
#include "polyAddCell.H"
#include "polyModifyFace.H"
#include "syncTools.H"
#include "faceSet.H"
#include "cellSet.H"
#include "pointSet.H"
#include "OFstream.H"
#include "Time.H"
#include "FaceCellWave.H"
#include "mapDistributePolyMesh.H"
#include "refinementData.H"
#include "refinementDistanceData.H"
#include "degenerateMatcher.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hexRef8, 0);

    //- Reduction class. If x and y are not equal assign value.
    template<int value>
    class ifEqEqOp
    {
        public:
        void operator()(label& x, const label y) const
        {
            x = (x == y) ? x : value;
        }
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::hexRef8::reorder
(
    const labelList& map,
    const label len,
    const label null,
    labelList& elems
)
{
    labelList newElems(len, null);

    forAll(elems, i)
    {
        label newI = map[i];

        if (newI >= len)
        {
            FatalErrorIn("hexRef8::reorder(..)") << abort(FatalError);
        }

        if (newI >= 0)
        {
            newElems[newI] = elems[i];
        }
    }

    elems.transfer(newElems);
}


void Foam::hexRef8::getFaceInfo
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


// Adds a face on top of existing faceI.
Foam::label Foam::hexRef8::addFace
(
    polyTopoChange& meshMod,
    const label faceI,
    const face& newFace,
    const label own,
    const label nei
) const
{
    label patchID, zoneID, zoneFlip;

    getFaceInfo(faceI, patchID, zoneID, zoneFlip);

    label newFaceI = -1;

    if ((nei == -1) || (own < nei))
    {
        // Ordering ok.
        newFaceI = meshMod.setAction
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
        newFaceI = meshMod.setAction
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
    return newFaceI;
}


// Adds an internal face from an edge. Assumes orientation correct.
// Problem is that the face is between four new vertices. So what do we provide
// as master? The only existing mesh item we have is the edge we have split.
// Have to be careful in only using it if it has internal faces since otherwise
// polyMeshMorph will complain (because it cannot generate a sensible mapping
// for the face)
Foam::label Foam::hexRef8::addInternalFace
(
    polyTopoChange& meshMod,
    const label meshFaceI,
    const label meshPointI,
    const face& newFace,
    const label own,
    const label nei
) const
{
    if (mesh_.isInternalFace(meshFaceI))
    {
        return meshMod.setAction
        (
            polyAddFace
            (
                newFace,                    // face
                own,                        // owner
                nei,                        // neighbour
                -1,                         // master point
                -1,                         // master edge
                meshFaceI,                  // master face for addition
                false,                      // flux flip
                -1,                         // patch for face
                -1,                         // zone for face
                false                       // face zone flip
            )
        );
    }
    else
    {
        // Two choices:
        // - append (i.e. create out of nothing - will not be mapped)
        //   problem: field does not get mapped.
        // - inflate from point.
        //   problem: does interpolative mapping which constructs full
        //   volPointInterpolation!

        // For now create out of nothing

        return meshMod.setAction
        (
            polyAddFace
            (
                newFace,                    // face
                own,                        // owner
                nei,                        // neighbour
                -1,                         // master point
                -1,                         // master edge
                -1,                         // master face for addition
                false,                      // flux flip
                -1,                         // patch for face
                -1,                         // zone for face
                false                       // face zone flip
            )
        );


        ////- Inflate-from-point:
        //// Check if point has any internal faces we can use.
        //label masterPointI = -1;
        //
        //const labelList& pFaces = mesh_.pointFaces()[meshPointI];
        //
        //forAll(pFaces, i)
        //{
        //    if (mesh_.isInternalFace(pFaces[i]))
        //    {
        //        // meshPoint uses internal faces so ok to inflate from it
        //        masterPointI = meshPointI;
        //
        //        break;
        //    }
        //}
        //
        //return meshMod.setAction
        //(
        //    polyAddFace
        //    (
        //        newFace,                    // face
        //        own,                        // owner
        //        nei,                        // neighbour
        //        masterPointI,               // master point
        //        -1,                         // master edge
        //        -1,                         // master face for addition
        //        false,                      // flux flip
        //        -1,                         // patch for face
        //        -1,                         // zone for face
        //        false                       // face zone flip
        //    )
        //);
    }
}


// Modifies existing faceI for either new owner/neighbour or new face points.
void Foam::hexRef8::modFace
(
    polyTopoChange& meshMod,
    const label faceI,
    const face& newFace,
    const label own,
    const label nei
) const
{
    label patchID, zoneID, zoneFlip;

    getFaceInfo(faceI, patchID, zoneID, zoneFlip);

    if
    (
        (own != mesh_.faceOwner()[faceI])
     || (
            mesh_.isInternalFace(faceI)
         && (nei != mesh_.faceNeighbour()[faceI])
        )
     || (newFace != mesh_.faces()[faceI])
    )
    {
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


// Bit complex way to determine the unrefined edge length.
Foam::scalar Foam::hexRef8::getLevel0EdgeLength() const
{
    if (cellLevel_.size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "hexRef8::getLevel0EdgeLength() const"
        )   << "Number of cells in mesh:" << mesh_.nCells()
            << " does not equal size of cellLevel:" << cellLevel_.size()
            << endl
            << "This might be because of a restart with inconsistent cellLevel."
            << abort(FatalError);
    }

    // Determine minimum edge length per refinement level
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const scalar GREAT2 = sqr(GREAT);

    label nLevels = gMax(cellLevel_)+1;

    scalarField typEdgeLenSqr(nLevels, GREAT2);


    // 1. Look only at edges surrounded by cellLevel cells only.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        // Per edge the cellLevel of connected cells. -1 if not set,
        // labelMax if different levels, otherwise levels of connected cells.
        labelList edgeLevel(mesh_.nEdges(), -1);

        forAll(cellLevel_, cellI)
        {
            const label cLevel = cellLevel_[cellI];

            const labelList& cEdges = mesh_.cellEdges(cellI);

            forAll(cEdges, i)
            {
                label edgeI = cEdges[i];

                if (edgeLevel[edgeI] == -1)
                {
                    edgeLevel[edgeI] = cLevel;
                }
                else if (edgeLevel[edgeI] == labelMax)
                {
                    // Already marked as on different cellLevels
                }
                else if (edgeLevel[edgeI] != cLevel)
                {
                    edgeLevel[edgeI] = labelMax;
                }
            }
        }

        // Make sure that edges with different levels on different processors
        // are also marked. Do the same test (edgeLevel != cLevel) on coupled
        // edges.
        syncTools::syncEdgeList
        (
            mesh_,
            edgeLevel,
            ifEqEqOp<labelMax>(),
            labelMin
        );

        // Now use the edgeLevel with a valid value to determine the
        // length per level.
        forAll(edgeLevel, edgeI)
        {
            const label eLevel = edgeLevel[edgeI];

            if (eLevel >= 0 && eLevel < labelMax)
            {
                const edge& e = mesh_.edges()[edgeI];

                scalar edgeLenSqr = magSqr(e.vec(mesh_.points()));

                typEdgeLenSqr[eLevel] = min(typEdgeLenSqr[eLevel], edgeLenSqr);
            }
        }
    }

    // Get the minimum per level over all processors. Note minimum so if
    // cells are not cubic we use the smallest edge side.
    Pstream::listCombineGather(typEdgeLenSqr, minEqOp<scalar>());
    Pstream::listCombineScatter(typEdgeLenSqr);

    if (debug)
    {
        Pout<< "hexRef8::getLevel0EdgeLength() :"
            << " After phase1: Edgelengths (squared) per refinementlevel:"
            << typEdgeLenSqr << endl;
    }


    // 2. For any levels where we haven't determined a valid length yet
    //    use any surrounding cell level. Here we use the max so we don't
    //    pick up levels between celllevel and higher celllevel (will have
    //    edges sized according to highest celllevel)
    //    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    scalarField maxEdgeLenSqr(nLevels, -GREAT2);

    forAll(cellLevel_, cellI)
    {
        const label cLevel = cellLevel_[cellI];

        const labelList& cEdges = mesh_.cellEdges(cellI);

        forAll(cEdges, i)
        {
            const edge& e = mesh_.edges()[cEdges[i]];

            scalar edgeLenSqr = magSqr(e.vec(mesh_.points()));

            maxEdgeLenSqr[cLevel] = max(maxEdgeLenSqr[cLevel], edgeLenSqr);
        }
    }

    Pstream::listCombineGather(maxEdgeLenSqr, maxEqOp<scalar>());
    Pstream::listCombineScatter(maxEdgeLenSqr);

    if (debug)
    {
        Pout<< "hexRef8::getLevel0EdgeLength() :"
            << " Crappy Edgelengths (squared) per refinementlevel:"
            << maxEdgeLenSqr << endl;
    }


    // 3. Combine the two sets of lengths
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(typEdgeLenSqr, levelI)
    {
        if (typEdgeLenSqr[levelI] == GREAT2 && maxEdgeLenSqr[levelI] >= 0)
        {
            typEdgeLenSqr[levelI] = maxEdgeLenSqr[levelI];
        }
    }

    if (debug)
    {
        Pout<< "hexRef8::getLevel0EdgeLength() :"
            << " Final Edgelengths (squared) per refinementlevel:"
            << typEdgeLenSqr << endl;
    }

    // Find lowest level present
    scalar level0Size = -1;

    forAll(typEdgeLenSqr, levelI)
    {
        scalar lenSqr = typEdgeLenSqr[levelI];

        if (lenSqr < GREAT2)
        {
            level0Size = Foam::sqrt(lenSqr)*(1<<levelI);

            if (debug)
            {
                Pout<< "hexRef8::getLevel0EdgeLength() :"
                    << " For level:" << levelI
                    << " have edgeLen:" << Foam::sqrt(lenSqr)
                    << " with equivalent level0 len:" << level0Size
                    << endl;
            }
            break;
        }
    }

    if (level0Size == -1)
    {
        FatalErrorIn("hexRef8::getLevel0EdgeLength()")
            << "Problem : typEdgeLenSqr:" << typEdgeLenSqr << abort(FatalError);
    }

    return level0Size;
}


// Check whether pointI is an anchor on cellI.
// If it is not check whether any other point on the face is an anchor cell.
Foam::label Foam::hexRef8::getAnchorCell
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const label cellI,
    const label faceI,
    const label pointI
) const
{
    if (cellAnchorPoints[cellI].size())
    {
        label index = findIndex(cellAnchorPoints[cellI], pointI);

        if (index != -1)
        {
            return cellAddedCells[cellI][index];
        }


        // pointI is not an anchor cell.
        // Maybe we are already a refined face so check all the face
        // vertices.
        const face& f = mesh_.faces()[faceI];

        forAll(f, fp)
        {
            label index = findIndex(cellAnchorPoints[cellI], f[fp]);

            if (index != -1)
            {
                return cellAddedCells[cellI][index];
            }
        }

        // Problem.
        dumpCell(cellI);
        Perr<< "cell:" << cellI << " anchorPoints:" << cellAnchorPoints[cellI]
            << endl;

        FatalErrorIn("hexRef8::getAnchorCell(..)")
            << "Could not find point " << pointI
            << " in the anchorPoints for cell " << cellI << endl
            << "Does your original mesh obey the 2:1 constraint and"
            << " did you use consistentRefinement to make your cells to refine"
            << " obey this constraint as well?"
            << abort(FatalError);

        return -1;
    }
    else
    {
        return cellI;
    }
}


// Get new owner and neighbour
void Foam::hexRef8::getFaceNeighbours
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const label faceI,
    const label pointI,

    label& own,
    label& nei
) const
{
    // Is owner split?
    own = getAnchorCell
    (
        cellAnchorPoints,
        cellAddedCells,
        mesh_.faceOwner()[faceI],
        faceI,
        pointI
    );

    if (mesh_.isInternalFace(faceI))
    {
        nei = getAnchorCell
        (
            cellAnchorPoints,
            cellAddedCells,
            mesh_.faceNeighbour()[faceI],
            faceI,
            pointI
        );
    }
    else
    {
        nei = -1;
    }
}


// Get point with the lowest pointLevel
Foam::label Foam::hexRef8::findMinLevel(const labelList& f) const
{
    label minLevel = labelMax;
    label minFp = -1;

    forAll(f, fp)
    {
        label level = pointLevel_[f[fp]];

        if (level < minLevel)
        {
            minLevel = level;
            minFp = fp;
        }
    }

    return minFp;
}


// Get point with the highest pointLevel
Foam::label Foam::hexRef8::findMaxLevel(const labelList& f) const
{
    label maxLevel = labelMin;
    label maxFp = -1;

    forAll(f, fp)
    {
        label level = pointLevel_[f[fp]];

        if (level > maxLevel)
        {
            maxLevel = level;
            maxFp = fp;
        }
    }

    return maxFp;
}


Foam::label Foam::hexRef8::countAnchors
(
    const labelList& f,
    const label anchorLevel
) const
{
    label nAnchors = 0;

    forAll(f, fp)
    {
        if (pointLevel_[f[fp]] <= anchorLevel)
        {
            nAnchors++;
        }
    }
    return nAnchors;
}


void Foam::hexRef8::dumpCell(const label cellI) const
{
    OFstream str(mesh_.time().path()/"cell_" + Foam::name(cellI) + ".obj");
    Pout<< "hexRef8 : Dumping cell as obj to " << str.name() << endl;

    const cell& cFaces = mesh_.cells()[cellI];

    Map<label> pointToObjVert;
    label objVertI = 0;

    forAll(cFaces, i)
    {
        const face& f = mesh_.faces()[cFaces[i]];

        forAll(f, fp)
        {
            if (pointToObjVert.insert(f[fp], objVertI))
            {
                meshTools::writeOBJ(str, mesh_.points()[f[fp]]);
                objVertI++;
            }
        }
    }

    forAll(cFaces, i)
    {
        const face& f = mesh_.faces()[cFaces[i]];

        forAll(f, fp)
        {
            label pointI = f[fp];
            label nexPointI = f[f.fcIndex(fp)];

            str << "l " << pointToObjVert[pointI]+1
                << ' ' << pointToObjVert[nexPointI]+1 << nl;
        }
    }
}


// Find point with certain pointLevel. Skip any higher levels.
Foam::label Foam::hexRef8::findLevel
(
    const label faceI,
    const face& f,
    const label startFp,
    const bool searchForward,
    const label wantedLevel
) const
{
    label fp = startFp;

    forAll(f, i)
    {
        label pointI = f[fp];

        if (pointLevel_[pointI] < wantedLevel)
        {
            dumpCell(mesh_.faceOwner()[faceI]);
            if (mesh_.isInternalFace(faceI))
            {
                dumpCell(mesh_.faceNeighbour()[faceI]);
            }

            FatalErrorIn("hexRef8::findLevel(..)")
                << "face:" << f
                << " level:" << UIndirectList<label>(pointLevel_, f)()
                << " startFp:" << startFp
                << " wantedLevel:" << wantedLevel
                << abort(FatalError);
        }
        else if (pointLevel_[pointI] == wantedLevel)
        {
            return fp;
        }

        if (searchForward)
        {
            fp = f.fcIndex(fp);
        }
        else
        {
            fp = f.rcIndex(fp);
        }
    }

    dumpCell(mesh_.faceOwner()[faceI]);
    if (mesh_.isInternalFace(faceI))
    {
        dumpCell(mesh_.faceNeighbour()[faceI]);
    }

    FatalErrorIn("hexRef8::findLevel(..)")
        << "face:" << f
        << " level:" << UIndirectList<label>(pointLevel_, f)()
        << " startFp:" << startFp
        << " wantedLevel:" << wantedLevel
        << abort(FatalError);

    return -1;
}


// Gets cell level such that the face has four points <= level.
Foam::label Foam::hexRef8::getAnchorLevel(const label faceI) const
{
    const face& f = mesh_.faces()[faceI];

    if (f.size() <= 4)
    {
        return pointLevel_[f[findMaxLevel(f)]];
    }
    else
    {
        label ownLevel = cellLevel_[mesh_.faceOwner()[faceI]];

        if (countAnchors(f, ownLevel) == 4)
        {
            return ownLevel;
        }
        else if (countAnchors(f, ownLevel+1) == 4)
        {
            return ownLevel+1;
        }
        else
        {
            return -1;
        }
    }
}


void Foam::hexRef8::checkInternalOrientation
(
    polyTopoChange& meshMod,
    const label cellI,
    const label faceI,
    const point& ownPt,
    const point& neiPt,
    const face& newFace
)
{
    face compactFace(identity(newFace.size()));
    pointField compactPoints(meshMod.points(), newFace);

    vector n(compactFace.normal(compactPoints));

    vector dir(neiPt - ownPt);

    if ((dir & n) < 0)
    {
        FatalErrorIn("checkInternalOrientation(..)")
            << "cell:" << cellI << " old face:" << faceI
            << " newFace:" << newFace << endl
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " neiPt:" << neiPt
            << abort(FatalError);
    }

    vector fcToOwn(compactFace.centre(compactPoints) - ownPt);

    scalar s = (fcToOwn&n) / (dir&n);

    if (s < 0.1 || s > 0.9)
    {
        FatalErrorIn("checkInternalOrientation(..)")
            << "cell:" << cellI << " old face:" << faceI
            << " newFace:" << newFace << endl
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " neiPt:" << neiPt
            << " s:" << s
            << abort(FatalError);
    }
}


void Foam::hexRef8::checkBoundaryOrientation
(
    polyTopoChange& meshMod,
    const label cellI,
    const label faceI,
    const point& ownPt,
    const point& boundaryPt,
    const face& newFace
)
{
    face compactFace(identity(newFace.size()));
    pointField compactPoints(meshMod.points(), newFace);

    vector n(compactFace.normal(compactPoints));

    vector dir(boundaryPt - ownPt);

    if ((dir & n) < 0)
    {
        FatalErrorIn("checkBoundaryOrientation(..)")
            << "cell:" << cellI << " old face:" << faceI
            << " newFace:" << newFace
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " boundaryPt:" << boundaryPt
            << abort(FatalError);
    }

    vector fcToOwn(compactFace.centre(compactPoints) - ownPt);

    scalar s = (fcToOwn&dir) / magSqr(dir);

    if (s < 0.7 || s > 1.3)
    {
        WarningIn("checkBoundaryOrientation(..)")
            << "cell:" << cellI << " old face:" << faceI
            << " newFace:" << newFace
            << " coords:" << compactPoints
            << " ownPt:" << ownPt
            << " boundaryPt:" << boundaryPt
            << " s:" << s
            << endl;
            //<< abort(FatalError);
    }
}


// If p0 and p1 are existing vertices check if edge is split and insert
// splitPoint.
void Foam::hexRef8::insertEdgeSplit
(
    const labelList& edgeMidPoint,
    const label p0,
    const label p1,
    DynamicList<label>& verts
) const
{
    if (p0 < mesh_.nPoints() && p1 < mesh_.nPoints())
    {
        label edgeI = meshTools::findEdge(mesh_, p0, p1);

        if (edgeI != -1 && edgeMidPoint[edgeI] != -1)
        {
            verts.append(edgeMidPoint[edgeI]);
        }
    }
}


// Internal faces are one per edge between anchor points. So one per midPoint
// between the anchor points. Here we store the information on the midPoint
// and if we have enough information:
// - two anchors
// - two face mid points
// we add the face. Note that this routine can get called anywhere from
// two times (two unrefined faces) to four times (two refined faces) so
// the first call that adds the information creates the face.
Foam::label Foam::hexRef8::storeMidPointInfo
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const labelList& cellMidPoint,
    const labelList& edgeMidPoint,
    const label cellI,
    const label faceI,
    const bool faceOrder,
    const label edgeMidPointI,
    const label anchorPointI,
    const label faceMidPointI,

    Map<edge>& midPointToAnchors,
    Map<edge>& midPointToFaceMids,
    polyTopoChange& meshMod
) const
{
    // See if need to store anchors.

    bool changed = false;
    bool haveTwoAnchors = false;

    Map<edge>::iterator edgeMidFnd = midPointToAnchors.find(edgeMidPointI);

    if (edgeMidFnd == midPointToAnchors.end())
    {
        midPointToAnchors.insert(edgeMidPointI, edge(anchorPointI, -1));
    }
    else
    {
        edge& e = edgeMidFnd();

        if (anchorPointI != e[0])
        {
            if (e[1] == -1)
            {
                e[1] = anchorPointI;
                changed = true;
            }
        }

        if (e[0] != -1 && e[1] != -1)
        {
            haveTwoAnchors = true;
        }
    }

    bool haveTwoFaceMids = false;

    Map<edge>::iterator faceMidFnd = midPointToFaceMids.find(edgeMidPointI);

    if (faceMidFnd == midPointToFaceMids.end())
    {
        midPointToFaceMids.insert(edgeMidPointI, edge(faceMidPointI, -1));
    }
    else
    {
        edge& e = faceMidFnd();

        if (faceMidPointI != e[0])
        {
            if (e[1] == -1)
            {
                e[1] = faceMidPointI;
                changed = true;
            }
        }

        if (e[0] != -1 && e[1] != -1)
        {
            haveTwoFaceMids = true;
        }
    }

    // Check if this call of storeMidPointInfo is the one that completed all
    // the necessary information.

    if (changed && haveTwoAnchors && haveTwoFaceMids)
    {
        const edge& anchors = midPointToAnchors[edgeMidPointI];
        const edge& faceMids = midPointToFaceMids[edgeMidPointI];

        label otherFaceMidPointI = faceMids.otherVertex(faceMidPointI);

        // Create face consistent with anchorI being the owner.
        // Note that the edges between the edge mid point and the face mids
        // might be marked for splitting. Note that these edge splits cannot
        // be between cellMid and face mids.

        DynamicList<label> newFaceVerts(4);
        if (faceOrder == (mesh_.faceOwner()[faceI] == cellI))
        {
            newFaceVerts.append(faceMidPointI);

            // Check & insert edge split if any
            insertEdgeSplit
            (
                edgeMidPoint,
                faceMidPointI,  // edge between faceMid and
                edgeMidPointI,  // edgeMid
                newFaceVerts
            );

            newFaceVerts.append(edgeMidPointI);

            insertEdgeSplit
            (
                edgeMidPoint,
                edgeMidPointI,
                otherFaceMidPointI,
                newFaceVerts
            );

            newFaceVerts.append(otherFaceMidPointI);
            newFaceVerts.append(cellMidPoint[cellI]);
        }
        else
        {
            newFaceVerts.append(otherFaceMidPointI);

            insertEdgeSplit
            (
                edgeMidPoint,
                otherFaceMidPointI,
                edgeMidPointI,
                newFaceVerts
            );

            newFaceVerts.append(edgeMidPointI);

            insertEdgeSplit
            (
                edgeMidPoint,
                edgeMidPointI,
                faceMidPointI,
                newFaceVerts
            );

            newFaceVerts.append(faceMidPointI);
            newFaceVerts.append(cellMidPoint[cellI]);
        }

        face newFace;
        newFace.transfer(newFaceVerts);

        label anchorCell0 = getAnchorCell
        (
            cellAnchorPoints,
            cellAddedCells,
            cellI,
            faceI,
            anchorPointI
        );
        label anchorCell1 = getAnchorCell
        (
            cellAnchorPoints,
            cellAddedCells,
            cellI,
            faceI,
            anchors.otherVertex(anchorPointI)
        );


        label own, nei;
        point ownPt, neiPt;

        if (anchorCell0 < anchorCell1)
        {
            own = anchorCell0;
            nei = anchorCell1;

            ownPt = mesh_.points()[anchorPointI];
            neiPt = mesh_.points()[anchors.otherVertex(anchorPointI)];

        }
        else
        {
            own = anchorCell1;
            nei = anchorCell0;
            newFace.flip();

            ownPt = mesh_.points()[anchors.otherVertex(anchorPointI)];
            neiPt = mesh_.points()[anchorPointI];
        }

        if (debug)
        {
            point ownPt, neiPt;

            if (anchorCell0 < anchorCell1)
            {
                ownPt = mesh_.points()[anchorPointI];
                neiPt = mesh_.points()[anchors.otherVertex(anchorPointI)];
            }
            else
            {
                ownPt = mesh_.points()[anchors.otherVertex(anchorPointI)];
                neiPt = mesh_.points()[anchorPointI];
            }

            checkInternalOrientation
            (
                meshMod,
                cellI,
                faceI,
                ownPt,
                neiPt,
                newFace
            );
        }

        return addInternalFace
        (
            meshMod,
            faceI,
            anchorPointI,
            newFace,
            own,
            nei
        );
    }
    else
    {
        return -1;
    }
}


// Creates all the 12 internal faces for cellI.
void Foam::hexRef8::createInternalFaces
(
    const labelListList& cellAnchorPoints,
    const labelListList& cellAddedCells,
    const labelList& cellMidPoint,
    const labelList& faceMidPoint,
    const labelList& faceAnchorLevel,
    const labelList& edgeMidPoint,
    const label cellI,

    polyTopoChange& meshMod
) const
{
    // Find in every face the cellLevel+1 points (from edge subdivision)
    // and the anchor points.

    const cell& cFaces = mesh_.cells()[cellI];
    const label cLevel = cellLevel_[cellI];

    // From edge mid to anchor points
    Map<edge> midPointToAnchors(24);
    // From edge mid to face mids
    Map<edge> midPointToFaceMids(24);

    // Storage for on-the-fly addressing
    DynamicList<label> storage;


    // Running count of number of internal faces added so far.
    label nFacesAdded = 0;

    forAll(cFaces, i)
    {
        label faceI = cFaces[i];

        const face& f = mesh_.faces()[faceI];
        const labelList& fEdges = mesh_.faceEdges(faceI, storage);

        // We are on the cellI side of face f. The face will have 1 or 4
        // cLevel points and lots of higher numbered ones.

        label faceMidPointI = -1;

        label nAnchors = countAnchors(f, cLevel);

        if (nAnchors == 1)
        {
            // Only one anchor point. So the other side of the face has already
            // been split using cLevel+1 and cLevel+2 points.

            // Find the one anchor.
            label anchorFp = -1;

            forAll(f, fp)
            {
                if (pointLevel_[f[fp]] <= cLevel)
                {
                    anchorFp = fp;
                    break;
                }
            }

            // Now the face mid point is the second cLevel+1 point
            label edgeMid = findLevel
            (
                faceI,
                f,
                f.fcIndex(anchorFp),
                true,
                cLevel+1
            );
            label faceMid = findLevel
            (
                faceI,
                f,
                f.fcIndex(edgeMid),
                true,
                cLevel+1
            );

            faceMidPointI = f[faceMid];
        }
        else if (nAnchors == 4)
        {
            // There is no face middle yet but the face will be marked for
            // splitting.

            faceMidPointI = faceMidPoint[faceI];
        }
        else
        {
            dumpCell(mesh_.faceOwner()[faceI]);
            if (mesh_.isInternalFace(faceI))
            {
                dumpCell(mesh_.faceNeighbour()[faceI]);
            }

            FatalErrorIn("createInternalFaces(..)")
                << "nAnchors:" << nAnchors
                << " faceI:" << faceI
                << abort(FatalError);
        }



        // Now loop over all the anchors (might be just one) and store
        // the edge mids connected to it. storeMidPointInfo will collect
        // all the info and combine it all.

        forAll(f, fp0)
        {
            label point0 = f[fp0];

            if (pointLevel_[point0] <= cLevel)
            {
                // Anchor.

                // Walk forward
                // ~~~~~~~~~~~~
                // to cLevel+1 or edgeMidPoint of this level.


                label edgeMidPointI = -1;

                label fp1 = f.fcIndex(fp0);

                if (pointLevel_[f[fp1]] <= cLevel)
                {
                    // Anchor. Edge will be split.
                    label edgeI = fEdges[fp0];

                    edgeMidPointI = edgeMidPoint[edgeI];

                    if (edgeMidPointI == -1)
                    {
                        dumpCell(cellI);

                        const labelList& cPoints = mesh_.cellPoints(cellI);

                        FatalErrorIn("createInternalFaces(..)")
                            << "cell:" << cellI << " cLevel:" << cLevel
                            << " cell points:" << cPoints
                            << " pointLevel:"
                            << UIndirectList<label>(pointLevel_, cPoints)()
                            << " face:" << faceI
                            << " f:" << f
                            << " pointLevel:"
                            << UIndirectList<label>(pointLevel_, f)()
                            << " faceAnchorLevel:" << faceAnchorLevel[faceI]
                            << " faceMidPoint:" << faceMidPoint[faceI]
                            << " faceMidPointI:" << faceMidPointI
                            << " fp:" << fp0
                            << abort(FatalError);
                    }
                }
                else
                {
                    // Search forward in face to clevel+1
                    label edgeMid = findLevel(faceI, f, fp1, true, cLevel+1);

                    edgeMidPointI = f[edgeMid];
                }

                label newFaceI = storeMidPointInfo
                (
                    cellAnchorPoints,
                    cellAddedCells,
                    cellMidPoint,
                    edgeMidPoint,

                    cellI,
                    faceI,
                    true,                   // mid point after anchor
                    edgeMidPointI,          // edgemid
                    point0,                 // anchor
                    faceMidPointI,

                    midPointToAnchors,
                    midPointToFaceMids,
                    meshMod
                );

                if (newFaceI != -1)
                {
                    nFacesAdded++;

                    if (nFacesAdded == 12)
                    {
                        break;
                    }
                }



                // Walk backward
                // ~~~~~~~~~~~~~

                label fpMin1 = f.rcIndex(fp0);

                if (pointLevel_[f[fpMin1]] <= cLevel)
                {
                    // Anchor. Edge will be split.
                    label edgeI = fEdges[fpMin1];

                    edgeMidPointI = edgeMidPoint[edgeI];

                    if (edgeMidPointI == -1)
                    {
                        dumpCell(cellI);

                        const labelList& cPoints = mesh_.cellPoints(cellI);

                        FatalErrorIn("createInternalFaces(..)")
                            << "cell:" << cellI << " cLevel:" << cLevel
                            << " cell points:" << cPoints
                            << " pointLevel:"
                            << UIndirectList<label>(pointLevel_, cPoints)()
                            << " face:" << faceI
                            << " f:" << f
                            << " pointLevel:"
                            << UIndirectList<label>(pointLevel_, f)()
                            << " faceAnchorLevel:" << faceAnchorLevel[faceI]
                            << " faceMidPoint:" << faceMidPoint[faceI]
                            << " faceMidPointI:" << faceMidPointI
                            << " fp:" << fp0
                            << abort(FatalError);
                    }
                }
                else
                {
                    // Search back to clevel+1
                    label edgeMid = findLevel
                    (
                        faceI,
                        f,
                        fpMin1,
                        false,
                        cLevel+1
                    );

                    edgeMidPointI = f[edgeMid];
                }

                newFaceI = storeMidPointInfo
                (
                    cellAnchorPoints,
                    cellAddedCells,
                    cellMidPoint,
                    edgeMidPoint,

                    cellI,
                    faceI,
                    false,                  // mid point before anchor
                    edgeMidPointI,          // edgemid
                    point0,                 // anchor
                    faceMidPointI,

                    midPointToAnchors,
                    midPointToFaceMids,
                    meshMod
                );

                if (newFaceI != -1)
                {
                    nFacesAdded++;

                    if (nFacesAdded == 12)
                    {
                        break;
                    }
                }
            }   // done anchor
        }   // done face

        if (nFacesAdded == 12)
        {
            break;
        }
    }
}


void Foam::hexRef8::walkFaceToMid
(
    const labelList& edgeMidPoint,
    const label cLevel,
    const label faceI,
    const label startFp,
    DynamicList<label>& faceVerts
) const
{
    const face& f = mesh_.faces()[faceI];
    const labelList& fEdges = mesh_.faceEdges(faceI);

    label fp = startFp;

    // Starting from fp store all (1 or 2) vertices until where the face
    // gets split

    while (true)
    {
        if (edgeMidPoint[fEdges[fp]] >= 0)
        {
            faceVerts.append(edgeMidPoint[fEdges[fp]]);
        }

        fp = f.fcIndex(fp);

        if (pointLevel_[f[fp]] <= cLevel)
        {
            // Next anchor. Have already append split point on edge in code
            // above.
            return;
        }
        else if (pointLevel_[f[fp]] == cLevel+1)
        {
            // Mid level
            faceVerts.append(f[fp]);

            return;
        }
        else if (pointLevel_[f[fp]] == cLevel+2)
        {
            // Store and continue to cLevel+1.
            faceVerts.append(f[fp]);
        }
    }
}


// Same as walkFaceToMid but now walk back.
void Foam::hexRef8::walkFaceFromMid
(
    const labelList& edgeMidPoint,
    const label cLevel,
    const label faceI,
    const label startFp,
    DynamicList<label>& faceVerts
) const
{
    const face& f = mesh_.faces()[faceI];
    const labelList& fEdges = mesh_.faceEdges(faceI);

    label fp = f.rcIndex(startFp);

    while (true)
    {
        if (pointLevel_[f[fp]] <= cLevel)
        {
            // anchor.
            break;
        }
        else if (pointLevel_[f[fp]] == cLevel+1)
        {
            // Mid level
            faceVerts.append(f[fp]);
            break;
        }
        else if (pointLevel_[f[fp]] == cLevel+2)
        {
            // Continue to cLevel+1.
        }
        fp = f.rcIndex(fp);
    }

    // Store
    while (true)
    {
        if (edgeMidPoint[fEdges[fp]] >= 0)
        {
            faceVerts.append(edgeMidPoint[fEdges[fp]]);
        }

        fp = f.fcIndex(fp);

        if (fp == startFp)
        {
            break;
        }
        faceVerts.append(f[fp]);
    }
}


// Updates refineCell (cells marked for refinement) so across all faces
// there will be 2:1 consistency after refinement.
Foam::label Foam::hexRef8::faceConsistentRefinement
(
    const bool maxSet,
    PackedBoolList& refineCell
) const
{
    label nChanged = 0;

    // Internal faces.
    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label own = mesh_.faceOwner()[faceI];
        label ownLevel = cellLevel_[own] + refineCell.get(own);

        label nei = mesh_.faceNeighbour()[faceI];
        label neiLevel = cellLevel_[nei] + refineCell.get(nei);

        if (ownLevel > (neiLevel+1))
        {
            if (maxSet)
            {
                refineCell.set(nei);
            }
            else
            {
                refineCell.unset(own);
            }
            nChanged++;
        }
        else if (neiLevel > (ownLevel+1))
        {
            if (maxSet)
            {
                refineCell.set(own);
            }
            else
            {
                refineCell.unset(nei);
            }
            nChanged++;
        }
    }


    // Coupled faces. Swap owner level to get neighbouring cell level.
    // (only boundary faces of neiLevel used)
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

    forAll(neiLevel, i)
    {
        label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];

        neiLevel[i] = cellLevel_[own] + refineCell.get(own);
    }

    // Swap to neighbour
    syncTools::swapBoundaryFaceList(mesh_, neiLevel);

    // Now we have neighbour value see which cells need refinement
    forAll(neiLevel, i)
    {
        label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];
        label ownLevel = cellLevel_[own] + refineCell.get(own);

        if (ownLevel > (neiLevel[i]+1))
        {
            if (!maxSet)
            {
                refineCell.unset(own);
                nChanged++;
            }
        }
        else if (neiLevel[i] > (ownLevel+1))
        {
            if (maxSet)
            {
                refineCell.set(own);
                nChanged++;
            }
        }
    }

    return nChanged;
}


// Debug: check if wanted refinement is compatible with 2:1
void Foam::hexRef8::checkWantedRefinementLevels
(
    const labelList& cellsToRefine
) const
{
    PackedBoolList refineCell(mesh_.nCells());
    forAll(cellsToRefine, i)
    {
        refineCell.set(cellsToRefine[i]);
    }

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label own = mesh_.faceOwner()[faceI];
        label ownLevel = cellLevel_[own] + refineCell.get(own);

        label nei = mesh_.faceNeighbour()[faceI];
        label neiLevel = cellLevel_[nei] + refineCell.get(nei);

        if (mag(ownLevel-neiLevel) > 1)
        {
            dumpCell(own);
            dumpCell(nei);
            FatalErrorIn
            (
                "hexRef8::checkWantedRefinementLevels(const labelList&)"
            )   << "cell:" << own
                << " current level:" << cellLevel_[own]
                << " level after refinement:" << ownLevel
                << nl
                << "neighbour cell:" << nei
                << " current level:" << cellLevel_[nei]
                << " level after refinement:" << neiLevel
                << nl
                << "which does not satisfy 2:1 constraints anymore."
                << abort(FatalError);
        }
    }

    // Coupled faces. Swap owner level to get neighbouring cell level.
    // (only boundary faces of neiLevel used)
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

    forAll(neiLevel, i)
    {
        label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];

        neiLevel[i] = cellLevel_[own] + refineCell.get(own);
    }

    // Swap to neighbour
    syncTools::swapBoundaryFaceList(mesh_, neiLevel);

    // Now we have neighbour value see which cells need refinement
    forAll(neiLevel, i)
    {
        label faceI = i + mesh_.nInternalFaces();

        label own = mesh_.faceOwner()[faceI];
        label ownLevel = cellLevel_[own] + refineCell.get(own);

        if (mag(ownLevel - neiLevel[i]) > 1)
        {
            label patchI = mesh_.boundaryMesh().whichPatch(faceI);

            dumpCell(own);
            FatalErrorIn
            (
                "hexRef8::checkWantedRefinementLevels(const labelList&)"
            )   << "Celllevel does not satisfy 2:1 constraint."
                << " On coupled face "
                << faceI
                << " on patch " << patchI << " "
                << mesh_.boundaryMesh()[patchI].name()
                << " owner cell " << own
                << " current level:" << cellLevel_[own]
                << " level after refinement:" << ownLevel
                << nl
                << " (coupled) neighbour cell will get refinement "
                << neiLevel[i]
                << abort(FatalError);
        }
    }
}


// Set instance for mesh files
void Foam::hexRef8::setInstance(const fileName& inst)
{
    if (debug)
    {
        Pout<< "hexRef8::setInstance(const fileName& inst) : "
            << "Resetting file instance to " << inst << endl;
    }

    cellLevel_.instance() = inst;
    pointLevel_.instance() = inst;
    level0Edge_.instance() = inst;
    history_.instance() = inst;
}


void Foam::hexRef8::collectLevelPoints
(
    const labelList& f,
    const label level,
    DynamicList<label>& points
) const
{
    forAll(f, fp)
    {
        if (pointLevel_[f[fp]] <= level)
        {
            points.append(f[fp]);
        }
    }
}


void Foam::hexRef8::collectLevelPoints
(
    const labelList& meshPoints,
    const labelList& f,
    const label level,
    DynamicList<label>& points
) const
{
    forAll(f, fp)
    {
        label pointI = meshPoints[f[fp]];
        if (pointLevel_[pointI] <= level)
        {
            points.append(pointI);
        }
    }
}


// Return true if we've found 6 quads. faces guaranteed to be outwards pointing.
bool Foam::hexRef8::matchHexShape
(
    const label cellI,
    const label cellLevel,
    DynamicList<face>& quads
) const
{
    const cell& cFaces = mesh_.cells()[cellI];

    // Work arrays
    DynamicList<label> verts(4);
    quads.clear();


    // 1. pick up any faces with four cellLevel points

    forAll(cFaces, i)
    {
        label faceI = cFaces[i];
        const face& f = mesh_.faces()[faceI];

        verts.clear();
        collectLevelPoints(f, cellLevel, verts);
        if (verts.size() == 4)
        {
            if (mesh_.faceOwner()[faceI] != cellI)
            {
                reverse(verts);
            }
            quads.append(face(0));
            labelList& quadVerts = quads.last();
            quadVerts.transfer(verts);
        }
    }


    if (quads.size() < 6)
    {
        Map<labelList> pointFaces(2*cFaces.size());

        forAll(cFaces, i)
        {
            label faceI = cFaces[i];
            const face& f = mesh_.faces()[faceI];

            // Pick up any faces with only one level point.
            // See if there are four of these where the commont point
            // is a level+1 point. This common point is then the mid of
            // a split face.

            verts.clear();
            collectLevelPoints(f, cellLevel, verts);
            if (verts.size() == 1)
            {
                // Add to pointFaces for any level+1 point (this might be
                // a midpoint of a split face)
                forAll(f, fp)
                {
                    label pointI = f[fp];
                    if (pointLevel_[pointI] == cellLevel+1)
                    {
                        Map<labelList>::iterator iter =
                            pointFaces.find(pointI);
                        if (iter != pointFaces.end())
                        {
                            labelList& pFaces = iter();
                            if (findIndex(pFaces, faceI) == -1)
                            {
                                pFaces.append(faceI);
                            }
                        }
                        else
                        {
                            pointFaces.insert
                            (
                                pointI,
                                labelList(1, faceI)
                            );
                        }
                    }
                }
            }
        }

        // 2. Check if we've collected any midPoints.
        forAllConstIter(Map<labelList>, pointFaces, iter)
        {
            const labelList& pFaces = iter();

            if (pFaces.size() == 4)
            {
                // Collect and orient.
                faceList fourFaces(pFaces.size());
                forAll(pFaces, pFaceI)
                {
                    label faceI = pFaces[pFaceI];
                    const face& f = mesh_.faces()[faceI];
                    if (mesh_.faceOwner()[faceI] == cellI)
                    {
                        fourFaces[pFaceI] = f;
                    }
                    else
                    {
                        fourFaces[pFaceI] = f.reverseFace();
                    }
                }

                primitivePatch bigFace
                (
                    SubList<face>(fourFaces, fourFaces.size()),
                    mesh_.points()
                );
                const labelListList& edgeLoops = bigFace.edgeLoops();

                if (edgeLoops.size() == 1)
                {
                    // Collect the 4 cellLevel points
                    verts.clear();
                    collectLevelPoints
                    (
                        bigFace.meshPoints(),
                        bigFace.edgeLoops()[0],
                        cellLevel,
                        verts
                    );

                    if (verts.size() == 4)
                    {
                        quads.append(face(0));
                        labelList& quadVerts = quads.last();
                        quadVerts.transfer(verts);
                    }
                }
            }
        }
    }

    return (quads.size() == 6);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh, read refinement data
Foam::hexRef8::hexRef8(const polyMesh& mesh, const bool readHistory)
:
    mesh_(mesh),
    cellLevel_
    (
        IOobject
        (
            "cellLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        labelList(mesh_.nCells(), 0)
    ),
    pointLevel_
    (
        IOobject
        (
            "pointLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        labelList(mesh_.nPoints(), 0)
    ),
    level0Edge_
    (
        IOobject
        (
            "level0Edge",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        dimensionedScalar("level0Edge", dimLength, getLevel0EdgeLength())
    ),
    history_
    (
        IOobject
        (
            "refinementHistory",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (readHistory ? mesh_.nCells() : 0)  // All cells visible if not be read
    ),
    faceRemover_(mesh_, GREAT),     // merge boundary faces wherever possible
    savedPointLevel_(0),
    savedCellLevel_(0)
{
    if (readHistory)
    {
        // Make sure we don't use the master-only reading. Bit of a hack for
        // now.
        regIOobject::fileCheckTypes oldType =
            regIOobject::fileModificationChecking;
        regIOobject::fileModificationChecking = regIOobject::timeStamp;
        history_.readOpt() = IOobject::READ_IF_PRESENT;
        if (history_.headerOk())
        {
            history_.read();
        }
        regIOobject::fileModificationChecking = oldType;
    }

    if (history_.active() && history_.visibleCells().size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "hexRef8::hexRef8(const polyMesh&)"
        )   << "History enabled but number of visible cells "
            << history_.visibleCells().size() << " in "
            << history_.objectPath()
            << " is not equal to the number of cells in the mesh "
            << mesh_.nCells()
            << abort(FatalError);
    }

    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorIn
        (
            "hexRef8::hexRef8(const polyMesh&)"
        )   << "Restarted from inconsistent cellLevel or pointLevel files."
            << endl
            << "cellLevel file " << cellLevel_.objectPath() << endl
            << "pointLevel file " << pointLevel_.objectPath() << endl
            << "Number of cells in mesh:" << mesh_.nCells()
            << " does not equal size of cellLevel:" << cellLevel_.size() << endl
            << "Number of points in mesh:" << mesh_.nPoints()
            << " does not equal size of pointLevel:" << pointLevel_.size()
            << abort(FatalError);
    }


    // Check refinement levels for consistency
    checkRefinementLevels(-1, labelList(0));


    // Check initial mesh for consistency

    //if (debug)
    {
        checkMesh();
    }
}


// Construct from components
Foam::hexRef8::hexRef8
(
    const polyMesh& mesh,
    const labelList& cellLevel,
    const labelList& pointLevel,
    const refinementHistory& history,
    const scalar level0Edge
)
:
    mesh_(mesh),
    cellLevel_
    (
        IOobject
        (
            "cellLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        cellLevel
    ),
    pointLevel_
    (
        IOobject
        (
            "pointLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pointLevel
    ),
    level0Edge_
    (
        IOobject
        (
            "level0Edge",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dimensionedScalar
        (
            "level0Edge",
            dimLength,
            (level0Edge >= 0 ? level0Edge : getLevel0EdgeLength())
        )
    ),
    history_
    (
        IOobject
        (
            "refinementHistory",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        history
    ),
    faceRemover_(mesh_, GREAT),     // merge boundary faces wherever possible
    savedPointLevel_(0),
    savedCellLevel_(0)
{
    if (history_.active() && history_.visibleCells().size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "hexRef8::hexRef8(const polyMesh&, const labelList&"
            ", const labelList&, const refinementHistory&)"
        )   << "History enabled but number of visible cells in it "
            << history_.visibleCells().size()
            << " is not equal to the number of cells in the mesh "
            << mesh_.nCells() << abort(FatalError);
    }

    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorIn
        (
            "hexRef8::hexRef8(const polyMesh&, const labelList&"
            ", const labelList&, const refinementHistory&)"
        )   << "Incorrect cellLevel or pointLevel size." << endl
            << "Number of cells in mesh:" << mesh_.nCells()
            << " does not equal size of cellLevel:" << cellLevel_.size() << endl
            << "Number of points in mesh:" << mesh_.nPoints()
            << " does not equal size of pointLevel:" << pointLevel_.size()
            << abort(FatalError);
    }

    // Check refinement levels for consistency
    checkRefinementLevels(-1, labelList(0));


    // Check initial mesh for consistency

    //if (debug)
    {
        checkMesh();
    }
}


// Construct from components
Foam::hexRef8::hexRef8
(
    const polyMesh& mesh,
    const labelList& cellLevel,
    const labelList& pointLevel,
    const scalar level0Edge
)
:
    mesh_(mesh),
    cellLevel_
    (
        IOobject
        (
            "cellLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        cellLevel
    ),
    pointLevel_
    (
        IOobject
        (
            "pointLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pointLevel
    ),
    level0Edge_
    (
        IOobject
        (
            "level0Edge",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dimensionedScalar
        (
            "level0Edge",
            dimLength,
            (level0Edge >= 0 ? level0Edge : getLevel0EdgeLength())
        )
    ),
    history_
    (
        IOobject
        (
            "refinementHistory",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        List<refinementHistory::splitCell8>(0),
        labelList(0)
    ),
    faceRemover_(mesh_, GREAT),     // merge boundary faces wherever possible
    savedPointLevel_(0),
    savedCellLevel_(0)
{
    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorIn
        (
            "hexRef8::hexRef8(const polyMesh&, const labelList&"
            ", const labelList&)"
        )   << "Incorrect cellLevel or pointLevel size." << endl
            << "Number of cells in mesh:" << mesh_.nCells()
            << " does not equal size of cellLevel:" << cellLevel_.size() << endl
            << "Number of points in mesh:" << mesh_.nPoints()
            << " does not equal size of pointLevel:" << pointLevel_.size()
            << abort(FatalError);
    }

    // Check refinement levels for consistency
    checkRefinementLevels(-1, labelList(0));

    // Check initial mesh for consistency

    //if (debug)
    {
        checkMesh();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::hexRef8::consistentRefinement
(
    const labelList& cellsToRefine,
    const bool maxSet
) const
{
    // Loop, modifying cellsToRefine, until no more changes to due to 2:1
    // conflicts.
    // maxSet = false : unselect cells to refine
    // maxSet = true  : select cells to refine

    // Go to straight boolList.
    PackedBoolList refineCell(mesh_.nCells());
    forAll(cellsToRefine, i)
    {
        refineCell.set(cellsToRefine[i]);
    }

    while (true)
    {
        label nChanged = faceConsistentRefinement(maxSet, refineCell);

        reduce(nChanged, sumOp<label>());

        if (debug)
        {
            Pout<< "hexRef8::consistentRefinement : Changed " << nChanged
                << " refinement levels due to 2:1 conflicts."
                << endl;
        }

        if (nChanged == 0)
        {
            break;
        }
    }


    // Convert back to labelList.
    label nRefined = 0;

    forAll(refineCell, cellI)
    {
        if (refineCell.get(cellI))
        {
            nRefined++;
        }
    }

    labelList newCellsToRefine(nRefined);
    nRefined = 0;

    forAll(refineCell, cellI)
    {
        if (refineCell.get(cellI))
        {
            newCellsToRefine[nRefined++] = cellI;
        }
    }

    if (debug)
    {
        checkWantedRefinementLevels(newCellsToRefine);
    }

    return newCellsToRefine;
}


// Given a list of cells to refine determine additional cells to refine
// such that the overall refinement:
// - satisfies maxFaceDiff (e.g. 2:1) across neighbouring faces
// - satisfies maxPointDiff (e.g. 4:1) across selected point connected
//   cells. This is used to ensure that e.g. cells on the surface are not
//   point connected to cells which are 8 times smaller.
Foam::labelList Foam::hexRef8::consistentSlowRefinement
(
    const label maxFaceDiff,
    const labelList& cellsToRefine,
    const labelList& facesToCheck,
    const label maxPointDiff,
    const labelList& pointsToCheck
) const
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();


    if (maxFaceDiff <= 0)
    {
        FatalErrorIn
        (
            "hexRef8::consistentSlowRefinement"
            "(const label, const labelList&, const labelList&"
            ", const label, const labelList&)"
        )   << "Illegal maxFaceDiff " << maxFaceDiff << nl
            << "Value should be >= 1" << exit(FatalError);
    }


    // Bit tricky. Say we want a distance of three cells between two
    // consecutive refinement levels. This is done by using FaceCellWave to
    // transport out the new refinement level. It gets decremented by one
    // every cell it crosses so if we initialize it to maxFaceDiff
    // we will get a field everywhere that tells us whether an unselected cell
    // needs refining as well.


    // Initial information about (distance to) cellLevel on all cells
    List<refinementData> allCellInfo(mesh_.nCells());

    // Initial information about (distance to) cellLevel on all faces
    List<refinementData> allFaceInfo(mesh_.nFaces());

    forAll(allCellInfo, cellI)
    {
        // maxFaceDiff since refinementData counts both
        // faces and cells.
        allCellInfo[cellI] = refinementData
        (
            maxFaceDiff*(cellLevel_[cellI]+1),// when cell is to be refined
            maxFaceDiff*cellLevel_[cellI]     // current level
        );
    }

    // Cells to be refined will have cellLevel+1
    forAll(cellsToRefine, i)
    {
        label cellI = cellsToRefine[i];

        allCellInfo[cellI].count() = allCellInfo[cellI].refinementCount();
    }


    // Labels of seed faces
    DynamicList<label> seedFaces(mesh_.nFaces()/100);
    // refinementLevel data on seed faces
    DynamicList<refinementData> seedFacesInfo(mesh_.nFaces()/100);

    // Dummy additional info for FaceCellWave
    int dummyTrackData = 0;


    // Additional buffer layer thickness by changing initial count. Usually
    // this happens on boundary faces. Bit tricky. Use allFaceInfo to mark
    // off thus marked faces so they're skipped in the next loop.
    forAll(facesToCheck, i)
    {
        label faceI = facesToCheck[i];

        if (allFaceInfo[faceI].valid(dummyTrackData))
        {
            // Can only occur if face has already gone through loop below.
            FatalErrorIn
            (
                "hexRef8::consistentSlowRefinement"
                "(const label, const labelList&, const labelList&"
                ", const label, const labelList&)"
            )   << "Argument facesToCheck seems to have duplicate entries!"
                << endl
                << "face:" << faceI << " occurs at positions "
                << findIndices(facesToCheck, faceI)
                << abort(FatalError);
        }


        const refinementData& ownData = allCellInfo[faceOwner[faceI]];

        if (mesh_.isInternalFace(faceI))
        {
            // Seed face if neighbouring cell (after possible refinement)
            // will be refined one more than the current owner or neighbour.

            const refinementData& neiData = allCellInfo[faceNeighbour[faceI]];

            label faceCount;
            label faceRefineCount;
            if (neiData.count() > ownData.count())
            {
                faceCount = neiData.count() + maxFaceDiff;
                faceRefineCount = faceCount + maxFaceDiff;
            }
            else
            {
                faceCount = ownData.count() + maxFaceDiff;
                faceRefineCount = faceCount + maxFaceDiff;
            }

            seedFaces.append(faceI);
            seedFacesInfo.append
            (
                refinementData
                (
                    faceRefineCount,
                    faceCount
                )
            );
            allFaceInfo[faceI] = seedFacesInfo.last();
        }
        else
        {
            label faceCount = ownData.count() + maxFaceDiff;
            label faceRefineCount = faceCount + maxFaceDiff;

            seedFaces.append(faceI);
            seedFacesInfo.append
            (
                refinementData
                (
                    faceRefineCount,
                    faceCount
                )
            );
            allFaceInfo[faceI] = seedFacesInfo.last();
        }
    }


    // Just seed with all faces inbetween different refinement levels for now
    // (alternatively only seed faces on cellsToRefine but that gives problems
    //  if no cells to refine)
    forAll(faceNeighbour, faceI)
    {
        // Check if face already handled in loop above
        if (!allFaceInfo[faceI].valid(dummyTrackData))
        {
            label own = faceOwner[faceI];
            label nei = faceNeighbour[faceI];

            // Seed face with transported data from highest cell.

            if (allCellInfo[own].count() > allCellInfo[nei].count())
            {
                allFaceInfo[faceI].updateFace
                (
                    mesh_,
                    faceI,
                    own,
                    allCellInfo[own],
                    FaceCellWave<refinementData, int>::propagationTol(),
                    dummyTrackData
                );
                seedFaces.append(faceI);
                seedFacesInfo.append(allFaceInfo[faceI]);
            }
            else if (allCellInfo[own].count() < allCellInfo[nei].count())
            {
                allFaceInfo[faceI].updateFace
                (
                    mesh_,
                    faceI,
                    nei,
                    allCellInfo[nei],
                    FaceCellWave<refinementData, int>::propagationTol(),
                    dummyTrackData
                );
                seedFaces.append(faceI);
                seedFacesInfo.append(allFaceInfo[faceI]);
            }
        }
    }

    // Seed all boundary faces with owner value. This is to make sure that
    // they are visited (probably only important for coupled faces since
    // these need to be visited from both sides)
    for (label faceI = mesh_.nInternalFaces(); faceI < mesh_.nFaces(); faceI++)
    {
        // Check if face already handled in loop above
        if (!allFaceInfo[faceI].valid(dummyTrackData))
        {
            label own = faceOwner[faceI];

            // Seed face with transported data from owner.
            refinementData faceData;
            faceData.updateFace
            (
                mesh_,
                faceI,
                own,
                allCellInfo[own],
                FaceCellWave<refinementData, int>::propagationTol(),
                dummyTrackData
            );
            seedFaces.append(faceI);
            seedFacesInfo.append(faceData);
        }
    }


    // face-cell-face transport engine
    FaceCellWave<refinementData, int> levelCalc
    (
        mesh_,
        allFaceInfo,
        allCellInfo,
        dummyTrackData
    );

    while (true)
    {
        if (debug)
        {
            Pout<< "hexRef8::consistentSlowRefinement : Seeded "
                << seedFaces.size() << " faces between cells with different"
                << " refinement level." << endl;
        }

        // Set seed faces
        levelCalc.setFaceInfo(seedFaces.shrink(), seedFacesInfo.shrink());
        seedFaces.clear();
        seedFacesInfo.clear();

        // Iterate until no change. Now 2:1 face difference should be satisfied
        levelCalc.iterate(mesh_.globalData().nTotalFaces()+1);


        // Now check point-connected cells (face-connected cells already ok):
        // - get per point max of connected cells
        // - sync across coupled points
        // - check cells against above point max

        if (maxPointDiff == -1)
        {
            // No need to do any point checking.
            break;
        }

        // Determine per point the max cell level. (done as count, not
        // as cell level purely for ease)
        labelList maxPointCount(mesh_.nPoints(), 0);

        forAll(maxPointCount, pointI)
        {
            label& pLevel = maxPointCount[pointI];

            const labelList& pCells = mesh_.pointCells(pointI);

            forAll(pCells, i)
            {
                pLevel = max(pLevel, allCellInfo[pCells[i]].count());
            }
        }

        // Sync maxPointCount to neighbour
        syncTools::syncPointList
        (
            mesh_,
            maxPointCount,
            maxEqOp<label>(),
            labelMin            // null value
        );

        // Update allFaceInfo from maxPointCount for all points to check
        // (usually on boundary faces)

        // Per face the new refinement data
        Map<refinementData> changedFacesInfo(pointsToCheck.size());

        forAll(pointsToCheck, i)
        {
            label pointI = pointsToCheck[i];

            // Loop over all cells using the point and check whether their
            // refinement level is much less than the maximum.

            const labelList& pCells = mesh_.pointCells(pointI);

            forAll(pCells, pCellI)
            {
                label cellI = pCells[pCellI];

                refinementData& cellInfo = allCellInfo[cellI];

                if
                (
                   !cellInfo.isRefined()
                 && (
                        maxPointCount[pointI]
                      > cellInfo.count() + maxFaceDiff*maxPointDiff
                    )
                )
                {
                    // Mark cell for refinement
                    cellInfo.count() = cellInfo.refinementCount();

                    // Insert faces of cell as seed faces.
                    const cell& cFaces = mesh_.cells()[cellI];

                    forAll(cFaces, cFaceI)
                    {
                        label faceI = cFaces[cFaceI];

                        refinementData faceData;
                        faceData.updateFace
                        (
                            mesh_,
                            faceI,
                            cellI,
                            cellInfo,
                            FaceCellWave<refinementData, int>::propagationTol(),
                            dummyTrackData
                        );

                        if (faceData.count() > allFaceInfo[faceI].count())
                        {
                            changedFacesInfo.insert(faceI, faceData);
                        }
                    }
                }
            }
        }

        label nChanged = changedFacesInfo.size();
        reduce(nChanged, sumOp<label>());

        if (nChanged == 0)
        {
            break;
        }


        // Transfer into seedFaces, seedFacesInfo
        seedFaces.setCapacity(changedFacesInfo.size());
        seedFacesInfo.setCapacity(changedFacesInfo.size());

        forAllConstIter(Map<refinementData>, changedFacesInfo, iter)
        {
            seedFaces.append(iter.key());
            seedFacesInfo.append(iter());
        }
    }


    if (debug)
    {
        for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
        {
            label own = mesh_.faceOwner()[faceI];
            label ownLevel =
                cellLevel_[own]
              + (allCellInfo[own].isRefined() ? 1 : 0);

            label nei = mesh_.faceNeighbour()[faceI];
            label neiLevel =
                cellLevel_[nei]
              + (allCellInfo[nei].isRefined() ? 1 : 0);

            if (mag(ownLevel-neiLevel) > 1)
            {
                dumpCell(own);
                dumpCell(nei);
                FatalErrorIn
                (
                    "hexRef8::consistentSlowRefinement"
                    "(const label, const labelList&, const labelList&"
                    ", const label, const labelList&)"
                )   << "cell:" << own
                    << " current level:" << cellLevel_[own]
                    << " current refData:" << allCellInfo[own]
                    << " level after refinement:" << ownLevel
                    << nl
                    << "neighbour cell:" << nei
                    << " current level:" << cellLevel_[nei]
                    << " current refData:" << allCellInfo[nei]
                    << " level after refinement:" << neiLevel
                    << nl
                    << "which does not satisfy 2:1 constraints anymore." << nl
                    << "face:" << faceI << " faceRefData:" << allFaceInfo[faceI]
                    << abort(FatalError);
            }
        }


        // Coupled faces. Swap owner level to get neighbouring cell level.
        // (only boundary faces of neiLevel used)

        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
        labelList neiCount(mesh_.nFaces()-mesh_.nInternalFaces());
        labelList neiRefCount(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(neiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];
            neiLevel[i] = cellLevel_[own];
            neiCount[i] = allCellInfo[own].count();
            neiRefCount[i] = allCellInfo[own].refinementCount();
        }

        // Swap to neighbour
        syncTools::swapBoundaryFaceList(mesh_, neiLevel);
        syncTools::swapBoundaryFaceList(mesh_, neiCount);
        syncTools::swapBoundaryFaceList(mesh_, neiRefCount);

        // Now we have neighbour value see which cells need refinement
        forAll(neiLevel, i)
        {
            label faceI = i+mesh_.nInternalFaces();

            label own = mesh_.faceOwner()[faceI];
            label ownLevel =
                cellLevel_[own]
              + (allCellInfo[own].isRefined() ? 1 : 0);

            label nbrLevel =
                neiLevel[i]
              + ((neiCount[i] >= neiRefCount[i]) ? 1 : 0);

            if (mag(ownLevel - nbrLevel) > 1)
            {
                dumpCell(own);
                label patchI = mesh_.boundaryMesh().whichPatch(faceI);

                FatalErrorIn
                (
                    "hexRef8::consistentSlowRefinement"
                    "(const label, const labelList&, const labelList&"
                    ", const label, const labelList&)"
                )   << "Celllevel does not satisfy 2:1 constraint."
                    << " On coupled face "
                    << faceI
                    << " refData:" << allFaceInfo[faceI]
                    << " on patch " << patchI << " "
                    << mesh_.boundaryMesh()[patchI].name() << nl
                    << "owner cell " << own
                    << " current level:" << cellLevel_[own]
                    << " current count:" << allCellInfo[own].count()
                    << " current refCount:"
                    << allCellInfo[own].refinementCount()
                    << " level after refinement:" << ownLevel
                    << nl
                    << "(coupled) neighbour cell"
                    << " has current level:" << neiLevel[i]
                    << " current count:" << neiCount[i]
                    << " current refCount:" << neiRefCount[i]
                    << " level after refinement:" << nbrLevel
                    << abort(FatalError);
            }
        }
    }

    // Convert back to labelList of cells to refine.

    label nRefined = 0;

    forAll(allCellInfo, cellI)
    {
        if (allCellInfo[cellI].isRefined())
        {
            nRefined++;
        }
    }

    // Updated list of cells to refine
    labelList newCellsToRefine(nRefined);
    nRefined = 0;

    forAll(allCellInfo, cellI)
    {
        if (allCellInfo[cellI].isRefined())
        {
            newCellsToRefine[nRefined++] = cellI;
        }
    }

    if (debug)
    {
        Pout<< "hexRef8::consistentSlowRefinement : From "
            << cellsToRefine.size() << " to " << newCellsToRefine.size()
            << " cells to refine." << endl;
    }

    return newCellsToRefine;
}


Foam::labelList Foam::hexRef8::consistentSlowRefinement2
(
    const label maxFaceDiff,
    const labelList& cellsToRefine,
    const labelList& facesToCheck
) const
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    if (maxFaceDiff <= 0)
    {
        FatalErrorIn
        (
            "hexRef8::consistentSlowRefinement2"
            "(const label, const labelList&, const labelList&)"
        )   << "Illegal maxFaceDiff " << maxFaceDiff << nl
            << "Value should be >= 1" << exit(FatalError);
    }

    const scalar level0Size = 2*maxFaceDiff*level0EdgeLength();


    // Bit tricky. Say we want a distance of three cells between two
    // consecutive refinement levels. This is done by using FaceCellWave to
    // transport out the 'refinement shell'. Anything inside the refinement
    // shell (given by a distance) gets marked for refinement.

    // Initial information about (distance to) cellLevel on all cells
    List<refinementDistanceData> allCellInfo(mesh_.nCells());

    // Initial information about (distance to) cellLevel on all faces
    List<refinementDistanceData> allFaceInfo(mesh_.nFaces());

    // Dummy additional info for FaceCellWave
    int dummyTrackData = 0;


    // Mark cells with wanted refinement level
    forAll(cellsToRefine, i)
    {
        label cellI = cellsToRefine[i];

        allCellInfo[cellI] = refinementDistanceData
        (
            level0Size,
            mesh_.cellCentres()[cellI],
            cellLevel_[cellI]+1             // wanted refinement
        );
    }
    // Mark all others with existing refinement level
    forAll(allCellInfo, cellI)
    {
        if (!allCellInfo[cellI].valid(dummyTrackData))
        {
            allCellInfo[cellI] = refinementDistanceData
            (
                level0Size,
                mesh_.cellCentres()[cellI],
                cellLevel_[cellI]           // wanted refinement
            );
        }
    }


    // Labels of seed faces
    DynamicList<label> seedFaces(mesh_.nFaces()/100);
    // refinementLevel data on seed faces
    DynamicList<refinementDistanceData> seedFacesInfo(mesh_.nFaces()/100);

    const pointField& cc = mesh_.cellCentres();

    forAll(facesToCheck, i)
    {
        label faceI = facesToCheck[i];

        if (allFaceInfo[faceI].valid(dummyTrackData))
        {
            // Can only occur if face has already gone through loop below.
            FatalErrorIn
            (
                "hexRef8::consistentSlowRefinement2"
                "(const label, const labelList&, const labelList&)"
            )   << "Argument facesToCheck seems to have duplicate entries!"
                << endl
                << "face:" << faceI << " occurs at positions "
                << findIndices(facesToCheck, faceI)
                << abort(FatalError);
        }

        label own = faceOwner[faceI];

        label ownLevel =
        (
            allCellInfo[own].valid(dummyTrackData)
          ? allCellInfo[own].originLevel()
          : cellLevel_[own]
        );

        if (!mesh_.isInternalFace(faceI))
        {
            // Do as if boundary face would have neighbour with one higher
            // refinement level.
            const point& fc = mesh_.faceCentres()[faceI];

            refinementDistanceData neiData
            (
                level0Size,
                2*fc - cc[own],    // est'd cell centre
                ownLevel+1
            );

            allFaceInfo[faceI].updateFace
            (
                mesh_,
                faceI,
                own,        // not used (should be nei)
                neiData,
                FaceCellWave<refinementDistanceData, int>::propagationTol(),
                dummyTrackData
            );
        }
        else
        {
            label nei = faceNeighbour[faceI];

            label neiLevel =
            (
                allCellInfo[nei].valid(dummyTrackData)
              ? allCellInfo[nei].originLevel()
              : cellLevel_[nei]
            );

            if (ownLevel == neiLevel)
            {
                // Fake as if nei>own or own>nei (whichever one 'wins')
                allFaceInfo[faceI].updateFace
                (
                    mesh_,
                    faceI,
                    nei,
                    refinementDistanceData(level0Size, cc[nei], neiLevel+1),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
                allFaceInfo[faceI].updateFace
                (
                    mesh_,
                    faceI,
                    own,
                    refinementDistanceData(level0Size, cc[own], ownLevel+1),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
            }
            else
            {
                // Difference in level anyway.
                allFaceInfo[faceI].updateFace
                (
                    mesh_,
                    faceI,
                    nei,
                    refinementDistanceData(level0Size, cc[nei], neiLevel),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
                allFaceInfo[faceI].updateFace
                (
                    mesh_,
                    faceI,
                    own,
                    refinementDistanceData(level0Size, cc[own], ownLevel),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
            }
        }
        seedFaces.append(faceI);
        seedFacesInfo.append(allFaceInfo[faceI]);
    }


    // Create some initial seeds to start walking from. This is only if there
    // are no facesToCheck.
    // Just seed with all faces inbetween different refinement levels for now
    forAll(faceNeighbour, faceI)
    {
        // Check if face already handled in loop above
        if (!allFaceInfo[faceI].valid(dummyTrackData))
        {
            label own = faceOwner[faceI];

            label ownLevel =
            (
                allCellInfo[own].valid(dummyTrackData)
              ? allCellInfo[own].originLevel()
              : cellLevel_[own]
            );

            label nei = faceNeighbour[faceI];

            label neiLevel =
            (
                allCellInfo[nei].valid(dummyTrackData)
              ? allCellInfo[nei].originLevel()
              : cellLevel_[nei]
            );

            if (ownLevel > neiLevel)
            {
                // Set face to owner data. (since face not yet would be copy)
                seedFaces.append(faceI);
                allFaceInfo[faceI].updateFace
                (
                    mesh_,
                    faceI,
                    own,
                    refinementDistanceData(level0Size, cc[own], ownLevel),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
                seedFacesInfo.append(allFaceInfo[faceI]);
            }
            else if (neiLevel > ownLevel)
            {
                seedFaces.append(faceI);
                allFaceInfo[faceI].updateFace
                (
                    mesh_,
                    faceI,
                    nei,
                    refinementDistanceData(level0Size, cc[nei], neiLevel),
                    FaceCellWave<refinementDistanceData, int>::propagationTol(),
                    dummyTrackData
                );
                seedFacesInfo.append(allFaceInfo[faceI]);
            }
        }
    }

    seedFaces.shrink();
    seedFacesInfo.shrink();

    // face-cell-face transport engine
    FaceCellWave<refinementDistanceData, int> levelCalc
    (
        mesh_,
        seedFaces,
        seedFacesInfo,
        allFaceInfo,
        allCellInfo,
        mesh_.globalData().nTotalCells()+1,
        dummyTrackData
    );


    //if (debug)
    //{
    //    // Dump wanted level
    //    volScalarField wantedLevel
    //    (
    //        IOobject
    //        (
    //            "wantedLevel",
    //            fMesh.time().timeName(),
    //            fMesh,
    //            IOobject::NO_READ,
    //            IOobject::AUTO_WRITE,
    //            false
    //        ),
    //        fMesh,
    //        dimensionedScalar("zero", dimless, 0)
    //    );
    //
    //    forAll(wantedLevel, cellI)
    //    {
    //        wantedLevel[cellI] = allCellInfo[cellI].wantedLevel(cc[cellI]);
    //    }
    //
    //    Pout<< "Writing " << wantedLevel.objectPath() << endl;
    //    wantedLevel.write();
    //}


    // Convert back to labelList of cells to refine.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // 1. Force original refinement cells to be picked up by setting the
    // originLevel of input cells to be a very large level (but within range
    // of 1<< shift inside refinementDistanceData::wantedLevel)
    forAll(cellsToRefine, i)
    {
        label cellI = cellsToRefine[i];

        allCellInfo[cellI].originLevel() = sizeof(label)*8-2;
        allCellInfo[cellI].origin() = cc[cellI];
    }

    // 2. Extend to 2:1. I don't understand yet why this is not done
    // 2. Extend to 2:1. For non-cube cells the scalar distance does not work
    // so make sure it at least provides 2:1.
    PackedBoolList refineCell(mesh_.nCells());
    forAll(allCellInfo, cellI)
    {
        label wanted = allCellInfo[cellI].wantedLevel(cc[cellI]);

        if (wanted > cellLevel_[cellI]+1)
        {
            refineCell.set(cellI);
        }
    }
    faceConsistentRefinement(true, refineCell);

    while (true)
    {
        label nChanged = faceConsistentRefinement(true, refineCell);

        reduce(nChanged, sumOp<label>());

        if (debug)
        {
            Pout<< "hexRef8::consistentSlowRefinement2 : Changed " << nChanged
                << " refinement levels due to 2:1 conflicts."
                << endl;
        }

        if (nChanged == 0)
        {
            break;
        }
    }

    // 3. Convert back to labelList.
    label nRefined = 0;

    forAll(refineCell, cellI)
    {
//        if (refineCell.get(cellI))
        if (refineCell[cellI])
        {
            nRefined++;
        }
    }

    labelList newCellsToRefine(nRefined);
    nRefined = 0;

    forAll(refineCell, cellI)
    {
//        if (refineCell.get(cellI))
        if (refineCell[cellI])
        {
            newCellsToRefine[nRefined++] = cellI;
        }
    }

    if (debug)
    {
        Pout<< "hexRef8::consistentSlowRefinement2 : From "
            << cellsToRefine.size() << " to " << newCellsToRefine.size()
            << " cells to refine." << endl;

        // Check that newCellsToRefine obeys at least 2:1.

        {
            cellSet cellsIn(mesh_, "cellsToRefineIn", cellsToRefine);
            Pout<< "hexRef8::consistentSlowRefinement2 : writing "
                << cellsIn.size() << " to cellSet "
                << cellsIn.objectPath() << endl;
            cellsIn.write();
        }
        {
            cellSet cellsOut(mesh_, "cellsToRefineOut", newCellsToRefine);
            Pout<< "hexRef8::consistentSlowRefinement2 : writing "
                << cellsOut.size() << " to cellSet "
                << cellsOut.objectPath() << endl;
            cellsOut.write();
        }

        // Extend to 2:1
        PackedBoolList refineCell(mesh_.nCells());
        forAll(newCellsToRefine, i)
        {
            refineCell.set(newCellsToRefine[i]);
        }
        const PackedBoolList savedRefineCell(refineCell);

        label nChanged = faceConsistentRefinement(true, refineCell);

        {
            cellSet cellsOut2
            (
                mesh_, "cellsToRefineOut2", newCellsToRefine.size()
            );
            forAll(refineCell, cellI)
            {
                if (refineCell.get(cellI))
                {
                    cellsOut2.insert(cellI);
                }
            }
            Pout<< "hexRef8::consistentSlowRefinement2 : writing "
                << cellsOut2.size() << " to cellSet "
                << cellsOut2.objectPath() << endl;
            cellsOut2.write();
        }

        if (nChanged > 0)
        {
            forAll(refineCell, cellI)
            {
                if (refineCell.get(cellI) && !savedRefineCell.get(cellI))
                {
                    dumpCell(cellI);
                    FatalErrorIn
                    (
                        "hexRef8::consistentSlowRefinement2"
                        "(const label, const labelList&, const labelList&)"
                    )   << "Cell:" << cellI << " cc:"
                        << mesh_.cellCentres()[cellI]
                        << " was not marked for refinement but does not obey"
                        << " 2:1 constraints."
                        << abort(FatalError);
                }
            }
        }
    }

    return newCellsToRefine;
}


// Top level driver to insert topo changes to do all refinement.
Foam::labelListList Foam::hexRef8::setRefinement
(
    const labelList& cellLabels,
    polyTopoChange& meshMod
)
{
    if (debug)
    {
        Pout<< "hexRef8::setRefinement :"
            << " Checking initial mesh just to make sure" << endl;

        checkMesh();
        // Cannot call checkRefinementlevels since hanging points might
        // get triggered by the mesher after subsetting.
        //checkRefinementLevels(-1, labelList(0));
    }

    // Clear any saved point/cell data.
    savedPointLevel_.clear();
    savedCellLevel_.clear();


    // New point/cell level. Copy of pointLevel for existing points.
    DynamicList<label> newCellLevel(cellLevel_.size());
    forAll(cellLevel_, cellI)
    {
        newCellLevel.append(cellLevel_[cellI]);
    }
    DynamicList<label> newPointLevel(pointLevel_.size());
    forAll(pointLevel_, pointI)
    {
        newPointLevel.append(pointLevel_[pointI]);
    }


    if (debug)
    {
        Pout<< "hexRef8::setRefinement :"
            << " Allocating " << cellLabels.size() << " cell midpoints."
            << endl;
    }


    // Mid point per refined cell.
    // -1 : not refined
    // >=0: label of mid point.
    labelList cellMidPoint(mesh_.nCells(), -1);

    forAll(cellLabels, i)
    {
        label cellI = cellLabels[i];

        label anchorPointI = mesh_.faces()[mesh_.cells()[cellI][0]][0];

        cellMidPoint[cellI] = meshMod.setAction
        (
            polyAddPoint
            (
                mesh_.cellCentres()[cellI],     // point
                anchorPointI,                   // master point
                -1,                             // zone for point
                true                            // supports a cell
            )
        );

        newPointLevel(cellMidPoint[cellI]) = cellLevel_[cellI]+1;
    }


    if (debug)
    {
        cellSet splitCells(mesh_, "splitCells", cellLabels.size());

        forAll(cellMidPoint, cellI)
        {
            if (cellMidPoint[cellI] >= 0)
            {
                splitCells.insert(cellI);
            }
        }

        Pout<< "hexRef8::setRefinement : Dumping " << splitCells.size()
            << " cells to split to cellSet " << splitCells.objectPath()
            << endl;

        splitCells.write();
    }



    // Split edges
    // ~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef8::setRefinement :"
            << " Allocating edge midpoints."
            << endl;
    }

    // Unrefined edges are ones between cellLevel or lower points.
    // If any cell using this edge gets split then the edge needs to be split.

    // -1  : no need to split edge
    // >=0 : label of introduced mid point
    labelList edgeMidPoint(mesh_.nEdges(), -1);

    // Note: Loop over cells to be refined or edges?
    forAll(cellMidPoint, cellI)
    {
        if (cellMidPoint[cellI] >= 0)
        {
            const labelList& cEdges = mesh_.cellEdges(cellI);

            forAll(cEdges, i)
            {
                label edgeI = cEdges[i];

                const edge& e = mesh_.edges()[edgeI];

                if
                (
                    pointLevel_[e[0]] <= cellLevel_[cellI]
                 && pointLevel_[e[1]] <= cellLevel_[cellI]
                )
                {
                    edgeMidPoint[edgeI] = 12345;    // mark need for splitting
                }
            }
        }
    }

    // Synchronize edgeMidPoint across coupled patches. Take max so that
    // any split takes precedence.
    syncTools::syncEdgeList
    (
        mesh_,
        edgeMidPoint,
        maxEqOp<label>(),
        labelMin
    );


    // Introduce edge points
    // ~~~~~~~~~~~~~~~~~~~~~

    {
        // Phase 1: calculate midpoints and sync.
        // This needs doing for if people do not write binary and we slowly
        // get differences.

        pointField edgeMids(mesh_.nEdges(), point(-GREAT, -GREAT, -GREAT));

        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                // Edge marked to be split.
                edgeMids[edgeI] = mesh_.edges()[edgeI].centre(mesh_.points());
            }
        }
        syncTools::syncEdgePositions
        (
            mesh_,
            edgeMids,
            maxEqOp<vector>(),
            point(-GREAT, -GREAT, -GREAT)
        );


        // Phase 2: introduce points at the synced locations.
        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                // Edge marked to be split. Replace edgeMidPoint with actual
                // point label.

                const edge& e = mesh_.edges()[edgeI];

                edgeMidPoint[edgeI] = meshMod.setAction
                (
                    polyAddPoint
                    (
                        edgeMids[edgeI],            // point
                        e[0],                       // master point
                        -1,                         // zone for point
                        true                        // supports a cell
                    )
                );

                newPointLevel(edgeMidPoint[edgeI]) =
                    max
                    (
                        pointLevel_[e[0]],
                        pointLevel_[e[1]]
                    )
                  + 1;
            }
        }
    }

    if (debug)
    {
        OFstream str(mesh_.time().path()/"edgeMidPoint.obj");

        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                const edge& e = mesh_.edges()[edgeI];

                meshTools::writeOBJ(str, e.centre(mesh_.points()));
            }
        }

        Pout<< "hexRef8::setRefinement :"
            << " Dumping edge centres to split to file " << str.name() << endl;
    }


    // Calculate face level
    // ~~~~~~~~~~~~~~~~~~~~
    // (after splitting)

    if (debug)
    {
        Pout<< "hexRef8::setRefinement :"
            << " Allocating face midpoints."
            << endl;
    }

    // Face anchor level. There are guaranteed 4 points with level
    // <= anchorLevel. These are the corner points.
    labelList faceAnchorLevel(mesh_.nFaces());

    for (label faceI = 0; faceI < mesh_.nFaces(); faceI++)
    {
        faceAnchorLevel[faceI] = getAnchorLevel(faceI);
    }

    // -1  : no need to split face
    // >=0 : label of introduced mid point
    labelList faceMidPoint(mesh_.nFaces(), -1);


    // Internal faces: look at cells on both sides. Uniquely determined since
    // face itself guaranteed to be same level as most refined neighbour.
    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        if (faceAnchorLevel[faceI] >= 0)
        {
            label own = mesh_.faceOwner()[faceI];
            label ownLevel = cellLevel_[own];
            label newOwnLevel = ownLevel + (cellMidPoint[own] >= 0 ? 1 : 0);

            label nei = mesh_.faceNeighbour()[faceI];
            label neiLevel = cellLevel_[nei];
            label newNeiLevel = neiLevel + (cellMidPoint[nei] >= 0 ? 1 : 0);

            if
            (
                newOwnLevel > faceAnchorLevel[faceI]
             || newNeiLevel > faceAnchorLevel[faceI]
            )
            {
                faceMidPoint[faceI] = 12345;    // mark to be split
            }
        }
    }

    // Coupled patches handled like internal faces except now all information
    // from neighbour comes from across processor.
    // Boundary faces are more complicated since the boundary face can
    // be more refined than its owner (or neighbour for coupled patches)
    // (does not happen if refining/unrefining only, but does e.g. when
    //  refinining and subsetting)

    {
        labelList newNeiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(newNeiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];
            label ownLevel = cellLevel_[own];
            label newOwnLevel = ownLevel + (cellMidPoint[own] >= 0 ? 1 : 0);

            newNeiLevel[i] = newOwnLevel;
        }

        // Swap.
        syncTools::swapBoundaryFaceList(mesh_, newNeiLevel);

        // So now we have information on the neighbour.

        forAll(newNeiLevel, i)
        {
            label faceI = i+mesh_.nInternalFaces();

            if (faceAnchorLevel[faceI] >= 0)
            {
                label own = mesh_.faceOwner()[faceI];
                label ownLevel = cellLevel_[own];
                label newOwnLevel = ownLevel + (cellMidPoint[own] >= 0 ? 1 : 0);

                if
                (
                    newOwnLevel > faceAnchorLevel[faceI]
                 || newNeiLevel[i] > faceAnchorLevel[faceI]
                )
                {
                    faceMidPoint[faceI] = 12345;    // mark to be split
                }
            }
        }
    }


    // Synchronize faceMidPoint across coupled patches. (logical or)
    syncTools::syncFaceList
    (
        mesh_,
        faceMidPoint,
        maxEqOp<label>()
    );



    // Introduce face points
    // ~~~~~~~~~~~~~~~~~~~~~

    {
        // Phase 1: determine mid points and sync. See comment for edgeMids
        // above
        pointField bFaceMids
        (
            mesh_.nFaces()-mesh_.nInternalFaces(),
            point(-GREAT, -GREAT, -GREAT)
        );

        forAll(bFaceMids, i)
        {
            label faceI = i+mesh_.nInternalFaces();

            if (faceMidPoint[faceI] >= 0)
            {
                bFaceMids[i] = mesh_.faceCentres()[faceI];
            }
        }
        syncTools::syncBoundaryFacePositions
        (
            mesh_,
            bFaceMids,
            maxEqOp<vector>()
        );

        forAll(faceMidPoint, faceI)
        {
            if (faceMidPoint[faceI] >= 0)
            {
                // Face marked to be split. Replace faceMidPoint with actual
                // point label.

                const face& f = mesh_.faces()[faceI];

                faceMidPoint[faceI] = meshMod.setAction
                (
                    polyAddPoint
                    (
                        (
                            faceI < mesh_.nInternalFaces()
                          ? mesh_.faceCentres()[faceI]
                          : bFaceMids[faceI-mesh_.nInternalFaces()]
                        ),                          // point
                        f[0],                       // master point
                        -1,                         // zone for point
                        true                        // supports a cell
                    )
                );

                // Determine the level of the corner points and midpoint will
                // be one higher.
                newPointLevel(faceMidPoint[faceI]) = faceAnchorLevel[faceI]+1;
            }
        }
    }

    if (debug)
    {
        faceSet splitFaces(mesh_, "splitFaces", cellLabels.size());

        forAll(faceMidPoint, faceI)
        {
            if (faceMidPoint[faceI] >= 0)
            {
                splitFaces.insert(faceI);
            }
        }

        Pout<< "hexRef8::setRefinement : Dumping " << splitFaces.size()
            << " faces to split to faceSet " << splitFaces.objectPath() << endl;

        splitFaces.write();
    }


    // Information complete
    // ~~~~~~~~~~~~~~~~~~~~
    // At this point we have all the information we need. We should no
    // longer reference the cellLabels to refine. All the information is:
    // - cellMidPoint >= 0 : cell needs to be split
    // - faceMidPoint >= 0 : face needs to be split
    // - edgeMidPoint >= 0 : edge needs to be split



    // Get the corner/anchor points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef8::setRefinement :"
            << " Finding cell anchorPoints (8 per cell)"
            << endl;
    }

    // There will always be 8 points on the hex that have were introduced
    // with the hex and will have the same or lower refinement level.

    // Per cell the 8 corner points.
    labelListList cellAnchorPoints(mesh_.nCells());

    {
        labelList nAnchorPoints(mesh_.nCells(), 0);

        forAll(cellMidPoint, cellI)
        {
            if (cellMidPoint[cellI] >= 0)
            {
                cellAnchorPoints[cellI].setSize(8);
            }
        }

        forAll(pointLevel_, pointI)
        {
            const labelList& pCells = mesh_.pointCells(pointI);

            forAll(pCells, pCellI)
            {
                label cellI = pCells[pCellI];

                if
                (
                    cellMidPoint[cellI] >= 0
                 && pointLevel_[pointI] <= cellLevel_[cellI]
                )
                {
                    if (nAnchorPoints[cellI] == 8)
                    {
                        dumpCell(cellI);
                        FatalErrorIn
                        (
                            "hexRef8::setRefinement(const labelList&"
                            ", polyTopoChange&)"
                        )   << "cell " << cellI
                            << " of level " << cellLevel_[cellI]
                            << " uses more than 8 points of equal or"
                            << " lower level" << nl
                            << "Points so far:" << cellAnchorPoints[cellI]
                            << abort(FatalError);
                    }

                    cellAnchorPoints[cellI][nAnchorPoints[cellI]++] = pointI;
                }
            }
        }

        forAll(cellMidPoint, cellI)
        {
            if (cellMidPoint[cellI] >= 0)
            {
                if (nAnchorPoints[cellI] != 8)
                {
                    dumpCell(cellI);

                    const labelList& cPoints = mesh_.cellPoints(cellI);

                    FatalErrorIn
                    (
                        "hexRef8::setRefinement(const labelList&"
                        ", polyTopoChange&)"
                    )   << "cell " << cellI
                        << " of level " << cellLevel_[cellI]
                        << " does not seem to have 8 points of equal or"
                        << " lower level" << endl
                        << "cellPoints:" << cPoints << endl
                        << "pointLevels:"
                        << UIndirectList<label>(pointLevel_, cPoints)() << endl
                        << abort(FatalError);
                }
            }
        }
    }


    // Add the cells
    // ~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef8::setRefinement :"
            << " Adding cells (1 per anchorPoint)"
            << endl;
    }

    // Per cell the 7 added cells (+ original cell)
    labelListList cellAddedCells(mesh_.nCells());

    forAll(cellAnchorPoints, cellI)
    {
        const labelList& cAnchors = cellAnchorPoints[cellI];

        if (cAnchors.size() == 8)
        {
            labelList& cAdded = cellAddedCells[cellI];
            cAdded.setSize(8);

            // Original cell at 0
            cAdded[0] = cellI;

            // Update cell level
            newCellLevel[cellI] = cellLevel_[cellI]+1;


            for (label i = 1; i < 8; i++)
            {
                cAdded[i] = meshMod.setAction
                (
                    polyAddCell
                    (
                        -1,                                 // master point
                        -1,                                 // master edge
                        -1,                                 // master face
                        cellI,                              // master cell
                        mesh_.cellZones().whichZone(cellI)  // zone for cell
                    )
                );

                newCellLevel(cAdded[i]) = cellLevel_[cellI]+1;
            }
        }
    }


    // Faces
    // ~~~~~
    // 1. existing faces that get split (into four always)
    // 2. existing faces that do not get split but only edges get split
    // 3. existing faces that do not get split but get new owner/neighbour
    // 4. new internal faces inside split cells.

    if (debug)
    {
        Pout<< "hexRef8::setRefinement :"
            << " Marking faces to be handled"
            << endl;
    }

    // Get all affected faces.
    PackedBoolList affectedFace(mesh_.nFaces());

    {
        forAll(cellMidPoint, cellI)
        {
            if (cellMidPoint[cellI] >= 0)
            {
                const cell& cFaces = mesh_.cells()[cellI];

                forAll(cFaces, i)
                {
                    affectedFace.set(cFaces[i]);
                }
            }
        }

        forAll(faceMidPoint, faceI)
        {
            if (faceMidPoint[faceI] >= 0)
            {
                affectedFace.set(faceI);
            }
        }

        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                const labelList& eFaces = mesh_.edgeFaces(edgeI);

                forAll(eFaces, i)
                {
                    affectedFace.set(eFaces[i]);
                }
            }
        }
    }


    // 1. Faces that get split
    // ~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef8::setRefinement : Splitting faces" << endl;
    }

    forAll(faceMidPoint, faceI)
    {
        if (faceMidPoint[faceI] >= 0 && affectedFace.get(faceI))
        {
            // Face needs to be split and hasn't yet been done in some way
            // (affectedFace - is impossible since this is first change but
            //  just for completeness)

            const face& f = mesh_.faces()[faceI];

            // Has original faceI been used (three faces added, original gets
            // modified)
            bool modifiedFace = false;

            label anchorLevel = faceAnchorLevel[faceI];

            face newFace(4);

            forAll(f, fp)
            {
                label pointI = f[fp];

                if (pointLevel_[pointI] <= anchorLevel)
                {
                    // point is anchor. Start collecting face.

                    DynamicList<label> faceVerts(4);

                    faceVerts.append(pointI);

                    // Walk forward to mid point.
                    // - if next is +2 midpoint is +1
                    // - if next is +1 it is midpoint
                    // - if next is +0 there has to be edgeMidPoint

                    walkFaceToMid
                    (
                        edgeMidPoint,
                        anchorLevel,
                        faceI,
                        fp,
                        faceVerts
                    );

                    faceVerts.append(faceMidPoint[faceI]);

                    walkFaceFromMid
                    (
                        edgeMidPoint,
                        anchorLevel,
                        faceI,
                        fp,
                        faceVerts
                    );

                    // Convert dynamiclist to face.
                    newFace.transfer(faceVerts);

                    //Pout<< "Split face:" << faceI << " verts:" << f
                    //    << " into quad:" << newFace << endl;

                    // Get new owner/neighbour
                    label own, nei;
                    getFaceNeighbours
                    (
                        cellAnchorPoints,
                        cellAddedCells,
                        faceI,
                        pointI,          // Anchor point

                        own,
                        nei
                    );


                    if (debug)
                    {
                        if (mesh_.isInternalFace(faceI))
                        {
                            label oldOwn = mesh_.faceOwner()[faceI];
                            label oldNei = mesh_.faceNeighbour()[faceI];

                            checkInternalOrientation
                            (
                                meshMod,
                                oldOwn,
                                faceI,
                                mesh_.cellCentres()[oldOwn],
                                mesh_.cellCentres()[oldNei],
                                newFace
                            );
                        }
                        else
                        {
                            label oldOwn = mesh_.faceOwner()[faceI];

                            checkBoundaryOrientation
                            (
                                meshMod,
                                oldOwn,
                                faceI,
                                mesh_.cellCentres()[oldOwn],
                                mesh_.faceCentres()[faceI],
                                newFace
                            );
                        }
                    }


                    if (!modifiedFace)
                    {
                        modifiedFace = true;

                        modFace(meshMod, faceI, newFace, own, nei);
                    }
                    else
                    {
                        addFace(meshMod, faceI, newFace, own, nei);
                    }
                }
            }

            // Mark face as having been handled
            affectedFace.unset(faceI);
        }
    }


    // 2. faces that do not get split but use edges that get split
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef8::setRefinement :"
            << " Adding edge splits to unsplit faces"
            << endl;
    }

    DynamicList<label> eFacesStorage;
    DynamicList<label> fEdgesStorage;

    forAll(edgeMidPoint, edgeI)
    {
        if (edgeMidPoint[edgeI] >= 0)
        {
            // Split edge. Check that face not already handled above.

            const labelList& eFaces = mesh_.edgeFaces(edgeI, eFacesStorage);

            forAll(eFaces, i)
            {
                label faceI = eFaces[i];

                if (faceMidPoint[faceI] < 0 && affectedFace.get(faceI))
                {
                    // Unsplit face. Add edge splits to face.

                    const face& f = mesh_.faces()[faceI];
                    const labelList& fEdges = mesh_.faceEdges
                    (
                        faceI,
                        fEdgesStorage
                    );

                    DynamicList<label> newFaceVerts(f.size());

                    forAll(f, fp)
                    {
                        newFaceVerts.append(f[fp]);

                        label edgeI = fEdges[fp];

                        if (edgeMidPoint[edgeI] >= 0)
                        {
                            newFaceVerts.append(edgeMidPoint[edgeI]);
                        }
                    }

                    face newFace;
                    newFace.transfer(newFaceVerts);

                    // The point with the lowest level should be an anchor
                    // point of the neighbouring cells.
                    label anchorFp = findMinLevel(f);

                    label own, nei;
                    getFaceNeighbours
                    (
                        cellAnchorPoints,
                        cellAddedCells,
                        faceI,
                        f[anchorFp],          // Anchor point

                        own,
                        nei
                    );


                    if (debug)
                    {
                        if (mesh_.isInternalFace(faceI))
                        {
                            label oldOwn = mesh_.faceOwner()[faceI];
                            label oldNei = mesh_.faceNeighbour()[faceI];

                            checkInternalOrientation
                            (
                                meshMod,
                                oldOwn,
                                faceI,
                                mesh_.cellCentres()[oldOwn],
                                mesh_.cellCentres()[oldNei],
                                newFace
                            );
                        }
                        else
                        {
                            label oldOwn = mesh_.faceOwner()[faceI];

                            checkBoundaryOrientation
                            (
                                meshMod,
                                oldOwn,
                                faceI,
                                mesh_.cellCentres()[oldOwn],
                                mesh_.faceCentres()[faceI],
                                newFace
                            );
                        }
                    }

                    modFace(meshMod, faceI, newFace, own, nei);

                    // Mark face as having been handled
                    affectedFace.unset(faceI);
                }
            }
        }
    }


    // 3. faces that do not get split but whose owner/neighbour change
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout<< "hexRef8::setRefinement :"
            << " Changing owner/neighbour for otherwise unaffected faces"
            << endl;
    }

    forAll(affectedFace, faceI)
    {
        if (affectedFace.get(faceI))
        {
            const face& f = mesh_.faces()[faceI];

            // The point with the lowest level should be an anchor
            // point of the neighbouring cells.
            label anchorFp = findMinLevel(f);

            label own, nei;
            getFaceNeighbours
            (
                cellAnchorPoints,
                cellAddedCells,
                faceI,
                f[anchorFp],          // Anchor point

                own,
                nei
            );

            modFace(meshMod, faceI, f, own, nei);

            // Mark face as having been handled
            affectedFace.unset(faceI);
        }
    }


    // 4. new internal faces inside split cells.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    // This is the hard one. We have to find the splitting points between
    // the anchor points. But the edges between the anchor points might have
    // been split (into two,three or four edges).

    if (debug)
    {
        Pout<< "hexRef8::setRefinement :"
            << " Create new internal faces for split cells"
            << endl;
    }

    forAll(cellMidPoint, cellI)
    {
        if (cellMidPoint[cellI] >= 0)
        {
            createInternalFaces
            (
                cellAnchorPoints,
                cellAddedCells,
                cellMidPoint,
                faceMidPoint,
                faceAnchorLevel,
                edgeMidPoint,
                cellI,
                meshMod
            );
        }
    }

    // Extend pointLevels and cellLevels for the new cells. Could also be done
    // in updateMesh but saves passing cellAddedCells out of this routine.

    // Check
    if (debug)
    {
        label minPointI = labelMax;
        label maxPointI = labelMin;

        forAll(cellMidPoint, cellI)
        {
            if (cellMidPoint[cellI] >= 0)
            {
                minPointI = min(minPointI, cellMidPoint[cellI]);
                maxPointI = max(maxPointI, cellMidPoint[cellI]);
            }
        }
        forAll(faceMidPoint, faceI)
        {
            if (faceMidPoint[faceI] >= 0)
            {
                minPointI = min(minPointI, faceMidPoint[faceI]);
                maxPointI = max(maxPointI, faceMidPoint[faceI]);
            }
        }
        forAll(edgeMidPoint, edgeI)
        {
            if (edgeMidPoint[edgeI] >= 0)
            {
                minPointI = min(minPointI, edgeMidPoint[edgeI]);
                maxPointI = max(maxPointI, edgeMidPoint[edgeI]);
            }
        }

        if (minPointI != labelMax && minPointI != mesh_.nPoints())
        {
            FatalErrorIn("hexRef8::setRefinement(..)")
                << "Added point labels not consecutive to existing mesh points."
                << nl
                << "mesh_.nPoints():" << mesh_.nPoints()
                << " minPointI:" << minPointI
                << " maxPointI:" << maxPointI
                << abort(FatalError);
        }
    }

    pointLevel_.transfer(newPointLevel);
    cellLevel_.transfer(newCellLevel);

    // Mark files as changed
    setInstance(mesh_.facesInstance());


    // Update the live split cells tree.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // New unrefinement structure
    if (history_.active())
    {
        if (debug)
        {
            Pout<< "hexRef8::setRefinement :"
                << " Updating refinement history to " << cellLevel_.size()
                << " cells" << endl;
        }

        // Extend refinement history for new cells
        history_.resize(cellLevel_.size());

        forAll(cellAddedCells, cellI)
        {
            const labelList& addedCells = cellAddedCells[cellI];

            if (addedCells.size())
            {
                // Cell was split.
                history_.storeSplit(cellI, addedCells);
            }
        }
    }

    // Compact cellAddedCells.

    labelListList refinedCells(cellLabels.size());

    forAll(cellLabels, i)
    {
        label cellI = cellLabels[i];

        refinedCells[i].transfer(cellAddedCells[cellI]);
    }

    return refinedCells;
}


void Foam::hexRef8::storeData
(
    const labelList& pointsToStore,
    const labelList& facesToStore,
    const labelList& cellsToStore
)
{
    savedPointLevel_.resize(2*pointsToStore.size());
    forAll(pointsToStore, i)
    {
        label pointI = pointsToStore[i];
        savedPointLevel_.insert(pointI, pointLevel_[pointI]);
    }

    savedCellLevel_.resize(2*cellsToStore.size());
    forAll(cellsToStore, i)
    {
        label cellI = cellsToStore[i];
        savedCellLevel_.insert(cellI, cellLevel_[cellI]);
    }
}


// Gets called after the mesh change. setRefinement will already have made
// sure the pointLevel_ and cellLevel_ are the size of the new mesh so we
// only need to account for reordering.
void Foam::hexRef8::updateMesh(const mapPolyMesh& map)
{
    Map<label> dummyMap(0);

    updateMesh(map, dummyMap, dummyMap, dummyMap);
}


// Gets called after the mesh change. setRefinement will already have made
// sure the pointLevel_ and cellLevel_ are the size of the new mesh so we
// only need to account for reordering.
void Foam::hexRef8::updateMesh
(
    const mapPolyMesh& map,
    const Map<label>& pointsToRestore,
    const Map<label>& facesToRestore,
    const Map<label>& cellsToRestore
)
{
    // Update celllevel
    if (debug)
    {
        Pout<< "hexRef8::updateMesh :"
            << " Updating various lists"
            << endl;
    }

    {
        const labelList& reverseCellMap = map.reverseCellMap();

        if (debug)
        {
            Pout<< "hexRef8::updateMesh :"
                << " reverseCellMap:" << map.reverseCellMap().size()
                << " cellMap:" << map.cellMap().size()
                << " nCells:" << mesh_.nCells()
                << " nOldCells:" << map.nOldCells()
                << " cellLevel_:" << cellLevel_.size()
                << " reversePointMap:" << map.reversePointMap().size()
                << " pointMap:" << map.pointMap().size()
                << " nPoints:" << mesh_.nPoints()
                << " nOldPoints:" << map.nOldPoints()
                << " pointLevel_:" << pointLevel_.size()
                << endl;
        }

        if (reverseCellMap.size() == cellLevel_.size())
        {
            // Assume it is after hexRef8 that this routine is called.
            // Just account for reordering. We cannot use cellMap since
            // then cells created from cells would get cellLevel_ of
            // cell they were created from.
            reorder(reverseCellMap, mesh_.nCells(), -1, cellLevel_);
        }
        else
        {
            // Map data
            const labelList& cellMap = map.cellMap();

            labelList newCellLevel(cellMap.size());
            forAll(cellMap, newCellI)
            {
                label oldCellI = cellMap[newCellI];

                if (oldCellI == -1)
                {
                    newCellLevel[newCellI] = -1;
                }
                else
                {
                    newCellLevel[newCellI] = cellLevel_[oldCellI];
                }
            }
            cellLevel_.transfer(newCellLevel);
        }

        // See if any cells to restore. This will be for some new cells
        // the corresponding old cell.
        forAllConstIter(Map<label>, cellsToRestore, iter)
        {
            label newCellI = iter.key();
            label storedCellI = iter();

            Map<label>::iterator fnd = savedCellLevel_.find(storedCellI);

            if (fnd == savedCellLevel_.end())
            {
                FatalErrorIn("hexRef8::updateMesh(const mapPolyMesh&)")
                    << "Problem : trying to restore old value for new cell "
                    << newCellI << " but cannot find old cell " << storedCellI
                    << " in map of stored values " << savedCellLevel_
                    << abort(FatalError);
            }
            cellLevel_[newCellI] = fnd();
        }

        //if (findIndex(cellLevel_, -1) != -1)
        //{
        //    WarningIn("hexRef8::updateMesh(const mapPolyMesh&)")
        //        << "Problem : "
        //        << "cellLevel_ contains illegal value -1 after mapping
        //        << " at cell " << findIndex(cellLevel_, -1) << endl
        //        << "This means that another program has inflated cells"
        //        << " (created cells out-of-nothing) and hence we don't know"
        //        << " their cell level. Continuing with illegal value."
        //        << abort(FatalError);
        //}
    }


    // Update pointlevel
    {
        const labelList& reversePointMap = map.reversePointMap();

        if (reversePointMap.size() == pointLevel_.size())
        {
            // Assume it is after hexRef8 that this routine is called.
            reorder(reversePointMap, mesh_.nPoints(), -1,  pointLevel_);
        }
        else
        {
            // Map data
            const labelList& pointMap = map.pointMap();

            labelList newPointLevel(pointMap.size());

            forAll(pointMap, newPointI)
            {
                label oldPointI = pointMap[newPointI];

                if (oldPointI == -1)
                {
                    //FatalErrorIn("hexRef8::updateMesh(const mapPolyMesh&)")
                    //    << "Problem : point " << newPointI
                    //    << " at " << mesh_.points()[newPointI]
                    //    << " does not originate from another point"
                    //    << " (i.e. is inflated)." << nl
                    //    << "Hence we cannot determine the new pointLevel"
                    //    << " for it." << abort(FatalError);
                    newPointLevel[newPointI] = -1;
                }
                else
                {
                    newPointLevel[newPointI] = pointLevel_[oldPointI];
                }
            }
            pointLevel_.transfer(newPointLevel);
        }

        // See if any points to restore. This will be for some new points
        // the corresponding old point (the one from the call to storeData)
        forAllConstIter(Map<label>, pointsToRestore, iter)
        {
            label newPointI = iter.key();
            label storedPointI = iter();

            Map<label>::iterator fnd = savedPointLevel_.find(storedPointI);

            if (fnd == savedPointLevel_.end())
            {
                FatalErrorIn("hexRef8::updateMesh(const mapPolyMesh&)")
                    << "Problem : trying to restore old value for new point "
                    << newPointI << " but cannot find old point "
                    << storedPointI
                    << " in map of stored values " << savedPointLevel_
                    << abort(FatalError);
            }
            pointLevel_[newPointI] = fnd();
        }

        //if (findIndex(pointLevel_, -1) != -1)
        //{
        //    WarningIn("hexRef8::updateMesh(const mapPolyMesh&)")
        //        << "Problem : "
        //        << "pointLevel_ contains illegal value -1 after mapping"
        //        << " at point" << findIndex(pointLevel_, -1) << endl
        //        << "This means that another program has inflated points"
        //        << " (created points out-of-nothing) and hence we don't know"
        //        << " their point level. Continuing with illegal value."
        //        //<< abort(FatalError);
        //}
    }

    // Update refinement tree
    if (history_.active())
    {
        history_.updateMesh(map);
    }

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // Update face removal engine
    faceRemover_.updateMesh(map);

    // Clear cell shapes
    cellShapesPtr_.clear();
}


// Gets called after mesh subsetting. Maps are from new back to old.
void Foam::hexRef8::subset
(
    const labelList& pointMap,
    const labelList& faceMap,
    const labelList& cellMap
)
{
    // Update celllevel
    if (debug)
    {
        Pout<< "hexRef8::subset :"
            << " Updating various lists"
            << endl;
    }

    if (history_.active())
    {
        WarningIn
        (
            "hexRef8::subset(const labelList&, const labelList&"
            ", const labelList&)"
        )   << "Subsetting will not work in combination with unrefinement."
            << nl
            << "Proceed at your own risk." << endl;
    }


    // Update celllevel
    {
        labelList newCellLevel(cellMap.size());

        forAll(cellMap, newCellI)
        {
            newCellLevel[newCellI] = cellLevel_[cellMap[newCellI]];
        }

        cellLevel_.transfer(newCellLevel);

        if (findIndex(cellLevel_, -1) != -1)
        {
            FatalErrorIn("hexRef8::subset(..)")
                << "Problem : "
                << "cellLevel_ contains illegal value -1 after mapping:"
                << cellLevel_
                << abort(FatalError);
        }
    }

    // Update pointlevel
    {
        labelList newPointLevel(pointMap.size());

        forAll(pointMap, newPointI)
        {
            newPointLevel[newPointI] = pointLevel_[pointMap[newPointI]];
        }

        pointLevel_.transfer(newPointLevel);

        if (findIndex(pointLevel_, -1) != -1)
        {
            FatalErrorIn("hexRef8::subset(..)")
                << "Problem : "
                << "pointLevel_ contains illegal value -1 after mapping:"
                << pointLevel_
                << abort(FatalError);
        }
    }

    // Update refinement tree
    if (history_.active())
    {
        history_.subset(pointMap, faceMap, cellMap);
    }

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // Nothing needs doing to faceRemover.
    //faceRemover_.subset(pointMap, faceMap, cellMap);

    // Clear cell shapes
    cellShapesPtr_.clear();
}


// Gets called after the mesh distribution
void Foam::hexRef8::distribute(const mapDistributePolyMesh& map)
{
    if (debug)
    {
        Pout<< "hexRef8::distribute :"
            << " Updating various lists"
            << endl;
    }

    // Update celllevel
    map.distributeCellData(cellLevel_);
    // Update pointlevel
    map.distributePointData(pointLevel_);

    // Update refinement tree
    if (history_.active())
    {
        history_.distribute(map);
    }

    // Update face removal engine
    faceRemover_.distribute(map);

    // Clear cell shapes
    cellShapesPtr_.clear();
}


void Foam::hexRef8::checkMesh() const
{
    const scalar smallDim = 1e-6 * mesh_.bounds().mag();

    if (debug)
    {
        Pout<< "hexRef8::checkMesh : Using matching tolerance smallDim:"
            << smallDim << endl;
    }

    // Check owner on coupled faces.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // There should be only one coupled face between two cells. Why? Since
    // otherwise mesh redistribution might cause multiple faces between two
    // cells
    {
        labelList nei(mesh_.nFaces()-mesh_.nInternalFaces());
        forAll(nei, i)
        {
            nei[i] = mesh_.faceOwner()[i+mesh_.nInternalFaces()];
        }

        // Replace data on coupled patches with their neighbour ones.
        syncTools::swapBoundaryFaceList(mesh_, nei);

        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (pp.coupled())
            {
                // Check how many faces between owner and neighbour. Should
                // be only one.
                HashTable<label, labelPair, labelPair::Hash<> >
                    cellToFace(2*pp.size());

                label faceI = pp.start();

                forAll(pp, i)
                {
                    label own = mesh_.faceOwner()[faceI];
                    label bFaceI = faceI-mesh_.nInternalFaces();

                    if (!cellToFace.insert(labelPair(own, nei[bFaceI]), faceI))
                    {
                        dumpCell(own);
                        FatalErrorIn("hexRef8::checkMesh()")
                            << "Faces do not seem to be correct across coupled"
                            << " boundaries" << endl
                            << "Coupled face " << faceI
                            << " between owner " << own
                            << " on patch " << pp.name()
                            << " and coupled neighbour " << nei[bFaceI]
                            << " has two faces connected to it:"
                            << faceI << " and "
                            << cellToFace[labelPair(own, nei[bFaceI])]
                            << abort(FatalError);
                    }

                    faceI++;
                }
            }
        }
    }

    // Check face areas.
    // ~~~~~~~~~~~~~~~~~

    {
        scalarField neiFaceAreas(mesh_.nFaces()-mesh_.nInternalFaces());
        forAll(neiFaceAreas, i)
        {
            neiFaceAreas[i] = mag(mesh_.faceAreas()[i+mesh_.nInternalFaces()]);
        }

        // Replace data on coupled patches with their neighbour ones.
        syncTools::swapBoundaryFaceList(mesh_, neiFaceAreas);

        forAll(neiFaceAreas, i)
        {
            label faceI = i+mesh_.nInternalFaces();

            const scalar magArea = mag(mesh_.faceAreas()[faceI]);

            if (mag(magArea - neiFaceAreas[i]) > smallDim)
            {
                const face& f = mesh_.faces()[faceI];
                label patchI = mesh_.boundaryMesh().whichPatch(faceI);

                dumpCell(mesh_.faceOwner()[faceI]);

                FatalErrorIn("hexRef8::checkMesh()")
                    << "Faces do not seem to be correct across coupled"
                    << " boundaries" << endl
                    << "Coupled face " << faceI
                    << " on patch " << patchI
                    << " " << mesh_.boundaryMesh()[patchI].name()
                    << " coords:" << UIndirectList<point>(mesh_.points(), f)()
                    << " has face area:" << magArea
                    << " (coupled) neighbour face area differs:"
                    << neiFaceAreas[i]
                    << " to within tolerance " << smallDim
                    << abort(FatalError);
            }
        }
    }


    // Check number of points on faces.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    {
        labelList nVerts(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(nVerts, i)
        {
            nVerts[i] = mesh_.faces()[i+mesh_.nInternalFaces()].size();
        }

        // Replace data on coupled patches with their neighbour ones.
        syncTools::swapBoundaryFaceList(mesh_, nVerts);

        forAll(nVerts, i)
        {
            label faceI = i+mesh_.nInternalFaces();

            const face& f = mesh_.faces()[faceI];

            if (f.size() != nVerts[i])
            {
                dumpCell(mesh_.faceOwner()[faceI]);

                label patchI = mesh_.boundaryMesh().whichPatch(faceI);

                FatalErrorIn("hexRef8::checkMesh()")
                    << "Faces do not seem to be correct across coupled"
                    << " boundaries" << endl
                    << "Coupled face " << faceI
                    << " on patch " << patchI
                    << " " << mesh_.boundaryMesh()[patchI].name()
                    << " coords:" << UIndirectList<point>(mesh_.points(), f)()
                    << " has size:" << f.size()
                    << " (coupled) neighbour face has size:"
                    << nVerts[i]
                    << abort(FatalError);
            }
        }
    }


    // Check points of face
    // ~~~~~~~~~~~~~~~~~~~~
    {
        // Anchor points.
        pointField anchorPoints(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(anchorPoints, i)
        {
            label faceI = i+mesh_.nInternalFaces();
            const point& fc = mesh_.faceCentres()[faceI];
            const face& f = mesh_.faces()[faceI];
            const vector anchorVec(mesh_.points()[f[0]] - fc);

            anchorPoints[i] = anchorVec;
        }

        // Replace data on coupled patches with their neighbour ones. Apply
        // rotation transformation (but not separation since is relative vector
        // to point on same face.
        syncTools::swapBoundaryFaceList(mesh_, anchorPoints);

        forAll(anchorPoints, i)
        {
            label faceI = i+mesh_.nInternalFaces();
            const point& fc = mesh_.faceCentres()[faceI];
            const face& f = mesh_.faces()[faceI];
            const vector anchorVec(mesh_.points()[f[0]] - fc);

            if (mag(anchorVec - anchorPoints[i]) > smallDim)
            {
                dumpCell(mesh_.faceOwner()[faceI]);

                label patchI = mesh_.boundaryMesh().whichPatch(faceI);

                FatalErrorIn("hexRef8::checkMesh()")
                    << "Faces do not seem to be correct across coupled"
                    << " boundaries" << endl
                    << "Coupled face " << faceI
                    << " on patch " << patchI
                    << " " << mesh_.boundaryMesh()[patchI].name()
                    << " coords:" << UIndirectList<point>(mesh_.points(), f)()
                    << " has anchor vector:" << anchorVec
                    << " (coupled) neighbour face anchor vector differs:"
                    << anchorPoints[i]
                    << " to within tolerance " << smallDim
                    << abort(FatalError);
            }
        }
    }

    if (debug)
    {
        Pout<< "hexRef8::checkMesh : Returning" << endl;
    }
}


void Foam::hexRef8::checkRefinementLevels
(
    const label maxPointDiff,
    const labelList& pointsToCheck
) const
{
    if (debug)
    {
        Pout<< "hexRef8::checkRefinementLevels :"
            << " Checking 2:1 refinement level" << endl;
    }

    if
    (
        cellLevel_.size() != mesh_.nCells()
     || pointLevel_.size() != mesh_.nPoints()
    )
    {
        FatalErrorIn("hexRef8::checkRefinementLevels(const label)")
            << "cellLevel size should be number of cells"
            << " and pointLevel size should be number of points."<< nl
            << "cellLevel:" << cellLevel_.size()
            << " mesh.nCells():" << mesh_.nCells() << nl
            << "pointLevel:" << pointLevel_.size()
            << " mesh.nPoints():" << mesh_.nPoints()
            << abort(FatalError);
    }


    // Check 2:1 consistency.
    // ~~~~~~~~~~~~~~~~~~~~~~

    {
        // Internal faces.
        for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
        {
            label own = mesh_.faceOwner()[faceI];
            label nei = mesh_.faceNeighbour()[faceI];

            if (mag(cellLevel_[own] - cellLevel_[nei]) > 1)
            {
                dumpCell(own);
                dumpCell(nei);

                FatalErrorIn
                (
                    "hexRef8::checkRefinementLevels(const label)"
                )   << "Celllevel does not satisfy 2:1 constraint." << nl
                    << "On face " << faceI << " owner cell " << own
                    << " has refinement " << cellLevel_[own]
                    << " neighbour cell " << nei << " has refinement "
                    << cellLevel_[nei]
                    << abort(FatalError);
            }
        }

        // Coupled faces. Get neighbouring value
        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(neiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];

            neiLevel[i] = cellLevel_[own];
        }

        // No separation
        syncTools::swapBoundaryFaceList(mesh_, neiLevel);

        forAll(neiLevel, i)
        {
            label faceI = i+mesh_.nInternalFaces();

            label own = mesh_.faceOwner()[faceI];

            if (mag(cellLevel_[own] - neiLevel[i]) > 1)
            {
                dumpCell(own);

                label patchI = mesh_.boundaryMesh().whichPatch(faceI);

                FatalErrorIn
                (
                    "hexRef8::checkRefinementLevels(const label)"
                )   << "Celllevel does not satisfy 2:1 constraint."
                    << " On coupled face " << faceI
                    << " on patch " << patchI << " "
                    << mesh_.boundaryMesh()[patchI].name()
                    << " owner cell " << own << " has refinement "
                    << cellLevel_[own]
                    << " (coupled) neighbour cell has refinement "
                    << neiLevel[i]
                    << abort(FatalError);
            }
        }
    }


    // Check pointLevel is synchronized
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    {
        labelList syncPointLevel(pointLevel_);

        // Get min level
        syncTools::syncPointList
        (
            mesh_,
            syncPointLevel,
            minEqOp<label>(),
            labelMax
        );


        forAll(syncPointLevel, pointI)
        {
            if (pointLevel_[pointI] != syncPointLevel[pointI])
            {
                FatalErrorIn
                (
                    "hexRef8::checkRefinementLevels(const label)"
                )   << "PointLevel is not consistent across coupled patches."
                    << endl
                    << "point:" << pointI << " coord:" << mesh_.points()[pointI]
                    << " has level " << pointLevel_[pointI]
                    << " whereas the coupled point has level "
                    << syncPointLevel[pointI]
                    << abort(FatalError);
            }
        }
    }


    // Check 2:1 across points (instead of faces)
    if (maxPointDiff != -1)
    {
        // Determine per point the max cell level.
        labelList maxPointLevel(mesh_.nPoints(), 0);

        forAll(maxPointLevel, pointI)
        {
            const labelList& pCells = mesh_.pointCells(pointI);

            label& pLevel = maxPointLevel[pointI];

            forAll(pCells, i)
            {
                pLevel = max(pLevel, cellLevel_[pCells[i]]);
            }
        }

        // Sync maxPointLevel to neighbour
        syncTools::syncPointList
        (
            mesh_,
            maxPointLevel,
            maxEqOp<label>(),
            labelMin            // null value
        );

        // Check 2:1 across boundary points
        forAll(pointsToCheck, i)
        {
            label pointI = pointsToCheck[i];

            const labelList& pCells = mesh_.pointCells(pointI);

            forAll(pCells, i)
            {
                label cellI = pCells[i];

                if
                (
                    mag(cellLevel_[cellI]-maxPointLevel[pointI])
                  > maxPointDiff
                )
                {
                    dumpCell(cellI);

                    FatalErrorIn
                    (
                        "hexRef8::checkRefinementLevels(const label)"
                    )   << "Too big a difference between"
                        << " point-connected cells." << nl
                        << "cell:" << cellI
                        << " cellLevel:" << cellLevel_[cellI]
                        << " uses point:" << pointI
                        << " coord:" << mesh_.points()[pointI]
                        << " which is also used by a cell with level:"
                        << maxPointLevel[pointI]
                        << abort(FatalError);
                }
            }
        }
    }


    //- Gives problems after first splitting off inside mesher.
    //// Hanging points
    //{
    //    // Any patches with points having only two edges.
    //
    //    boolList isHangingPoint(mesh_.nPoints(), false);
    //
    //    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    //
    //    forAll(patches, patchI)
    //    {
    //        const polyPatch& pp = patches[patchI];
    //
    //        const labelList& meshPoints = pp.meshPoints();
    //
    //        forAll(meshPoints, i)
    //        {
    //            label pointI = meshPoints[i];
    //
    //            const labelList& pEdges = mesh_.pointEdges()[pointI];
    //
    //            if (pEdges.size() == 2)
    //            {
    //                isHangingPoint[pointI] = true;
    //            }
    //        }
    //    }
    //
    //    syncTools::syncPointList
    //    (
    //        mesh_,
    //        isHangingPoint,
    //        andEqOp<bool>(),        // only if all decide it is hanging point
    //        true,                   // null
    //        false                   // no separation
    //    );
    //
    //    //OFstream str(mesh_.time().path()/"hangingPoints.obj");
    //
    //    label nHanging = 0;
    //
    //    forAll(isHangingPoint, pointI)
    //    {
    //        if (isHangingPoint[pointI])
    //        {
    //            nHanging++;
    //
    //            Pout<< "Hanging boundary point " << pointI
    //                << " at " << mesh_.points()[pointI]
    //                << endl;
    //            //meshTools::writeOBJ(str, mesh_.points()[pointI]);
    //        }
    //    }
    //
    //    if (returnReduce(nHanging, sumOp<label>()) > 0)
    //    {
    //        FatalErrorIn
    //        (
    //            "hexRef8::checkRefinementLevels(const label)"
    //        )   << "Detected a point used by two edges only (hanging point)"
    //            << nl << "This is not allowed"
    //            << abort(FatalError);
    //    }
    //}
}


const Foam::cellShapeList& Foam::hexRef8::cellShapes() const
{
    if (cellShapesPtr_.empty())
    {
        if (debug)
        {
            Pout<< "hexRef8::cellShapes() : calculating splitHex cellShapes."
                << " cellLevel:" << cellLevel_.size()
                << " pointLevel:" << pointLevel_.size()
                << endl;
        }

        const cellShapeList& meshShapes = mesh_.cellShapes();
        cellShapesPtr_.reset(new cellShapeList(meshShapes));

        label nSplitHex = 0;
        label nUnrecognised = 0;

        forAll(cellLevel_, cellI)
        {
            if (meshShapes[cellI].model().index() == 0)
            {
                label level = cellLevel_[cellI];

                // Return true if we've found 6 quads
                DynamicList<face> quads;
                bool haveQuads = matchHexShape
                (
                    cellI,
                    level,
                    quads
                );

                if (haveQuads)
                {
                    faceList faces(quads.xfer());
                    cellShapesPtr_()[cellI] = degenerateMatcher::match(faces);
                    nSplitHex++;
                }
                else
                {
                    nUnrecognised++;
                }
            }
        }
        if (debug)
        {
            Pout<< "hexRef8::cellShapes() :"
                << " nCells:" << mesh_.nCells() << " of which" << nl
                << "    primitive:" << (mesh_.nCells()-nSplitHex-nUnrecognised)
                << nl
                << "    split-hex:" << nSplitHex << nl
                << "    poly     :" << nUnrecognised << nl
                << endl;
        }
    }
    return cellShapesPtr_();
}



//
// Unrefinement
// ~~~~~~~~~~~~
//


Foam::labelList Foam::hexRef8::getSplitPoints() const
{
    if (debug)
    {
        checkRefinementLevels(-1, labelList(0));
    }

    if (debug)
    {
        Pout<< "hexRef8::getSplitPoints :"
            << " Calculating unrefineable points" << endl;
    }


    if (!history_.active())
    {
        FatalErrorIn("hexRef8::getSplitPoints()")
            << "Only call if constructed with history capability"
            << abort(FatalError);
    }

    // Master cell
    // -1 undetermined
    // -2 certainly not split point
    // >= label of master cell
    labelList splitMaster(mesh_.nPoints(), -1);
    labelList splitMasterLevel(mesh_.nPoints(), 0);

    // Unmark all with not 8 cells
    //const labelListList& pointCells = mesh_.pointCells();

    for (label pointI = 0; pointI < mesh_.nPoints(); pointI++)
    {
        const labelList& pCells = mesh_.pointCells(pointI);

        if (pCells.size() != 8)
        {
            splitMaster[pointI] = -2;
        }
    }

    // Unmark all with different master cells
    const labelList& visibleCells = history_.visibleCells();

    forAll(visibleCells, cellI)
    {
        const labelList& cPoints = mesh_.cellPoints(cellI);

        if (visibleCells[cellI] != -1 && history_.parentIndex(cellI) >= 0)
        {
            label parentIndex = history_.parentIndex(cellI);

            // Check same master.
            forAll(cPoints, i)
            {
                label pointI = cPoints[i];

                label masterCellI = splitMaster[pointI];

                if (masterCellI == -1)
                {
                    // First time visit of point. Store parent cell and
                    // level of the parent cell (with respect to cellI). This
                    // is additional guarantee that we're referring to the
                    // same master at the same refinement level.

                    splitMaster[pointI] = parentIndex;
                    splitMasterLevel[pointI] = cellLevel_[cellI] - 1;
                }
                else if (masterCellI == -2)
                {
                    // Already decided that point is not splitPoint
                }
                else if
                (
                    (masterCellI != parentIndex)
                 || (splitMasterLevel[pointI] != cellLevel_[cellI] - 1)
                )
                {
                    // Different masters so point is on two refinement
                    // patterns
                    splitMaster[pointI] = -2;
                }
            }
        }
        else
        {
            // Either not visible or is unrefined cell
            forAll(cPoints, i)
            {
                label pointI = cPoints[i];

                splitMaster[pointI] = -2;
            }
        }
    }

    // Unmark boundary faces
    for
    (
        label faceI = mesh_.nInternalFaces();
        faceI < mesh_.nFaces();
        faceI++
    )
    {
        const face& f = mesh_.faces()[faceI];

        forAll(f, fp)
        {
            splitMaster[f[fp]] = -2;
        }
    }


    // Collect into labelList

    label nSplitPoints = 0;

    forAll(splitMaster, pointI)
    {
        if (splitMaster[pointI] >= 0)
        {
            nSplitPoints++;
        }
    }

    labelList splitPoints(nSplitPoints);
    nSplitPoints = 0;

    forAll(splitMaster, pointI)
    {
        if (splitMaster[pointI] >= 0)
        {
            splitPoints[nSplitPoints++] = pointI;
        }
    }

    return splitPoints;
}


//void Foam::hexRef8::markIndex
//(
//    const label maxLevel,
//    const label level,
//    const label index,
//    const label markValue,
//    labelList& indexValues
//) const
//{
//    if (level < maxLevel && indexValues[index] == -1)
//    {
//        // Mark
//        indexValues[index] = markValue;
//
//        // Mark parent
//        const splitCell8& split = history_.splitCells()[index];
//
//        if (split.parent_ >= 0)
//        {
//            markIndex
//            (
//              maxLevel, level+1, split.parent_, markValue, indexValues);
//            )
//        }
//    }
//}
//
//
//// Get all cells which (down to level) originate from the same cell.
//// level=0 returns cell only, level=1 returns the 8 cells this cell
//// originates from, level=2 returns 64 cells etc.
//// If the cell does not originate from refinement returns just itself.
//void Foam::hexRef8::markCellClusters
//(
//    const label maxLevel,
//    labelList& cluster
//) const
//{
//    cluster.setSize(mesh_.nCells());
//    cluster = -1;
//
//    const DynamicList<splitCell8>& splitCells = history_.splitCells();
//
//    // Mark all splitCells down to level maxLevel with a cell originating from
//    // it.
//
//    labelList indexLevel(splitCells.size(), -1);
//
//    forAll(visibleCells, cellI)
//    {
//        label index = visibleCells[cellI];
//
//        if (index >= 0)
//        {
//            markIndex(maxLevel, 0, index, cellI, indexLevel);
//        }
//    }
//
//    // Mark cells with splitCell
//}


Foam::labelList Foam::hexRef8::consistentUnrefinement
(
    const labelList& pointsToUnrefine,
    const bool maxSet
) const
{
    if (debug)
    {
        Pout<< "hexRef8::consistentUnrefinement :"
            << " Determining 2:1 consistent unrefinement" << endl;
    }

    if (maxSet)
    {
        FatalErrorIn
        (
            "hexRef8::consistentUnrefinement(const labelList&, const bool"
        )   << "maxSet not implemented yet."
            << abort(FatalError);
    }

    // Loop, modifying pointsToUnrefine, until no more changes to due to 2:1
    // conflicts.
    // maxSet = false : unselect points to refine
    // maxSet = true: select points to refine

    // Maintain boolList for pointsToUnrefine and cellsToUnrefine
    PackedBoolList unrefinePoint(mesh_.nPoints());

    forAll(pointsToUnrefine, i)
    {
        label pointI = pointsToUnrefine[i];

        unrefinePoint.set(pointI);
    }


    while (true)
    {
        // Construct cells to unrefine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        PackedBoolList unrefineCell(mesh_.nCells());

        forAll(unrefinePoint, pointI)
        {
            if (unrefinePoint.get(pointI))
            {
                const labelList& pCells = mesh_.pointCells(pointI);

                forAll(pCells, j)
                {
                    unrefineCell.set(pCells[j]);
                }
            }
        }


        label nChanged = 0;


        // Check 2:1 consistency taking refinement into account
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Internal faces.
        for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
        {
            label own = mesh_.faceOwner()[faceI];
            label ownLevel = cellLevel_[own] - unrefineCell.get(own);

            label nei = mesh_.faceNeighbour()[faceI];
            label neiLevel = cellLevel_[nei] - unrefineCell.get(nei);

            if (ownLevel < (neiLevel-1))
            {
                // Since was 2:1 this can only occur if own is marked for
                // unrefinement.

                if (maxSet)
                {
                    unrefineCell.set(nei);
                }
                else
                {
                    // could also combine with unset:
                    // if (!unrefineCell.unset(own))
                    // {
                    //     FatalErrorIn("hexRef8::consistentUnrefinement(..)")
                    //         << "problem cell already unset"
                    //         << abort(FatalError);
                    // }
                    if (unrefineCell.get(own) == 0)
                    {
                        FatalErrorIn("hexRef8::consistentUnrefinement(..)")
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.unset(own);
                }
                nChanged++;
            }
            else if (neiLevel < (ownLevel-1))
            {
                if (maxSet)
                {
                    unrefineCell.set(own);
                }
                else
                {
                    if (unrefineCell.get(nei) == 0)
                    {
                        FatalErrorIn("hexRef8::consistentUnrefinement(..)")
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.unset(nei);
                }
                nChanged++;
            }
        }


        // Coupled faces. Swap owner level to get neighbouring cell level.
        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

        forAll(neiLevel, i)
        {
            label own = mesh_.faceOwner()[i+mesh_.nInternalFaces()];

            neiLevel[i] = cellLevel_[own] - unrefineCell.get(own);
        }

        // Swap to neighbour
        syncTools::swapBoundaryFaceList(mesh_, neiLevel);

        forAll(neiLevel, i)
        {
            label faceI = i+mesh_.nInternalFaces();
            label own = mesh_.faceOwner()[faceI];
            label ownLevel = cellLevel_[own] - unrefineCell.get(own);

            if (ownLevel < (neiLevel[i]-1))
            {
                if (!maxSet)
                {
                    if (unrefineCell.get(own) == 0)
                    {
                        FatalErrorIn("hexRef8::consistentUnrefinement(..)")
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.unset(own);
                    nChanged++;
                }
            }
            else if (neiLevel[i] < (ownLevel-1))
            {
                if (maxSet)
                {
                    if (unrefineCell.get(own) == 1)
                    {
                        FatalErrorIn("hexRef8::consistentUnrefinement(..)")
                            << "problem" << abort(FatalError);
                    }

                    unrefineCell.set(own);
                    nChanged++;
                }
            }
        }

        reduce(nChanged, sumOp<label>());

        if (debug)
        {
            Pout<< "hexRef8::consistentUnrefinement :"
                << " Changed " << nChanged
                << " refinement levels due to 2:1 conflicts."
                << endl;
        }

        if (nChanged == 0)
        {
            break;
        }


        // Convert cellsToUnrefine back into points to unrefine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Knock out any point whose cell neighbour cannot be unrefined.
        forAll(unrefinePoint, pointI)
        {
            if (unrefinePoint.get(pointI))
            {
                const labelList& pCells = mesh_.pointCells(pointI);

                forAll(pCells, j)
                {
                    if (!unrefineCell.get(pCells[j]))
                    {
                        unrefinePoint.unset(pointI);
                        break;
                    }
                }
            }
        }
    }


    // Convert back to labelList.
    label nSet = 0;

    forAll(unrefinePoint, pointI)
    {
        if (unrefinePoint.get(pointI))
        {
            nSet++;
        }
    }

    labelList newPointsToUnrefine(nSet);
    nSet = 0;

    forAll(unrefinePoint, pointI)
    {
        if (unrefinePoint.get(pointI))
        {
            newPointsToUnrefine[nSet++] = pointI;
        }
    }

    return newPointsToUnrefine;
}


void Foam::hexRef8::setUnrefinement
(
    const labelList& splitPointLabels,
    polyTopoChange& meshMod
)
{
    if (!history_.active())
    {
        FatalErrorIn
        (
            "hexRef8::setUnrefinement(const labelList&, polyTopoChange&)"
        )   << "Only call if constructed with history capability"
            << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "hexRef8::setUnrefinement :"
            << " Checking initial mesh just to make sure" << endl;

        checkMesh();

        forAll(cellLevel_, cellI)
        {
            if (cellLevel_[cellI] < 0)
            {
                FatalErrorIn
                (
                    "hexRef8::setUnrefinement"
                    "("
                        "const labelList&, "
                        "polyTopoChange&"
                    ")"
                )   << "Illegal cell level " << cellLevel_[cellI]
                    << " for cell " << cellI
                    << abort(FatalError);
            }
        }


        // Write to sets.
        pointSet pSet(mesh_, "splitPoints", splitPointLabels);
        pSet.write();

        cellSet cSet(mesh_, "splitPointCells", splitPointLabels.size());

        forAll(splitPointLabels, i)
        {
            const labelList& pCells = mesh_.pointCells(splitPointLabels[i]);

            forAll(pCells, j)
            {
                cSet.insert(pCells[j]);
            }
        }
        cSet.write();

        Pout<< "hexRef8::setRefinement : Dumping " << pSet.size()
            << " points and "
            << cSet.size() << " cells for unrefinement to" << nl
            << "    pointSet " << pSet.objectPath() << nl
            << "    cellSet " << cSet.objectPath()
            << endl;
    }


    labelList cellRegion;
    labelList cellRegionMaster;
    labelList facesToRemove;

    {
        labelHashSet splitFaces(12*splitPointLabels.size());

        forAll(splitPointLabels, i)
        {
            const labelList& pFaces = mesh_.pointFaces()[splitPointLabels[i]];

            forAll(pFaces, j)
            {
                splitFaces.insert(pFaces[j]);
            }
        }

        // Check with faceRemover what faces will get removed. Note that this
        // can be more (but never less) than splitFaces provided.
        faceRemover_.compatibleRemoves
        (
            splitFaces.toc(),   // pierced faces
            cellRegion,         // per cell -1 or region it is merged into
            cellRegionMaster,   // per region the master cell
            facesToRemove       // new faces to be removed.
        );

        if (facesToRemove.size() != splitFaces.size())
        {
            FatalErrorIn
            (
                "hexRef8::setUnrefinement(const labelList&, polyTopoChange&)"
            )   << "Ininitial set of split points to unrefine does not"
                << " seem to be consistent or not mid points of refined cells"
                << abort(FatalError);
        }
    }

    // Redo the region master so it is consistent with our master.
    // This will guarantee that the new cell (for which faceRemover uses
    // the region master) is already compatible with our refinement structure.

    forAll(splitPointLabels, i)
    {
        label pointI = splitPointLabels[i];

        // Get original cell label

        const labelList& pCells = mesh_.pointCells(pointI);

        // Check
        if (pCells.size() != 8)
        {
            FatalErrorIn
            (
                "hexRef8::setUnrefinement(const labelList&, polyTopoChange&)"
            )   << "splitPoint " << pointI
                << " should have 8 cells using it. It has " << pCells
                << abort(FatalError);
        }


        // Check that the lowest numbered pCells is the master of the region
        // (should be guaranteed by directRemoveFaces)
        //if (debug)
        {
            label masterCellI = min(pCells);

            forAll(pCells, j)
            {
                label cellI = pCells[j];

                label region = cellRegion[cellI];

                if (region == -1)
                {
                    FatalErrorIn("hexRef8::setUnrefinement(..)")
                        << "Ininitial set of split points to unrefine does not"
                        << " seem to be consistent or not mid points"
                        << " of refined cells" << nl
                        << "cell:" << cellI << " on splitPoint " << pointI
                        << " has no region to be merged into"
                        << abort(FatalError);
                }

                if (masterCellI != cellRegionMaster[region])
                {
                    FatalErrorIn("hexRef8::setUnrefinement(..)")
                        << "cell:" << cellI << " on splitPoint:" << pointI
                        << " in region " << region
                        << " has master:" << cellRegionMaster[region]
                        << " which is not the lowest numbered cell"
                        << " among the pointCells:" << pCells
                        << abort(FatalError);
                }
            }
        }
    }

    // Insert all commands to combine cells. Never fails so don't have to
    // test for success.
    faceRemover_.setRefinement
    (
        facesToRemove,
        cellRegion,
        cellRegionMaster,
        meshMod
    );

    // Remove the 8 cells that originated from merging around the split point
    // and adapt cell levels (not that pointLevels stay the same since points
    // either get removed or stay at the same position.
    forAll(splitPointLabels, i)
    {
        label pointI = splitPointLabels[i];

        const labelList& pCells = mesh_.pointCells(pointI);

        label masterCellI = min(pCells);

        forAll(pCells, j)
        {
            cellLevel_[pCells[j]]--;
        }

        history_.combineCells(masterCellI, pCells);
    }

    // Mark files as changed
    setInstance(mesh_.facesInstance());

    // history_.updateMesh will take care of truncating.
}


// Write refinement to polyMesh directory.
bool Foam::hexRef8::write() const
{
    bool writeOk =
        cellLevel_.write()
     && pointLevel_.write()
     && level0Edge_.write();

    if (history_.active())
    {
        writeOk = writeOk && history_.write();
    }

    return writeOk;
}


// ************************************************************************* //
