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

#include "cellSplitter.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "polyAddCell.H"
#include "polyAddFace.H"
#include "polyAddPoint.H"
#include "polyModifyFace.H"
#include "mapPolyMesh.H"
#include "meshTools.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(cellSplitter, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cellSplitter::getFaceInfo
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


// Find the new owner of faceI (since the original cell has been split into
// newCells
Foam::label Foam::cellSplitter::newOwner
(
    const label faceI,
    const Map<labelList>& cellToCells
) const
{
    label oldOwn = mesh_.faceOwner()[faceI];

    Map<labelList>::const_iterator fnd = cellToCells.find(oldOwn);

    if (fnd == cellToCells.end())
    {
        // Unsplit cell
        return oldOwn;
    }
    else
    {
        // Look up index of face in the cells' faces.

        const labelList& newCells = fnd();

        const cell& cFaces = mesh_.cells()[oldOwn];

        return newCells[findIndex(cFaces, faceI)];
    }
}


Foam::label Foam::cellSplitter::newNeighbour
(
    const label faceI,
    const Map<labelList>& cellToCells
) const
{
    label oldNbr = mesh_.faceNeighbour()[faceI];

    Map<labelList>::const_iterator fnd = cellToCells.find(oldNbr);

    if (fnd == cellToCells.end())
    {
        // Unsplit cell
        return oldNbr;
    }
    else
    {
        // Look up index of face in the cells' faces.

        const labelList& newCells = fnd();

        const cell& cFaces = mesh_.cells()[oldNbr];

        return newCells[findIndex(cFaces, faceI)];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::cellSplitter::cellSplitter(const polyMesh& mesh)
:
    mesh_(mesh),
    addedPoints_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellSplitter::~cellSplitter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cellSplitter::setRefinement
(
    const Map<point>& cellToMidPoint,
    polyTopoChange& meshMod
)
{
    addedPoints_.clear();
    addedPoints_.resize(cellToMidPoint.size());


    //
    // Introduce cellToMidPoints.
    //

    forAllConstIter(Map<point>, cellToMidPoint, iter)
    {
        label cellI = iter.key();

        label anchorPoint = mesh_.cellPoints()[cellI][0];

        label addedPointI =
            meshMod.setAction
            (
                polyAddPoint
                (
                    iter(),         // point
                    anchorPoint,    // master point
                    -1,             // zone for point
                    true            // supports a cell
                )
            );
        addedPoints_.insert(cellI, addedPointI);

        //Pout<< "Added point " << addedPointI
        //    << iter() << " in cell " << cellI << " with centre "
        //    << mesh_.cellCentres()[cellI] << endl;
    }


    //
    // Add cells (first one is modified original cell)
    //

    Map<labelList> cellToCells(cellToMidPoint.size());

    forAllConstIter(Map<point>, cellToMidPoint, iter)
    {
        label cellI = iter.key();

        const cell& cFaces = mesh_.cells()[cellI];

        // Cells created for this cell.
        labelList newCells(cFaces.size());

        // First pyramid is the original cell
        newCells[0] = cellI;

        // Add other pyramids
        for (label i = 1; i < cFaces.size(); i++)
        {
            label addedCellI =
                meshMod.setAction
                (
                    polyAddCell
                    (
                        -1,     // master point
                        -1,     // master edge
                        -1,     // master face
                        cellI,  // master cell
                        -1      // zone
                    )
                );

            newCells[i] = addedCellI;
        }

        cellToCells.insert(cellI, newCells);

        //Pout<< "Split cell " << cellI
        //    << " with centre " << mesh_.cellCentres()[cellI] << nl
        //    << " faces:" << cFaces << nl
        //    << " into :" << newCells << endl;
    }


    //
    // Introduce internal faces. These go from edges of the cell to the mid
    // point.
    //

    forAllConstIter(Map<point>, cellToMidPoint, iter)
    {
        label cellI = iter.key();

        label midPointI = addedPoints_[cellI];

        const cell& cFaces = mesh_.cells()[cellI];

        const labelList& cEdges = mesh_.cellEdges()[cellI];

        forAll(cEdges, i)
        {
            label edgeI = cEdges[i];
            const edge& e = mesh_.edges()[edgeI];

            // Get the faces on the cell using the edge
            label face0, face1;
            meshTools::getEdgeFaces(mesh_, cellI, edgeI, face0, face1);

            // Get the cells on both sides of the face by indexing into cFaces.
            // (since newly created cells are stored in cFaces order)
            const labelList& newCells = cellToCells[cellI];

            label cell0 = newCells[findIndex(cFaces, face0)];
            label cell1 = newCells[findIndex(cFaces, face1)];

            if (cell0 < cell1)
            {
                // Construct face to midpoint that is pointing away from
                // (pyramid split off from) cellI

                const face& f0 = mesh_.faces()[face0];

                label index = findIndex(f0, e[0]);

                bool edgeInFaceOrder = (f0[f0.fcIndex(index)] == e[1]);

                // Check if cellI is the face owner

                face newF(3);
                if (edgeInFaceOrder == (mesh_.faceOwner()[face0] == cellI))
                {
                    // edge used in face order.
                    newF[0] = e[1];
                    newF[1] = e[0];
                    newF[2] = midPointI;
                }
                else
                {
                    newF[0] = e[0];
                    newF[1] = e[1];
                    newF[2] = midPointI;
                }

                // Now newF points away from cell0
                meshMod.setAction
                (
                    polyAddFace
                    (
                        newF,                       // face
                        cell0,                      // owner
                        cell1,                      // neighbour
                        -1,                         // master point
                        -1,                         // master edge
                        face0,                      // master face for addition
                        false,                      // flux flip
                        -1,                         // patch for face
                        -1,                         // zone for face
                        false                       // face zone flip
                    )
                );
            }
            else
            {
                // Construct face to midpoint that is pointing away from
                // (pyramid split off from) cellI

                const face& f1 = mesh_.faces()[face1];

                label index = findIndex(f1, e[0]);

                bool edgeInFaceOrder = (f1[f1.fcIndex(index)] == e[1]);

                // Check if cellI is the face owner

                face newF(3);
                if (edgeInFaceOrder == (mesh_.faceOwner()[face1] == cellI))
                {
                    // edge used in face order.
                    newF[0] = e[1];
                    newF[1] = e[0];
                    newF[2] = midPointI;
                }
                else
                {
                    newF[0] = e[0];
                    newF[1] = e[1];
                    newF[2] = midPointI;
                }

                // Now newF points away from cell1
                meshMod.setAction
                (
                    polyAddFace
                    (
                        newF,                       // face
                        cell1,                      // owner
                        cell0,                      // neighbour
                        -1,                         // master point
                        -1,                         // master edge
                        face0,                      // master face for addition
                        false,                      // flux flip
                        -1,                         // patch for face
                        -1,                         // zone for face
                        false                       // face zone flip
                    )
                );
            }
        }
    }


    //
    // Update all existing faces for split owner or neighbour.
    //


    // Mark off affected face.
    boolList faceUpToDate(mesh_.nFaces(), true);

    forAllConstIter(Map<point>, cellToMidPoint, iter)
    {
        label cellI = iter.key();

        const cell& cFaces = mesh_.cells()[cellI];

        forAll(cFaces, i)
        {
            label faceI = cFaces[i];

            faceUpToDate[faceI] = false;
        }
    }

    forAll(faceUpToDate, faceI)
    {
        if (!faceUpToDate[faceI])
        {
            const face& f = mesh_.faces()[faceI];

            if (mesh_.isInternalFace(faceI))
            {
                label newOwn = newOwner(faceI, cellToCells);
                label newNbr = newNeighbour(faceI, cellToCells);

                if (newOwn < newNbr)
                {
                    meshMod.setAction
                    (
                        polyModifyFace
                        (
                            f,
                            faceI,
                            newOwn,         // owner
                            newNbr,         // neighbour
                            false,          // flux flip
                            -1,             // patch for face
                            false,          // remove from zone
                            -1,             // zone for face
                            false           // face zone flip
                        )
                    );
                }
                else
                {
                    meshMod.setAction
                    (
                        polyModifyFace
                        (
                            f.reverseFace(),
                            faceI,
                            newNbr,         // owner
                            newOwn,         // neighbour
                            false,          // flux flip
                            -1,             // patch for face
                            false,          // remove from zone
                            -1,             // zone for face
                            false           // face zone flip
                        )
                    );
                }

            }
            else
            {
                label newOwn = newOwner(faceI, cellToCells);

                label patchID, zoneID, zoneFlip;
                getFaceInfo(faceI, patchID, zoneID, zoneFlip);

                meshMod.setAction
                (
                    polyModifyFace
                    (
                        mesh_.faces()[faceI],
                        faceI,
                        newOwn,         // owner
                        -1,             // neighbour
                        false,          // flux flip
                        patchID,        // patch for face
                        false,          // remove from zone
                        zoneID,         // zone for face
                        zoneFlip        // face zone flip
                    )
                );
            }

            faceUpToDate[faceI] = true;
        }
    }
}


void Foam::cellSplitter::updateMesh(const mapPolyMesh& morphMap)
{
    // Create copy since we're deleting entries. Only if both cell and added
    // point get mapped do they get inserted.
    Map<label> newAddedPoints(addedPoints_.size());

    forAllConstIter(Map<label>, addedPoints_, iter)
    {
        label oldCellI = iter.key();

        label newCellI = morphMap.reverseCellMap()[oldCellI];

        label oldPointI = iter();

        label newPointI = morphMap.reversePointMap()[oldPointI];

        if (newCellI >= 0 && newPointI >= 0)
        {
            newAddedPoints.insert(newCellI, newPointI);
        }
    }

    // Copy
    addedPoints_.transfer(newAddedPoints);
}


// ************************************************************************* //
