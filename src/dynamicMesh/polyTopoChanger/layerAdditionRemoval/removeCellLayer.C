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

#include "layerAdditionRemoval.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "polyTopoChange.H"
#include "oppositeFace.H"
#include "polyTopoChanger.H"
#include "polyRemoveCell.H"
#include "polyRemoveFace.H"
#include "polyRemovePoint.H"
#include "polyModifyFace.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::layerAdditionRemoval::validCollapse() const
{
    // Check for valid layer collapse
    // - no boundary-to-boundary collapse

    if (debug)
    {
        Pout<< "Checking layer collapse for object " << name() << endl;
    }

    // Grab the face collapse mapping
    const polyMesh& mesh = topoChanger().mesh();

    const labelList& ftc = facesPairing();
    const labelList& mf = mesh.faceZones()[faceZoneID_.index()];

    label nBoundaryHits = 0;

    forAll(mf, facei)
    {
        if
        (
            !mesh.isInternalFace(mf[facei])
         && !mesh.isInternalFace(ftc[facei])
        )
        {
            nBoundaryHits++;
        }
    }


    if (debug)
    {
        Pout<< "Finished checking layer collapse for object "
            << name() <<".  Number of boundary-on-boundary hits: "
            << nBoundaryHits << endl;
    }

    if (returnReduce(nBoundaryHits, sumOp<label>()) > 0)
    {
        return false;
    }
    else
    {
        return true;
    }
}


void Foam::layerAdditionRemoval::removeCellLayer
(
    polyTopoChange& ref
) const
{
    // Algorithm for layer removal.  Second phase: topological change
    // 2)  Go through all the faces of the master cell layer and remove
    //     the ones that are not in the master face zone.
    //
    // 3)  Grab all the faces coming out of points that are collapsed
    //     and select the ones that are not marked for removal.  Go
    //     through the remaining faces and replace the point to remove by
    //     the equivalent point in the master face zone.
    if (debug)
    {
        Pout<< "Removing the cell layer for object " << name() << endl;
    }

    const polyMesh& mesh = topoChanger().mesh();

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    const labelList& ptc = pointsPairing();
    const labelList& ftc = facesPairing();

    // Remove all the cells from the master layer
    const labelList& mc =
        mesh.faceZones()[faceZoneID_.index()].masterCells();

    forAll(mc, facei)
    {
        label slaveSideCell = own[ftc[facei]];

        if (mesh.isInternalFace(ftc[facei]) && slaveSideCell == mc[facei])
        {
            // Owner cell of the face is being removed.
            // Grab the neighbour instead
            slaveSideCell = nei[ftc[facei]];
        }

        ref.setAction(polyRemoveCell(mc[facei], slaveSideCell));
    }

    // Remove all the faces from the master layer cells which are not in
    // the master face layer
    labelHashSet facesToRemoveMap(mc.size()*primitiveMesh::facesPerCell_);

    const cellList& cells = mesh.cells();

    forAll(mc, celli)
    {
        const cell& curCell = cells[mc[celli]];

        forAll(curCell, facei)
        {
            // Check if the face is in the master zone.  If not, remove it
            if
            (
                mesh.faceZones().whichZone(curCell[facei])
             != faceZoneID_.index()
            )
            {
                facesToRemoveMap.insert(curCell[facei]);
            }
        }
    }

    forAllConstIter(labelHashSet, facesToRemoveMap, iter)
    {
        ref.setAction(polyRemoveFace(iter.key()));
    }

    // Remove all points that will be collapsed
    forAll(ptc, pointi)
    {
        ref.setAction(polyRemovePoint(ptc[pointi]));
    }

    // Grab all faces coming off points to be deleted.  If the face
    // has not been deleted, replace the removed point with the
    // equivalent from the master layer.

    // Make a map of all point to be removed, giving the master point label
    // for each of them

    Map<label> removedPointMap(2*ptc.size());

    const labelList& meshPoints =
        mesh.faceZones()[faceZoneID_.index()]().meshPoints();

    forAll(ptc, pointi)
    {
        removedPointMap.insert(ptc[pointi], meshPoints[pointi]);
    }

    const labelListList& pf = mesh.pointFaces();

    const faceList& faces = mesh.faces();

    // Make a list of faces to be modified using the map to avoid duplicates
    labelHashSet facesToModify(ptc.size()*primitiveMesh::facesPerPoint_);

    forAll(ptc, pointi)
    {
        const labelList& curFaces = pf[ptc[pointi]];

        forAll(curFaces, facei)
        {
            if (!facesToRemoveMap.found(curFaces[facei]))
            {
                facesToModify.insert(curFaces[facei]);
            }
        }
    }

    labelList ftm = facesToModify.toc();

    if (debug > 1)
    {
        Pout<< "faces to modify: " << ftm << endl;
    }

    forAll(ftm, facei)
    {
        // For every face to modify, copy the face and re-map the vertices.
        // It is known all the faces will be changed since they hang off
        // re-mapped vertices
        label curFaceID = ftm[facei];

        face newFace(faces[curFaceID]);

        forAll(newFace, pointi)
        {
            Map<label>::iterator rpmIter =
                removedPointMap.find(newFace[pointi]);

            if (rpmIter != removedPointMap.end())
            {
                // Point mapped. Replace it
                newFace[pointi] = rpmIter();
            }
        }

        if (debug > 1)
        {
            Pout<< "face label: " << curFaceID
                << " old face: " << faces[curFaceID]
                << " new face: " << newFace << endl;
        }

        // Get face zone and its flip
        label modifiedFaceZone = mesh.faceZones().whichZone(curFaceID);
        bool modifiedFaceZoneFlip = false;

        if (modifiedFaceZone >= 0)
        {
            const faceZone& fz = mesh.faceZones()[modifiedFaceZone];
            modifiedFaceZoneFlip = fz.flipMap()[fz.whichFace(curFaceID)];
        }

        label newNeighbour = -1;

        if (curFaceID < mesh.nInternalFaces())
        {
            newNeighbour = nei[curFaceID];
        }

        // Modify the face
        ref.setAction
        (
            polyModifyFace
            (
                newFace,                // modified face
                curFaceID,              // label of face being modified
                own[curFaceID],         // owner
                newNeighbour,           // neighbour
                false,                  // face flip
                mesh.boundaryMesh().whichPatch(curFaceID),// patch for face
                false,                  // remove from zone
                modifiedFaceZone,       // zone for face
                modifiedFaceZoneFlip    // face flip in zone
            )
        );
    }

    // Modify the faces in the master layer to point past the removed cells

    const labelList& mf = mesh.faceZones()[faceZoneID_.index()];
    const boolList& mfFlip = mesh.faceZones()[faceZoneID_.index()].flipMap();

    forAll(mf, facei)
    {
        // Grab the owner and neighbour of the faces to be collapsed and get rid
        // of the cell to be removed
        label masterSideCell = own[mf[facei]];

        if (masterSideCell == mc[facei])
        {
            if (mesh.isInternalFace(mf[facei]))
            {
                // Owner cell of the face is being removed.
                // Grab the neighbour instead
                masterSideCell = nei[mf[facei]];
            }
            else
            {
                masterSideCell = -1;
            }
        }

        label slaveSideCell = own[ftc[facei]];

        if (slaveSideCell == mc[facei])
        {
            if (mesh.isInternalFace(ftc[facei]))
            {
                // Owner cell of the face is being removed.
                // Grab the neighbour instead
                slaveSideCell = nei[ftc[facei]];
            }
            else
            {
                slaveSideCell = -1;
            }
        }

        // Find out if the face needs to be flipped
        label newOwner = -1;
        label newNeighbour = -1;
        bool flipFace = false;
        label newPatchID = -1;
        label newZoneID = -1;

        // Is any of the faces a boundary face?  If so, grab the patch
        // A boundary-to-boundary collapse is checked for in validCollapse()
        // and cannot happen here.

        if (!mesh.isInternalFace(mf[facei]))
        {
            // Master is the boundary face: it gets a new owner but no flip
            newOwner = slaveSideCell;
            newNeighbour = -1;
            flipFace = false;
            newPatchID = mesh.boundaryMesh().whichPatch(mf[facei]);
            newZoneID = mesh.faceZones().whichZone(mf[facei]);
        }
        else if (!mesh.isInternalFace(ftc[facei]))
        {
            // Slave is the boundary face: grab its patch
            newOwner = slaveSideCell;
            newNeighbour = -1;

            // Find out if the face flip is necessary
            if (own[mf[facei]] == slaveSideCell)
            {
                flipFace = false;
            }
            else
            {
                flipFace = true;
            }

            newPatchID = mesh.boundaryMesh().whichPatch(ftc[facei]);

            // The zone of the master face is preserved
            newZoneID = mesh.faceZones().whichZone(mf[facei]);
        }
        else
        {
            // Both faces are internal: flip is decided based on which of the
            // new cells around it has a lower label
            newOwner = min(masterSideCell, slaveSideCell);
            newNeighbour = max(masterSideCell, slaveSideCell);

            if (newOwner == own[mf[facei]] || newNeighbour == nei[mf[facei]])
            {
                flipFace = false;
            }
            else
            {
                flipFace = true;
            }

            newPatchID = -1;

            // The zone of the master face is preserved
            newZoneID = mesh.faceZones().whichZone(mf[facei]);
        }

        // Modify the face and flip if necessary
        face newFace = faces[mf[facei]];
        bool zoneFlip = mfFlip[facei];

        if (flipFace)
        {
            newFace.flip();
            zoneFlip = !zoneFlip;
        }

        if (debug > 1)
        {
            Pout<< "Modifying face " << mf[facei]
                << " newFace: " << newFace << nl
                << " newOwner: " << newOwner
                << " newNeighbour: " << newNeighbour
                << " flipFace: " << flipFace
                << " newPatchID: " << newPatchID
                << " newZoneID: " << newZoneID << nl
                << " oldOwn: " << own[mf[facei]]
                << " oldNei: " << nei[mf[facei]] << endl;
        }

        ref.setAction
        (
            polyModifyFace
            (
                newFace,       // modified face
                mf[facei],     // label of face being modified
                newOwner,      // owner
                newNeighbour,  // neighbour
                flipFace,      // flip
                newPatchID,    // patch for face
                false,         // remove from zone
                newZoneID,     // new zone
                zoneFlip       // face flip in zone
            )
        );
    }
}


// ************************************************************************* //
