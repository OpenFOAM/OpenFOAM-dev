/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "slidingInterface.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "polyTopoChange.H"
#include "polyTopoChanger.H"
#include "polyModifyFace.H"
#include "polyModifyPoint.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::slidingInterface::decoupleInterface
(
    polyTopoChange& ref
) const
{
    if (debug)
    {
        Pout<< "void slidingInterface::decoupleInterface("
            << "polyTopoChange& ref) const : "
            << "Decoupling sliding interface " << name() << endl;
    }

    if (!attached_)
    {
        if (debug)
        {
            Pout<< "void slidingInterface::decoupleInterface("
                << "polyTopoChange& ref) const : "
                << "Interface already decoupled." << endl;
        }

        return;
    }

    // Clear previous couple
    clearCouple(ref);

    const polyMesh& mesh = topoChanger().mesh();
    const faceList& faces = mesh.faces();
    const cellList& cells = mesh.cells();

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    // Master side

    const primitiveFacePatch& masterPatch =
        mesh.faceZones()[masterFaceZoneID_.index()]();

    const labelList& masterPatchAddr =
        mesh.faceZones()[masterFaceZoneID_.index()];

    const boolList& masterPatchFlip =
        mesh.faceZones()[masterFaceZoneID_.index()].flipMap();

    const labelList& masterFc = masterFaceCells();

    // Recover faces in master patch

    forAll(masterPatchAddr, faceI)
    {
        // Make a copy of the face and turn it if necessary
        face newFace = faces[masterPatchAddr[faceI]];

        if (masterPatchFlip[faceI])
        {
            newFace.flip();
        }

        ref.setAction
        (
            polyModifyFace
            (
                newFace,                         // new face
                masterPatchAddr[faceI],          // master face index
                masterFc[faceI],                 // owner
                -1,                              // neighbour
                false,                           // flux flip
                masterPatchID_.index(),          // patch ID
                false,                           // remove from zone
                masterFaceZoneID_.index(),       // zone ID
                false                            // zone flip.  Face corrected
            )
        );

        // Pout<< "Modifying master patch face no "
        //     << masterPatchAddr[faceI]
        //     << " face: " << faces[masterPatchAddr[faceI]]
        //     << " old owner: " << own[masterPatchAddr[faceI]]
        //     << " new owner: " << masterFc[faceI]
        //     << endl;
    }

    // Slave side

    const primitiveFacePatch& slavePatch =
        mesh.faceZones()[slaveFaceZoneID_.index()]();

    const labelList& slavePatchAddr =
        mesh.faceZones()[slaveFaceZoneID_.index()];

    const boolList& slavePatchFlip =
        mesh.faceZones()[slaveFaceZoneID_.index()].flipMap();

    const labelList& slaveFc = slaveFaceCells();

    // Grab retired point mapping
    const Map<label>& rpm = retiredPointMap();

    // Recover faces in slave patch

    forAll(slavePatchAddr, faceI)
    {
        // Make a copy of face and turn it if necessary
        face newFace = faces[slavePatchAddr[faceI]];

        if (slavePatchFlip[faceI])
        {
            newFace.flip();
        }

        // Recover retired points on the slave side
        forAll(newFace, pointI)
        {
            Map<label>::const_iterator rpmIter = rpm.find(newFace[pointI]);
            if (rpmIter != rpm.end())
            {
                // Master of retired point; grab its original
                // Pout<< "Reinstating retired point: " << newFace[pointI]
                //     << " with old: " << rpm.find(newFace[pointI])()
                //     << endl;

                newFace[pointI] = rpmIter();
            }
        }

        ref.setAction
        (
            polyModifyFace
            (
                newFace,                         // new face
                slavePatchAddr[faceI],           // master face index
                slaveFc[faceI],                  // owner
                -1,                              // neighbour
                false,                           // flux flip
                slavePatchID_.index(),           // patch ID
                false,                           // remove from zone
                slaveFaceZoneID_.index(),        // zone ID
                false                            // zone flip.  Face corrected
            )
        );
    }

    // Re-create the master stick-out faces

    // Grab the list of faces in the layer
    const labelList& masterStickOuts = masterStickOutFaces();

    forAll(masterStickOuts, faceI)
    {
        // Renumber the face and remove additional points

        const label curFaceID = masterStickOuts[faceI];

        const face& oldFace = faces[curFaceID];

        DynamicList<label> newFaceLabels(oldFace.size());

        bool changed = false;

        forAll(oldFace, pointI)
        {
            // Check if the point is removed
            if (ref.pointRemoved(oldFace[pointI]))
            {
                // Point removed; skip it
                changed = true;
            }
            else
            {
                newFaceLabels.append(oldFace[pointI]);
            }
        }

        if (changed)
        {
            if (newFaceLabels.size() < 3)
            {
                FatalErrorIn
                (
                    "void slidingInterface::decoupleInterface("
                    "polyTopoChange& ref) const"
                )   << "Face " << curFaceID << " reduced to less than "
                    << "3 points.  Topological/cutting error." << nl
                    << "Old face: " << oldFace << " new face: " << newFaceLabels
                    << abort(FatalError);
            }

            // Get face zone and its flip
            label modifiedFaceZone = mesh.faceZones().whichZone(curFaceID);
            bool modifiedFaceZoneFlip = false;

            if (modifiedFaceZone >= 0)
            {
                modifiedFaceZoneFlip =
                    mesh.faceZones()[modifiedFaceZone].flipMap()
                    [
                        mesh.faceZones()[modifiedFaceZone].whichFace(curFaceID)
                    ];
            }

            face newFace;
            newFace.transfer(newFaceLabels);

            // Pout<< "Modifying master stick-out face " << curFaceID
            //     << " old face: " << oldFace
            //     << " new face: " << newFace
            //     << endl;

            // Modify the face
            ref.setAction
            (
                polyModifyFace
                (
                    newFace,                // modified face
                    curFaceID,              // label of face being modified
                    own[curFaceID],         // owner
                    nei[curFaceID],         // neighbour
                    false,                  // face flip
                    mesh.boundaryMesh().whichPatch(curFaceID), // patch for face
                    false,                  // remove from zone
                    modifiedFaceZone,       // zone for face
                    modifiedFaceZoneFlip    // face flip in zone
                )
            );
        }
    }

    // Re-create the slave stick-out faces

    labelHashSet slaveLayerCellFaceMap
    (
        primitiveMesh::facesPerCell_*(masterPatch.size() + slavePatch.size())
    );

    forAll(slaveFc, faceI)
    {
        const labelList& curFaces = cells[slaveFc[faceI]];

        forAll(curFaces, faceI)
        {
            // Check if the face belongs to the slave face zone; and
            // if it has been removed; if not add it
            if
            (
                mesh.faceZones().whichZone(curFaces[faceI])
             != slaveFaceZoneID_.index()
             && !ref.faceRemoved(curFaces[faceI])

            )
            {
                slaveLayerCellFaceMap.insert(curFaces[faceI]);
            }
        }
    }

    // Grab the list of faces in the layer
    const labelList& slaveStickOuts = slaveStickOutFaces();

    // Grab master point mapping
    const Map<label>& masterPm = masterPatch.meshPointMap();

    forAll(slaveStickOuts, faceI)
    {
        // Renumber the face and remove additional points

        const label curFaceID = slaveStickOuts[faceI];

        const face& oldFace = faces[curFaceID];

        DynamicList<label> newFaceLabels(oldFace.size());

        bool changed = false;

        forAll(oldFace, pointI)
        {
            // Check if the point is removed or retired
            if (rpm.found(oldFace[pointI]))
            {
                // Master of retired point; grab its original
                changed = true;

                // Pout<< "Reinstating retired point: " << oldFace[pointI]
                //     << " with old: " << rpm.find(oldFace[pointI])()
                //     << endl;

                newFaceLabels.append(rpm.find(oldFace[pointI])());
            }
            else if (ref.pointRemoved(oldFace[pointI]))
            {
                // Point removed; skip it
                changed = true;
            }
            else if (masterPm.found(oldFace[pointI]))
            {
                // Point from master patch only; skip it
                changed = true;
            }
            else
            {
                newFaceLabels.append(oldFace[pointI]);
            }
        }

        if (changed)
        {
            if (newFaceLabels.size() < 3)
            {
                FatalErrorIn
                (
                    "void slidingInterface::decoupleInterface("
                    "polyTopoChange& ref) const"
                )   << "Face " << curFaceID << " reduced to less than "
                    << "3 points.  Topological/cutting error." << nl
                    << "Old face: " << oldFace << " new face: " << newFaceLabels
                    << abort(FatalError);
            }

            // Get face zone and its flip
            label modifiedFaceZone = mesh.faceZones().whichZone(curFaceID);
            bool modifiedFaceZoneFlip = false;

            if (modifiedFaceZone >= 0)
            {
                modifiedFaceZoneFlip =
                    mesh.faceZones()[modifiedFaceZone].flipMap()
                    [
                        mesh.faceZones()[modifiedFaceZone].whichFace(curFaceID)
                    ];
            }

            face newFace;
            newFace.transfer(newFaceLabels);

            // Pout<< "Modifying slave stick-out face " << curFaceID
            //     << " old face: " << oldFace
            //     << " new face: " << newFace
            //     << endl;

            // Modify the face
            ref.setAction
            (
                polyModifyFace
                (
                    newFace,                // modified face
                    curFaceID,              // label of face being modified
                    own[curFaceID],         // owner
                    nei[curFaceID],         // neighbour
                    false,                  // face flip
                    mesh.boundaryMesh().whichPatch(curFaceID), // patch for face
                    false,                  // remove from zone
                    modifiedFaceZone,       // zone for face
                    modifiedFaceZoneFlip    // face flip in zone
                )
            );
        }
    }

    // Bring all slave patch points back to life
    const pointField& points = mesh.points();

    const labelList& slaveMeshPoints =
        mesh.faceZones()[slaveFaceZoneID_.index()]().meshPoints();

    forAll(slaveMeshPoints, pointI)
    {
        ref.setAction
        (
            polyModifyPoint
            (
                slaveMeshPoints[pointI],             // point ID
                points[slaveMeshPoints[pointI]],     // point
                false,                               // remove from zone
                mesh.pointZones().whichZone(slaveMeshPoints[pointI]), // zone
                true                                // in a cell
            )
        );
    }

    // Clear the retired point numbering
    retiredPointMapPtr_->clear();

    // Finished decoupling
    attached_ = false;

    if (debug)
    {
        Pout<< "void slidingInterface::coupleInterface("
            << "polyTopoChange& ref) const : "
            << "Finished decoupling sliding interface " << name() << endl;
    }
}


// ************************************************************************* //
