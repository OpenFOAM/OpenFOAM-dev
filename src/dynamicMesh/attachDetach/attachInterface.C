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

#include "attachDetach.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "polyTopoChange.H"
#include "polyTopoChanger.H"
#include "polyRemovePoint.H"
#include "polyRemoveFace.H"
#include "polyModifyFace.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::attachDetach::positionDifference_ = 1e-8;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::attachDetach::attachInterface
(
    polyTopoChange& ref
) const
{
    // Algorithm:
    // 1. Create the reverse patch out of the slave faces.
    // 2. Go through all the mesh points from the master and slave patch.
    //    If the point labels are different, insert them into the point
    //    renumbering list and remove them from the mesh.
    // 3. Remove all faces from the slave patch
    // 4. Modify all the faces from the master patch by making them internal
    //    between the faceCell cells for the two patches.  If the master owner
    //    is higher than the slave owner, turn the face around
    // 5. Get all the faces attached to the slave patch points.
    //    If they have not been removed, renumber them using the
    //    point renumbering list.

    if (debug)
    {
        Pout<< "void attachDetach::attachInterface("
            << "polyTopoChange& ref) const "
            << " for object " << name() << " : "
            << "Attaching interface" << endl;
    }

    const polyMesh& mesh = topoChanger().mesh();
    const faceList& faces = mesh.faces();
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    const polyPatch& masterPatch = mesh.boundaryMesh()[masterPatchID_.index()];
    const polyPatch& slavePatch = mesh.boundaryMesh()[slavePatchID_.index()];

    const label masterPatchStart = masterPatch.start();
    const label slavePatchStart = slavePatch.start();

    const labelList& slaveMeshPoints = slavePatch.meshPoints();

    const Map<label>& removedPointMap = pointMatchMap();

    const labelList removedPoints = removedPointMap.toc();

    forAll(removedPoints, pointI)
    {
        //Pout<< "Removing point:" << removedPoints[pointI]
        //    << " currently at:" << ref.points()[removedPoints[pointI]]
        //    << endl;

        ref.setAction(polyRemovePoint(removedPoints[pointI]));
    }

// Pout<< "Points to be mapped: " << removedPoints << endl;
    // Remove all faces from the slave patch
    forAll(slavePatch, i)
    {
        //Pout<< "Removing face " << i + slavePatchStart
        //    << " with verts:" << ref.faces()[i + slavePatchStart]
        //    << " at:"
        //    <<  UIndirectList<point>
        //        (
        //            ref.points(),
        //            ref.faces()[i + slavePatchStart]
        //        )
        //    << endl;
        ref.setAction(polyRemoveFace(i + slavePatchStart));
    }

    // Modify the faces from the master patch
    const labelList& masterFaceCells = masterPatch.faceCells();
    const labelList& slaveFaceCells = slavePatch.faceCells();

    const boolList& mfFlip = mesh.faceZones()[faceZoneID_.index()].flipMap();

    forAll(masterFaceCells, faceI)
    {
        // If slave neighbour is greater than master, face does not need
        // turning.  Modify it to become internal
        if (masterFaceCells[faceI] < slaveFaceCells[faceI])
        {
            ref.setAction
            (
                polyModifyFace
                (
                    faces[masterPatchStart + faceI], // modified face
                    masterPatchStart + faceI,    // label of face being modified
                    masterFaceCells[faceI],          // owner
                    slaveFaceCells[faceI],           // neighbour
                    false,                           // face flip
                    -1,                              // patch for face
                    false,                           // remove from zone
                    faceZoneID_.index(),             // zone for face
                    mfFlip[faceI]                    // face flip in zone
                )
            );
        }
        else
        {
            // Flip required
            ref.setAction
            (
                polyModifyFace
                (
                    faces[masterPatchStart + faceI].reverseFace(), // mod face
                    masterPatchStart + faceI,    // label of face being modified
                    slaveFaceCells[faceI],        // owner
                    masterFaceCells[faceI],       // neighbour
                    true,                         // face flip
                    -1,                           // patch for face
                    false,                        // remove from zone
                    faceZoneID_.index(),          // zone for face
                    !mfFlip[faceI]                // face flip in zone
                )
            );
        }
    }

    // Renumber faces affected by point removal
// Pout<< "slaveMeshPoints: " << slaveMeshPoints << endl;
    // Make a map of faces that need to be renumbered
    labelHashSet facesToModifyMap
    (
        slaveMeshPoints.size()*primitiveMesh::facesPerPoint_
    );

    const labelListList& pf = mesh.pointFaces();

    // Grab all the faces off the points in the slave patch.  If the face has
    //  not been removed, add it to the map of faces to renumber
    forAll(slaveMeshPoints, pointI)
    {
        const labelList& curFaces = pf[slaveMeshPoints[pointI]];

        forAll(curFaces, faceI)
        {
            if (!ref.faceRemoved(curFaces[faceI]))
            {
                facesToModifyMap.insert(curFaces[faceI]);
            }
        }
    }

    // Grab the faces to be renumbered
    const labelList ftm = facesToModifyMap.toc();

    forAll(ftm, faceI)
    {
        // For every face to modify, copy the face and re-map the vertices.
        // It is known all the faces will be changed since they hang off
        // re-mapped vertices
        label curFaceID = ftm[faceI];

        face newFace(faces[curFaceID]);

        forAll(newFace, pointI)
        {
            Map<label>::const_iterator rpmIter =
                removedPointMap.find(newFace[pointI]);

            if (rpmIter != removedPointMap.end())
            {
                // Point mapped. Replace it
                newFace[pointI] = rpmIter();
            }
        }

        // Pout<< "face label: " << curFaceID
        //     << " old face: " << faces[curFaceID]
        //     << " new face: " << newFace
        //     << endl;

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


        label patchID = mesh.boundaryMesh().whichPatch(curFaceID);
        label neiCell;
        if (patchID == -1)
        {
            neiCell = nei[curFaceID];
        }
        else
        {
            neiCell = -1;
        }


        // Modify the face
        ref.setAction
        (
            polyModifyFace
            (
                newFace,                // modified face
                curFaceID,              // label of face being modified
                own[curFaceID],         // owner
                neiCell,                // neighbour
                false,                  // face flip
                patchID,                // patch for face
                false,                  // remove from zone
                modifiedFaceZone,       // zone for face
                modifiedFaceZoneFlip    // face flip in zone
            )
        );
    }

    if (debug)
    {
        Pout<< "void attachDetach::attachInterface("
            << "polyTopoChange& ref) const "
            << " for object " << name() << " : "
            << "Finished attaching interface" << endl;
    }
}


void Foam::attachDetach::modifyMotionPoints
(
    pointField& motionPoints
) const
{
    const Map<label>& removedPointMap = pointMatchMap();

    const labelList removedPoints = removedPointMap.toc();

    if (debug)
    {
        Pout<< "void attachDetach::modifyMotionPoints("
            << "pointField& motionPoints) const "
            << " for object " << name() << " : "
            << "Adjusting motion points." << endl;

        // Calculate the difference in motion point positions
        scalar pointDiff = 0;

        forAll(removedPoints, pointI)
        {
            pointDiff +=
                mag
                (
                    motionPoints[removedPoints[pointI]]
                  - motionPoints[removedPointMap.find(removedPoints[pointI])()]
                );
        }

        if (pointDiff > removedPoints.size()*positionDifference_)
        {
            Pout<< "Point motion difference = " << pointDiff << endl;
        }
    }

    // Put the slave point on top of the master point
    forAll(removedPoints, pointI)
    {
        motionPoints[removedPoints[pointI]] =
            motionPoints[removedPointMap.find(removedPoints[pointI])()];
    }

}


// ************************************************************************* //
