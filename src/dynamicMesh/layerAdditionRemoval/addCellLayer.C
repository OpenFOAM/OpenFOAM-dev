/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
#include "polyTopoChanger.H"
#include "polyAddPoint.H"
#include "polyAddCell.H"
#include "polyAddFace.H"
#include "polyModifyFace.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::layerAdditionRemoval::extrusionDir() const
{
    const polyMesh& mesh = topoChanger().mesh();
    const primitiveFacePatch& masterFaceLayer =
        mesh.faceZones()[faceZoneID_.index()]();

    const pointField& points = mesh.points();
    const labelList& mp = masterFaceLayer.meshPoints();

    tmp<vectorField> textrusionDir(new vectorField(mp.size()));
    vectorField& extrusionDir = textrusionDir();

    if (setLayerPairing())
    {
        if (debug)
        {
            Pout<< "void layerAdditionRemoval::extrusionDir() const "
                << " for object " << name() << " : "
                << "Using edges for point insertion" << endl;
        }

        // Detected a valid layer.  Grab the point and face collapse mapping
        const labelList& ptc = pointsPairing();

        forAll(extrusionDir, mpI)
        {
            extrusionDir[mpI] = points[ptc[mpI]] - points[mp[mpI]];
        }
    }
    else
    {
        if (debug)
        {
            Pout<< "void layerAdditionRemoval::extrusionDir() const "
                << " for object " << name() << " : "
                << "A valid layer could not be found in front of "
                << "the addition face layer.  Using face-based "
                << "point normals for point addition"
                << endl;
        }

        extrusionDir = minLayerThickness_*masterFaceLayer.pointNormals();
    }

    return textrusionDir;
}


void Foam::layerAdditionRemoval::addCellLayer
(
    polyTopoChange& ref
) const
{
    // Insert the layer addition instructions into the topological change

    // Algorithm:
    // 1.  For every point in the master zone add a new point
    // 2.  For every face in the master zone add a cell
    // 3.  Go through all the edges of the master zone.  For all internal edges,
    //     add a face with the correct orientation and make it internal.
    //     For all boundary edges, find the patch it belongs to and add the face
    //     to this patch
    // 4.  Create a set of new faces on the patch image and assign them to be
    //     between the old master cells and the newly created cells
    // 5.  Modify all the faces in the patch such that they are located
    //     between the old slave cells and newly created cells

    if (debug)
    {
        Pout<< "void layerAdditionRemoval::addCellLayer("
            << "polyTopoChange& ref) const for object " << name() << " : "
            << "Adding cell layer" << endl;
    }

    // Create the points

    const polyMesh& mesh = topoChanger().mesh();
    const primitiveFacePatch& masterFaceLayer =
        mesh.faceZones()[faceZoneID_.index()]();

    const pointField& points = mesh.points();
    const labelList& mp = masterFaceLayer.meshPoints();

    // Get the extrusion direction for the added points

    const vectorField pointOffsets(extrusionDir());

    // Add the new points
    labelList addedPoints(mp.size());

    forAll(mp, pointI)
    {
        // Add the point nominal thickness away from the master zone point
        // and grab the label
        addedPoints[pointI] =
            ref.setAction
            (
                polyAddPoint
                (
                    points[mp[pointI]]                  // point
                  + addDelta_*pointOffsets[pointI],
                    mp[pointI],                         // master point
                    -1,                                 // zone for point
                    true                                // supports a cell
                )
            );
    }

    if (debug > 1)
    {
        Pout<< "mp: " << mp << " addedPoints: " << addedPoints << endl;
    }

    // Create the cells

    const labelList& mc =
        mesh.faceZones()[faceZoneID_.index()].masterCells();
    const labelList& sc =
        mesh.faceZones()[faceZoneID_.index()].slaveCells();

    const labelList& mf = mesh.faceZones()[faceZoneID_.index()];
    const boolList& mfFlip = mesh.faceZones()[faceZoneID_.index()].flipMap();

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    labelList addedCells(mf.size());

    forAll(mf, faceI)
    {
        label cellI = mc[faceI];
        label zoneI =  mesh.cellZones().whichZone(cellI);

        addedCells[faceI] =
            ref.setAction
            (
                polyAddCell
                (
                    -1,           // master point
                    -1,           // master edge
                    mf[faceI],    // master face
                    -1,           // master cell
                    zoneI         // zone for cell
                )
            );
    }

    // Create the new faces

    // Grab the local faces from the master zone
    const faceList& zoneFaces = masterFaceLayer.localFaces();

    // Create the new faces by copying the master zone.  All the added
    // faces need to point into the newly created cells, which means
    // all the faces from the master layer are flipped.  The flux flip
    // is determined from the original master layer face and the face
    // owner: if the master cell is equal to the face owner the flux
    // remains the same; otherwise it is flipped

    forAll(zoneFaces, faceI)
    {
        const face oldFace = zoneFaces[faceI].reverseFace();

        face newFace(oldFace.size());

        forAll(oldFace, pointI)
        {
            newFace[pointI] = addedPoints[oldFace[pointI]];
        }

        bool flipFaceFlux = false;

        // Flip the face as necessary
        if
        (
           !mesh.isInternalFace(mf[faceI])
         || mc[faceI] == nei[mf[faceI]]
        )
        {
            flipFaceFlux = true;
        }

        // Add the face
        ref.setAction
        (
            polyAddFace
            (
                newFace,               // face
                mc[faceI],             // owner
                addedCells[faceI],     // neighbour
                -1,                    // master point
                -1,                    // master edge
                mf[faceI],             // master face for addition
                flipFaceFlux,          // flux flip
                -1,                    // patch for face
                -1,                    // zone for face
                false                  // face zone flip
            )
        );

        if (debug > 1)
        {
            Pout<< "adding face: " << newFace
                << " own: " << mc[faceI]
                << " nei: " << addedCells[faceI]
                << endl;
        }
    }

    // Modify the faces from the master zone for the new neighbour

    const faceList& faces = mesh.faces();

    // Pout<< "mfFlip: " << mfFlip << endl;

    forAll(mf, faceI)
    {
        const label curfaceID = mf[faceI];

        // If the face is internal, modify its owner to be the newly
        // created cell.  No flip is necessary
        if (!mesh.isInternalFace(curfaceID))
        {
            ref.setAction
            (
                polyModifyFace
                (
                    faces[curfaceID],            // modified face
                    curfaceID,                   // label of face being modified
                    addedCells[faceI],           // owner
                    -1,                          // neighbour
                    false,                       // face flip
                    mesh.boundaryMesh().whichPatch(curfaceID),// patch for face
                    false,                       // remove from zone
                    faceZoneID_.index(),         // zone for face
                    mfFlip[faceI]                // face flip in zone
                )
            );

            if (debug > 1)
            {
                Pout<< "Modifying a boundary face. Face: " << curfaceID
                    << " flip: " << mfFlip[faceI]
                    << endl;
            }
        }

        // If slave cell is owner, the face remains the same (but with
        // a new neighbour - the newly created cell).  Otherwise, the
        // face is flipped.
        else if (sc[faceI] == own[curfaceID])
        {
            // Orientation is good, just change neighbour
            ref.setAction
            (
                polyModifyFace
                (
                    faces[curfaceID],            // modified face
                    curfaceID,                   // label of face being modified
                    own[curfaceID],              // owner
                    addedCells[faceI],           // neighbour
                    false,                       // face flip
                    mesh.boundaryMesh().whichPatch(curfaceID),// patch for face
                    false,                       // remove from zone
                    faceZoneID_.index(),         // zone for face
                    mfFlip[faceI]                // face flip in zone
                )
            );

            if (debug > 1)
            {
                Pout<< "modify face, no flip " << curfaceID
                    << " own: " << own[curfaceID]
                    << " nei: " << addedCells[faceI]
                    << endl;
            }
        }
        else
        {
            // Change in orientation; flip face
            ref.setAction
            (
                polyModifyFace
                (
                    faces[curfaceID].reverseFace(), // modified face
                    curfaceID,                   // label of face being modified
                    nei[curfaceID],                 // owner
                    addedCells[faceI],              // neighbour
                    true,                           // face flip
                    mesh.boundaryMesh().whichPatch(curfaceID), // patch for face
                    false,                          // remove from zone
                    faceZoneID_.index(),            // zone for face
                    !mfFlip[faceI]                  // face flip in zone
                )
            );

            if (debug > 1)
            {
                Pout<< "modify face, with flip " << curfaceID
                    << " own: " << own[curfaceID]
                    << " nei: " << addedCells[faceI]
                    << endl;
            }
        }
    }

    // Create new faces from edges

    const edgeList& zoneLocalEdges = masterFaceLayer.edges();

    const labelListList& edgeFaces = masterFaceLayer.edgeFaces();

    label nInternalEdges = masterFaceLayer.nInternalEdges();

    const labelList& meshEdges =
        mesh.faceZones()[faceZoneID_.index()].meshEdges();

    // Do all internal edges
    for (label curEdgeID = 0; curEdgeID < nInternalEdges; curEdgeID++)
    {
        face newFace(4);

        newFace[0] = mp[zoneLocalEdges[curEdgeID].start()];
        newFace[1] = mp[zoneLocalEdges[curEdgeID].end()];
        newFace[2] = addedPoints[zoneLocalEdges[curEdgeID].end()];
        newFace[3] = addedPoints[zoneLocalEdges[curEdgeID].start()];

        ref.setAction
        (
            polyAddFace
            (
                newFace,                               // face
                addedCells[edgeFaces[curEdgeID][0]],   // owner
                addedCells[edgeFaces[curEdgeID][1]],   // neighbour
                -1,                                    // master point
                meshEdges[curEdgeID],                  // master edge
                -1,                                    // master face
                false,                                 // flip flux
                -1,                                    // patch for face
                -1,                                    // zone for face
                false                                  // face zone flip
            )
        );

        if (debug > 1)
        {
            Pout<< "Add internal face off edge: " << newFace
                << " own: " << addedCells[edgeFaces[curEdgeID][0]]
                << " nei: " << addedCells[edgeFaces[curEdgeID][1]]
                << endl;
        }
    }

    // Prepare creation of faces from boundary edges.
    // Note:
    // A tricky part of the algorithm is finding the patch into which the
    // newly created face will be added.  For this purpose, take the edge
    // and grab all the faces coming from it.  From the set of faces
    // eliminate the internal faces and find the boundary face from the rest.
    //  If there is more than one boundary face (excluding the ones in
    // the master zone), the patch with the lower label is selected.

    const labelListList& meshEdgeFaces = mesh.edgeFaces();

    const faceZoneMesh& zoneMesh = mesh.faceZones();

    // Do all boundary edges

    for
    (
        label curEdgeID = nInternalEdges;
        curEdgeID < zoneLocalEdges.size();
        curEdgeID++
    )
    {
        face newFace(4);
        newFace[0] = mp[zoneLocalEdges[curEdgeID].start()];
        newFace[1] = mp[zoneLocalEdges[curEdgeID].end()];
        newFace[2] = addedPoints[zoneLocalEdges[curEdgeID].end()];
        newFace[3] = addedPoints[zoneLocalEdges[curEdgeID].start()];

        // Determine the patch for insertion
        const labelList& curFaces = meshEdgeFaces[meshEdges[curEdgeID]];

        label patchID = -1;
        label zoneID = -1;

        forAll(curFaces, faceI)
        {
            const label cf = curFaces[faceI];

            if (!mesh.isInternalFace(cf))
            {
                // Face not internal. Check if it is in the zone
                if (zoneMesh.whichZone(cf) != faceZoneID_.index())
                {
                    // Found the face in a boundary patch which is not in zone
                    patchID = mesh.boundaryMesh().whichPatch(cf);
                    zoneID = mesh.faceZones().whichZone(cf);

                    break;
                }
            }
        }

        if (patchID < 0)
        {
            FatalErrorIn
            (
                "void Foam::layerAdditionRemoval::setRefinement("
                "polyTopoChange& ref)"
            )   << "Cannot find patch for edge " << meshEdges[curEdgeID]
                << ". Edge: " << mesh.edges()[meshEdges[curEdgeID]]
                << abort(FatalError);
        }

        ref.setAction
        (
            polyAddFace
            (
                newFace,                               // face
                addedCells[edgeFaces[curEdgeID][0]],   // owner
                -1,                                    // neighbour
                -1,                                    // master point
                meshEdges[curEdgeID],                  // master edge
                -1,                                    // master face
                false,                                 // flip flux
                patchID,                               // patch for face
                zoneID,                                // zone
                false                                  // zone face flip
            )
        );

        if (debug > 1)
        {
            Pout<< "add boundary face: " << newFace
                << " into patch " << patchID
                << " own: " << addedCells[edgeFaces[curEdgeID][0]]
                << endl;
        }
    }

    // Modify the remaining faces of the master cells to reconnect to the new
    // layer of faces.
    // Algorithm: Go through all the cells of the master zone and make
    // a map of faces to avoid duplicates.  Do not insert the faces in
    // the master patch (as they have already been dealt with).  Make
    // a master layer point renumbering map, which for every point in
    // the master layer gives its new label. Loop through all faces in
    // the map and attempt to renumber them using the master layer
    // point renumbering map.  Once the face is renumbered, compare it
    // with the original face; if they are the same, the face has not
    // changed; if not, modify the face but keep all of its old
    // attributes (apart from the vertex numbers).

    // Create the map of faces in the master cell layer
    labelHashSet masterCellFaceMap(primitiveMesh::facesPerCell_*mc.size());

    const cellList& cells = mesh.cells();

    forAll(mc, cellI)
    {
        const labelList& curFaces = cells[mc[cellI]];

        forAll(curFaces, faceI)
        {
            // Check if the face belongs to the master zone; if not add it
            if (zoneMesh.whichZone(curFaces[faceI]) != faceZoneID_.index())
            {
                masterCellFaceMap.insert(curFaces[faceI]);
            }
        }
    }

    // Create the master layer point map
    Map<label> masterLayerPointMap(2*mp.size());

    forAll(mp, pointI)
    {
        masterLayerPointMap.insert
        (
            mp[pointI],
            addedPoints[pointI]
        );
    }

    // Grab the list of faces of the master layer
    const labelList masterCellFaces = masterCellFaceMap.toc();

    forAll(masterCellFaces, faceI)
    {
        // Attempt to renumber the face using the masterLayerPointMap.
        // Missing point remain the same

        const label curFaceID = masterCellFaces[faceI];

        const face& oldFace = faces[curFaceID];

        face newFace(oldFace.size());

        bool changed = false;

        forAll(oldFace, pointI)
        {
            if (masterLayerPointMap.found(oldFace[pointI]))
            {
                changed = true;

                newFace[pointI] = masterLayerPointMap.find(oldFace[pointI])();
            }
            else
            {
                newFace[pointI] = oldFace[pointI];
            }
        }

        // If the face has changed, create a modification entry
        if (changed)
        {
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

            if (mesh.isInternalFace(curFaceID))
            {
                ref.setAction
                (
                    polyModifyFace
                    (
                        newFace,                // modified face
                        curFaceID,              // label of face being modified
                        own[curFaceID],         // owner
                        nei[curFaceID],         // neighbour
                        false,                  // face flip
                        -1,                     // patch for face
                        false,                  // remove from zone
                        modifiedFaceZone,       // zone for face
                        modifiedFaceZoneFlip    // face flip in zone
                    )
                );

                if (debug > 1)
                {
                    Pout<< "modifying stick-out face. Internal Old face: "
                        << oldFace
                        << " new face: " << newFace
                        << " own: " << own[curFaceID]
                        << " nei: " << nei[curFaceID]
                        << endl;
                }
            }
            else
            {
                ref.setAction
                (
                    polyModifyFace
                    (
                        newFace,                // modified face
                        curFaceID,              // label of face being modified
                        own[curFaceID],         // owner
                        -1,                     // neighbour
                        false,                  // face flip
                        mesh.boundaryMesh().whichPatch(curFaceID),
                                                // patch for face
                        false,                  // remove from zone
                        modifiedFaceZone,       // zone for face
                        modifiedFaceZoneFlip    // face flip in zone
                    )
                );

                if (debug > 1)
                {
                    Pout<< "modifying stick-out face. Boundary Old face: "
                        << oldFace
                        << " new face: " << newFace
                        << " own: " << own[curFaceID]
                        << " patch: "
                        << mesh.boundaryMesh().whichPatch(curFaceID)
                        << endl;
                }
            }
        }
    }

    if (debug)
    {
        Pout<< "void layerAdditionRemoval::addCellLayer(polyTopoChange&) const "
            << " for object " << name() << ": "
            << "Finished adding cell layer" << endl;
    }
}


// ************************************************************************* //
