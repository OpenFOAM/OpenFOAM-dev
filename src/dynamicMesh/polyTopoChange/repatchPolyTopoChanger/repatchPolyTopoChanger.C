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

Description
    A mesh which allows changes in the patch distribution of the
    faces.  The change in patching is set using changePatchID.  For a
    boundary face, a new patch ID is given.  If the face is internal,
    it will be added to the first patch and its opposite to the second
    patch (take care with face orientation!).

\*---------------------------------------------------------------------------*/

#include "repatchPolyTopoChanger.H"
#include "polyTopoChanger.H"
#include "mapPolyMesh.H"
#include "polyModifyFace.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::polyTopoChange& Foam::repatchPolyTopoChanger::meshMod()
{
    if (meshModPtr_.empty())
    {
        meshModPtr_.reset(new polyTopoChange(mesh_));
    }
    return meshModPtr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::repatchPolyTopoChanger::repatchPolyTopoChanger(polyMesh& mesh)
:
    mesh_(mesh),
    meshModPtr_(NULL)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::repatchPolyTopoChanger::changePatches
(
    const List<polyPatch*>& patches
)
{
    if (meshModPtr_.valid())
    {
        FatalErrorIn
        (
            "repatchPolyTopoChanger::changePatches(const List<polyPatch*>&)"
        )   << "Cannot change patches after having changed faces. " << nl
            << "Please call changePatches first."
            << exit(FatalError);
    }
    meshModPtr_.clear();
    mesh_.removeBoundary();
    mesh_.addPatches(patches);
}


void Foam::repatchPolyTopoChanger::changePatchID
(
    const label faceID,
    const label patchID
)
{
    if (polyTopoChanger::debug)
    {
        // Check that the request is possible
        if
        (
            faceID >= mesh_.faces().size()
         || patchID >= mesh_.boundaryMesh().size()
         || mesh_.isInternalFace(faceID)
        )
        {
            FatalErrorIn
            (
                "void Foam::repatchPolyTopoChanger::changePatchID\n"
                "(\n"
                "    const label faceID,\n"
                "    const label patchID\n"
                ")\n"
            )   << "Error in definition.  faceID: " << faceID
                << " patchID: " << patchID << ".  "
                << "Labels out of range or internal face."
                << abort(FatalError);
        }
    }

    const label zoneID = mesh_.faceZones().whichZone(faceID);

    bool zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh_.faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(faceID)];
    }

    meshMod().setAction
    (
        polyModifyFace
        (
            mesh_.faces()[faceID],              // face
            faceID,                             // face ID
            mesh_.faceOwner()[faceID],          // owner
            -1,                                 // neighbour
            false,                              // flip flux
            patchID,                            // patch ID
            false,                              // remove from zone
            zoneID,                             // zone ID
            zoneFlip                            // zone flip
        )
    );
}


void Foam::repatchPolyTopoChanger::setFaceZone
(
    const label faceID,
    const label zoneID,
    const bool zoneFlip
)
{
    if (polyTopoChanger::debug)
    {
        // Check that the request is possible
        if (faceID > mesh_.faces().size())
        {
            FatalErrorIn
            (
                "void Foam::repatchPolyTopoChanger::setFaceZone"
                "(\n"
                "    const label faceID,\n"
                "    const label zoneID,\n"
                "    const bool flip\n"
                ")\n"
            )   << "Error in definition.  faceID: " << faceID
                << "out of range."
                << abort(FatalError);
        }
    }

    meshMod().setAction
    (
        polyModifyFace
        (
            mesh_.faces()[faceID],              // face
            faceID,                             // face ID
            mesh_.faceOwner()[faceID],          // owner
            mesh_.faceNeighbour()[faceID],      // neighbour
            false,                              // flip flux
            mesh_.boundaryMesh().whichPatch(faceID),  // patch ID
            true,                               // remove from zone
            zoneID,                             // zone ID
            zoneFlip                            // zone flip
        )
    );
}


void Foam::repatchPolyTopoChanger::changeAnchorPoint
(
    const label faceID,
    const label fp
)
{
    if (polyTopoChanger::debug)
    {
        // Check that the request is possible
        if (faceID > mesh_.faces().size())
        {
            FatalErrorIn
            (
                "void Foam::repatchPolyTopoChanger::setFaceZone"
                "(\n"
                "    const label faceID,\n"
                "    const label zoneID,\n"
                "    const bool flip\n"
                ")\n"
            )   << "Error in definition.  faceID: " << faceID
                << "out of range."
                << abort(FatalError);
        }
    }

    const face& f = mesh_.faces()[faceID];

    if ((fp < 0) || (fp >= f.size()))
    {
        FatalErrorIn
        (
            "void Foam::repatchPolyTopoChanger::changeAnchorPoint"
            "(\n"
            "    const label faceID,\n"
            "    const label fp\n"
            ")\n"
        )   << "Error in definition.  Face point: " << fp
            << "indexes out of face " << f
            << abort(FatalError);
    }

    label patchID = mesh_.boundaryMesh().whichPatch(faceID);

    const label zoneID = mesh_.faceZones().whichZone(faceID);

    bool zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh_.faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(faceID)];
    }

    if (fp == 0)
    {
        // Do dummy modify to keep patch ordering.
        meshMod().setAction
        (
            polyModifyFace
            (
                f,                                  // face
                faceID,                             // face ID
                mesh_.faceOwner()[faceID],          // owner
                -1,                                 // neighbour
                false,                              // flip flux
                patchID,                            // patch ID
                false,                              // remove from zone
                zoneID,                             // zone ID
                zoneFlip                            // zone flip
            )
        );
    }
    else
    {
        // Construct new face with fp as first point.

        face newFace(f.size());

        label fVert = fp;

        forAll(f, i)
        {
            newFace[i] = f[fVert++];

            if (fVert == f.size())
            {
                fVert = 0;
            }
        }


        meshMod().setAction
        (
            polyModifyFace
            (
                newFace,                            // face
                faceID,                             // face ID
                mesh_.faceOwner()[faceID],          // owner
                -1,                                 // neighbour
                false,                              // flip flux
                patchID,                            // patch ID
                false,                              // remove from zone
                zoneID,                             // zone ID
                zoneFlip                            // zone flip
            )
        );
    }
}


void Foam::repatchPolyTopoChanger::repatch()
{
    // Change mesh, no inflation
    meshMod().changeMesh(mesh_, false);

    // Clear topo change for the next operation
    meshModPtr_.clear();
}


// ************************************************************************* //
