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

Application
    mergeBaffles

Description
    Detects faces that share points (baffles) and merge them into internal
    faces.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "syncTools.H"
#include "faceSet.H"
#include "pointSet.H"
#include "meshTools.H"
#include "polyTopoChange.H"
#include "indirectPrimitivePatch.H"
#include "processorPolyPatch.H"
#include "localPointRegion.H"
#include "duplicatePoints.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void mergeDuplicateBoundaryFaces
(
    const polyMesh& mesh,
    polyTopoChange& meshMod
)
{
    // Get all duplicate face labels in the boundary
    labelList duplicates = localPointRegion::findDuplicateFaces
    (
        mesh,
        identityMap(mesh.nFaces() - mesh.nInternalFaces())
      + mesh.nInternalFaces()
    );

    // Check that none are on processor patches
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    forAll(duplicates, bFacei)
    {
        if (duplicates[bFacei] != -1)
        {
            label facei = mesh.nInternalFaces() + bFacei;
            label patchi = patches.whichPatch(facei);

            if (isA<processorPolyPatch>(patches[patchi]))
            {
                FatalErrorInFunction
                    << "Duplicate face " << facei
                    << " is on a processorPolyPatch."
                    << "This is not allowed." << nl
                    << "Face:" << facei
                    << " is on patch:" << patches[patchi].name()
                    << abort(FatalError);
            }
        }
    }

    const faceList& faces = mesh.faces();
    const labelList& faceOwner = mesh.faceOwner();

    forAll(duplicates, bFacei)
    {
        label otherFacei = duplicates[bFacei];

        if (otherFacei != -1 && otherFacei > bFacei)
        {
            // Two duplicate faces. Merge.

            label face0 = mesh.nInternalFaces() + bFacei;
            label face1 = mesh.nInternalFaces() + otherFacei;

            label own0 = faceOwner[face0];
            label own1 = faceOwner[face1];

            if (own0 < own1)
            {
                // Use face0 as the new internal face.

                meshMod.removeFace(face1, -1);
                meshMod.modifyFace
                (
                    faces[face0],           // modified face
                    face0,                  // label of face being modified
                    own0,                   // owner
                    own1,                   // neighbour
                    false,                  // face flip
                    -1                      // patch for face
                );
            }
            else
            {
                // Use face1 as the new internal face.

                meshMod.removeFace(face0, -1);
                meshMod.modifyFace
                (
                    faces[face1],           // modified face
                    face1,                  // label of face being modified
                    own1,                   // owner
                    own0,                   // neighbour
                    false,                  // face flip
                    -1                      // patch for face
                );
            }
        }
    }
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Detect faces that share points (baffles)\n"
        "and merge them into internal faces."
    );

    #include "addOverwriteOption.H"
    #include "addMeshOption.H"
    #include "addRegionOption.H"
    argList::addBoolOption
    (
        "fields",
        "update fields"
    );

    #include "setRootCase.H"
    #include "createTimeNoFunctionObjects.H"
    #include "createSpecifiedMeshNoChangers.H"

    const bool overwrite  = args.optionFound("overwrite");
    const bool fields     = args.optionFound("fields");

    const word oldInstance = mesh.pointsInstance();

    // Read objects in time directory
    IOobjectList objects(mesh, runTime.name());

    if (fields) Info<< "Reading geometric fields" << nl << endl;

    #include "readVolFields.H"
    #include "readSurfaceFields.H"
    #include "readPointFields.H"

    // Mesh change engine
    polyTopoChange meshMod(mesh);

    // Merge duplicate boundary faces into internal faces.
    mergeDuplicateBoundaryFaces(mesh, meshMod);

    if (!overwrite)
    {
        runTime++;
    }

    // Change the mesh
    autoPtr<polyTopoChangeMap> map = meshMod.changeMesh(mesh);

    // Update fields
    mesh.topoChange(map);

    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }

    Info<< "Writing mesh to time " << runTime.name() << endl;
    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
