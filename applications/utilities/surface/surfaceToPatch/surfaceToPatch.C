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

Application
    surfaceToPatch

Description
    Reads surface and applies surface regioning to a mesh. Uses boundaryMesh
    to do the hard work.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "boundaryMesh.H"
#include "polyMesh.H"
#include "faceSet.H"
#include "polyTopoChange.H"
#include "polyModifyFace.H"
#include "globalMeshData.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Adds empty patch if not yet there. Returns patchID.
label addPatch(polyMesh& mesh, const word& patchName)
{
    label patchI = mesh.boundaryMesh().findPatchID(patchName);

    if (patchI == -1)
    {
        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        List<polyPatch*> newPatches(patches.size() + 1);

        patchI = 0;

        // Copy all old patches
        forAll(patches, i)
        {
            const polyPatch& pp = patches[i];

            newPatches[patchI] =
                pp.clone
                (
                    patches,
                    patchI,
                    pp.size(),
                    pp.start()
                ).ptr();

            patchI++;
        }

        // Add zero-sized patch
        newPatches[patchI] =
            new polyPatch
            (
                patchName,
                0,
                mesh.nFaces(),
                patchI,
                patches,
                polyPatch::typeName
            );

        mesh.removeBoundary();
        mesh.addPatches(newPatches);

        Pout<< "Created patch " << patchName << " at " << patchI << endl;
    }
    else
    {
        Pout<< "Reusing patch " << patchName << " at " << patchI << endl;
    }

    return patchI;
}


// Repatch single face. Return true if patch changed.
bool repatchFace
(
    const polyMesh& mesh,
    const boundaryMesh& bMesh,
    const labelList& nearest,
    const labelList& surfToMeshPatch,
    const label faceI,
    polyTopoChange& meshMod
)
{
    bool changed = false;

    label bFaceI = faceI - mesh.nInternalFaces();

    if (nearest[bFaceI] != -1)
    {
        // Use boundary mesh one.
        label bMeshPatchID = bMesh.whichPatch(nearest[bFaceI]);

        label patchID = surfToMeshPatch[bMeshPatchID];

        if (patchID != mesh.boundaryMesh().whichPatch(faceI))
        {
            label own = mesh.faceOwner()[faceI];

            label zoneID = mesh.faceZones().whichZone(faceI);

            bool zoneFlip = false;

            if (zoneID >= 0)
            {
                const faceZone& fZone = mesh.faceZones()[zoneID];

                zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
            }

            meshMod.setAction
            (
                polyModifyFace
                (
                    mesh.faces()[faceI],// modified face
                    faceI,              // label of face being modified
                    own,                // owner
                    -1,                 // neighbour
                    false,              // face flip
                    patchID,            // patch for face
                    false,              // remove from zone
                    zoneID,             // zone for face
                    zoneFlip            // face flip in zone
                )
            );

            changed = true;
        }
    }
    else
    {
        changed = false;
    }
    return changed;
}



int main(int argc, char *argv[])
{
    argList::addNote
    (
        "reads surface and applies surface regioning to a mesh"
    );

    argList::noParallel();
    argList::validArgs.append("surfaceFile");
    argList::addOption
    (
        "faceSet",
        "name",
        "only repatch the faces in specified faceSet"
    );
    argList::addOption
    (
        "tol",
        "scalar",
        "search tolerance as fraction of mesh size (default 1e-3)"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    const fileName surfName = args[1];

    Info<< "Reading surface from " << surfName << " ..." << endl;

    word setName;
    const bool readSet = args.optionReadIfPresent("faceSet", setName);

    if (readSet)
    {
        Info<< "Repatching only the faces in faceSet " << setName
            << " according to nearest surface triangle ..." << endl;
    }
    else
    {
        Info<< "Patching all boundary faces according to nearest surface"
            << " triangle ..." << endl;
    }

    const scalar searchTol = args.optionLookupOrDefault("tol", 1e-3);

    // Get search box. Anything not within this box will not be considered.
    const boundBox& meshBb = mesh.bounds();
    const vector searchSpan = searchTol * meshBb.span();

    Info<< "All boundary faces further away than " << searchTol
        << " of mesh bounding box " << meshBb
        << " will keep their patch label ..." << endl;


    Info<< "Before patching:" << nl
        << "    patch\tsize" << endl;

    forAll(mesh.boundaryMesh(), patchI)
    {
        Info<< "    " << mesh.boundaryMesh()[patchI].name() << '\t'
            << mesh.boundaryMesh()[patchI].size() << nl;
    }
    Info<< endl;


    boundaryMesh bMesh;

    // Load in the surface.
    bMesh.readTriSurface(surfName);

    // Add all the boundaryMesh patches to the mesh.
    const PtrList<boundaryPatch>& bPatches = bMesh.patches();

    // Map from surface patch ( = boundaryMesh patch) to polyMesh patch
    labelList patchMap(bPatches.size());

    forAll(bPatches, i)
    {
        patchMap[i] = addPatch(mesh, bPatches[i].name());
    }

    // Obtain nearest face in bMesh for each boundary face in mesh that
    // is within search span.
    // Note: should only determine for faceSet if working with that.
    labelList nearest(bMesh.getNearest(mesh, searchSpan));

    {
        // Dump unmatched faces to faceSet for debugging.
        faceSet unmatchedFaces(mesh, "unmatchedFaces", nearest.size()/100);

        forAll(nearest, bFaceI)
        {
            if (nearest[bFaceI] == -1)
            {
                unmatchedFaces.insert(mesh.nInternalFaces() + bFaceI);
            }
        }

        Pout<< "Writing all " << unmatchedFaces.size()
            << " unmatched faces to faceSet "
            << unmatchedFaces.name()
            << endl;

        unmatchedFaces.write();
    }


    polyTopoChange meshMod(mesh);

    label nChanged = 0;

    if (readSet)
    {
        faceSet faceLabels(mesh, setName);
        Info<< "Read " << faceLabels.size() << " faces to repatch ..." << endl;

        forAllConstIter(faceSet, faceLabels, iter)
        {
            label faceI = iter.key();

            if (repatchFace(mesh, bMesh, nearest, patchMap, faceI, meshMod))
            {
                nChanged++;
            }
        }
    }
    else
    {
        forAll(nearest, bFaceI)
        {
            label faceI = mesh.nInternalFaces() + bFaceI;

            if (repatchFace(mesh, bMesh, nearest, patchMap, faceI, meshMod))
            {
                nChanged++;
            }
        }
    }

    Pout<< "Changed " << nChanged << " boundary faces." << nl << endl;

    if (nChanged > 0)
    {
        meshMod.changeMesh(mesh, false);

        Info<< "After patching:" << nl
            << "    patch\tsize" << endl;

        forAll(mesh.boundaryMesh(), patchI)
        {
            Info<< "    " << mesh.boundaryMesh()[patchI].name() << '\t'
                << mesh.boundaryMesh()[patchI].size() << endl;
        }
        Info<< endl;


        runTime++;

        // Write resulting mesh
        Info<< "Writing modified mesh to time " << runTime.value() << endl;
        mesh.write();
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
