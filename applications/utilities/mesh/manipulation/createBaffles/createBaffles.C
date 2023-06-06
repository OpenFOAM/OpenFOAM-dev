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

Application
    createBaffles

Description
    Makes faces into boundary faces. Does not duplicate points.

    Notes:
    - If any coupled patch face is selected for baffling the opposite member
      has to be selected for baffling as well.
    - If fields are being modified then boundary conditions must be specified
      within the patch's patchFields sub-dictionary.
    - Any patches left with no faces will be removed, except for patches of
      coupled, internal, of non-conformal type

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyTopoChange.H"
#include "polyModifyFace.H"
#include "polyAddFace.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "fvMeshMapper.H"
#include "faceSelection.H"
#include "searchableSurface.H"
#include "fvMeshTools.H"
#include "systemDict.H"
#include "processorPolyPatch.H"
#include "internalPolyPatch.H"
#include "nonConformalPolyPatch.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void filterPatches(fvMesh& mesh, const HashSet<word>& bafflePatches)
{
    // Remove any zero-sized patches, except for constraint types, and any
    // specified in the system dictionary

    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    label newPatchi = 0;

    // List of patch's new indices
    labelList oldToNew(bMesh.size(), -1);

    // Add all the kept non-processor patches to the list first
    forAll(bMesh, patchi)
    {
        const polyPatch& pp = bMesh[patchi];

        if (!isA<processorPolyPatch>(pp))
        {
            if
            (
                bafflePatches.found(pp.name())
             || fvPatch::constraintType(pp.type())
             || returnReduce(pp.size(), sumOp<label>())
            )
            {
                oldToNew[patchi] = newPatchi++;
            }
        }
    }

    // Now add all the processor patches
    forAll(bMesh, patchi)
    {
        const polyPatch& pp = bMesh[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            oldToNew[patchi] = newPatchi++;
        }
    }

    // Note how many patches are to be kept
    const label nKeepPatches = newPatchi;

    // If there are any to remove ...
    if (nKeepPatches != bMesh.size())
    {
        Info<< endl << "Removing zero-sized patches:" << endl << incrIndent;

        // Add all the removed patches to the list
        forAll(oldToNew, patchi)
        {
            if (oldToNew[patchi] == -1)
            {
                Info<< indent << bMesh[patchi].name()
                    << " type " << bMesh[patchi].type()
                    << " at position " << patchi << endl;

                oldToNew[patchi] = newPatchi++;
            }
        }

        Info<< decrIndent;

        // Permute the mesh, keeping only the ones not to be removed
        fvMeshTools::reorderPatches(mesh, oldToNew, nKeepPatches, true);

        Info<< endl;
    }
}


void modifyOrAddFace
(
    polyTopoChange& meshMod,
    const face& f,
    const label facei,
    const label own,
    const bool flipFaceFlux,
    const label newPatchi,
    const label zoneID,
    const bool zoneFlip,
    PackedBoolList& modifiedFace
)
{
    if (!modifiedFace[facei])
    {
        // First usage of face. Modify.
        meshMod.setAction
        (
            polyModifyFace
            (
                f,                          // modified face
                facei,                      // label of face
                own,                        // owner
                -1,                         // neighbour
                flipFaceFlux,               // face flip
                newPatchi,                  // patch for face
                false,                      // remove from zone
                zoneID,                     // zone for face
                zoneFlip                    // face flip in zone
            )
        );

        modifiedFace[facei] = 1;
    }
    else
    {
        // Second usage of face. Add.
        meshMod.setAction
        (
            polyAddFace
            (
                f,                          // modified face
                own,                        // owner
                -1,                         // neighbour
                -1,                         // master point
                -1,                         // master edge
                facei,                      // master face
                flipFaceFlux,               // face flip
                newPatchi,                  // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face flip in zone
            )
        );
    }
}


label createFaces
(
    const bool internalFacesOnly,
    const fvMesh& mesh,
    const faceZone& fZone,
    const label newOwnerPatchi,
    const label newNeighbourPatchi,
    polyTopoChange& meshMod,
    PackedBoolList& modifiedFace
)
{
    label nModified = 0;

    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    // Pass 1. Do selected side of zone
    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        const label zoneFacei = fZone.whichFace(facei);

        if (zoneFacei != -1)
        {
            if (!fZone.flipMap()[zoneFacei])
            {
                // Use owner side of face
                modifyOrAddFace
                (
                    meshMod,
                    mesh.faces()[facei],    // modified face
                    facei,                  // label of face
                    mesh.faceOwner()[facei],// owner
                    false,                  // face flip
                    newOwnerPatchi,         // patch for face
                    fZone.index(),          // zone for face
                    false,                  // face flip in zone
                    modifiedFace            // modify or add status
                );
            }
            else
            {
                // Use neighbour side of face.
                // To keep faceZone pointing out of original neighbour
                // we don't need to set faceFlip since that cell
                // now becomes the owner
                modifyOrAddFace
                (
                    meshMod,
                    mesh.faces()[facei].reverseFace(), // modified face
                    facei,                      // label of face
                    mesh.faceNeighbour()[facei],// owner
                    true,                       // face flip
                    newOwnerPatchi,             // patch for face
                    fZone.index(),              // zone for face
                    false,                      // face flip in zone
                    modifiedFace                // modify or add status
                );
            }

            nModified++;
        }
    }

    // Pass 2. Do other side of zone
    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        label zoneFacei = fZone.whichFace(facei);

        if (zoneFacei != -1)
        {
            if (!fZone.flipMap()[zoneFacei])
            {
                // Use neighbour side of face
                modifyOrAddFace
                (
                    meshMod,
                    mesh.faces()[facei].reverseFace(),  // modified face
                    facei,                          // label of face
                    mesh.faceNeighbour()[facei],    // owner
                    true,                           // face flip
                    newNeighbourPatchi,             // patch for face
                    fZone.index(),                  // zone for face
                    true,                           // face flip in zone
                    modifiedFace                    // modify or add
                );
            }
            else
            {
                // Use owner side of face
                modifyOrAddFace
                (
                    meshMod,
                    mesh.faces()[facei],    // modified face
                    facei,                  // label of face
                    mesh.faceOwner()[facei],// owner
                    false,                  // face flip
                    newNeighbourPatchi,     // patch for face
                    fZone.index(),          // zone for face
                    true,                   // face flip in zone
                    modifiedFace            // modify or add status
                );
            }
        }
    }

    // Modify any boundary faces
    forAll(bMesh, patchi)
    {
        const polyPatch& pp = bMesh[patchi];

        if (pp.coupled() || !internalFacesOnly)
        {
            forAll(pp, i)
            {
                const label facei = pp.start() + i;
                const label zoneFacei = fZone.whichFace(facei);

                if (zoneFacei != -1)
                {
                    const polyPatch& newPp =
                        fZone.flipMap()[zoneFacei]
                      ? bMesh[newNeighbourPatchi]
                      : bMesh[newOwnerPatchi];

                    // We cannot move coupled faces to different coupled
                    // faces. Generate an error if this is attempted.
                    if (pp.coupled() && newPp.coupled())
                    {
                        FatalErrorInFunction
                            << "Face on coupled patch \"" << pp.name()
                            << "\" selected for conversion to coupled "
                            << "patch \"" << newPp.name() << "\". "
                            << "Conversion from coupled patch to coupled "
                            << "patch is not allowed."
                            << exit(FatalError);
                    }

                    modifyOrAddFace
                    (
                        meshMod,
                        mesh.faces()[facei],        // modified face
                        facei,                      // label of face
                        mesh.faceOwner()[facei],    // owner
                        false,                      // face flip
                        newPp.index(),              // patch for face
                        fZone.index(),              // zone for face
                        fZone.flipMap()[zoneFacei], // face flip in zone
                        modifiedFace                // modify or add
                    );

                    nModified++;
                }
            }
        }
    }

    return returnReduce(nModified, sumOp<label>());
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Makes internal faces into boundary faces.\n"
        "Does not duplicate points."
    );
    #include "addDictOption.H"
    #include "addOverwriteOption.H"
    #include "addDictOption.H"
    #include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTimeNoFunctionObjects.H"
    #include "createNamedMesh.H"

    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    const bool overwrite = args.optionFound("overwrite");

    const word oldInstance = mesh.pointsInstance();

    // Read the system dictionary and global controls
    const dictionary dict(systemDict("createBafflesDict", args, mesh));
    Switch internalFacesOnly = dict.lookup<Switch>("internalFacesOnly");
    const Switch fields = dict.lookupOrDefault("fields", false);

    // Create face selections
    PtrList<faceSelection> selectors;
    {
        const dictionary& selectionsDict = dict.subDict("baffles");

        label n = 0;
        forAllConstIter(dictionary, selectionsDict, iter)
        {
            if (iter().isDict())
            {
                n++;
            }
        }
        selectors.setSize(n);

        n = 0;
        forAllConstIter(dictionary, selectionsDict, iter)
        {
            if (iter().isDict())
            {
                selectors.set
                (
                    n++,
                    faceSelection::New(iter().keyword(), mesh, iter().dict())
                );
            }
        }
    }

    if (internalFacesOnly)
    {
        Info<< "Not converting faces on non-coupled patches." << nl << endl;
    }

    // Read fields
    IOobjectList objects(mesh, runTime.name());
    if (fields) Info<< "Reading geometric fields" << nl << endl;
    #include "readVolFields.H"
    #include "readSurfaceFields.H"
    #include "readPointFields.H"
    if (fields) Info<< endl;

    // Creating faceZones for selectors that are not already faceZones
    forAll(selectors, selectorI)
    {
        const word& name = selectors[selectorI].name();

        if (mesh.faceZones().findZoneID(name) == -1)
        {
            mesh.faceZones().clearAddressing();

            const label sz = mesh.faceZones().size();

            labelList addr(0);
            boolList flip(0);

            mesh.faceZones().setSize(sz+1);
            mesh.faceZones().set
            (
                sz,
                new faceZone(name, addr, flip, sz, mesh.faceZones())
            );
        }
    }

    // Per face, its zone ID and flip status
    labelList faceToZoneID(mesh.nFaces(), -1);
    boolList faceToFlip(mesh.nFaces(), false);
    forAll(selectors, selectorI)
    {
        const word& name = selectors[selectorI].name();
        label zoneID = mesh.faceZones().findZoneID(name);

        selectors[selectorI].select(zoneID, faceToZoneID, faceToFlip);
    }

    // Add faces to faceZones
    labelList nFaces(mesh.faceZones().size(), 0);
    forAll(faceToZoneID, facei)
    {
        label zoneID = faceToZoneID[facei];
        if (zoneID != -1)
        {
            nFaces[zoneID]++;
        }
    }
    forAll(selectors, selectorI)
    {
        const word& name = selectors[selectorI].name();
        const label zoneID = mesh.faceZones().findZoneID(name);

        label& n = nFaces[zoneID];
        labelList addr(n);
        boolList flip(n);

        n = 0;
        forAll(faceToZoneID, facei)
        {
            label zone = faceToZoneID[facei];
            if (zone == zoneID)
            {
                addr[n] = facei;
                flip[n] = faceToFlip[facei];
                n++;
            }
        }

        Info<< "Created zone '" << name << "' with "
            << returnReduce(n, sumOp<label>()) << " faces" << endl;

        mesh.faceZones().set
        (
            zoneID,
            new faceZone(name, addr, flip, zoneID, mesh.faceZones())
        );
    }

    Info<< endl;

    // Keywords for the two patch entries in each selector dictionary. Note
    // "owner" and "neighbour" are preferred; "master" and "slave" are still
    // permitted for backwards compatibility.
    const FixedList<wordList, 2> patchEntryKeywords
    ({
        wordList({"owner", "master"}),
        wordList({"neighbour", "slave"})
    });

    // Create the baffles. Notes:
    //
    // - This is done in multiple steps
    //   - Patches are created with 'calculated' patchFields
    //   - Faces are moved into these patches
    //   - The patchFields are changed to the desired type
    // - This helps resolve ordering issues
    //   - patchFields cannot be created at the same time as patches since a
    //     coupled patchField constructor frequently needs access to the
    //     neighbouring patch
    //   - patchFields need to be created after the faces have been added to
    //     the patch so that they don't later have to be mapped from nothing
    //

    // Add patches
    HashSet<word> bafflePatches;
    forAll(selectors, selectorI)
    {
        const word& groupName = selectors[selectorI].name();
        const dictionary& dict =
            selectors[selectorI].dict().optionalSubDict("patches");

        forAll(patchEntryKeywords, i)
        {
            dictionary patchDict =
                dict.lookupEntryBackwardsCompatible
                (
                    patchEntryKeywords[i],
                    false,
                    false
                ).dict();

            const word patchName(patchDict.lookup<word>("name"));

            if (bMesh.findPatchID(patchName) == -1)
            {
                Info<< "Adding patch '" << patchName << "' to the mesh" << endl;

                // Create the patch
                patchDict.set("nFaces", 0);
                patchDict.set("startFace", 0);
                autoPtr<polyPatch> ppPtr
                (
                    polyPatch::New
                    (
                        patchName,
                        patchDict,
                        0,
                        bMesh
                    )
                );
                polyPatch& pp = ppPtr();

                // Add it to the group
                if (!groupName.empty() && !pp.inGroup(groupName))
                {
                    pp.inGroups().append(groupName);
                }

                // Add the patch to the mesh, creating calculated boundary
                // conditions everywhere. These will be changed later.
                fvMeshTools::addPatch
                (
                    mesh,
                    pp,
                    dictionary(),   // do not set specialised patchFields
                    calculatedFvPatchField<scalar>::typeName,
                    true            // parallel sync'ed addition
                );
            }
            else
            {
                Info<< "Patch '" << patchName
                    << "' already exists in the mesh" << endl;

                patchDict.remove("name");
                patchDict.remove("patchFields");

                if (!patchDict.empty())
                {
                    WarningInFunction
                        << "Patch settings found in " << patchDict.name()
                        << " for patch '" << patchName << "', but this patch "
                        << "already exists so these settings will not be used"
                        << endl;
                }
            }

            bafflePatches.insert(patchName);
        }
    }

    Info<< endl;

    // Make sure patches and zoneFaces are synchronised across couples
    mesh.boundaryMesh().checkParallelSync(true);
    mesh.faceZones().checkParallelSync(true);

    // Do the topology changes. Notes:
    //
    // - Loop in incrementing face order (not necessary if faceZone ordered).
    //   Preserves any existing ordering on patch faces.
    // - Two passes, do non-flip faces first and flip faces second. This
    //   guarantees that when e.g. creating a cyclic all faces from one
    //   side come first and faces from the other side next.
    //
    polyTopoChange meshMod(mesh);
    {
        PackedBoolList modifiedFace(mesh.nFaces());

        forAll(selectors, selectorI)
        {
            const word& groupName = selectors[selectorI].name();

            const faceZone& fZone = mesh.faceZones()[groupName];

            const dictionary& dict =
                selectors[selectorI].dict().optionalSubDict("patches");

            FixedList<label, 2> newPatchIDs;
            forAll(patchEntryKeywords, i)
            {
                const dictionary& patchDict =
                    dict.lookupEntryBackwardsCompatible
                    (
                        patchEntryKeywords[i],
                        false,
                        false
                    ).dict();

                const word patchName(patchDict.lookup<word>("name"));

                newPatchIDs[i] = bMesh.findPatchID(patchName);
            }

            const label nModified =
                createFaces
                (
                    internalFacesOnly,
                    mesh,
                    fZone,
                    newPatchIDs[0],
                    newPatchIDs[1],
                    meshMod,
                    modifiedFace
                );

            Info<< "Converted " << nModified
                << " faces into boundary faces on ";
            if (newPatchIDs[0] == newPatchIDs[1])
            {
                Info<< "patch '" << bMesh[newPatchIDs[0]].name() << "'";
            }
            else
            {
                Info<< "patches '" << bMesh[newPatchIDs[0]].name()
                    << "' and '" << bMesh[newPatchIDs[1]].name() << "'";
            }
            Info<< endl;
        }
    }

    Info<< endl;

    if (!overwrite)
    {
        runTime++;
    }

    // Change the mesh. Change points directly (no inflation).
    autoPtr<polyTopoChangeMap> map = meshMod.changeMesh(mesh, false);

    // Update fields
    mesh.topoChange(map);

    // Change the patch fields
    HashSet<word> bafflePatchFields;
    forAll(selectors, selectorI)
    {
        const dictionary& dict =
            selectors[selectorI].dict().optionalSubDict("patches");

        forAll(patchEntryKeywords, i)
        {
            const dictionary& patchDict =
                dict.lookupEntryBackwardsCompatible
                (
                    patchEntryKeywords[i],
                    false,
                    false
                ).dict();

            const word patchName(patchDict.lookup<word>("name"));

            if (!bafflePatchFields.found(patchName))
            {
                bafflePatchFields.insert(patchName);

                dictionary patchFieldsDict =
                    patchDict.subOrEmptyDict("patchFields");

                fvMeshTools::setPatchFields
                (
                    mesh,
                    bMesh.findPatchID(patchName),
                    patchFieldsDict
                );
            }
            else
            {
                if (patchDict.found("patchFields"))
                {
                    WarningInFunction
                        << "Patch field settings found in " << patchDict.name()
                        << " for patch '" << patchName << "', but fields have "
                        << "already been set for this patch so these settings "
                        << "will not be used" << endl;
                }
            }
        }
    }

    // Move mesh (since morphing might not do this)
    if (map().hasMotionPoints())
    {
        mesh.setPoints(map().preMotionPoints());
    }

    // Remove any now zero-sized patches
    filterPatches(mesh, bafflePatches);

    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }

    Info<< "Writing mesh to " << runTime.name() << endl;

    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
