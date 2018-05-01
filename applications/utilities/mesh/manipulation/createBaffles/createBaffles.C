/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Makes internal faces into boundary faces. Does not duplicate points, unlike
    mergeOrSplitBaffles.

    Note: if any coupled patch face is selected for baffling the opposite
    member has to be selected for baffling as well.

    - if the patch already exists will not override it nor its fields
    - if the patch does not exist it will be created together with 'calculated'
      patchfields unless the field is mentioned in the patchFields section.
    - any 0-sized patches (since faces have been moved out) will get removed

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyTopoChange.H"
#include "polyModifyFace.H"
#include "polyAddFace.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMeshMapper.H"
#include "faceSelection.H"

#include "fvMeshTools.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label addPatch
(
    fvMesh& mesh,
    const word& patchName,
    const word& groupName,
    const dictionary& patchDict
)
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    if (pbm.findPatchID(patchName) == -1)
    {
        autoPtr<polyPatch> ppPtr
        (
            polyPatch::New
            (
                patchName,
                patchDict,
                0,
                pbm
            )
        );
        polyPatch& pp = ppPtr();

        if (!groupName.empty() && !pp.inGroup(groupName))
        {
            pp.inGroups().append(groupName);
        }


        // Add patch, create calculated everywhere
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
            << "' already exists.  Only "
            << "moving patch faces - type will remain the same"
            << endl;
    }

    return pbm.findPatchID(patchName);
}


// Filter out the empty patches.
void filterPatches(fvMesh& mesh, const HashSet<word>& addedPatchNames)
{
    // Remove any zero-sized ones. Assumes
    // - processor patches are already only there if needed
    // - all other patches are available on all processors
    // - but coupled ones might still be needed, even if zero-size
    //   (e.g. processorCyclic)
    // See also logic in createPatch.
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    labelList oldToNew(pbm.size(), -1);
    label newPatchi = 0;
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (!isA<processorPolyPatch>(pp))
        {
            if
            (
                isA<coupledPolyPatch>(pp)
             || returnReduce(pp.size(), sumOp<label>())
             || addedPatchNames.found(pp.name())
            )
            {
                // Coupled (and unknown size) or uncoupled and used
                oldToNew[patchi] = newPatchi++;
            }
        }
    }

    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            oldToNew[patchi] = newPatchi++;
        }
    }


    const label nKeepPatches = newPatchi;

    // Shuffle unused ones to end
    if (nKeepPatches != pbm.size())
    {
        Info<< endl
            << "Removing zero-sized patches:" << endl << incrIndent;

        forAll(oldToNew, patchi)
        {
            if (oldToNew[patchi] == -1)
            {
                Info<< indent << pbm[patchi].name()
                    << " type " << pbm[patchi].type()
                    << " at position " << patchi << endl;
                oldToNew[patchi] = newPatchi++;
            }
        }
        Info<< decrIndent;

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
        // Second or more usage of face. Add.
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


// Create faces for fZone faces. Usually newMasterPatches, newSlavePatches
// only size one but can be more for duplicate baffle sets
void createFaces
(
    const bool internalFacesOnly,
    const fvMesh& mesh,
    const faceZone& fZone,
    const labelList& newMasterPatches,
    const labelList& newSlavePatches,
    polyTopoChange& meshMod,
    PackedBoolList& modifiedFace,
    label& nModified
)
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    forAll(newMasterPatches, i)
    {
        // Pass 1. Do selected side of zone
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            label zoneFacei = fZone.whichFace(facei);

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
                        newMasterPatches[i],    // patch for face
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
                        mesh.faces()[facei].reverseFace(),  // modified face
                        facei,                      // label of face
                        mesh.faceNeighbour()[facei],// owner
                        true,                       // face flip
                        newMasterPatches[i],        // patch for face
                        fZone.index(),              // zone for face
                        false,                      // face flip in zone
                        modifiedFace                // modify or add status
                    );
                }

                nModified++;
            }
        }


        // Pass 2. Do other side of zone
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
                        newSlavePatches[i],             // patch for face
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
                        newSlavePatches[i],     // patch for face
                        fZone.index(),          // zone for face
                        true,                   // face flip in zone
                        modifiedFace            // modify or add status
                    );
                }
            }
        }


        // Modify any boundary faces
        // ~~~~~~~~~~~~~~~~~~~~~~~~~

        // Normal boundary:
        // - move to new patch. Might already be back-to-back baffle
        // you want to add cyclic to. Do warn though.
        //
        // Processor boundary:
        // - do not move to cyclic
        // - add normal patches though.

        // For warning once per patch.
        labelHashSet patchWarned;

        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];

            const label newMasterPatchi = newMasterPatches[i];
            const label newSlavePatchi = newSlavePatches[i];

            if
            (
                pp.coupled()
             && (
                    pbm[newMasterPatchi].coupled()
                 || pbm[newSlavePatchi].coupled()
                )
            )
            {
                // Do not allow coupled faces to be moved to different
                // coupled patches.
            }
            else if (pp.coupled() || !internalFacesOnly)
            {
                forAll(pp, i)
                {
                    label facei = pp.start()+i;

                    label zoneFacei = fZone.whichFace(facei);

                    if (zoneFacei != -1)
                    {
                        if (patchWarned.insert(patchi))
                        {
                            WarningInFunction
                                << "Found boundary face (in patch "
                                << pp.name()
                                << ") in faceZone " << fZone.name()
                                << " to convert to baffle patches "
                                << pbm[newMasterPatchi].name() << "/"
                                << pbm[newSlavePatchi].name()
                                << endl
                                << "    Set internalFacesOnly to true in the"
                                << " createBaffles control dictionary if you"
                                << " don't wish to convert boundary faces."
                                << endl;
                        }

                        modifyOrAddFace
                        (
                            meshMod,
                            mesh.faces()[facei],        // modified face
                            facei,                      // label of face
                            mesh.faceOwner()[facei],    // owner
                            false,                      // face flip
                            fZone.flipMap()[zoneFacei]
                          ? newSlavePatchi
                          : newMasterPatchi,            // patch for face
                            fZone.index(),              // zone for face
                            fZone.flipMap()[zoneFacei], // face flip in zone
                            modifiedFace                // modify or add
                        );

                        nModified++;
                    }
                }
            }
        }
    }
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
    #include "createTime.H"
    runTime.functionObjects().off();
    #include "createNamedMesh.H"


    const bool overwrite = args.optionFound("overwrite");

    const word oldInstance = mesh.pointsInstance();

    const word dictName("createBafflesDict");
    #include "setSystemMeshDictionaryIO.H"

    Switch internalFacesOnly(false);

    Switch fields(false);

    PtrList<faceSelection> selectors;
    {
        Info<< "Reading baffle criteria from " << dictName << nl << endl;
        IOdictionary dict(dictIO);

        dict.lookup("internalFacesOnly") >> internalFacesOnly;
        fields = dict.lookupOrDefault("fields", false);

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


    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    if (fields) Info<< "Reading geometric fields" << nl << endl;

    PtrList<volScalarField> vsFlds;
    if (fields) ReadFields(mesh, objects, vsFlds);

    PtrList<volVectorField> vvFlds;
    if (fields) ReadFields(mesh, objects, vvFlds);

    PtrList<volSphericalTensorField> vstFlds;
    if (fields) ReadFields(mesh, objects, vstFlds);

    PtrList<volSymmTensorField> vsymtFlds;
    if (fields) ReadFields(mesh, objects, vsymtFlds);

    PtrList<volTensorField> vtFlds;
    if (fields) ReadFields(mesh, objects, vtFlds);

    PtrList<surfaceScalarField> ssFlds;
    if (fields) ReadFields(mesh, objects, ssFlds);

    PtrList<surfaceVectorField> svFlds;
    if (fields) ReadFields(mesh, objects, svFlds);

    PtrList<surfaceSphericalTensorField> sstFlds;
    if (fields) ReadFields(mesh, objects, sstFlds);

    PtrList<surfaceSymmTensorField> ssymtFlds;
    if (fields) ReadFields(mesh, objects, ssymtFlds);

    PtrList<surfaceTensorField> stFlds;
    if (fields) ReadFields(mesh, objects, stFlds);




    // Creating (if necessary) faceZones
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(selectors, selectorI)
    {
        const word& name = selectors[selectorI].name();

        if (mesh.faceZones().findZoneID(name) == -1)
        {
            mesh.faceZones().clearAddressing();
            label sz = mesh.faceZones().size();

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


    // Select faces
    // ~~~~~~~~~~~~

    //- Per face zoneID it is in and flip status.
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
        label zoneID = mesh.faceZones().findZoneID(name);

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

        Info<< "Created zone " << name
            << " at index " << zoneID
            << " with " << n << " faces" << endl;

        mesh.faceZones().set
        (
            zoneID,
            new faceZone(name, addr, flip, zoneID, mesh.faceZones())
        );
    }



    // Count patches to add
    // ~~~~~~~~~~~~~~~~~~~~
    HashSet<word> bafflePatches;
    {
        forAll(selectors, selectorI)
        {
            const dictionary& dict = selectors[selectorI].dict();

            if (dict.found("patches"))
            {
                const dictionary& patchSources = dict.subDict("patches");
                forAllConstIter(dictionary, patchSources, iter)
                {
                    const word patchName(iter().dict()["name"]);
                    bafflePatches.insert(patchName);
                }
            }
            else
            {
                const word masterName = selectors[selectorI].name() + "_master";
                bafflePatches.insert(masterName);
                const word slaveName = selectors[selectorI].name() + "_slave";
                bafflePatches.insert(slaveName);
            }
        }
    }


    // Create baffles
    // ~~~~~~~~~~~~~~
    // Is done in multiple steps
    // - create patches with 'calculated' patchFields
    // - move faces into these patches
    // - change the patchFields to the wanted type
    // This order is done so e.g. fixedJump works:
    // - you cannot create patchfields at the same time as patches since
    //   they do an evaluate upon construction
    // - you want to create the patchField only after you have faces
    //   so you don't get the 'create-from-nothing' mapping problem.


    // Pass 1: add patches
    // ~~~~~~~~~~~~~~~~~~~

    // HashSet<word> addedPatches;
    {
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();
        forAll(selectors, selectorI)
        {
            const dictionary& dict = selectors[selectorI].dict();
            const word& groupName = selectors[selectorI].name();

            if (dict.found("patches"))
            {
                const dictionary& patchSources = dict.subDict("patches");
                forAllConstIter(dictionary, patchSources, iter)
                {
                    const word patchName(iter().dict()["name"]);

                    if (pbm.findPatchID(patchName) == -1)
                    {
                        dictionary patchDict = iter().dict();
                        patchDict.set("nFaces", 0);
                        patchDict.set("startFace", 0);

                        // Note: do not set coupleGroup if constructed from
                        //       baffles so you have freedom specifying it
                        //       yourself.
                        // patchDict.set("coupleGroup", groupName);

                        addPatch(mesh, patchName, groupName, patchDict);
                    }
                    else
                    {
                        Info<< "Patch '" << patchName
                            << "' already exists.  Only "
                            << "moving patch faces - type will remain the same"
                            << endl;
                    }
                }
            }
            else
            {
                const dictionary& patchSource = dict.subDict("patchPairs");
                const word masterName = groupName + "_master";
                const word slaveName = groupName + "_slave";

                word groupNameMaster = groupName;
                word groupNameSlave = groupName;


                dictionary patchDictMaster(patchSource);
                patchDictMaster.set("nFaces", 0);
                patchDictMaster.set("startFace", 0);
                patchDictMaster.set("coupleGroup", groupName);

                dictionary patchDictSlave(patchDictMaster);

                // Note: This is added for the particular case where we want
                // master and slave in different groupNames
                // (ie 3D thermal baffles)

                Switch sameGroup
                (
                    patchSource.lookupOrDefault("sameGroup", true)
                );
                if (!sameGroup)
                {
                    groupNameMaster = groupName + "Group_master";
                    groupNameSlave = groupName + "Group_slave";
                    patchDictMaster.set("coupleGroup", groupNameMaster);
                    patchDictSlave.set("coupleGroup", groupNameSlave);
                }

                addPatch(mesh, masterName, groupNameMaster, patchDictMaster);
                addPatch(mesh, slaveName, groupNameSlave, patchDictSlave);
            }
        }
    }


    // Make sure patches and zoneFaces are synchronised across couples
    mesh.boundaryMesh().checkParallelSync(true);
    mesh.faceZones().checkParallelSync(true);



    // Mesh change container
    polyTopoChange meshMod(mesh);

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();


    // Do the actual changes. Note:
    // - loop in incrementing face order (not necessary if faceZone ordered).
    //   Preserves any existing ordering on patch faces.
    // - two passes, do non-flip faces first and flip faces second. This
    //   guarantees that when e.g. creating a cyclic all faces from one
    //   side come first and faces from the other side next.

    // Whether first use of face (modify) or consecutive (add)
    PackedBoolList modifiedFace(mesh.nFaces());
    label nModified = 0;

    forAll(selectors, selectorI)
    {
        const word& name = selectors[selectorI].name();
        label zoneID = mesh.faceZones().findZoneID(name);
        const faceZone& fZone = mesh.faceZones()[zoneID];

        const dictionary& dict = selectors[selectorI].dict();

        DynamicList<label> newMasterPatches;
        DynamicList<label> newSlavePatches;

        if (dict.found("patches"))
        {
            const dictionary& patchSources = dict.subDict("patches");

            bool master = true;
            forAllConstIter(dictionary, patchSources, iter)
            {
                const word patchName(iter().dict()["name"]);
                label patchi = pbm.findPatchID(patchName);
                if (master)
                {
                    newMasterPatches.append(patchi);
                }
                else
                {
                    newSlavePatches.append(patchi);
                }
                master = !master;
            }
        }
        else
        {
            const word masterName = selectors[selectorI].name() + "_master";
            newMasterPatches.append(pbm.findPatchID(masterName));

            const word slaveName = selectors[selectorI].name() + "_slave";
            newSlavePatches.append(pbm.findPatchID(slaveName));
        }



        createFaces
        (
            internalFacesOnly,
            mesh,
            fZone,
            newMasterPatches,
            newSlavePatches,
            meshMod,
            modifiedFace,
            nModified
        );
    }


    Info<< "Converted " << returnReduce(nModified, sumOp<label>())
        << " faces into boundary faces in patches "
        << bafflePatches.sortedToc() << nl << endl;

    if (!overwrite)
    {
        runTime++;
    }

    // Change the mesh. Change points directly (no inflation).
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

    // Update fields
    mesh.updateMesh(map);



    // Correct boundary faces mapped-out-of-nothing.
    // This is just a hack to correct the value field.
    {
        fvMeshMapper mapper(mesh, map);
        bool hasWarned = false;

        forAllConstIter(HashSet<word>, bafflePatches, iter)
        {
            label patchi = mesh.boundaryMesh().findPatchID(iter.key());

            const fvPatchMapper& pm = mapper.boundaryMap()[patchi];

            if (pm.sizeBeforeMapping() == 0)
            {
                if (!hasWarned)
                {
                    hasWarned = true;
                    WarningInFunction
                        << "Setting field on boundary faces to zero." << endl
                        << "You might have to edit these fields." << endl;
                }

                fvMeshTools::zeroPatchFields(mesh, patchi);
            }
        }
    }


    // Pass 2: change patchFields
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();
        forAll(selectors, selectorI)
        {
            const dictionary& dict = selectors[selectorI].dict();
            if (dict.found("patches"))
            {
                const dictionary& patchSources = dict.subDict("patches");

                forAllConstIter(dictionary, patchSources, iter)
                {
                    const word patchName(iter().dict()["name"]);
                    label patchi = pbm.findPatchID(patchName);

                    if (iter().dict().found("patchFields"))
                    {
                        const dictionary& patchFieldsDict =
                            iter().dict().subDict
                            (
                                "patchFields"
                            );

                        fvMeshTools::setPatchFields
                        (
                            mesh,
                            patchi,
                            patchFieldsDict
                        );
                    }
                }
            }
            else
            {
                const dictionary& patchSource = dict.subDict("patchPairs");

                Switch sameGroup
                (
                    patchSource.lookupOrDefault("sameGroup", true)
                );

                const word& groupName = selectors[selectorI].name();

                if (patchSource.found("patchFields"))
                {
                    dictionary patchFieldsDict = patchSource.subDict
                    (
                        "patchFields"
                    );

                    if (sameGroup)
                    {
                        // Add coupleGroup to all entries
                        forAllIter(dictionary, patchFieldsDict, iter)
                        {
                            if (iter().isDict())
                            {
                                dictionary& dict = iter().dict();
                                dict.set("coupleGroup", groupName);
                            }
                        }

                        const labelList& patchIDs =
                            pbm.groupPatchIDs()[groupName];

                        forAll(patchIDs, i)
                        {
                            fvMeshTools::setPatchFields
                            (
                                mesh,
                                patchIDs[i],
                                patchFieldsDict
                            );
                        }
                    }
                    else
                    {
                        const word masterPatchName(groupName + "_master");
                        const word slavePatchName(groupName + "_slave");

                        label patchiMaster = pbm.findPatchID(masterPatchName);
                        label patchiSlave = pbm.findPatchID(slavePatchName);

                        fvMeshTools::setPatchFields
                        (
                            mesh,
                            patchiMaster,
                            patchFieldsDict
                        );

                        fvMeshTools::setPatchFields
                        (
                            mesh,
                            patchiSlave,
                            patchFieldsDict
                        );
                    }
                }
            }
        }
    }


    // Move mesh (since morphing might not do this)
    if (map().hasMotionPoints())
    {
        mesh.movePoints(map().preMotionPoints());
    }


    // Remove any now zero-sized patches
    filterPatches(mesh, bafflePatches);


    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }

    Info<< "Writing mesh to " << runTime.timeName() << endl;

    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
