/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    PDRMesh

Description
    Mesh and field preparation utility for PDR type simulations.

    Reads
    - cellSet giving blockedCells
    - faceSets giving blockedFaces and the patch they should go into

    NOTE: To avoid exposing wrong fields values faceSets should include
    faces contained in the blockedCells cellset.

    - coupledFaces reads coupledFacesSet to introduces mixe-coupled baffles

    Subsets out the blocked cells and splits the blockedFaces and updates
    fields.

    The face splitting is done by duplicating the faces. No duplication of
    points for now (so checkMesh will show a lot of 'duplicate face' messages)

\*---------------------------------------------------------------------------*/

#include "fvMeshSubset.H"
#include "argList.H"
#include "cellSet.H"
#include "IOobjectList.H"
#include "volFields.H"
#include "mapPolyMesh.H"
#include "faceSet.H"
#include "cellSet.H"
#include "syncTools.H"
#include "polyTopoChange.H"
#include "polyModifyFace.H"
#include "polyAddFace.H"
#include "regionSplit.H"
#include "Tuple2.H"
#include "cyclicFvPatch.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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


template<class Type>
void subsetVolFields
(
    const fvMeshSubset& subsetter,
    const IOobjectList& objectsList,
    const label patchi,
    const Type& exposedValue,
    const word GeomVolType,
    PtrList<GeometricField<Type, fvPatchField, volMesh>>& subFields
)
{
    const fvMesh& baseMesh = subsetter.baseMesh();

    label i = 0;

    forAllConstIter(IOobjectList , objectsList, iter)
    {
        if (iter()->headerClassName() == GeomVolType)
        {
            const word fieldName = iter()->name();

            Info<< "Subsetting field " << fieldName << endl;

            GeometricField<Type, fvPatchField, volMesh> volField
            (
                *iter(),
                baseMesh
            );

            subFields.set(i, subsetter.interpolate(volField));

            // Explicitly set exposed faces (in patchi) to exposedValue.
            if (patchi >= 0)
            {
                fvPatchField<Type>& fld =
                    subFields[i++].boundaryFieldRef()[patchi];

                label newStart = fld.patch().patch().start();

                label oldPatchi = subsetter.patchMap()[patchi];

                if (oldPatchi == -1)
                {
                    // New patch. Reset whole value.
                    fld = exposedValue;
                }
                else
                {
                    // Reset those faces that originate from different patch
                    // or internal faces.
                    label oldSize = volField.boundaryField()[oldPatchi].size();
                    label oldStart = volField.boundaryField()
                    [
                        oldPatchi
                    ].patch().patch().start();

                    forAll(fld, j)
                    {
                        label oldFacei = subsetter.faceMap()[newStart+j];

                        if (oldFacei < oldStart || oldFacei >= oldStart+oldSize)
                        {
                            fld[j] = exposedValue;
                        }
                    }
                }
            }
        }
    }
}


template<class Type>
void subsetSurfaceFields
(
    const fvMeshSubset& subsetter,
    const IOobjectList& objectsList,
    const label patchi,
    const Type& exposedValue,
    const word GeomSurfType,
    PtrList<GeometricField<Type, fvsPatchField, surfaceMesh>>& subFields
)
{
    const fvMesh& baseMesh = subsetter.baseMesh();

    label i(0);

    forAllConstIter(IOobjectList , objectsList, iter)
    {
        if (iter()->headerClassName() == GeomSurfType)
        {
            const word& fieldName = iter.key();

            Info<< "Subsetting field " << fieldName << endl;

            GeometricField<Type, fvsPatchField, surfaceMesh> volField
            (
                *iter(),
                baseMesh
            );

            subFields.set(i, subsetter.interpolate(volField));


            // Explicitly set exposed faces (in patchi) to exposedValue.
            if (patchi >= 0)
            {
                fvsPatchField<Type>& fld =
                    subFields[i++].boundaryFieldRef()[patchi];

                label newStart = fld.patch().patch().start();

                label oldPatchi = subsetter.patchMap()[patchi];

                if (oldPatchi == -1)
                {
                    // New patch. Reset whole value.
                    fld = exposedValue;
                }
                else
                {
                    // Reset those faces that originate from different patch
                    // or internal faces.
                    label oldSize = volField.boundaryField()[oldPatchi].size();
                    label oldStart = volField.boundaryField()
                    [
                        oldPatchi
                    ].patch().patch().start();

                    forAll(fld, j)
                    {
                        label oldFacei = subsetter.faceMap()[newStart+j];

                        if (oldFacei < oldStart || oldFacei >= oldStart+oldSize)
                        {
                            fld[j] = exposedValue;
                        }
                    }
                }
            }
        }
    }
}


// Faces introduced into zero-sized patches don't get a value at all.
// This is hack to set them to an initial value.
template<class GeoField>
void initCreatedPatches
(
    const fvMesh& mesh,
    const mapPolyMesh& map,
    const typename GeoField::value_type initValue
)
{
    HashTable<const GeoField*> fields
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );

    for
    (
        typename HashTable<const GeoField*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        GeoField& field = const_cast<GeoField&>(*fieldIter());

        typename GeoField::Boundary& fieldBf =
            field.boundaryFieldRef();

        forAll(fieldBf, patchi)
        {
            if (map.oldPatchSizes()[patchi] == 0)
            {
                // Not mapped.
                fieldBf[patchi] = initValue;

                if (fieldBf[patchi].fixesValue())
                {
                    fieldBf[patchi] == initValue;
                }
            }
        }
    }
}


void createCoupledBaffles
(
    fvMesh& mesh,
    const labelList& coupledWantedPatch,
    polyTopoChange& meshMod,
    PackedBoolList&  modifiedFace
)
{
    const meshFaceZones& faceZones = mesh.faceZones();

    forAll(coupledWantedPatch, facei)
    {
        if (coupledWantedPatch[facei] != -1)
        {
            const face& f = mesh.faces()[facei];
            label zoneID = faceZones.whichZone(facei);
            bool zoneFlip = false;

            if (zoneID >= 0)
            {
                const faceZone& fZone = faceZones[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
            }

            // Use owner side of face
            modifyOrAddFace
            (
                meshMod,
                f,                          // modified face
                facei,                      // label of face
                mesh.faceOwner()[facei],    // owner
                false,                      // face flip
                coupledWantedPatch[facei],  // patch for face
                zoneID,                     // zone for face
                zoneFlip,                   // face flip in zone
                modifiedFace                // modify or add status
            );

            if (mesh.isInternalFace(facei))
            {
                label zoneID = faceZones.whichZone(facei);
                bool zoneFlip = false;

                if (zoneID >= 0)
                {
                    const faceZone& fZone = faceZones[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
                }
                // Use neighbour side of face
                modifyOrAddFace
                (
                    meshMod,
                    f.reverseFace(),            // modified face
                    facei,                      // label of face
                    mesh.faceNeighbour()[facei],// owner
                    false,                      // face flip
                    coupledWantedPatch[facei],  // patch for face
                    zoneID,                     // zone for face
                    zoneFlip,                   // face flip in zone
                    modifiedFace                // modify or add status
                );
            }
        }
    }
}


void createCyclicCoupledBaffles
(
    fvMesh& mesh,
    const labelList& cyclicMasterPatch,
    const labelList& cyclicSlavePatch,
    polyTopoChange& meshMod,
    PackedBoolList&  modifiedFace
)
{
    const meshFaceZones& faceZones = mesh.faceZones();

    forAll(cyclicMasterPatch, facei)
    {
        if (cyclicMasterPatch[facei] != -1)
        {
            const face& f = mesh.faces()[facei];

            label zoneID = faceZones.whichZone(facei);
            bool zoneFlip = false;

            if (zoneID >= 0)
            {
                const faceZone& fZone = faceZones[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
            }

            modifyOrAddFace
            (
                meshMod,
                f.reverseFace(),                    // modified face
                facei,                              // label of face
                mesh.faceNeighbour()[facei],        // owner
                false,                              // face flip
                cyclicMasterPatch[facei],           // patch for face
                zoneID,                             // zone for face
                zoneFlip,                           // face flip in zone
                modifiedFace                        // modify or add
            );
        }
    }

    forAll(cyclicSlavePatch, facei)
    {
        if (cyclicSlavePatch[facei] != -1)
        {
            const face& f = mesh.faces()[facei];
            if (mesh.isInternalFace(facei))
            {
                label zoneID = faceZones.whichZone(facei);
                bool zoneFlip = false;

                if (zoneID >= 0)
                {
                    const faceZone& fZone = faceZones[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
                }
            // Use owner side of face
                modifyOrAddFace
                (
                    meshMod,
                    f,                          // modified face
                    facei,                      // label of face
                    mesh.faceOwner()[facei],    // owner
                    false,                      // face flip
                    cyclicSlavePatch[facei],    // patch for face
                    zoneID,                     // zone for face
                    zoneFlip,                   // face flip in zone
                    modifiedFace                // modify or add status
                );
            }
        }
    }
}


void createBaffles
(
    fvMesh& mesh,
    const labelList& wantedPatch,
    polyTopoChange& meshMod
)
{
    const meshFaceZones& faceZones = mesh.faceZones();
    Info << "faceZone:createBaffle " << faceZones << endl;
    forAll(wantedPatch, facei)
    {
        if (wantedPatch[facei] != -1)
        {
            const face& f = mesh.faces()[facei];

            label zoneID = faceZones.whichZone(facei);
            bool zoneFlip = false;

            if (zoneID >= 0)
            {
                const faceZone& fZone = faceZones[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
            }

            meshMod.setAction
            (
                polyModifyFace
                (
                    f,                          // modified face
                    facei,                      // label of face
                    mesh.faceOwner()[facei],    // owner
                    -1,                         // neighbour
                    false,                      // face flip
                    wantedPatch[facei],         // patch for face
                    false,                      // remove from zone
                    zoneID,                     // zone for face
                    zoneFlip                    // face flip in zone
                )
            );

            if (mesh.isInternalFace(facei))
            {
                label zoneID = faceZones.whichZone(facei);
                bool zoneFlip = false;

                if (zoneID >= 0)
                {
                    const faceZone& fZone = faceZones[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
                }

                meshMod.setAction
                (
                    polyAddFace
                    (
                        f.reverseFace(),            // modified face
                        mesh.faceNeighbour()[facei],// owner
                        -1,                         // neighbour
                        -1,                         // masterPointID
                        -1,                         // masterEdgeID
                        facei,                      // masterFaceID,
                        false,                      // face flip
                        wantedPatch[facei],         // patch for face
                        zoneID,                     // zone for face
                        zoneFlip                    // face flip in zone
                    )
                );
            }
        }
    }
}


// Wrapper around find patch. Also makes sure same patch in parallel.
label findPatch(const polyBoundaryMesh& patches, const word& patchName)
{
    label patchi = patches.findPatchID(patchName);

    if (patchi == -1)
    {
        FatalErrorInFunction
            << "Illegal patch " << patchName
            << nl << "Valid patches are " << patches.names()
            << exit(FatalError);
    }

    // Check same patch for all procs
    {
        label newPatch = patchi;
        reduce(newPatch, minOp<label>());

        if (newPatch != patchi)
        {
            FatalErrorInFunction
                << "Patch " << patchName
                << " should have the same patch index on all processors." << nl
                << "On my processor it has index " << patchi
                << " ; on some other processor it has index " << newPatch
                << exit(FatalError);
        }
    }
    return patchi;
}



int main(int argc, char *argv[])
{
    #include "addOverwriteOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();
    #include "createMesh.H"

    // Read control dictionary
    // ~~~~~~~~~~~~~~~~~~~~~~~

    IOdictionary dict
    (
        IOobject
        (
            "PDRMeshDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    // Per faceSet the patch to put the baffles into
    const List<Pair<word>> setsAndPatches(dict.lookup("blockedFaces"));

    // Per faceSet the patch to put the coupled baffles into
    DynamicList<FixedList<word, 3>> coupledAndPatches(10);
    const dictionary& functionDicts = dict.subDict("coupledFaces");
    forAllConstIter(dictionary, functionDicts, iter)
    {
        // safety:
        if (!iter().isDict())
        {
            continue;
        }
        const word& key = iter().keyword();

        const dictionary& dict = iter().dict();
        const word cyclicName = dict.lookup("cyclicMasterPatchName");
        const word wallName = dict.lookup("wallPatchName");
        FixedList<word, 3> nameAndType;
        nameAndType[0] = key;
        nameAndType[1] = wallName;
        nameAndType[2] = cyclicName;
        coupledAndPatches.append(nameAndType);
    }

    forAll(setsAndPatches, setI)
    {
        Info<< "Faces in faceSet " << setsAndPatches[setI][0]
            << " become baffles in patch " << setsAndPatches[setI][1]
            << endl;
    }

    forAll(coupledAndPatches, setI)
    {
        Info<< "Faces in faceSet " << coupledAndPatches[setI][0]
            << " become coupled baffles in patch " << coupledAndPatches[setI][1]
            << endl;
    }

    // All exposed faces that are not explicitly marked to be put into a patch
    const word defaultPatch(dict.lookup("defaultPatch"));

    Info<< "Faces that get exposed become boundary faces in patch "
        << defaultPatch << endl;

    const word blockedSetName(dict.lookup("blockedCells"));

    Info<< "Reading blocked cells from cellSet " << blockedSetName
        << endl;

    const bool overwrite = args.optionFound("overwrite");


    // Read faceSets, lookup patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Check that face sets don't have coincident faces
    labelList wantedPatch(mesh.nFaces(), -1);
    forAll(setsAndPatches, setI)
    {
        faceSet fSet(mesh, setsAndPatches[setI][0]);

        label patchi = findPatch
        (
            mesh.boundaryMesh(),
            setsAndPatches[setI][1]
        );

        forAllConstIter(faceSet, fSet, iter)
        {
            if (wantedPatch[iter.key()] != -1)
            {
                FatalErrorInFunction
                    << "Face " << iter.key()
                    << " is in faceSet " << setsAndPatches[setI][0]
                    << " destined for patch " << setsAndPatches[setI][1]
                    << " but also in patch " << wantedPatch[iter.key()]
                    << exit(FatalError);
            }
            wantedPatch[iter.key()] = patchi;
        }
    }

    // Per face the patch for coupled baffle or -1.
    labelList coupledWantedPatch(mesh.nFaces(), -1);
    labelList cyclicWantedPatch_half0(mesh.nFaces(), -1);
    labelList cyclicWantedPatch_half1(mesh.nFaces(), -1);

    forAll(coupledAndPatches, setI)
    {
        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        const label cyclicId =
            findPatch(patches, coupledAndPatches[setI][2]);

        const label cyclicSlaveId = findPatch
        (
            patches,
            refCast<const cyclicFvPatch>
            (
                mesh.boundary()[cyclicId]
            ).neighbFvPatch().name()
        );

        faceSet fSet(mesh, coupledAndPatches[setI][0]);
        label patchi = findPatch(patches, coupledAndPatches[setI][1]);

        forAllConstIter(faceSet, fSet, iter)
        {
            if (coupledWantedPatch[iter.key()] != -1)
            {
                FatalErrorInFunction
                    << "Face " << iter.key()
                    << " is in faceSet " << coupledAndPatches[setI][0]
                    << " destined for patch " << coupledAndPatches[setI][1]
                    << " but also in patch " << coupledWantedPatch[iter.key()]
                    << exit(FatalError);
            }
                coupledWantedPatch[iter.key()] = patchi;
                cyclicWantedPatch_half0[iter.key()] = cyclicId;
                cyclicWantedPatch_half1[iter.key()] = cyclicSlaveId;
        }
    }

    // Exposed faces patch
    label defaultPatchi = findPatch(mesh.boundaryMesh(), defaultPatch);


    //
    // Removing blockedCells (like subsetMesh)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //

    // Create mesh subsetting engine
    fvMeshSubset subsetter(mesh);

    {

        cellSet blockedCells(mesh, blockedSetName);

        // invert
        blockedCells.invert(mesh.nCells());

        // Create subsetted mesh.
        subsetter.setLargeCellSubset(blockedCells, defaultPatchi, true);
    }


    // Subset wantedPatch. Note that might also include boundary faces
    // that have been exposed by subsetting.
    wantedPatch = IndirectList<label>(wantedPatch, subsetter.faceMap())();

    coupledWantedPatch = IndirectList<label>
    (
        coupledWantedPatch,
        subsetter.faceMap()
    )();

    cyclicWantedPatch_half0 = IndirectList<label>
    (
        cyclicWantedPatch_half0,
        subsetter.faceMap()
    )();

    cyclicWantedPatch_half1 = IndirectList<label>
    (
        cyclicWantedPatch_half1,
        subsetter.faceMap()
    )();

    // Read all fields in time and constant directories
    IOobjectList objects(mesh, runTime.timeName());
    IOobjectList timeObjects(IOobjectList(mesh, mesh.facesInstance()));
    forAllConstIter(IOobjectList, timeObjects, iter)
    {
        if
        (
            iter()->headerClassName() == volScalarField::typeName
         || iter()->headerClassName() == volVectorField::typeName
         || iter()->headerClassName() == volSphericalTensorField::typeName
         || iter()->headerClassName() == volTensorField::typeName
         || iter()->headerClassName() == volSymmTensorField::typeName
         || iter()->headerClassName() == surfaceScalarField::typeName
         || iter()->headerClassName() == surfaceVectorField::typeName
         || iter()->headerClassName()
            == surfaceSphericalTensorField::typeName
         || iter()->headerClassName() == surfaceSymmTensorField::typeName
         || iter()->headerClassName() == surfaceTensorField::typeName
        )
        {
            objects.add(*iter());
        }
    }

    // Read vol fields and subset.

    wordList scalarNames(objects.names(volScalarField::typeName));
    PtrList<volScalarField> scalarFlds(scalarNames.size());
    subsetVolFields
    (
        subsetter,
        objects,
        defaultPatchi,
        scalar(Zero),
        volScalarField::typeName,
        scalarFlds
    );

    wordList vectorNames(objects.names(volVectorField::typeName));
    PtrList<volVectorField> vectorFlds(vectorNames.size());
    subsetVolFields
    (
        subsetter,
        objects,
        defaultPatchi,
        vector(Zero),
        volVectorField::typeName,
        vectorFlds
    );

    wordList sphericalTensorNames
    (
        objects.names(volSphericalTensorField::typeName)
    );
    PtrList<volSphericalTensorField> sphericalTensorFlds
    (
        sphericalTensorNames.size()
    );
    subsetVolFields
    (
        subsetter,
        objects,
        defaultPatchi,
        sphericalTensor(Zero),
        volSphericalTensorField::typeName,
        sphericalTensorFlds
    );

    wordList symmTensorNames(objects.names(volSymmTensorField::typeName));
    PtrList<volSymmTensorField> symmTensorFlds(symmTensorNames.size());
    subsetVolFields
    (
        subsetter,
        objects,
        defaultPatchi,
        symmTensor(Zero),
        volSymmTensorField::typeName,
        symmTensorFlds
    );

    wordList tensorNames(objects.names(volTensorField::typeName));
    PtrList<volTensorField> tensorFlds(tensorNames.size());
    subsetVolFields
    (
        subsetter,
        objects,
        defaultPatchi,
        tensor(Zero),
        volTensorField::typeName,
        tensorFlds
    );

    // Read surface fields and subset.

    wordList surfScalarNames(objects.names(surfaceScalarField::typeName));
    PtrList<surfaceScalarField> surfScalarFlds(surfScalarNames.size());
    subsetSurfaceFields
    (
        subsetter,
        objects,
        defaultPatchi,
        scalar(Zero),
        surfaceScalarField::typeName,
        surfScalarFlds
    );

    wordList surfVectorNames(objects.names(surfaceVectorField::typeName));
    PtrList<surfaceVectorField> surfVectorFlds(surfVectorNames.size());
    subsetSurfaceFields
    (
        subsetter,
        objects,
        defaultPatchi,
        vector(Zero),
        surfaceVectorField::typeName,
        surfVectorFlds
    );

    wordList surfSphericalTensorNames
    (
        objects.names(surfaceSphericalTensorField::typeName)
    );
    PtrList<surfaceSphericalTensorField> surfSphericalTensorFlds
    (
        surfSphericalTensorNames.size()
    );
    subsetSurfaceFields
    (
        subsetter,
        objects,
        defaultPatchi,
        sphericalTensor(Zero),
        surfaceSphericalTensorField::typeName,
        surfSphericalTensorFlds
    );

    wordList surfSymmTensorNames
    (
        objects.names(surfaceSymmTensorField::typeName)
    );

    PtrList<surfaceSymmTensorField> surfSymmTensorFlds
    (
        surfSymmTensorNames.size()
    );

    subsetSurfaceFields
    (
        subsetter,
        objects,
        defaultPatchi,
        symmTensor(Zero),
        surfaceSymmTensorField::typeName,
        surfSymmTensorFlds
    );

    wordList surfTensorNames(objects.names(surfaceTensorField::typeName));
    PtrList<surfaceTensorField> surfTensorFlds(surfTensorNames.size());
    subsetSurfaceFields
    (
        subsetter,
        objects,
        defaultPatchi,
        tensor(Zero),
        surfaceTensorField::typeName,
        surfTensorFlds
    );

    if (!overwrite)
    {
        runTime++;
    }

    Info<< "Writing mesh without blockedCells to time " << runTime.value()
        << endl;

    // Subsetting adds 'subset' prefix. Rename field to be like original.
    forAll(scalarFlds, i)
    {
        scalarFlds[i].rename(scalarNames[i]);
        scalarFlds[i].writeOpt() = IOobject::AUTO_WRITE;
        scalarFlds[i].checkIn();
    }
    forAll(vectorFlds, i)
    {
        vectorFlds[i].rename(vectorNames[i]);
        vectorFlds[i].writeOpt() = IOobject::AUTO_WRITE;
        vectorFlds[i].checkIn();
    }
    forAll(sphericalTensorFlds, i)
    {
        sphericalTensorFlds[i].rename(sphericalTensorNames[i]);
        sphericalTensorFlds[i].writeOpt() = IOobject::AUTO_WRITE;
        sphericalTensorFlds[i].checkIn();
    }
    forAll(symmTensorFlds, i)
    {
        symmTensorFlds[i].rename(symmTensorNames[i]);
        symmTensorFlds[i].writeOpt() = IOobject::AUTO_WRITE;
        symmTensorFlds[i].checkIn();
    }
    forAll(tensorFlds, i)
    {
        tensorFlds[i].rename(tensorNames[i]);
        tensorFlds[i].writeOpt() = IOobject::AUTO_WRITE;
        tensorFlds[i].checkIn();
    }

    // Surface ones.
    forAll(surfScalarFlds, i)
    {
        surfScalarFlds[i].rename(surfScalarNames[i]);
        surfScalarFlds[i].writeOpt() = IOobject::AUTO_WRITE;
        surfScalarFlds[i].checkIn();
    }
    forAll(surfVectorFlds, i)
    {
        surfVectorFlds[i].rename(surfVectorNames[i]);
        surfVectorFlds[i].writeOpt() = IOobject::AUTO_WRITE;
        surfVectorFlds[i].checkIn();
    }
    forAll(surfSphericalTensorFlds, i)
    {
        surfSphericalTensorFlds[i].rename(surfSphericalTensorNames[i]);
        surfSphericalTensorFlds[i].writeOpt() = IOobject::AUTO_WRITE;
        surfSphericalTensorFlds[i].checkIn();
    }
    forAll(surfSymmTensorFlds, i)
    {
        surfSymmTensorFlds[i].rename(surfSymmTensorNames[i]);
        surfSymmTensorFlds[i].writeOpt() = IOobject::AUTO_WRITE;
        surfSymmTensorFlds[i].checkIn();
    }
    forAll(surfTensorNames, i)
    {
        surfTensorFlds[i].rename(surfTensorNames[i]);
        surfTensorFlds[i].writeOpt() = IOobject::AUTO_WRITE;
        surfTensorFlds[i].checkIn();
    }

    subsetter.subMesh().write();


    //
    // Splitting blockedFaces
    // ~~~~~~~~~~~~~~~~~~~~~~
    //

    // Synchronise wantedPatch across coupled patches.
    syncTools::syncFaceList
    (
        subsetter.subMesh(),
        wantedPatch,
        maxEqOp<label>()
    );

    // Synchronise coupledWantedPatch across coupled patches.
    syncTools::syncFaceList
    (
        subsetter.subMesh(),
        coupledWantedPatch,
        maxEqOp<label>()
    );

    // Synchronise cyclicWantedPatch across coupled patches.
    syncTools::syncFaceList
    (
        subsetter.subMesh(),
        cyclicWantedPatch_half0,
        maxEqOp<label>()
    );

    // Synchronise cyclicWantedPatch across coupled patches.
    syncTools::syncFaceList
    (
        subsetter.subMesh(),
        cyclicWantedPatch_half1,
        maxEqOp<label>()
    );

    // Topochange container
    polyTopoChange meshMod(subsetter.subMesh());


    // Whether first use of face (modify) or consecutive (add)
    PackedBoolList modifiedFace(mesh.nFaces());

    // Create coupled wall-side baffles
    createCoupledBaffles
    (
        subsetter.subMesh(),
        coupledWantedPatch,
        meshMod,
        modifiedFace
    );

    // Create coupled master/slave cyclic baffles
    createCyclicCoupledBaffles
    (
        subsetter.subMesh(),
        cyclicWantedPatch_half0,
        cyclicWantedPatch_half1,
        meshMod,
        modifiedFace
    );

    // Create wall baffles
    createBaffles
    (
        subsetter.subMesh(),
        wantedPatch,
        meshMod
    );

    if (!overwrite)
    {
        runTime++;
    }

    // Change the mesh. Change points directly (no inflation).
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(subsetter.subMesh(), false);

    // Update fields
    subsetter.subMesh().updateMesh(map);

    // Fix faces that get mapped to zero-sized patches (these don't get any
    // value)
    initCreatedPatches<volScalarField>
    (
        subsetter.subMesh(),
        map,
        0.0
    );
    initCreatedPatches<volVectorField>
    (
        subsetter.subMesh(),
        map,
        Zero
    );
    initCreatedPatches<volSphericalTensorField>
    (
        subsetter.subMesh(),
        map,
        Zero
    );
    initCreatedPatches<volSymmTensorField>
    (
        subsetter.subMesh(),
        map,
        Zero
    );
    initCreatedPatches<volTensorField>
    (
        subsetter.subMesh(),
        map,
        Zero
    );

    initCreatedPatches<surfaceScalarField>
    (
        subsetter.subMesh(),
        map,
        0.0
    );
    initCreatedPatches<surfaceVectorField>
    (
        subsetter.subMesh(),
        map,
        Zero
    );
    initCreatedPatches<surfaceSphericalTensorField>
    (
        subsetter.subMesh(),
        map,
        Zero
    );
    initCreatedPatches<surfaceSymmTensorField>
    (
        subsetter.subMesh(),
        map,
        Zero
    );
    initCreatedPatches<surfaceTensorField>
    (
        subsetter.subMesh(),
        map,
        Zero
    );


    // Move mesh (since morphing might not do this)
    if (map().hasMotionPoints())
    {
        subsetter.subMesh().movePoints(map().preMotionPoints());
    }

    Info<< "Writing mesh with split blockedFaces to time " << runTime.value()
        << endl;

    subsetter.subMesh().write();


    //
    // Removing inaccessible regions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //

    // Determine connected regions. regionSplit is the labelList with the
    // region per cell.
    regionSplit cellRegion(subsetter.subMesh());

    if (cellRegion.nRegions() > 1)
    {
        WarningInFunction
            << "Removing blocked faces and cells created "
            << cellRegion.nRegions()
            << " regions that are not connected via a face." << nl
            << "    This is not supported in solvers." << nl
            << "    Use" << nl << nl
            << "    splitMeshRegions <root> <case> -largestOnly" << nl << nl
            << "    to extract a single region of the mesh." << nl
            << "    This mesh will be written to a new timedirectory"
            << " so might have to be moved back to constant/" << nl
            << endl;

        word startFrom(runTime.controlDict().lookup("startFrom"));

        if (startFrom != "latestTime")
        {
            WarningInFunction
                << "To run splitMeshRegions please set your"
                << " startFrom entry to latestTime" << endl;
        }
    }

    Info << nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
