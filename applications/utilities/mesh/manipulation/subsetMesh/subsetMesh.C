/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    subsetMesh

Description
    Selects a section of mesh based on a cellSet.

    The utility sub-sets the mesh to choose only a part of interest. Check the
    cellSet/topoSet utilities to see how to select cells based on various
    shapes.

    The mesh will subset all points, faces and cells needed to make a sub-mesh
    but will not preserve attached boundary types.

\*---------------------------------------------------------------------------*/

#include "fvMeshSubset.H"
#include "argList.H"
#include "cellSet.H"
#include "IOobjectList.H"
#include "volFields.H"
#include "polyTopoChangeMap.H"
#include "hexRef8Data.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
void subsetFields
(
    const fvMeshSubset& subsetter,
    const typename GeoField::Mesh& mesh,
    const wordList& fieldNames,
    PtrList<GeoField>& subFields
)
{
    forAll(fieldNames, i)
    {
        Info<< "Subsetting " << GeoField::typeName
            << " " << fieldNames[i] << endl;

        GeoField fld
        (
            IOobject
            (
                fieldNames[i],
                subsetter.baseMesh().time().name(),
                subsetter.baseMesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        subFields.set(i, subsetter.interpolate(fld));
    }
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "select a mesh subset based on a cellSet"
    );

    #include "addOverwriteOption.H"
    #include "addMeshOption.H"
    #include "addRegionOption.H"
    argList::validArgs.append("cellSet");
    argList::addOption
    (
        "patch",
        "name",
        "add exposed internal faces to specified patch instead of to "
        "'oldInternalFaces'"
    );
    argList::addOption
    (
        "resultTime",
        "time",
        "specify a time for the resulting mesh"
    );
    argList::addBoolOption
    (
        "noFields",
        "do not update fields"
    );

    #include "setRootCase.H"
    #include "createTimeNoFunctionObjects.H"

    Foam::word meshRegionName = polyMesh::defaultRegion;
    args.optionReadIfPresent("region", meshRegionName);

    #include "createSpecifiedMeshNoChangers.H"


    const word setName = args[1];

    word meshInstance = mesh.pointsInstance();
    word fieldsInstance = runTime.name();

    const bool overwrite = args.optionFound("overwrite");
    const bool specifiedInstance = args.optionReadIfPresent
    (
        "resultTime",
        fieldsInstance
    );
    if (specifiedInstance)
    {
        // Set both mesh and field to this time
        meshInstance = fieldsInstance;
    }
    const bool fields = !args.optionFound("noFields");


    Info<< "Reading cell set from " << setName << endl << endl;

    // Create mesh subsetting engine
    fvMeshSubset subsetter(mesh);

    label patchi = -1;

    if (args.optionFound("patch"))
    {
        const word patchName = args["patch"];

        patchi = mesh.boundaryMesh().findIndex(patchName);

        if (patchi == -1)
        {
            FatalErrorInFunction
                << nl << "Valid patches are " << mesh.boundaryMesh().names()
                << exit(FatalError);
        }

        Info<< "Adding exposed internal faces to patch " << patchName << endl
            << endl;
    }
    else
    {
        Info<< "Adding exposed internal faces to a patch called"
            << " \"oldInternalFaces\" (created if necessary)" << endl
            << endl;
    }


    // Read hexRef8 data, if any
    hexRef8Data refData
    (
        IOobject
        (
            "dummy",
            mesh.facesInstance(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );


    cellSet currentSet(mesh, setName);
    subsetter.setLargeCellSubset(currentSet, patchi, true);


    if (fields)
    {
        const pointMesh& pMesh = pointMesh::New(mesh);

        IOobjectList objects(mesh, runTime.name());

        // Subset all the fields
        #define SubsetTypeGeoFields(Type, GeoField, mesh)                      \
            wordList Type##GeoField##Names                                     \
            (                                                                  \
                objects.names(GeoField<Type>::typeName)                        \
            );                                                                 \
            PtrList<GeoField<Type>> Type##GeoField##s                          \
            (                                                                  \
                Type##GeoField##Names.size()                                   \
            );                                                                 \
            subsetFields                                                       \
            (                                                                  \
                subsetter,                                                     \
                mesh,                                                          \
                Type##GeoField##Names,                                         \
                Type##GeoField##s                                              \
            );
        FOR_ALL_FIELD_TYPES(SubsetTypeGeoFields, VolField, mesh);
        FOR_ALL_FIELD_TYPES(SubsetTypeGeoFields, VolInternalField, mesh);
        FOR_ALL_FIELD_TYPES(SubsetTypeGeoFields, SurfaceField, mesh);
        FOR_ALL_FIELD_TYPES(SubsetTypeGeoFields, PointField, pMesh);

        // Set the instance
        if (overwrite || specifiedInstance)
        {
            runTime.setTime(instant(fieldsInstance), 0);
            subsetter.subMesh().setInstance(meshInstance);
        }
        else
        {
            runTime++;
        }

        Info<< "Writing subsetted mesh and fields to time "
            << runTime.name() << endl;

        subsetter.subMesh().write();

        // Subsetting adds a 'subset' prefix. Rename the fields to remove this.
        #define RenameTypeGeoFields(Type, GeoField)                            \
            forAll(Type##GeoField##s, i)                                       \
            {                                                                  \
                Type##GeoField##s[i].rename(Type##GeoField##Names[i]);         \
                Type##GeoField##s[i].write();                                  \
            }
        FOR_ALL_FIELD_TYPES(RenameTypeGeoFields, VolField);
        FOR_ALL_FIELD_TYPES(RenameTypeGeoFields, VolInternalField);
        FOR_ALL_FIELD_TYPES(RenameTypeGeoFields, SurfaceField);
        FOR_ALL_FIELD_TYPES(RenameTypeGeoFields, PointField);
    }
    else
    {
        // Set the instance
        if (overwrite || specifiedInstance)
        {
            runTime.setTime(instant(fieldsInstance), 0);
            subsetter.subMesh().setInstance(meshInstance);
        }
        else
        {
            runTime++;
        }

        Info<< "Writing subsetted mesh to time "
            << runTime.name() << endl;

        subsetter.subMesh().write();
    }


    // Map the hexRef8 data
    polyTopoChangeMap map
    (
        mesh,
        mesh.nPoints(),                 // nOldPoints,
        mesh.nFaces(),                  // nOldFaces,
        mesh.nCells(),                  // nOldCells,
        labelList(subsetter.pointMap()),// pointMap,
        List<objectMap>(0),             // pointsFromPoints,
        labelList(0),                   // faceMap,
        List<objectMap>(0),             // facesFromFaces,
        labelList(subsetter.cellMap()), // cellMap,
        List<objectMap>(0),             // cellsFromCells,
        labelList(0),                   // reversePointMap,
        labelList(0),                   // reverseFaceMap,
        labelList(0),                   // reverseCellMap,
        labelHashSet(0),                // flipFaceFlux,
        labelListList(0),               // patchPointMap,
        labelList(0),                   // oldPatchSizes,
        labelList(0),                   // oldPatchStarts,
        labelList(0),                   // oldPatchNMeshPoints,
        autoPtr<scalarField>()          // oldCellVolumesPtr
    );
    refData.topoChange(map);
    refData.write();


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
