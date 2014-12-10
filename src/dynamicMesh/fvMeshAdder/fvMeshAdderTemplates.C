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

#include "volFields.H"
#include "surfaceFields.H"
#include "emptyFvPatchField.H"
#include "directFvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvMeshAdder::MapVolField
(
    const mapAddedPolyMesh& meshMap,

    GeometricField<Type, fvPatchField, volMesh>& fld,
    const GeometricField<Type, fvPatchField, volMesh>& fldToAdd
)
{
    const fvMesh& mesh = fld.mesh();

    // Internal field
    // ~~~~~~~~~~~~~~

    {
        // Store old internal field
        Field<Type> oldInternalField(fld.internalField());

        // Modify internal field
        Field<Type>& intFld = fld.internalField();

        intFld.setSize(mesh.nCells());

        intFld.rmap(oldInternalField, meshMap.oldCellMap());
        intFld.rmap(fldToAdd.internalField(), meshMap.addedCellMap());
    }


    // Patch fields from old mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    typename GeometricField<Type, fvPatchField, volMesh>::
    GeometricBoundaryField& bfld = fld.boundaryField();

    {
        const labelList& oldPatchMap = meshMap.oldPatchMap();
        const labelList& oldPatchStarts = meshMap.oldPatchStarts();
        const labelList& oldPatchSizes = meshMap.oldPatchSizes();

        // Reorder old patches in order of new ones. Put removed patches at end.

        label unusedPatchI = 0;

        forAll(oldPatchMap, patchI)
        {
            label newPatchI = oldPatchMap[patchI];

            if (newPatchI != -1)
            {
                unusedPatchI++;
            }
        }

        label nUsedPatches = unusedPatchI;

        // Reorder list for patchFields
        labelList oldToNew(oldPatchMap.size());

        forAll(oldPatchMap, patchI)
        {
            label newPatchI = oldPatchMap[patchI];

            if (newPatchI != -1)
            {
                oldToNew[patchI] = newPatchI;
            }
            else
            {
                oldToNew[patchI] = unusedPatchI++;
            }
        }


        // Sort deleted ones last so is now in newPatch ordering
        bfld.reorder(oldToNew);
        // Extend to covers all patches
        bfld.setSize(mesh.boundaryMesh().size());
        // Delete unused patches
        for
        (
            label newPatchI = nUsedPatches;
            newPatchI < bfld.size();
            newPatchI++
        )
        {
            bfld.set(newPatchI, NULL);
        }


        // Map old values
        // ~~~~~~~~~~~~~~

        forAll(oldPatchMap, patchI)
        {
            label newPatchI = oldPatchMap[patchI];

            if (newPatchI != -1)
            {
                labelList newToOld
                (
                    calcPatchMap
                    (
                        oldPatchStarts[patchI],
                        oldPatchSizes[patchI],
                        meshMap.oldFaceMap(),
                        mesh.boundaryMesh()[newPatchI],
                        -1              // unmapped value
                    )
                );

                directFvPatchFieldMapper patchMapper(newToOld);


                // Create new patchField with same type as existing one.
                // Note:
                // - boundaryField already in new order so access with newPatchI
                // - fld.boundaryField()[newPatchI] both used for type and old
                //   value
                // - hope that field mapping allows aliasing since old and new
                //   are same memory!
                bfld.set
                (
                    newPatchI,
                    fvPatchField<Type>::New
                    (
                        bfld[newPatchI],                // old field
                        mesh.boundary()[newPatchI],     // new fvPatch
                        fld.dimensionedInternalField(), // new internal field
                        patchMapper                     // mapper (new to old)
                    )
                );
            }
        }
    }



    // Patch fields from added mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        const labelList& addedPatchMap = meshMap.addedPatchMap();

        // Add addedMesh patches
        forAll(addedPatchMap, patchI)
        {
            label newPatchI = addedPatchMap[patchI];

            if (newPatchI != -1)
            {
                const polyPatch& newPatch = mesh.boundaryMesh()[newPatchI];
                const polyPatch& oldPatch =
                    fldToAdd.mesh().boundaryMesh()[patchI];

                if (!bfld(newPatchI))
                {
                    // First occurrence of newPatchI. Map from existing
                    // patchField

                    // From new patch faces to patch faces on added mesh.
                    labelList newToAdded
                    (
                        calcPatchMap
                        (
                            oldPatch.start(),
                            oldPatch.size(),
                            meshMap.addedFaceMap(),
                            newPatch,
                            -1          // unmapped values
                        )
                    );

                    directFvPatchFieldMapper patchMapper(newToAdded);

                    bfld.set
                    (
                        newPatchI,
                        fvPatchField<Type>::New
                        (
                            fldToAdd.boundaryField()[patchI], // added field
                            mesh.boundary()[newPatchI],       // new fvPatch
                            fld.dimensionedInternalField(),   // new int. field
                            patchMapper                       // mapper
                        )
                    );
                }
                else
                {
                    // PatchField will have correct size already. Just slot in
                    // my elements.

                    labelList addedToNew(oldPatch.size(), -1);
                    forAll(addedToNew, i)
                    {
                        label addedFaceI = oldPatch.start()+i;
                        label newFaceI = meshMap.addedFaceMap()[addedFaceI];
                        label patchFaceI = newFaceI-newPatch.start();
                        if (patchFaceI >= 0 && patchFaceI < newPatch.size())
                        {
                            addedToNew[i] = patchFaceI;
                        }
                    }

                    bfld[newPatchI].rmap
                    (
                        fldToAdd.boundaryField()[patchI],
                        addedToNew
                    );
                }
            }
        }
    }
}


template<class Type>
void Foam::fvMeshAdder::MapVolFields
(
    const mapAddedPolyMesh& meshMap,
    const fvMesh& mesh,
    const fvMesh& meshToAdd
)
{
    HashTable<const GeometricField<Type, fvPatchField, volMesh>*> fields
    (
        mesh.objectRegistry::lookupClass
        <GeometricField<Type, fvPatchField, volMesh> >
        ()
    );

    HashTable<const GeometricField<Type, fvPatchField, volMesh>*> fieldsToAdd
    (
        meshToAdd.objectRegistry::lookupClass
        <GeometricField<Type, fvPatchField, volMesh> >
        ()
    );

    // It is necessary to enforce that all old-time fields are stored
    // before the mapping is performed.  Otherwise, if the
    // old-time-level field is mapped before the field itself, sizes
    // will not match.

    for
    (
        typename HashTable<const GeometricField<Type, fvPatchField, volMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        const_cast<GeometricField<Type, fvPatchField, volMesh>*>(fieldIter())
            ->storeOldTimes();
    }


    for
    (
        typename HashTable<const GeometricField<Type, fvPatchField, volMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        GeometricField<Type, fvPatchField, volMesh>& fld =
            const_cast<GeometricField<Type, fvPatchField, volMesh>&>
            (
                *fieldIter()
            );

        if (fieldsToAdd.found(fld.name()))
        {
            const GeometricField<Type, fvPatchField, volMesh>& fldToAdd =
                *fieldsToAdd[fld.name()];

            MapVolField<Type>(meshMap, fld, fldToAdd);
        }
        else
        {
            WarningIn("fvMeshAdder::MapVolFields(..)")
                << "Not mapping field " << fld.name()
                << " since not present on mesh to add"
                << endl;
        }
    }
}


template<class Type>
void Foam::fvMeshAdder::MapSurfaceField
(
    const mapAddedPolyMesh& meshMap,

    GeometricField<Type, fvsPatchField, surfaceMesh>& fld,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& fldToAdd
)
{
    const fvMesh& mesh = fld.mesh();
    const labelList& oldPatchStarts = meshMap.oldPatchStarts();

    typename GeometricField<Type, fvsPatchField, surfaceMesh>::
    GeometricBoundaryField& bfld = fld.boundaryField();

    // Internal field
    // ~~~~~~~~~~~~~~

    // Store old internal field
    {
        Field<Type> oldField(fld);

        // Modify internal field
        Field<Type>& intFld = fld.internalField();

        intFld.setSize(mesh.nInternalFaces());

        intFld.rmap(oldField, meshMap.oldFaceMap());
        intFld.rmap(fldToAdd, meshMap.addedFaceMap());


        // Faces that were boundary faces but are not anymore.
        // Use owner value (so lowest numbered cell, i.e. from 'old' not 'added'
        // mesh)
        forAll(bfld, patchI)
        {
            const fvsPatchField<Type>& pf = bfld[patchI];

            label start = oldPatchStarts[patchI];

            forAll(pf, i)
            {
                label newFaceI = meshMap.oldFaceMap()[start + i];

                if (newFaceI >= 0 && newFaceI < mesh.nInternalFaces())
                {
                    intFld[newFaceI] = pf[i];
                }
            }
        }
    }


    // Patch fields from old mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        const labelList& oldPatchMap = meshMap.oldPatchMap();
        const labelList& oldPatchSizes = meshMap.oldPatchSizes();

        // Reorder old patches in order of new ones. Put removed patches at end.

        label unusedPatchI = 0;

        forAll(oldPatchMap, patchI)
        {
            label newPatchI = oldPatchMap[patchI];

            if (newPatchI != -1)
            {
                unusedPatchI++;
            }
        }

        label nUsedPatches = unusedPatchI;

        // Reorder list for patchFields
        labelList oldToNew(oldPatchMap.size());

        forAll(oldPatchMap, patchI)
        {
            label newPatchI = oldPatchMap[patchI];

            if (newPatchI != -1)
            {
                oldToNew[patchI] = newPatchI;
            }
            else
            {
                oldToNew[patchI] = unusedPatchI++;
            }
        }


        // Sort deleted ones last so is now in newPatch ordering
        bfld.reorder(oldToNew);
        // Extend to covers all patches
        bfld.setSize(mesh.boundaryMesh().size());
        // Delete unused patches
        for
        (
            label newPatchI = nUsedPatches;
            newPatchI < bfld.size();
            newPatchI++
        )
        {
            bfld.set(newPatchI, NULL);
        }


        // Map old values
        // ~~~~~~~~~~~~~~

        forAll(oldPatchMap, patchI)
        {
            label newPatchI = oldPatchMap[patchI];

            if (newPatchI != -1)
            {
                labelList newToOld
                (
                    calcPatchMap
                    (
                        oldPatchStarts[patchI],
                        oldPatchSizes[patchI],
                        meshMap.oldFaceMap(),
                        mesh.boundaryMesh()[newPatchI],
                        -1      // unmapped value
                    )
                );

                directFvPatchFieldMapper patchMapper(newToOld);

                // Create new patchField with same type as existing one.
                // Note:
                // - boundaryField already in new order so access with newPatchI
                // - bfld[newPatchI] both used for type and old
                //   value
                // - hope that field mapping allows aliasing since old and new
                //   are same memory!
                bfld.set
                (
                    newPatchI,
                    fvsPatchField<Type>::New
                    (
                        bfld[newPatchI],                // old field
                        mesh.boundary()[newPatchI],     // new fvPatch
                        fld.dimensionedInternalField(), // new internal field
                        patchMapper                     // mapper (new to old)
                    )
                );
            }
        }
    }



    // Patch fields from added mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        const labelList& addedPatchMap = meshMap.addedPatchMap();

        // Add addedMesh patches
        forAll(addedPatchMap, patchI)
        {
            label newPatchI = addedPatchMap[patchI];

            if (newPatchI != -1)
            {
                const polyPatch& newPatch = mesh.boundaryMesh()[newPatchI];
                const polyPatch& oldPatch =
                    fldToAdd.mesh().boundaryMesh()[patchI];

                if (!bfld(newPatchI))
                {
                    // First occurrence of newPatchI. Map from existing
                    // patchField

                    // From new patch faces to patch faces on added mesh.
                    labelList newToAdded
                    (
                        calcPatchMap
                        (
                            oldPatch.start(),
                            oldPatch.size(),
                            meshMap.addedFaceMap(),
                            newPatch,
                            -1                  // unmapped values
                        )
                    );

                    directFvPatchFieldMapper patchMapper(newToAdded);

                    bfld.set
                    (
                        newPatchI,
                        fvsPatchField<Type>::New
                        (
                            fldToAdd.boundaryField()[patchI],// added field
                            mesh.boundary()[newPatchI],      // new fvPatch
                            fld.dimensionedInternalField(),  // new int. field
                            patchMapper                      // mapper
                        )
                    );
                }
                else
                {
                    // PatchField will have correct size already. Just slot in
                    // my elements.

                    labelList addedToNew(oldPatch.size(), -1);
                    forAll(addedToNew, i)
                    {
                        label addedFaceI = oldPatch.start()+i;
                        label newFaceI = meshMap.addedFaceMap()[addedFaceI];
                        label patchFaceI = newFaceI-newPatch.start();
                        if (patchFaceI >= 0 && patchFaceI < newPatch.size())
                        {
                            addedToNew[i] = patchFaceI;
                        }
                    }

                    bfld[newPatchI].rmap
                    (
                        fldToAdd.boundaryField()[patchI],
                        addedToNew
                    );
                }
            }
        }
    }
}


template<class Type>
void Foam::fvMeshAdder::MapSurfaceFields
(
    const mapAddedPolyMesh& meshMap,
    const fvMesh& mesh,
    const fvMesh& meshToAdd
)
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fldType;

    HashTable<const fldType*> fields
    (
        mesh.objectRegistry::lookupClass<fldType>()
    );

    HashTable<const fldType*> fieldsToAdd
    (
        meshToAdd.objectRegistry::lookupClass<fldType>()
    );

    // It is necessary to enforce that all old-time fields are stored
    // before the mapping is performed.  Otherwise, if the
    // old-time-level field is mapped before the field itself, sizes
    // will not match.

    for
    (
        typename HashTable<const fldType*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        const_cast<fldType*>(fieldIter())
            ->storeOldTimes();
    }


    for
    (
        typename HashTable<const fldType*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        fldType& fld = const_cast<fldType&>(*fieldIter());

        if (fieldsToAdd.found(fld.name()))
        {
            const fldType& fldToAdd = *fieldsToAdd[fld.name()];

            MapSurfaceField<Type>(meshMap, fld, fldToAdd);
        }
        else
        {
            WarningIn("fvMeshAdder::MapSurfaceFields(..)")
                << "Not mapping field " << fld.name()
                << " since not present on mesh to add"
                << endl;
        }
    }
}


// ************************************************************************* //
