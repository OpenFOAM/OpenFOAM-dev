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

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "directFvPatchFieldMapper.H"
#include "directPointPatchFieldMapper.H"

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
        Field<Type> oldInternalField(fld.primitiveField());

        // Modify internal field
        Field<Type>& intFld = fld.primitiveFieldRef();

        intFld.setSize(mesh.nCells());

        intFld.rmap(oldInternalField, meshMap.oldCellMap());
        intFld.rmap(fldToAdd.primitiveField(), meshMap.addedCellMap());
    }


    // Patch fields from old mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    typename GeometricField<Type, fvPatchField, volMesh>::
    Boundary& bfld = fld.boundaryFieldRef();

    {
        const labelList& oldPatchMap = meshMap.oldPatchMap();
        const labelList& oldPatchStarts = meshMap.oldPatchStarts();
        const labelList& oldPatchSizes = meshMap.oldPatchSizes();

        // Reorder old patches in order of new ones. Put removed patches at end.

        label unusedPatchi = 0;

        forAll(oldPatchMap, patchi)
        {
            const label newPatchi = oldPatchMap[patchi];

            if (newPatchi != -1)
            {
                unusedPatchi++;
            }
        }

        label nUsedPatches = unusedPatchi;

        // Reorder list for patchFields
        labelList oldToNew(oldPatchMap.size());

        forAll(oldPatchMap, patchi)
        {
            const label newPatchi = oldPatchMap[patchi];

            if (newPatchi != -1)
            {
                oldToNew[patchi] = newPatchi;
            }
            else
            {
                oldToNew[patchi] = unusedPatchi++;
            }
        }


        // Sort deleted ones last so is now in newPatch ordering
        bfld.reorder(oldToNew);
        // Extend to covers all patches
        bfld.setSize(mesh.boundaryMesh().size());
        // Delete unused patches
        for
        (
            label newPatchi = nUsedPatches;
            newPatchi < bfld.size();
            newPatchi++
        )
        {
            bfld.set(newPatchi, nullptr);
        }


        // Map old values
        // ~~~~~~~~~~~~~~

        forAll(oldPatchMap, patchi)
        {
            const label newPatchi = oldPatchMap[patchi];

            if (newPatchi != -1)
            {
                labelList newToOld
                (
                    calcPatchMap
                    (
                        oldPatchStarts[patchi],
                        oldPatchSizes[patchi],
                        meshMap.oldFaceMap(),
                        mesh.boundaryMesh()[newPatchi],
                        -1              // unmapped value
                    )
                );

                directFvPatchFieldMapper patchMapper(newToOld);


                // Create new patchField with same type as existing one.
                // Note:
                // - boundaryField already in new order so access with newPatchi
                // - fld.boundaryField()[newPatchi] both used for type and old
                //   value
                // - hope that field mapping allows aliasing since old and new
                //   are same memory!
                bfld.set
                (
                    newPatchi,
                    fvPatchField<Type>::New
                    (
                        bfld[newPatchi],                // old field
                        mesh.boundary()[newPatchi],     // new fvPatch
                        fld(), // new internal field
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
        forAll(addedPatchMap, patchi)
        {
            const label newPatchi = addedPatchMap[patchi];

            if (newPatchi != -1)
            {
                const polyPatch& newPatch = mesh.boundaryMesh()[newPatchi];
                const polyPatch& oldPatch =
                    fldToAdd.mesh().boundaryMesh()[patchi];

                if (!bfld(newPatchi))
                {
                    // First occurrence of newPatchi. Map from existing
                    // patchField

                    // From new patch faces to patch faces on added mesh.
                    const labelList newToAdded
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
                        newPatchi,
                        fvPatchField<Type>::New
                        (
                            fldToAdd.boundaryField()[patchi], // added field
                            mesh.boundary()[newPatchi],       // new fvPatch
                            fld(),   // new int. field
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
                        const label addedFacei = oldPatch.start() + i;
                        const label newFacei =
                            meshMap.addedFaceMap()[addedFacei];
                        const label patchFacei = newFacei-newPatch.start();

                        if (patchFacei >= 0 && patchFacei < newPatch.size())
                        {
                            addedToNew[i] = patchFacei;
                        }
                    }

                    bfld[newPatchi].rmap
                    (
                        fldToAdd.boundaryField()[patchi],
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
        <GeometricField<Type, fvPatchField, volMesh>>
        ()
    );

    HashTable<const GeometricField<Type, fvPatchField, volMesh>*> fieldsToAdd
    (
        meshToAdd.objectRegistry::lookupClass
        <GeometricField<Type, fvPatchField, volMesh>>
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
        if (debug)
        {
            Pout<< "MapVolFields : Storing old time for " << fieldIter()->name()
                << endl;
        }

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

            if (debug)
            {
                Pout<< "MapVolFields : mapping " << fld.name()
                    << " and " << fldToAdd.name() << endl;
            }

            MapVolField<Type>(meshMap, fld, fldToAdd);
        }
        else
        {
            WarningInFunction
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
    Boundary& bfld = fld.boundaryFieldRef();

    // Internal field
    // ~~~~~~~~~~~~~~

    // Store old internal field
    {
        Field<Type> oldField(fld);

        // Modify internal field
        Field<Type>& intFld = fld.primitiveFieldRef();

        intFld.setSize(mesh.nInternalFaces());

        intFld.rmap(oldField, meshMap.oldFaceMap());
        intFld.rmap(fldToAdd, meshMap.addedFaceMap());


        // Faces that were boundary faces but are not anymore.
        // Use owner value (so lowest numbered cell, i.e. from 'old' not 'added'
        // mesh)
        forAll(bfld, patchi)
        {
            const fvsPatchField<Type>& pf = bfld[patchi];

            const label start = oldPatchStarts[patchi];

            forAll(pf, i)
            {
                const label newFacei = meshMap.oldFaceMap()[start + i];

                if (newFacei >= 0 && newFacei < mesh.nInternalFaces())
                {
                    intFld[newFacei] = pf[i];
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

        label unusedPatchi = 0;

        forAll(oldPatchMap, patchi)
        {
            const label newPatchi = oldPatchMap[patchi];

            if (newPatchi != -1)
            {
                unusedPatchi++;
            }
        }

        label nUsedPatches = unusedPatchi;

        // Reorder list for patchFields
        labelList oldToNew(oldPatchMap.size());

        forAll(oldPatchMap, patchi)
        {
            const label newPatchi = oldPatchMap[patchi];

            if (newPatchi != -1)
            {
                oldToNew[patchi] = newPatchi;
            }
            else
            {
                oldToNew[patchi] = unusedPatchi++;
            }
        }


        // Sort deleted ones last so is now in newPatch ordering
        bfld.reorder(oldToNew);
        // Extend to covers all patches
        bfld.setSize(mesh.boundaryMesh().size());
        // Delete unused patches
        for
        (
            label newPatchi = nUsedPatches;
            newPatchi < bfld.size();
            newPatchi++
        )
        {
            bfld.set(newPatchi, nullptr);
        }


        // Map old values
        // ~~~~~~~~~~~~~~

        forAll(oldPatchMap, patchi)
        {
            const label newPatchi = oldPatchMap[patchi];

            if (newPatchi != -1)
            {
                const labelList newToOld
                (
                    calcPatchMap
                    (
                        oldPatchStarts[patchi],
                        oldPatchSizes[patchi],
                        meshMap.oldFaceMap(),
                        mesh.boundaryMesh()[newPatchi],
                        -1      // unmapped value
                    )
                );

                directFvPatchFieldMapper patchMapper(newToOld);

                // Create new patchField with same type as existing one.
                // Note:
                // - boundaryField already in new order so access with newPatchi
                // - bfld[newPatchi] both used for type and old
                //   value
                // - hope that field mapping allows aliasing since old and new
                //   are same memory!
                bfld.set
                (
                    newPatchi,
                    fvsPatchField<Type>::New
                    (
                        bfld[newPatchi],                // old field
                        mesh.boundary()[newPatchi],     // new fvPatch
                        fld(), // new internal field
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
        forAll(addedPatchMap, patchi)
        {
            const label newPatchi = addedPatchMap[patchi];

            if (newPatchi != -1)
            {
                const polyPatch& newPatch = mesh.boundaryMesh()[newPatchi];
                const polyPatch& oldPatch =
                    fldToAdd.mesh().boundaryMesh()[patchi];

                if (!bfld(newPatchi))
                {
                    // First occurrence of newPatchi. Map from existing
                    // patchField

                    // From new patch faces to patch faces on added mesh.
                    const labelList newToAdded
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
                        newPatchi,
                        fvsPatchField<Type>::New
                        (
                            fldToAdd.boundaryField()[patchi],// added field
                            mesh.boundary()[newPatchi],      // new fvPatch
                            fld(),  // new int. field
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
                        const label addedFacei = oldPatch.start() + i;
                        const label newFacei =
                            meshMap.addedFaceMap()[addedFacei];
                        const label patchFacei = newFacei-newPatch.start();

                        if (patchFacei >= 0 && patchFacei < newPatch.size())
                        {
                            addedToNew[i] = patchFacei;
                        }
                    }

                    bfld[newPatchi].rmap
                    (
                        fldToAdd.boundaryField()[patchi],
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
        if (debug)
        {
            Pout<< "MapSurfaceFields : Storing old time for "
                << fieldIter()->name() << endl;
        }

        const_cast<fldType*>(fieldIter())->storeOldTimes();
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

            if (debug)
            {
                Pout<< "MapSurfaceFields : mapping " << fld.name()
                    << " and " << fldToAdd.name() << endl;
            }

            MapSurfaceField<Type>(meshMap, fld, fldToAdd);
        }
        else
        {
            WarningInFunction
                << "Not mapping field " << fld.name()
                << " since not present on mesh to add"
                << endl;
        }
    }
}


template<class Type>
void Foam::fvMeshAdder::MapPointField
(
    const pointMesh& mesh,
    const mapAddedPolyMesh& meshMap,
    const labelListList& oldMeshPoints,

    GeometricField<Type, pointPatchField, pointMesh>& fld,
    const GeometricField<Type, pointPatchField, pointMesh>& fldToAdd
)
{
    // This is a bit tricky:
    // - mesh pointed to by fld is invalid
    // - pointPatches pointed to be fld are invalid

    typename GeometricField<Type, pointPatchField, pointMesh>::
    Boundary& bfld = fld.boundaryFieldRef();

    // Internal field
    // ~~~~~~~~~~~~~~

    // Store old internal field
    {
        Field<Type> oldField(fld);

        // Modify internal field
        Field<Type>& intFld = fld.primitiveFieldRef();

        intFld.setSize(mesh.size());

        intFld.rmap(oldField, meshMap.oldPointMap());
        intFld.rmap(fldToAdd, meshMap.addedPointMap());
    }


    // Patch fields from old mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        const labelList& oldPatchMap = meshMap.oldPatchMap();

        // Reorder old patches in order of new ones. Put removed patches at end.

        label unusedPatchi = 0;

        forAll(oldPatchMap, patchi)
        {
            const label newPatchi = oldPatchMap[patchi];

            if (newPatchi != -1)
            {
                unusedPatchi++;
            }
        }

        label nUsedPatches = unusedPatchi;

        // Reorder list for patchFields
        labelList oldToNew(oldPatchMap.size());

        forAll(oldPatchMap, patchi)
        {
            const label newPatchi = oldPatchMap[patchi];

            if (newPatchi != -1)
            {
                oldToNew[patchi] = newPatchi;
            }
            else
            {
                oldToNew[patchi] = unusedPatchi++;
            }
        }


        // Sort deleted ones last so is now in newPatch ordering
        bfld.reorder(oldToNew);

        // Extend to covers all patches
        bfld.setSize(mesh.boundary().size());

        // Delete unused patches
        for
        (
            label newPatchi = nUsedPatches;
            newPatchi < bfld.size();
            newPatchi++
        )
        {
            bfld.set(newPatchi, nullptr);
        }


        // Map old values
        // ~~~~~~~~~~~~~~

        forAll(oldPatchMap, patchi)
        {
            const label newPatchi = oldPatchMap[patchi];

            if (newPatchi != -1)
            {
                const labelList& oldMp = oldMeshPoints[patchi];
                const pointPatch& newPp = mesh.boundary()[newPatchi];
                const labelList& newMeshPoints = newPp.meshPoints();
                const labelList& oldPointMap = meshMap.oldPointMap();

                Map<label> newMeshPointMap(2*newMeshPoints.size());
                forAll(newMeshPoints, ppi)
                {
                    newMeshPointMap.insert(newMeshPoints[ppi], ppi);
                }

                labelList newToOld(newPp.size(), -1);
                forAll(oldMp, oldPointi)
                {
                    const label newPointi = oldPointMap[oldMp[oldPointi]];

                    Map<label>::const_iterator fnd =
                        newMeshPointMap.find(newPointi);

                    if (fnd == newMeshPointMap.end())
                    {
                        // Possibly an internal point
                    }
                    else
                    {
                        // Part of new patch
                        newToOld[fnd()] = oldPointi;
                    }
                }

                directPointPatchFieldMapper patchMapper(newToOld);

                // Create new patchField with same type as existing one.
                // Note:
                // - boundaryField already in new order so access with newPatchi
                // - bfld[newPatchi] both used for type and old
                //   value
                // - hope that field mapping allows aliasing since old and new
                //   are same memory!
                bfld.set
                (
                    newPatchi,
                    pointPatchField<Type>::New
                    (
                        bfld[newPatchi],                // old field
                        mesh.boundary()[newPatchi],     // new pointPatch
                        fld(),                          // new internal field
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
        forAll(addedPatchMap, patchi)
        {
            const label newPatchi = addedPatchMap[patchi];

            if (newPatchi != -1)
            {
                const pointPatch& oldPatch = fldToAdd.mesh().boundary()[patchi];
                const labelList& oldMp = oldPatch.meshPoints();
                const pointPatch& newPp = mesh.boundary()[newPatchi];
                const labelList& newMeshPoints = newPp.meshPoints();
                const labelList& addedPointMap = meshMap.addedPointMap();

                Map<label> newMpm(2*newMeshPoints.size());
                forAll(newMeshPoints, ppi)
                {
                    newMpm.insert(newMeshPoints[ppi], ppi);
                }

                if (!bfld(newPatchi))
                {
                    // First occurrence of newPatchi. Map from existing
                    // patchField

                    labelList newToAdded(newPp.size(), -1);
                    forAll(oldMp, oldPointi)
                    {
                        const label newPointi = addedPointMap[oldMp[oldPointi]];

                        Map<label>::const_iterator fnd = newMpm.find(newPointi);
                        if (fnd == newMpm.end())
                        {
                            // Possibly an internal point
                        }
                        else
                        {
                            // Part of new patch
                            newToAdded[fnd()] = oldPointi;
                        }
                    }

                    directPointPatchFieldMapper patchMapper(newToAdded);

                    bfld.set
                    (
                        newPatchi,
                        pointPatchField<Type>::New
                        (
                            fldToAdd.boundaryField()[patchi],// added field
                            mesh.boundary()[newPatchi],      // new pointPatch
                            fld(),                           // new int. field
                            patchMapper                      // mapper
                        )
                    );
                }
                else
                {
                    // PatchField will have correct size already. Just slot in
                    // my elements.

                    labelList oldToNew(oldPatch.size(), -1);
                    forAll(oldMp, oldPointi)
                    {
                        const label newPointi = addedPointMap[oldMp[oldPointi]];

                        Map<label>::const_iterator fnd = newMpm.find(newPointi);
                        if (fnd != newMpm.end())
                        {
                            // Part of new patch
                            oldToNew[oldPointi] = fnd();
                        }
                    }

                    bfld[newPatchi].rmap
                    (
                        fldToAdd.boundaryField()[patchi],
                        oldToNew
                    );
                }
            }
        }
    }
}


template<class Type>
void Foam::fvMeshAdder::MapPointFields
(
    const mapAddedPolyMesh& meshMap,
    const pointMesh& mesh,
    const labelListList& oldMeshPoints,
    const objectRegistry& meshToAdd
)
{
    typedef GeometricField<Type, pointPatchField, pointMesh> fldType;

    HashTable<const fldType*> fields(mesh.thisDb().lookupClass<fldType>());
    HashTable<const fldType*> fieldsToAdd(meshToAdd.lookupClass<fldType>());

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
        if (debug)
        {
            Pout<< "MapPointFields : Storing old time for "
                << fieldIter()->name() << endl;
        }

        const_cast<fldType*>(fieldIter())->storeOldTimes();
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

            if (debug)
            {
                Pout<< "MapPointFields : mapping " << fld.name()
                    << " and " << fldToAdd.name() << endl;
            }

            MapPointField<Type>(mesh, meshMap, oldMeshPoints, fld, fldToAdd);
        }
        else
        {
            WarningInFunction
                << "Not mapping field " << fld.name()
                << " since not present on mesh to add"
                << endl;
        }
    }
}


template<class Type>
void Foam::fvMeshAdder::MapDimField
(
    const mapAddedPolyMesh& meshMap,

    DimensionedField<Type, volMesh>& fld,
    const DimensionedField<Type, volMesh>& fldToAdd
)
{
    const fvMesh& mesh = fld.mesh();

    // Store old field
    Field<Type> oldField(fld);

    fld.setSize(mesh.nCells());

    fld.rmap(oldField, meshMap.oldCellMap());
    fld.rmap(fldToAdd, meshMap.addedCellMap());
}


template<class Type>
void Foam::fvMeshAdder::MapDimFields
(
    const mapAddedPolyMesh& meshMap,
    const fvMesh& mesh,
    const fvMesh& meshToAdd
)
{
    typedef DimensionedField<Type, volMesh> fldType;

    // Note: use strict flag on lookupClass to avoid picking up
    //       volFields
    HashTable<const fldType*> fields
    (
        mesh.objectRegistry::lookupClass<fldType>(true)
    );

    HashTable<const fldType*> fieldsToAdd
    (
        meshToAdd.objectRegistry::lookupClass<fldType>(true)
    );

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

            if (debug)
            {
                Pout<< "MapDimFields : mapping " << fld.name()
                    << " and " << fldToAdd.name() << endl;
            }

            MapDimField<Type>(meshMap, fld, fldToAdd);
        }
        else
        {
            WarningIn("fvMeshAdder::MapDimFields(..)")
                << "Not mapping field " << fld.name()
                << " since not present on mesh to add"
                << endl;
        }
    }
}


// ************************************************************************* //
