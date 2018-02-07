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

\*---------------------------------------------------------------------------*/

#include "fvFieldReconstructor.H"
#include "Time.H"
#include "PtrList.H"
#include "fvPatchFields.H"
#include "emptyFvPatch.H"
#include "emptyFvPatchField.H"
#include "emptyFvsPatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::volMesh>>
Foam::fvFieldReconstructor::reconstructFvVolumeInternalField
(
    const IOobject& fieldIoObject,
    const PtrList<DimensionedField<Type, volMesh>>& procFields
) const
{
    // Create the internalField
    Field<Type> internalField(mesh_.nCells());

    forAll(procMeshes_, proci)
    {
        const DimensionedField<Type, volMesh>& procField = procFields[proci];

        // Set the cell values in the reconstructed field
        internalField.rmap
        (
            procField.field(),
            cellProcAddressing_[proci]
        );
    }

    return tmp<DimensionedField<Type, volMesh>>
    (
        new DimensionedField<Type, volMesh>
        (
            fieldIoObject,
            mesh_,
            procFields[0].dimensions(),
            internalField
        )
    );
}


template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::volMesh>>
Foam::fvFieldReconstructor::reconstructFvVolumeInternalField
(
    const IOobject& fieldIoObject
) const
{
    // Read the field for all the processors
    PtrList<DimensionedField<Type, volMesh>> procFields
    (
        procMeshes_.size()
    );

    forAll(procMeshes_, proci)
    {
        procFields.set
        (
            proci,
            new DimensionedField<Type, volMesh>
            (
                IOobject
                (
                    fieldIoObject.name(),
                    procMeshes_[proci].time().timeName(),
                    procMeshes_[proci],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                procMeshes_[proci]
            )
        );
    }


    return reconstructFvVolumeInternalField
    (
        IOobject
        (
            fieldIoObject.name(),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        procFields
    );
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::fvFieldReconstructor::reconstructFvVolumeField
(
    const IOobject& fieldIoObject,
    const PtrList<GeometricField<Type, fvPatchField, volMesh>>& procFields
) const
{
    // Create the internalField
    Field<Type> internalField(mesh_.nCells());

    // Create the patch fields
    PtrList<fvPatchField<Type>> patchFields(mesh_.boundary().size());

    forAll(procFields, proci)
    {
        const GeometricField<Type, fvPatchField, volMesh>& procField =
            procFields[proci];

        // Set the cell values in the reconstructed field
        internalField.rmap
        (
            procField.primitiveField(),
            cellProcAddressing_[proci]
        );

        // Set the boundary patch values in the reconstructed field
        forAll(boundaryProcAddressing_[proci], patchi)
        {
            // Get patch index of the original patch
            const label curBPatch = boundaryProcAddressing_[proci][patchi];

            // Get addressing slice for this patch
            const labelList::subList cp =
                procField.mesh().boundary()[patchi].patchSlice
                (
                    faceProcAddressing_[proci]
                );

            // check if the boundary patch is not a processor patch
            if (curBPatch >= 0)
            {
                // Regular patch. Fast looping

                if (!patchFields(curBPatch))
                {
                    patchFields.set
                    (
                        curBPatch,
                        fvPatchField<Type>::New
                        (
                            procField.boundaryField()[patchi],
                            mesh_.boundary()[curBPatch],
                            DimensionedField<Type, volMesh>::null(),
                            fvPatchFieldReconstructor
                            (
                                mesh_.boundary()[curBPatch].size()
                            )
                        )
                    );
                }

                const label curPatchStart =
                    mesh_.boundaryMesh()[curBPatch].start();

                labelList reverseAddressing(cp.size());

                forAll(cp, facei)
                {
                    // Check
                    if (cp[facei] <= 0)
                    {
                        FatalErrorInFunction
                            << "Processor " << proci
                            << " patch "
                            << procField.mesh().boundary()[patchi].name()
                            << " face " << facei
                            << " originates from reversed face since "
                            << cp[facei]
                            << exit(FatalError);
                    }

                    // Subtract one to take into account offsets for
                    // face direction.
                    reverseAddressing[facei] = cp[facei] - 1 - curPatchStart;
                }


                patchFields[curBPatch].rmap
                (
                    procField.boundaryField()[patchi],
                    reverseAddressing
                );
            }
            else
            {
                const Field<Type>& curProcPatch =
                    procField.boundaryField()[patchi];

                // In processor patches, there's a mix of internal faces (some
                // of them turned) and possible cyclics. Slow loop
                forAll(cp, facei)
                {
                    // Subtract one to take into account offsets for
                    // face direction.
                    label curF = cp[facei] - 1;

                    // Is the face on the boundary?
                    if (curF >= mesh_.nInternalFaces())
                    {
                        label curBPatch = mesh_.boundaryMesh().whichPatch(curF);

                        if (!patchFields(curBPatch))
                        {
                            patchFields.set
                            (
                                curBPatch,
                                fvPatchField<Type>::New
                                (
                                    mesh_.boundary()[curBPatch].type(),
                                    mesh_.boundary()[curBPatch],
                                    DimensionedField<Type, volMesh>::null()
                                )
                            );
                        }

                        // add the face
                        label curPatchFace =
                            mesh_.boundaryMesh()
                                [curBPatch].whichFace(curF);

                        patchFields[curBPatch][curPatchFace] =
                            curProcPatch[facei];
                    }
                }
            }
        }
    }

    forAll(mesh_.boundary(), patchi)
    {
        // add empty patches
        if
        (
            isType<emptyFvPatch>(mesh_.boundary()[patchi])
         && !patchFields(patchi)
        )
        {
            patchFields.set
            (
                patchi,
                fvPatchField<Type>::New
                (
                    emptyFvPatchField<Type>::typeName,
                    mesh_.boundary()[patchi],
                    DimensionedField<Type, volMesh>::null()
                )
            );
        }
    }


    // Now construct and write the field
    // setting the internalField and patchFields
    return tmp<GeometricField<Type, fvPatchField, volMesh>>
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            fieldIoObject,
            mesh_,
            procFields[0].dimensions(),
            internalField,
            patchFields
        )
    );
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::fvFieldReconstructor::reconstructFvVolumeField
(
    const IOobject& fieldIoObject
) const
{
    // Read the field for all the processors
    PtrList<GeometricField<Type, fvPatchField, volMesh>> procFields
    (
        procMeshes_.size()
    );

    forAll(procMeshes_, proci)
    {
        procFields.set
        (
            proci,
            new GeometricField<Type, fvPatchField, volMesh>
            (
                IOobject
                (
                    fieldIoObject.name(),
                    procMeshes_[proci].time().timeName(),
                    procMeshes_[proci],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                procMeshes_[proci]
            )
        );
    }

    return reconstructFvVolumeField
    (
        IOobject
        (
            fieldIoObject.name(),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        procFields
    );
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fvFieldReconstructor::reconstructFvSurfaceField
(
    const IOobject& fieldIoObject,
    const PtrList<GeometricField<Type, fvsPatchField, surfaceMesh>>& procFields
) const
{
    // Create the internalField
    Field<Type> internalField(mesh_.nInternalFaces());

    // Create the patch fields
    PtrList<fvsPatchField<Type>> patchFields(mesh_.boundary().size());


    forAll(procMeshes_, proci)
    {
        const GeometricField<Type, fvsPatchField, surfaceMesh>& procField =
            procFields[proci];

        // Set the face values in the reconstructed field

        // It is necessary to create a copy of the addressing array to
        // take care of the face direction offset trick.
        //
        {
            const labelList& faceMap = faceProcAddressing_[proci];

            // Correctly oriented copy of internal field
            Field<Type> procInternalField(procField.primitiveField());

            // Addressing into original field
            labelList curAddr(procInternalField.size());

            forAll(procInternalField, addrI)
            {
                curAddr[addrI] = mag(faceMap[addrI])-1;
                if (faceMap[addrI] < 0)
                {
                    procInternalField[addrI] = -procInternalField[addrI];
                }
            }

            // Map
            internalField.rmap(procInternalField, curAddr);
        }

        // Set the boundary patch values in the reconstructed field
        forAll(boundaryProcAddressing_[proci], patchi)
        {
            // Get patch index of the original patch
            const label curBPatch = boundaryProcAddressing_[proci][patchi];

            // Get addressing slice for this patch
            const labelList::subList cp =
                procMeshes_[proci].boundary()[patchi].patchSlice
                (
                    faceProcAddressing_[proci]
                );

            // check if the boundary patch is not a processor patch
            if (curBPatch >= 0)
            {
                // Regular patch. Fast looping

                if (!patchFields(curBPatch))
                {
                    patchFields.set
                    (
                        curBPatch,
                        fvsPatchField<Type>::New
                        (
                            procField.boundaryField()[patchi],
                            mesh_.boundary()[curBPatch],
                            DimensionedField<Type, surfaceMesh>::null(),
                            fvPatchFieldReconstructor
                            (
                                mesh_.boundary()[curBPatch].size()
                            )
                        )
                    );
                }

                const label curPatchStart =
                    mesh_.boundaryMesh()[curBPatch].start();

                labelList reverseAddressing(cp.size());

                forAll(cp, facei)
                {
                    // Subtract one to take into account offsets for
                    // face direction.
                    reverseAddressing[facei] = cp[facei] - 1 - curPatchStart;
                }

                patchFields[curBPatch].rmap
                (
                    procField.boundaryField()[patchi],
                    reverseAddressing
                );
            }
            else
            {
                const Field<Type>& curProcPatch =
                    procField.boundaryField()[patchi];

                // In processor patches, there's a mix of internal faces (some
                // of them turned) and possible cyclics. Slow loop
                forAll(cp, facei)
                {
                    label curF = cp[facei] - 1;

                    // Is the face turned the right side round
                    if (curF >= 0)
                    {
                        // Is the face on the boundary?
                        if (curF >= mesh_.nInternalFaces())
                        {
                            label curBPatch =
                                mesh_.boundaryMesh().whichPatch(curF);

                            if (!patchFields(curBPatch))
                            {
                                patchFields.set
                                (
                                    curBPatch,
                                    fvsPatchField<Type>::New
                                    (
                                        mesh_.boundary()[curBPatch].type(),
                                        mesh_.boundary()[curBPatch],
                                        DimensionedField<Type, surfaceMesh>
                                           ::null()
                                    )
                                );
                            }

                            // add the face
                            label curPatchFace =
                                mesh_.boundaryMesh()
                                [curBPatch].whichFace(curF);

                            patchFields[curBPatch][curPatchFace] =
                                curProcPatch[facei];
                        }
                        else
                        {
                            // Internal face
                            internalField[curF] = curProcPatch[facei];
                        }
                    }
                }
            }
        }
    }

    forAll(mesh_.boundary(), patchi)
    {
        // add empty patches
        if
        (
            isType<emptyFvPatch>(mesh_.boundary()[patchi])
         && !patchFields(patchi)
        )
        {
            patchFields.set
            (
                patchi,
                fvsPatchField<Type>::New
                (
                    emptyFvsPatchField<Type>::typeName,
                    mesh_.boundary()[patchi],
                    DimensionedField<Type, surfaceMesh>::null()
                )
            );
        }
    }


    // Now construct and write the field
    // setting the internalField and patchFields
    return tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            fieldIoObject,
            mesh_,
            procFields[0].dimensions(),
            internalField,
            patchFields
        )
    );
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fvFieldReconstructor::reconstructFvSurfaceField
(
    const IOobject& fieldIoObject
) const
{
    // Read the field for all the processors
    PtrList<GeometricField<Type, fvsPatchField, surfaceMesh>> procFields
    (
        procMeshes_.size()
    );

    forAll(procMeshes_, proci)
    {
        procFields.set
        (
            proci,
            new GeometricField<Type, fvsPatchField, surfaceMesh>
            (
                IOobject
                (
                    fieldIoObject.name(),
                    procMeshes_[proci].time().timeName(),
                    procMeshes_[proci],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                procMeshes_[proci]
            )
        );
    }

    return reconstructFvSurfaceField
    (
        IOobject
        (
            fieldIoObject.name(),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        procFields
    );
}


template<class Type>
void Foam::fvFieldReconstructor::reconstructFvVolumeInternalFields
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
)
{
    const word& fieldClassName = DimensionedField<Type, volMesh>::typeName;

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< "    Reconstructing " << fieldClassName << "s\n" << endl;

        forAllConstIter(IOobjectList, fields, fieldIter)
        {
            if
            (
                selectedFields.empty()
             || selectedFields.found(fieldIter()->name())
            )
            {
                Info<< "        " << fieldIter()->name() << endl;

                reconstructFvVolumeInternalField<Type>(*fieldIter())().write();

                nReconstructed_++;
            }
        }
        Info<< endl;
    }
}


template<class Type>
void Foam::fvFieldReconstructor::reconstructFvVolumeFields
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
)
{
    const word& fieldClassName =
        GeometricField<Type, fvPatchField, volMesh>::typeName;

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< "    Reconstructing " << fieldClassName << "s\n" << endl;

        forAllConstIter(IOobjectList, fields, fieldIter)
        {
            if
            (
                selectedFields.empty()
             || selectedFields.found(fieldIter()->name())
            )
            {
                Info<< "        " << fieldIter()->name() << endl;

                reconstructFvVolumeField<Type>(*fieldIter())().write();

                nReconstructed_++;
            }
        }
        Info<< endl;
    }
}


template<class Type>
void Foam::fvFieldReconstructor::reconstructFvSurfaceFields
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
)
{
    const word& fieldClassName =
        GeometricField<Type, fvsPatchField, surfaceMesh>::typeName;

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< "    Reconstructing " << fieldClassName << "s\n" << endl;

        forAllConstIter(IOobjectList, fields, fieldIter)
        {
            if
            (
                selectedFields.empty()
             || selectedFields.found(fieldIter()->name())
            )
            {
                Info<< "        " << fieldIter()->name() << endl;

                reconstructFvSurfaceField<Type>(*fieldIter())().write();

                nReconstructed_++;
            }
        }
        Info<< endl;
    }
}


// ************************************************************************* //
