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

\*---------------------------------------------------------------------------*/

#include "fvFieldReconstructor.H"
#include "Time.H"
#include "PtrList.H"
#include "fvPatchFields.H"
#include "emptyFvPatch.H"
#include "emptyFvPatchField.H"
#include "emptyFvsPatchField.H"
#include "processorCyclicFvPatch.H"
#include "reverseFvPatchFieldMapper.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvFieldReconstructor::rmapFaceToFace
(
    Field<Type>& toField,
    const Field<Type>& fromField,
    const labelUList& addressing,
    const bool isFlux
)
{
    forAll(addressing, i)
    {
        toField[mag(addressing[i]) - 1] =
            (isFlux && addressing[i] < 0 ? -1 : +1)*fromField[i];
    }
}


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
    Field<Type> internalField(completeMesh_.nCells());

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
            completeMesh_,
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
    PtrList<DimensionedField<Type, volMesh>>
        procFields(procMeshes_.size());

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
                    procMeshes_[proci].time().name(),
                    procMeshes_[proci],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
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
            completeMesh_.time().name(),
            completeMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        procFields
    );
}


template<class Type>
Foam::tmp<Foam::VolField<Type>>
Foam::fvFieldReconstructor::reconstructFvVolumeField
(
    const IOobject& fieldIoObject,
    const PtrList<VolField<Type>>& procFields
) const
{
    // Create the internalField
    Field<Type> internalField(completeMesh_.nCells());

    // Create the patch fields
    PtrList<fvPatchField<Type>> patchFields(completeMesh_.boundary().size());

    forAll(procFields, proci)
    {
        const VolField<Type>& procField =
            procFields[proci];

        // Set the cell values in the reconstructed field
        internalField.rmap
        (
            procField.primitiveField(),
            cellProcAddressing_[proci]
        );

        // Set the boundary patch values in the reconstructed field
        forAll(procField.boundaryField(), procPatchi)
        {
            const fvPatch& procPatch =
                procMeshes_[proci].boundary()[procPatchi];

            const label completePatchi = completePatchID(proci, procPatchi);

            if (completePatchi == procPatchi)
            {
                if (!patchFields(completePatchi))
                {
                    patchFields.set
                    (
                        completePatchi,
                        fvPatchField<Type>::New
                        (
                            procField.boundaryField()[procPatchi],
                            completeMesh_.boundary()[completePatchi],
                            DimensionedField<Type, volMesh>::null(),
                            setSizeFvPatchFieldMapper
                            (
                                completeMesh_.boundary()[completePatchi].size()
                            )
                        )
                    );
                }

                patchFields[completePatchi].map
                (
                    procField.boundaryField()[procPatchi],
                    reverseFvPatchFieldMapper
                    (
                        faceProcAddressingBf_[proci][procPatchi] - 1
                    )
                );
            }
            else if (isA<processorCyclicFvPatch>(procPatch))
            {
                if (!patchFields(completePatchi))
                {
                    patchFields.set
                    (
                        completePatchi,
                        fvPatchField<Type>::New
                        (
                            completeMesh_.boundary()[completePatchi].type(),
                            completeMesh_.boundary()[completePatchi],
                            DimensionedField<Type, volMesh>::null()
                        )
                    );
                }

                if (patchFields[completePatchi].overridesConstraint())
                {
                    OStringStream str;
                    str << "\nThe field \"" << procFields[0].name()
                        << "\" on cyclic patch \""
                        << patchFields[completePatchi].patch().name()
                        << "\" cannot be reconstructed as it is not a cyclic "
                        << "patch field. A \"patchType cyclic;\" setting has "
                        << "been used to override the cyclic patch type.\n\n"
                        << "Cyclic patches like this with non-cyclic boundary "
                        << "conditions should be confined to a single "
                        << "processor using decomposition constraints.";
                    FatalErrorInFunction
                        << stringOps::breakIntoIndentedLines(str.str()).c_str()
                        << exit(FatalError);
                }

                patchFields[completePatchi].map
                (
                    procField.boundaryField()[procPatchi],
                    reverseFvPatchFieldMapper
                    (
                        faceProcAddressingBf_[proci][procPatchi] - 1
                    )
                );
            }
        }
    }

    // Construct and return the field
    return tmp<VolField<Type>>
    (
        new VolField<Type>
        (
            fieldIoObject,
            completeMesh_,
            procFields[0].dimensions(),
            internalField,
            patchFields
        )
    );
}


template<class Type>
Foam::tmp<Foam::VolField<Type>>
Foam::fvFieldReconstructor::reconstructFvVolumeField
(
    const IOobject& fieldIoObject
) const
{
    PtrList<VolField<Type>>
        procFields(procMeshes_.size());

    forAll(procMeshes_, proci)
    {
        procFields.set
        (
            proci,
            new VolField<Type>
            (
                IOobject
                (
                    fieldIoObject.name(),
                    procMeshes_[proci].time().name(),
                    procMeshes_[proci],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
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
            completeMesh_.time().name(),
            completeMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        procFields
    );
}


template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::fvFieldReconstructor::reconstructFvSurfaceField
(
    const IOobject& fieldIoObject,
    const PtrList<SurfaceField<Type>>& procFields
) const
{
    // Create the internalField
    Field<Type> internalField(completeMesh_.nInternalFaces());

    // Create the patch fields
    PtrList<fvsPatchField<Type>> patchFields(completeMesh_.boundary().size());

    forAll(procMeshes_, proci)
    {
        const SurfaceField<Type>& procField =
            procFields[proci];

        // Set the internal face values in the reconstructed field
        rmapFaceToFace
        (
            internalField,
            procField.primitiveField(),
            SubList<label>
            (
                faceProcAddressing_[proci],
                procMeshes_[proci].nInternalFaces()
            ),
            isFlux(procFields[proci])
        );

        // Set the boundary patch values in the reconstructed field
        forAll(procField.boundaryField(), procPatchi)
        {
            const fvPatch& procPatch =
                procMeshes_[proci].boundary()[procPatchi];

            const label completePatchi = completePatchID(proci, procPatchi);

            if (completePatchi == procPatchi)
            {
                if (!patchFields(completePatchi))
                {
                    patchFields.set
                    (
                        completePatchi,
                        fvsPatchField<Type>::New
                        (
                            procField.boundaryField()[procPatchi],
                            completeMesh_.boundary()[completePatchi],
                            DimensionedField<Type, surfaceMesh>::null(),
                            setSizeFvPatchFieldMapper
                            (
                                completeMesh_.boundary()[completePatchi].size()
                            )
                        )
                    );
                }

                patchFields[completePatchi].map
                (
                    procField.boundaryField()[procPatchi],
                    reverseFvPatchFieldMapper
                    (
                        faceProcAddressingBf_[proci][procPatchi] - 1
                    )
                );
            }
            else if (isA<processorCyclicFvPatch>(procPatch))
            {
                if (!patchFields(completePatchi))
                {
                    patchFields.set
                    (
                        completePatchi,
                        fvsPatchField<Type>::New
                        (
                            completeMesh_.boundary()[completePatchi].type(),
                            completeMesh_.boundary()[completePatchi],
                            DimensionedField<Type, surfaceMesh>::null()
                        )
                    );
                }

                patchFields[completePatchi].map
                (
                    procField.boundaryField()[procPatchi],
                    reverseFvPatchFieldMapper
                    (
                        faceProcAddressingBf_[proci][procPatchi] - 1
                    )
                );
            }
            else if (isA<processorFvPatch>(procPatch))
            {
                rmapFaceToFace
                (
                    internalField,
                    procField.boundaryField()[procPatchi],
                    faceProcAddressingBf_[proci][procPatchi],
                    isFlux(procFields[proci])
                );
            }
        }
    }

    // Construct and return the field
    return tmp<SurfaceField<Type>>
    (
        new SurfaceField<Type>
        (
            fieldIoObject,
            completeMesh_,
            procFields[0].dimensions(),
            internalField,
            patchFields
        )
    );
}


template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::fvFieldReconstructor::reconstructFvSurfaceField
(
    const IOobject& fieldIoObject
) const
{
    PtrList<SurfaceField<Type>>
        procFields(procMeshes_.size());

    forAll(procMeshes_, proci)
    {
        procFields.set
        (
            proci,
            new SurfaceField<Type>
            (
                IOobject
                (
                    fieldIoObject.name(),
                    procMeshes_[proci].time().name(),
                    procMeshes_[proci],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
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
            completeMesh_.time().name(),
            completeMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
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
        VolField<Type>::typeName;

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
        SurfaceField<Type>::typeName;

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
