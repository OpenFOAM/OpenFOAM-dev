/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "pointFieldReconstructor.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh>>
Foam::pointFieldReconstructor::reconstructField(const IOobject& fieldIoObject)
{
    // Read the field for all the processors
    PtrList<GeometricField<Type, pointPatchField, pointMesh>> procFields
    (
        procMeshes_.size()
    );

    forAll(procMeshes_, proci)
    {
        procFields.set
        (
            proci,
            new GeometricField<Type, pointPatchField, pointMesh>
            (
                IOobject
                (
                    fieldIoObject.name(),
                    procMeshes_[proci]().time().timeName(),
                    procMeshes_[proci](),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                procMeshes_[proci]
            )
        );
    }


    // Create the internalField
    Field<Type> internalField(mesh_.size());

    // Create the patch fields
    PtrList<pointPatchField<Type>> patchFields(mesh_.boundary().size());


    forAll(procMeshes_, proci)
    {
        const GeometricField<Type, pointPatchField, pointMesh>&
            procField = procFields[proci];

        // Get processor-to-global addressing for use in rmap
        const labelList& procToGlobalAddr = pointProcAddressing_[proci];

        // Set the cell values in the reconstructed field
        internalField.rmap
        (
            procField.primitiveField(),
            procToGlobalAddr
        );

        // Set the boundary patch values in the reconstructed field
        forAll(procField.boundaryField(), patchi)
        {
            // Get patch index of the original patch
            const label curBPatch =
                patchi < mesh_.boundary().size() ? patchi : -1;

            // check if the boundary patch is not a processor patch
            if (curBPatch != -1)
            {
                if (!patchFields(curBPatch))
                {
                    patchFields.set(
                        curBPatch,
                        pointPatchField<Type>::New
                        (
                            procField.boundaryField()[patchi],
                            mesh_.boundary()[curBPatch],
                            DimensionedField<Type, pointMesh>::null(),
                            pointPatchFieldReconstructor
                            (
                                mesh_.boundary()[curBPatch].size()
                            )
                        )
                    );
                }

                patchFields[curBPatch].rmap
                (
                    procField.boundaryField()[patchi],
                    patchPointAddressing_[proci][patchi]
                );
            }
        }
    }

    // Construct and write the field
    // setting the internalField and patchFields
    return tmp<GeometricField<Type, pointPatchField, pointMesh>>
    (
        new GeometricField<Type, pointPatchField, pointMesh>
        (
            IOobject
            (
                fieldIoObject.name(),
                mesh_().time().timeName(),
                mesh_(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            procFields[0].dimensions(),
            internalField,
            patchFields
        )
    );
}


// Reconstruct and write all point fields
template<class Type>
void Foam::pointFieldReconstructor::reconstructFields
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
)
{
    word fieldClassName
    (
        GeometricField<Type, pointPatchField, pointMesh>::typeName
    );

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< "    Reconstructing " << fieldClassName << "s\n" << endl;

        forAllConstIter(IOobjectList, fields, fieldIter)
        {
            if
            (
                !selectedFields.size()
             || selectedFields.found(fieldIter()->name())
            )
            {
                Info<< "        " << fieldIter()->name() << endl;

                reconstructField<Type>(*fieldIter())().write();
            }
        }

        Info<< endl;
    }
}


// ************************************************************************* //
