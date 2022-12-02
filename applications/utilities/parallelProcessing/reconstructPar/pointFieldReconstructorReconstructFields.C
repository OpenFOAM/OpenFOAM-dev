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
Foam::tmp<Foam::PointField<Type>>
Foam::pointFieldReconstructor::reconstructField(const IOobject& fieldIoObject)
{
    // Read the field for all the processors
    PtrList<PointField<Type>> procFields
    (
        procMeshes_.size()
    );

    forAll(procMeshes_, proci)
    {
        procFields.set
        (
            proci,
            new PointField<Type>
            (
                IOobject
                (
                    fieldIoObject.name(),
                    procMeshes_[proci]().time().name(),
                    procMeshes_[proci](),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                procMeshes_[proci]
            )
        );
    }


    // Create the internalField
    Field<Type> internalField(completeMesh_.size());

    // Create the patch fields
    PtrList<pointPatchField<Type>> patchFields(completeMesh_.boundary().size());


    forAll(procMeshes_, proci)
    {
        const PointField<Type>&
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
                patchi < completeMesh_.boundary().size() ? patchi : -1;

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
                            completeMesh_.boundary()[curBPatch],
                            DimensionedField<Type, pointMesh>::null(),
                            pointPatchFieldReconstructor
                            (
                                completeMesh_.boundary()[curBPatch].size()
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
    return tmp<PointField<Type>>
    (
        new PointField<Type>
        (
            IOobject
            (
                fieldIoObject.name(),
                completeMesh_().time().name(),
                completeMesh_(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            completeMesh_,
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
        PointField<Type>::typeName
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

                nReconstructed_++;
            }
        }

        Info<< endl;
    }
}


// ************************************************************************* //
