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

#include "pointFieldDecomposer.H"
#include "fvMesh.H"
#include "processorPointPatchFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::PtrList<Foam::PointField<Type>>
Foam::pointFieldDecomposer::decomposeField
(
    const IOobject& fieldIoObject
) const
{
    // Read the field
    const PointField<Type> field
    (
        IOobject
        (
            fieldIoObject.name(),
            completeMesh_.db().time().name(),
            completeMesh_.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ),
        completeMesh_
    );

    // Construct the processor fields
    PtrList<PointField<Type>> procFields(procMeshes_.size());
    forAll(procMeshes_, proci)
    {
        const pointMesh& procMesh = pointMesh::New(procMeshes_[proci]);

        // Create and map the internal field values
        Field<Type> internalField
        (
            field.primitiveField(),
            pointProcAddressing_[proci]
        );

        // Create a list of pointers for the patchFields
        PtrList<pointPatchField<Type>> patchFields
        (
            procMesh.boundary().size()
        );

        // Create and map the patch field values
        forAll(procMesh.boundary(), patchi)
        {
            if (patchi < completeMesh_.boundary().size())
            {
                patchFields.set
                (
                    patchi,
                    pointPatchField<Type>::New
                    (
                        field.boundaryField()[patchi],
                        procMesh.boundary()[patchi],
                        DimensionedField<Type, pointMesh>::null(),
                        patchFieldDecomposers_[proci][patchi]
                    )
                );
            }
            else
            {
                patchFields.set
                (
                    patchi,
                    new processorPointPatchField<Type>
                    (
                        procMesh.boundary()[patchi],
                        DimensionedField<Type, pointMesh>::null()
                    )
                );
            }
        }

        // Create the field for the processor
        procFields.set
        (
            proci,
            new PointField<Type>
            (
                IOobject
                (
                    field.name(),
                    procMesh().time().name(),
                    procMesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                procMesh,
                field.dimensions(),
                internalField,
                patchFields
            )
        );
    }

    return procFields;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::pointFieldDecomposer::decomposeFields
(
    const IOobjectList& objects
)
{
    const word& fieldClassName = PointField<Type>::typeName;

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< nl << "    Decomposing " << fieldClassName << "s" << nl << endl;

        forAllConstIter(IOobjectList, fields, fieldIter)
        {
            Info<< "        " << fieldIter()->name() << endl;

            PtrList<PointField<Type>> procFields =
                decomposeField<Type>(*fieldIter());

            forAll(procFields, proci)
            {
                procFields[proci].write();
            }
        }
    }
}


// ************************************************************************* //
