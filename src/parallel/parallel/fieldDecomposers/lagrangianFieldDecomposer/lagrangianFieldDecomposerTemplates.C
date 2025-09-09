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

\*---------------------------------------------------------------------------*/

#include "lagrangianFieldDecomposer.H"
#include "cloud.H"
#include "IOobjectList.H"
#include "CompactIOField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template
<
    class Type,
    template<class> class ReadIOContainer,
    template<class> class IOContainer
>
void Foam::lagrangianFieldDecomposer::decomposeFields
(
    const IOobjectList& objects
) const
{
    const word& fieldClassName = ReadIOContainer<Type>::typeName;

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< nl << "    Decomposing " << fieldClassName << "s" << nl << endl;

        forAllConstIter(IOobjectList, fields, fieldIter)
        {
            Info<< "        " << fieldIter()->name() << endl;

            const PtrList<IOContainer<Type>> procFields =
                decomposeField<Type, ReadIOContainer, IOContainer>
                (
                    *fieldIter()
                );

            forAll(procFields, proci)
            {
                procFields[proci].write();
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    class Type,
    template<class> class ReadIOContainer,
    template<class> class IOContainer
>
Foam::PtrList<IOContainer<Type>>
Foam::lagrangianFieldDecomposer::decomposeField
(
    const IOobject& fieldIoObject
) const
{
    // Read the complete field
    const ReadIOContainer<Type> field(fieldIoObject);

    // Construct the processor fields
    PtrList<IOContainer<Type>> procFields(procMeshes_.size());
    forAll(procMeshes_, proci)
    {
        procFields.set
        (
            proci,
            new IOContainer<Type>
            (
                IOobject
                (
                    field.name(),
                    procMeshes_[proci].time().name(),
                    lagrangian::cloud::prefix/cloudName_,
                    procMeshes_[proci],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                Field<Type>(field, particleProcAddressing_[proci])
            )
        );
    }

    return procFields;
}


template<class Type>
void Foam::lagrangianFieldDecomposer::decomposeFields
(
    const IOobjectList& objects
) const
{
    decomposeFields<Type, IOField, IOField>(objects);
    decomposeFields<Field<Type>, IOField, CompactIOField>(objects);
    decomposeFields<Field<Type>, CompactIOField, CompactIOField>(objects);
}


// ************************************************************************* //
