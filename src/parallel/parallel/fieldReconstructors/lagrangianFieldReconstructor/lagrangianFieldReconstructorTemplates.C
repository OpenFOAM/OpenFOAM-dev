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

#include "lagrangianFieldReconstructor.H"
#include "IOobjectList.H"
#include "CompactIOField.H"
#include "cloud.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FieldType>
bool Foam::lagrangianFieldReconstructor::reconstructs
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
)
{
    IOobjectList fields = objects.lookupClass(FieldType::typeName);

    if (fields.size() && selectedFields.empty())
    {
        return true;
    }

    forAllConstIter(IOobjectList, fields, fieldIter)
    {
        if (selectedFields.found(fieldIter()->name()))
        {
            return true;
        }
    }

    return false;
}


template
<
    class Type,
    template<class> class ReadIOContainer,
    template<class> class IOContainer
>
void Foam::lagrangianFieldReconstructor::reconstructFields
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
) const
{
    const word& fieldClassName = ReadIOContainer<Type>::typeName;

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< nl << "    Reconstructing lagrangian "
            << fieldClassName << "s" << nl << endl;

        forAllConstIter(IOobjectList, fields, fieldIter)
        {
            if
            (
                selectedFields.empty()
             || selectedFields.found(fieldIter()->name())
            )
            {
                Info<< "        " << fieldIter()->name() << endl;

                reconstructField<Type, ReadIOContainer, IOContainer>
                (
                    *fieldIter()
                )().write();
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
Foam::tmp<IOContainer<Type>>
Foam::lagrangianFieldReconstructor::reconstructField
(
    const IOobject& fieldIoObject
) const
{
    // Construct the complete field
    tmp<IOContainer<Type>> tfield
    (
        new IOContainer<Type>
        (
            IOobject
            (
                fieldIoObject.name(),
                completeMesh_.time().name(),
                lagrangian::cloud::prefix/cloudName_,
                completeMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Field<Type>(0)
        )
    );
    Field<Type>& field = tfield.ref();

    // Combine the processor fields into the complete field
    forAll(procMeshes_, proci)
    {
        typeIOobject<ReadIOContainer<Type>> localIOobject
        (
            fieldIoObject.name(),
            procMeshes_[proci].time().name(),
            lagrangian::cloud::prefix/cloudName_,
            procMeshes_[proci],
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (localIOobject.headerOk())
        {
            IOContainer<Type> fieldi(localIOobject);

            label offset = field.size();
            field.setSize(offset + fieldi.size());

            forAll(fieldi, j)
            {
                field[offset + j] = fieldi[j];
            }
        }
    }

    return tfield;
}


template<class Type>
void Foam::lagrangianFieldReconstructor::reconstructFields
(
    const IOobjectList& o,
    const HashSet<word>& sf
) const
{
    reconstructFields<Type, IOField, IOField>(o, sf);
    reconstructFields<Field<Type>, IOField, CompactIOField>(o, sf);
    reconstructFields<Field<Type>, CompactIOField, CompactIOField>(o, sf);
}


// ************************************************************************* //
