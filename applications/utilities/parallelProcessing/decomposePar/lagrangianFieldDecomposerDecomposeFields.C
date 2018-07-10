/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

#include "lagrangianFieldDecomposer.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::lagrangianFieldDecomposer::readFields
(
    const label cloudI,
    const IOobjectList& lagrangianObjects,
    PtrList<PtrList<IOField<Type>>>& lagrangianFields
)
{
    // Search list of objects for lagrangian fields
    IOobjectList lagrangianTypeObjects
    (
        lagrangianObjects.lookupClass(IOField<Type>::typeName)
    );

    lagrangianFields.set
    (
        cloudI,
        new PtrList<IOField<Type>>
        (
            lagrangianTypeObjects.size()
        )
    );

    label lagrangianFieldi = 0;
    forAllIter(IOobjectList, lagrangianTypeObjects, iter)
    {
        lagrangianFields[cloudI].set
        (
            lagrangianFieldi++,
            new IOField<Type>(*iter())
        );
    }
}


template<class Type>
void Foam::lagrangianFieldDecomposer::readFieldFields
(
    const label cloudI,
    const IOobjectList& lagrangianObjects,
    PtrList<PtrList<CompactIOField<Field<Type>, Type>>>& lagrangianFields
)
{
    // Search list of objects for lagrangian fields
    IOobjectList lagrangianTypeObjectsA
    (
        lagrangianObjects.lookupClass(IOField<Field<Type>>::typeName)
    );

    IOobjectList lagrangianTypeObjectsB
    (
        lagrangianObjects.lookupClass
        (
            CompactIOField<Field<Type>,
            Type>::typeName
        )
    );

    lagrangianFields.set
    (
        cloudI,
        new PtrList<CompactIOField<Field<Type>, Type>>
        (
            lagrangianTypeObjectsA.size() + lagrangianTypeObjectsB.size()
        )
    );

    label lagrangianFieldi = 0;

    forAllIter(IOobjectList, lagrangianTypeObjectsA, iter)
    {
        lagrangianFields[cloudI].set
        (
            lagrangianFieldi++,
            new CompactIOField<Field<Type>, Type>(*iter())
        );
    }

    forAllIter(IOobjectList, lagrangianTypeObjectsB, iter)
    {
        lagrangianFields[cloudI].set
        (
            lagrangianFieldi++,
            new CompactIOField<Field<Type>, Type>(*iter())
        );
    }
}


template<class Type>
Foam::tmp<Foam::IOField<Type>>
Foam::lagrangianFieldDecomposer::decomposeField
(
    const word& cloudName,
    const IOField<Type>& field
) const
{
    // Create and map the internal field values
    Field<Type> procField(field, particleIndices_);

    // Create the field for the processor
    return tmp<IOField<Type>>
    (
        new IOField<Type>
        (
            IOobject
            (
                field.name(),
                procMesh_.time().timeName(),
                cloud::prefix/cloudName,
                procMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            procField
        )
    );
}


template<class Type>
Foam::tmp<Foam::CompactIOField<Foam::Field<Type>, Type>>
Foam::lagrangianFieldDecomposer::decomposeFieldField
(
    const word& cloudName,
    const CompactIOField<Field<Type>, Type>& field
) const
{
    // Create and map the internal field values
    Field<Field<Type>> procField(field, particleIndices_);

    // Create the field for the processor
    return tmp<CompactIOField<Field<Type>, Type>>
    (
        new CompactIOField<Field<Type>, Type>
        (
            IOobject
            (
                field.name(),
                procMesh_.time().timeName(),
                cloud::prefix/cloudName,
                procMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            procField
        )
    );
}


template<class GeoField>
void Foam::lagrangianFieldDecomposer::decomposeFields
(
    const word& cloudName,
    const PtrList<GeoField>& fields
) const
{
    // if (particleIndices_.size())
    {
        bool valid = particleIndices_.size() > 0;
        forAll(fields, fieldi)
        {
            decomposeField(cloudName, fields[fieldi])().write(valid);
        }
    }
}


template<class GeoField>
void Foam::lagrangianFieldDecomposer::decomposeFieldFields
(
    const word& cloudName,
    const PtrList<GeoField>& fields
) const
{
    // if (particleIndices_.size())
    {
        bool valid = particleIndices_.size() > 0;
        forAll(fields, fieldi)
        {
            decomposeFieldField(cloudName, fields[fieldi])().write(valid);
        }
    }
}


// ************************************************************************* //
