/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
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

#include "fvMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
typename Foam::pTraits<Type>::labelType Foam::fvMesh::validComponents() const
{
    return pow
    (
        this->solutionD(),
        pTraits
        <
            typename powProduct<Vector<label>,
            pTraits<Type>::rank>::type
        >::zero
    );
}


template<class GeoField>
Foam::UPtrList<GeoField> Foam::fvMesh::fields(const bool strict) const
{
    HashTable<GeoField*> fields
    (
        const_cast<fvMesh&>(*this).lookupClass<GeoField>(strict)
    );
    UPtrList<GeoField> curFields(fields.size());

    label i = 0;
    forAllIter(typename HashTable<GeoField*>, fields, iter)
    {
        if (!geometryFields.found(iter()->name()))
        {
            curFields.set(i++, iter());
        }
    }
    curFields.setSize(i);

    return curFields;
}


template<class GeoField>
Foam::UPtrList<GeoField> Foam::fvMesh::curFields() const
{
    HashTable<GeoField*> fields
    (
        const_cast<fvMesh&>(*this).lookupClass<GeoField>()
    );
    UPtrList<GeoField> curFields(fields.size());

    label i = 0;
    forAllIter(typename HashTable<GeoField*>, fields, iter)
    {
        if (!geometryFields.found(iter()->name()) && !iter()->isOldTime())
        {
            curFields.set(i++, iter());
        }
    }
    curFields.setSize(i);

    return curFields;
}


template<class Type, template<class> class GeoField>
void Foam::fvMesh::storeOldTimeFields()
{
    UPtrList<GeoField<Type>> curFields(this->curFields<GeoField<Type>>());

    forAll(curFields, i)
    {
        curFields[i].storeOldTimes();
    }
}


template<template<class> class GeoField>
void Foam::fvMesh::storeOldTimeFields()
{
    #define StoreOldTimeFields(Type, nullArg) \
        storeOldTimeFields<Type, GeoField>();
    FOR_ALL_FIELD_TYPES(StoreOldTimeFields);
    #undef StoreOldTimeFields
}


template<class Type, template<class> class GeoField>
void Foam::fvMesh::nullOldestTimeFields()
{
    UPtrList<GeoField<Type>> curFields(this->curFields<GeoField<Type>>());

    forAll(curFields, i)
    {
        curFields[i].nullOldestTime();
    }
}


template<template<class> class GeoField>
void Foam::fvMesh::nullOldestTimeFields()
{
    #define nullOldestTimeFields(Type, nullArg) \
        nullOldestTimeFields<Type, GeoField>();
    FOR_ALL_FIELD_TYPES(nullOldestTimeFields);
    #undef nullOldestTimeFields
}


// ************************************************************************* //
