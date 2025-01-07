/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "LagrangianMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, template<class> class GeoField>
void Foam::LagrangianMesh::appendSpecifiedField
(
    const LagrangianSubMesh& appendMesh,
    GeoField<Type>& geoField,
    const Field<Type>& field
) const
{
    geoField.resize(size());

    appendMesh.sub(geoField).ref().primitiveFieldRef() = field;
}


template<class Type, template<class> class GeoField>
bool Foam::LagrangianMesh::appendSpecifiedField
(
    const LagrangianSubMesh& appendMesh,
    const word& fieldName,
    const Field<Type>& field
) const
{
    if (foundObject<GeoField<Type>>(fieldName))
    {
        appendSpecifiedField<Type, GeoField>
        (
            appendMesh,
            lookupObjectRef<GeoField<Type>>(fieldName),
            field
        );

        return true;
    }
    else
    {
        return false;
    }
}


template<class Type, class ... FieldNamesAndFields>
Foam::wordHashSet Foam::LagrangianMesh::appendSpecifiedFields
(
    const LagrangianSubMesh& appendMesh,
    const word& fieldName,
    const Field<Type>& field,
    const FieldNamesAndFields& ... fieldNamesAndFields
) const
{
    const bool found =
        appendSpecifiedField<Type, LagrangianField>
        (
            appendMesh,
            fieldName,
            field
        )
     || appendSpecifiedField<Type, LagrangianDynamicField>
        (
            appendMesh,
            fieldName,
            field
        )
     || appendSpecifiedField<Type, LagrangianInternalField>
        (
            appendMesh,
            fieldName,
            field
        )
     || appendSpecifiedField<Type, LagrangianInternalDynamicField>
        (
            appendMesh,
            fieldName,
            field
        );

    if (!found)
    {
        static const word& TypeName = pTraits<Type>::typeName;

        wordList fieldNames;
        fieldNames.append(toc<LagrangianField<Type>>());
        fieldNames.append(toc<LagrangianDynamicField<Type>>());
        fieldNames.append(toc<LagrangianInternalField<Type>>());
        fieldNames.append(toc<LagrangianInternalDynamicField<Type>>());
        sort(fieldNames);

        FatalErrorInFunction
            << TypeName.capitalise() << " values were specified for "
            << fieldName << " but a " << TypeName << " Lagrangian field of "
            << "that name was not found." << nl << nl
            << "Available " << TypeName << " Lagrangian fields are:" << nl
            << fieldNames << exit(FatalError);
    }

    wordHashSet result
    (
        appendSpecifiedFields(appendMesh, fieldNamesAndFields ...)
    );

    result.insert(fieldName);

    return result;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class GeoField>
Foam::HashTable<const GeoField*>
Foam::LagrangianMesh::lookupCurrentFields(const bool strict) const
{
    HashTable<const GeoField*> fields(lookupClass<GeoField>(strict));

    forAllIter(typename HashTable<GeoField*>, fields, iter)
    {
        if (iter()->isOldTime())
        {
            fields.erase(iter);
        }
    }

    return fields;
}


template<class GeoField>
Foam::HashTable<GeoField*>
Foam::LagrangianMesh::lookupCurrentFields(const bool strict)
{
    HashTable<GeoField*> fields(lookupClass<GeoField>(strict));

    forAllIter(typename HashTable<GeoField*>, fields, iter)
    {
        if (iter()->isOldTime())
        {
            fields.erase(iter);
        }
    }

    return fields;
}


template<class ... FieldNamesAndFields>
Foam::LagrangianSubMesh Foam::LagrangianMesh::inject
(
    const barycentricField& coordinates,
    const labelField& celli,
    const labelField& facei,
    const labelField& faceTrii,
    const FieldNamesAndFields& ... fieldNamesAndFields
)
{
    const LagrangianSubMesh injectionMesh =
        append(coordinates, celli, facei, faceTrii);

    const wordHashSet specifiedFieldNames =
        appendSpecifiedFields(injectionMesh, fieldNamesAndFields ...);

    injectUnspecifiedFields(injectionMesh, specifiedFieldNames);

    return injectionMesh;
}


template<class ... FieldNamesAndFields>
Foam::LagrangianSubMesh Foam::LagrangianMesh::inject
(
    const LagrangianInjection& injection,
    const barycentricField& coordinates,
    const labelField& celli,
    const labelField& facei,
    const labelField& faceTrii,
    const FieldNamesAndFields& ... fieldNamesAndFields
)
{
    const LagrangianSubMesh injectionMesh =
        append(coordinates, celli, facei, faceTrii);

    const wordHashSet specifiedFieldNames =
        appendSpecifiedFields(injectionMesh, fieldNamesAndFields ...);

    injectUnspecifiedFields(injection, injectionMesh, specifiedFieldNames);

    return injectionMesh;
}


template<class ... FieldNamesAndFields>
Foam::LagrangianSubMesh Foam::LagrangianMesh::birth
(
    const labelList& parents,
    const FieldNamesAndFields& ... fieldNamesAndFields
)
{
    const LagrangianSubMesh birthMesh = append(parents);

    const wordHashSet specifiedFieldNames =
        appendSpecifiedFields(birthMesh, fieldNamesAndFields ...);

    birthUnspecifiedFields(parents, birthMesh, specifiedFieldNames);

    return birthMesh;
}


template<class Enumeration>
Foam::labelList Foam::LagrangianMesh::partition
(
    const label nGroups,
    const UList<elementGroup<Enumeration>>& elementsGroups
)
{
    return
        partition
        (
            nGroups,
            reinterpret_cast<const UList<labelPair>&>(elementsGroups)
        );
}


// ************************************************************************* //
