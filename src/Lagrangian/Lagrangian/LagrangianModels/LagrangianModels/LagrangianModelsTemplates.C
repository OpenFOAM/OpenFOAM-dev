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

#include "LagrangianModels.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template
<
    class Type,
    template<class> class PrimitiveField,
    class ... AlphaRhoFieldTypes
>
Foam::tmp<Foam::LagrangianEqn<Type>> Foam::LagrangianModels::sourceTerm
(
    const LagrangianSubField<Type, PrimitiveField>& eqnField,
    const LagrangianSubScalarField& deltaT,
    const AlphaRhoFieldTypes& ... alphaRhoFields
) const
{
    checkApplied();

    tmp<LagrangianEqn<Type>> tEqn
    (
        new LagrangianEqn<Type>
        (
            "S(" + LagrangianModel::fieldsName(alphaRhoFields ...) + ")",
            eqnField
        )
    );
    LagrangianEqn<Type>& eqn = tEqn.ref();

    const PtrListDictionary<LagrangianModel>& modelList(*this);

    const word fieldName = LagrangianModel::fieldName(alphaRhoFields ...);

    forAll(modelList, i)
    {
        const LagrangianModel& model = modelList[i];

        if (model.addsSupToField(fieldName))
        {
            addSupFields_[i].insert(fieldName);

            model.addSup(deltaT, toSubField(alphaRhoFields)() ..., eqn);
        }
    }

    return tEqn;
}


template
<
    class ModelType,
    class FieldSourceType,
    class ... ModelAndFieldSourceTypes
>
struct Foam::LagrangianModels::modelTypeFieldSourceType
<
    ModelType,
    FieldSourceType,
    ModelAndFieldSourceTypes ...
>
{
    static void insert(const LagrangianModel& model, HashTable<word>& result)
    {
        if (isA<ModelType>(model))
        {
            result.insert(model.name(), FieldSourceType::typeName);
        }
        else
        {
            modelTypeFieldSourceType<ModelAndFieldSourceTypes ...>::insert
            (
                model,
                result
            );
        }
    }
};


template<>
struct Foam::LagrangianModels::modelTypeFieldSourceType<>
{
    static void insert(const LagrangianModel& model, HashTable<word>& result)
    {}
};


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type, template<class> class PrimitiveField>
bool Foam::LagrangianModels::addsSupToField
(
    const LagrangianSubField<Type, PrimitiveField>& field
) const
{
    return addsSupToField(LagrangianModel::fieldName(field));
}


template<class ... ModelAndFieldSourceTypes>
Foam::HashTable<Foam::word>
Foam::LagrangianModels::modelTypeFieldSourceTypes() const
{
    const PtrListDictionary<LagrangianModel>& modelList(*this);

    HashTable<word> result;
    forAll(modelList, i)
    {
        modelTypeFieldSourceType<ModelAndFieldSourceTypes ...>::insert
        (
            modelList[i],
            result
        );
    }

    return result;
}


template<class Type, template<class> class PrimitiveField>
Foam::tmp<Foam::LagrangianEqn<Type>> Foam::LagrangianModels::source
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubField<Type, PrimitiveField>& field
) const
{
    return sourceTerm(field, deltaT, field);
}


template
<
    class Type,
    template<class> class PrimitiveField,
    template<class> class PrimitiveEqnField
>
Foam::tmp<Foam::LagrangianEqn<Type>> Foam::LagrangianModels::sourceProxy
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubField<Type, PrimitiveField>& field,
    const LagrangianSubField<Type, PrimitiveEqnField>& eqnField
) const
{
    return sourceTerm(eqnField, deltaT, field);
}


template<class Type, template<class> class PrimitiveField>
Foam::tmp<Foam::LagrangianEqn<Type>> Foam::LagrangianModels::source
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubField<scalar, PrimitiveField>& m,
    const LagrangianSubField<Type, PrimitiveField>& field
) const
{
    return sourceTerm(field, deltaT, m, field);
}


template
<
    class Type,
    template<class> class PrimitiveField,
    template<class> class PrimitiveEqnField
>
Foam::tmp<Foam::LagrangianEqn<Type>> Foam::LagrangianModels::sourceProxy
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubField<scalar, PrimitiveField>& m,
    const LagrangianSubField<Type, PrimitiveField>& field,
    const LagrangianSubField<Type, PrimitiveEqnField>& eqnField
) const
{
    return sourceTerm(eqnField, deltaT, m, field);
}


// ************************************************************************* //
