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

Namespace
    Foam

Description
    Functions for generating tables of interfacial models

\*---------------------------------------------------------------------------*/

#ifndef generateInterfacialModels_H
#define generateInterfacialModels_H

#include "phaseSystem.H"
#include "TypeSet.H"
#include "dispersedPhaseInterface.H"
#include "segregatedPhaseInterface.H"
#include "displacedPhaseInterface.H"
#include "sidedPhaseInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class ModelType>
word modelName()
{
    word name = ModelType::typeName;

    // Extract the innermost part of the template
    const word::size_type i0 = name.find_last_of('<');
    if (i0 != word::npos)
    {
        const word::size_type i1 = name.find_first_of('>', i0 + 1);
        if (i1 != word::npos)
        {
            name = name(i0 + 1, i1 - i0 - 1);
        }
    }

    // Strip "Model" off the end of the name
    if (name(name.size() - 5, 5) == "Model")
    {
        name = name(name.size() - 5);
    }

    return name;
}


template<class ModelType>
const dictionary& modelSubDict(const dictionary& dict)
{
    if (dict.size() != 1)
    {
        FatalErrorInFunction
            << "Too many matching entries for construction of a "
            << ModelType::typeName << nl << dict.toc()
            << exit(FatalError);
    }

    if (!dict.first()->isDict())
    {
        FatalErrorInFunction
            << "Non-sub-dictionary entries found for specification of a "
            << ModelType::typeName
            << exit(FatalError);
    }

    return dict.first()->dict();
}


template<class ModelType, typename = void>
struct ModelPhaseInterfaceTypes
:
    public std::false_type
{
    typedef TypeSet<> required;

    typedef TypeSet<> allowed;

    typedef TypeSet
    <
        dispersedPhaseInterface,
        segregatedPhaseInterface,
        displacedPhaseInterface,
        sidedPhaseInterface
    > prohibited;
};

template<class ModelType>
struct ModelPhaseInterfaceTypes<ModelType, VoidT<typename ModelType::modelType>>
:
    public std::true_type
{
    typedef typename
        TypeSetConcatenate
        <
            typename ModelPhaseInterfaceTypes
            <
                typename ModelType::modelType
            >::required,
            typename ModelType::requiredPhaseInterfaces
        >::type required;

    typedef typename
        TypeSetConcatenate
        <
            typename ModelPhaseInterfaceTypes
            <
                typename ModelType::modelType
            >::allowed,
            typename ModelType::allowedPhaseInterfaces
        >::type allowed;

    typedef typename
        TypeSetRemove
        <
            typename ModelPhaseInterfaceTypes<void>::prohibited,
            allowed
        >::type prohibited;
};


template<class ModelType>
bool isModelPhaseInterfaceType(const phaseInterface& interface)
{
    return
        ModelPhaseInterfaceTypes<ModelType>::required::isAll(interface)
     && !ModelPhaseInterfaceTypes<ModelType>::prohibited::isAny(interface);
}


template<class InterfaceTypeSet>
struct PhaseInterfaceInfo;

template<class InterfaceType, class ... InterfaceTypes>
struct PhaseInterfaceInfo<TypeSet<InterfaceType, InterfaceTypes ...>>
{
    string operator()(const char* andOr) const
    {
        const string subS =
            PhaseInterfaceInfo<TypeSet<InterfaceTypes ...>>()(andOr);

        string type(InterfaceType::typeName_());
        type.removeTrailing("PhaseInterface");

        const string s =
            type + " (i.e., contains '_" + InterfaceType::separator() + "_')";

        return
            sizeof ... (InterfaceTypes) == 0
          ? s.c_str()
          : sizeof ... (InterfaceTypes) == 1
          ? s + ' ' + andOr + ' ' + subS
          : s + ", " + subS;
    }
};

template<>
struct PhaseInterfaceInfo<TypeSet<>>
{
    inline string operator()(const char* andOr) const
    {
        return string();
    }
};


template<class ModelType>
void checkInterfacialModelsDict
(
    const phaseSystem& fluid,
    const dictionary& dict,
    const wordHashSet& ignoreKeys = wordHashSet()
)
{
    forAllConstIter(dictionary, dict, iter)
    {
        // Get the keyword name and skip if it is ignored
        const word& modelName = iter().keyword();
        if (ignoreKeys.found(modelName)) continue;

        // Construct the associated interface
        autoPtr<phaseInterface> modelInterfacePtr =
            phaseInterface::New(fluid, modelName);

        // Check the interface against the model's permitted interface types
        if (!isModelPhaseInterfaceType<ModelType>(modelInterfacePtr()))
        {
            FatalIOErrorInFunction(dict)
                << "The interface " << modelName
                << " is not of suitable type for construction of a "
                << ModelType::typeName << '.';

            const string requiredInfo =
                PhaseInterfaceInfo
                <
                    typename ModelPhaseInterfaceTypes<ModelType>::required
                >()("and");

            if (requiredInfo.size())
            {
                FatalIOError
                    << " The interface must be of type "
                    << requiredInfo.c_str() << '.';
            }

            const string prohibitedInfo =
                PhaseInterfaceInfo
                <
                    typename ModelPhaseInterfaceTypes<ModelType>::prohibited
                >()("or");

            if (prohibitedInfo.size())
            {
                FatalIOError
                    << " The interface must NOT be of type "
                    << prohibitedInfo.c_str() << '.';
            }

            FatalIOError
                << exit(FatalIOError);
        }
    }
}


template<class ModelType, class ... InterfaceTypes, class ... Args>
void generateInterfacialModels
(
    PtrList<phaseInterface>& interfaces,
    PtrList<ModelType>& models,
    const phaseSystem& fluid,
    const dictionary& dict,
    const wordHashSet& ignoreKeys,
    const phaseInterface& interface,
    const Args& ... args
)
{
    // Construct sub-dictionaries and associated interfaces
    hashedWordList names;
    PtrList<dictionary> dicts;
    forAllConstIter(dictionary, dict, iter)
    {
        // Get the keyword name and skip if it is ignored
        const word& modelName = iter().keyword();
        if (ignoreKeys.found(modelName)) continue;

        // Get the model sub dictionary and construct the associated interface
        const dictionary& modelDict = iter().dict();
        autoPtr<phaseInterface> modelInterfacePtr =
            phaseInterface::New(fluid, modelName);

        // Check the interface against the model's permitted interface types
        if
        (
            isNull(interface)
         && !isModelPhaseInterfaceType<ModelType>(modelInterfacePtr())
        ) continue;

        // Convert the interface to the first specified type possible
        autoPtr<phaseInterface> interfacePtr =
            TypeSet<InterfaceTypes ...>::clone(modelInterfacePtr());
        if (!interfacePtr.valid())
        {
            FatalIOErrorInFunction(dict)
                << "The interface " << modelName
                << " is not of suitable type for construction of a "
                << ModelType::typeName
                << exit(FatalIOError);
        }

        // If constructing for a specific interface then combine with this
        // interface. This ensures interface information propagates through
        // hierarchical model generation.
        if (notNull(interface))
        {
            interfacePtr = phaseInterface::New(interface, interfacePtr());
        }

        // Find an existing dictionary to add to or create a new one
        const word name = interfacePtr->name();
        if (!names.found(name))
        {
            names.append(name);
            dicts.append(new dictionary(dict.name()));
            interfaces.append(interfacePtr.ptr());
            models.append(nullptr);
        }

        // Add the model dictionary
        dicts[names[name]].add(modelName, modelDict);
    }

    // Construct the models
    forAll(interfaces, i)
    {
        models.set(i, ModelType::New(dicts[i], interfaces[i], args ...));
    }
}


template<class ModelType, class ... Args>
Foam::HashPtrTable
<
    ModelType,
    Foam::phaseInterfaceKey,
    Foam::phaseInterfaceKey::hash
> generateInterfacialModels
(
    const phaseSystem& fluid,
    const dictionary& dict,
    const wordHashSet& ignoreKeys,
    const bool ignoreNonModelPhaseInterfaceTypes,
    const Args& ... args
)
{
    // Check that the dictionary has the correct interface types in it
    if (!ignoreNonModelPhaseInterfaceTypes)
    {
        checkInterfacialModelsDict<ModelType>(fluid, dict, ignoreKeys);
    }

    // Construct lists of interfaces and models
    PtrList<phaseInterface> listInterfaces;
    PtrList<ModelType> listModels;
    generateInterfacialModels<ModelType, phaseInterface>
    (
        listInterfaces,
        listModels,
        fluid,
        dict,
        ignoreKeys,
        NullObjectRef<phaseInterface>(),
        args ...
    );

    // Transfer to a keyed table
    HashPtrTable<ModelType, phaseInterfaceKey, phaseInterfaceKey::hash> models;
    forAll(listInterfaces, i)
    {
        models.insert(listInterfaces[i], listModels.set(i, nullptr).ptr());
    }

    return models;
}


template<class ModelType>
Foam::HashPtrTable
<
    ModelType,
    Foam::phaseInterfaceKey,
    Foam::phaseInterfaceKey::hash
> generateInterfacialModels
(
    const phaseSystem& fluid,
    const dictionary& dict
)
{
    return
        generateInterfacialModels<ModelType>
        (
            fluid,
            dict,
            wordHashSet(),
            false
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
