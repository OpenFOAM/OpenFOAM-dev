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

#include "phaseSystem.H"
#include "sidedPhaseInterface.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::dictionary Foam::phaseSystem::interfacialDict(const word& name) const
{
    bool found = false;
    dictionary dict(name);

    // If it is a dictionary then merge it in
    if (this->isDict(name))
    {
        found = true;
        dict.merge(this->subDict(name));
    }

    // Function to add old-format list/table entries
    auto add = [&](const word& sidePhaseName)
    {
        const word nameSidePhaseName =
            IOobject::groupName(name, sidePhaseName);

        if (!this->found(nameSidePhaseName)) return;

        found = true;

        List<wordListAndType<Type>> wlats(this->lookup(nameSidePhaseName));

        forAll(wlats, i)
        {
            word keyword =
                phaseInterface::oldNamePartsToName(*this, wlats[i].wl);

            if (sidePhaseName != word::null)
            {
                keyword.append
                (
                    "_"
                  + sidedPhaseInterface::separator()
                  + "_"
                  + sidePhaseName
                );
            }

            dict.add(keyword, wlats[i].t);
        }
    };

    // If a dictionary entry was not found then try a list/table entry
    if (!found)
    {
        add(word::null);
    }

    // Add sided models
    forAll(phases(), phasei)
    {
        add(phases()[phasei].name());
    }

    // Barf if nothing was found
    if (!found)
    {
        return this->subDict(name);
    }

    return dict;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::phaseSystem::fillFields
(
    const word& name,
    const dimensionSet& dims,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& fieldList
) const
{
    forAll(this->phaseModels_, phasei)
    {
        if (fieldList.set(phasei))
        {
            continue;
        }

        const phaseModel& phase = this->phaseModels_[phasei];

        fieldList.set
        (
            phasei,
            new GeometricField<Type, PatchField, GeoMesh>
            (
                IOobject
                (
                    IOobject::groupName(name, phase.name()),
                    this->mesh_.time().name(),
                    this->mesh_
                ),
                this->mesh_,
                dimensioned<Type>("zero", dims, pTraits<Type>::zero)
            )
        );
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::phaseSystem::fillFields
(
    const word& name,
    const dimensionSet& dims,
    HashPtrTable<GeometricField<Type, PatchField, GeoMesh>>& fieldTable
) const
{
    forAll(this->phaseModels_, phasei)
    {
        const phaseModel& phase = this->phaseModels_[phasei];

        if (fieldTable.set(phase.name()))
        {
            continue;
        }

        fieldTable.set
        (
            phase.name(),
            new GeometricField<Type, PatchField, GeoMesh>
            (
                IOobject
                (
                    IOobject::groupName(name, phase.name()),
                    this->mesh_.time().name(),
                    this->mesh_
                ),
                this->mesh_,
                dimensioned<Type>("zero", dims, pTraits<Type>::zero)
            )
        );
    }
}


template<class ModelType>
Foam::word Foam::phaseSystem::modelName() const
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


template<class ModelType, class ... InterfaceTypes>
void Foam::phaseSystem::generateInterfacialModels
(
    const dictionary& dict,
    const phaseInterface& interface,
    PtrList<phaseInterface>& interfaces,
    PtrList<ModelType>& models
) const
{
    // Construct sub-dictionaries and associated interfaces
    hashedWordList names;
    PtrList<dictionary> dicts;
    forAllConstIter(dictionary, dict, iter)
    {
        // Get the model sub dictionary and its associated interface
        const dictionary& modelDict = iter().dict();
        autoPtr<phaseInterface> modelInterfacePtr =
            phaseInterface::New(*this, iter().keyword());

        // Cast the interface down to the first specified type possible
        autoPtr<phaseInterface> interfacePtr;
        List<bool>
        ({
            interfacePtr.empty()
         && isA<InterfaceTypes>(modelInterfacePtr())
         && (
                interfacePtr.set
                (
                    new InterfaceTypes
                    (
                        refCast<InterfaceTypes>(modelInterfacePtr())
                    )
                ),
                true
            )
            ...
        });
        if (!interfacePtr.valid())
        {
            FatalErrorInFunction
                << "Interface " << modelInterfacePtr->name()
                << " is not of suitable type for construction of a "
                << ModelType::typeName
                << exit(FatalError);
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
            dicts.append(new dictionary(name));
            interfaces.append(interfacePtr.ptr());
            models.append(nullptr);
        }

        // Add the model dictionary
        dicts[names[name]].add
        (
            modelInterfacePtr->name(),
            modelDict
        );
    }

    // Construct the models
    forAll(interfaces, i)
    {
        models.set(i, ModelType::New(dicts[i], interfaces[i]));
    }
}


template<class ModelType>
void Foam::phaseSystem::generateInterfacialModels
(
    const dictionary& dict,
    HashTable
    <
        autoPtr<ModelType>,
        phaseInterfaceKey,
        phaseInterfaceKey::hash
    >& models
) const
{
    // Construct lists of interfaces and models
    PtrList<phaseInterface> listInterfaces;
    PtrList<ModelType> listModels;
    generateInterfacialModels<ModelType, phaseInterface>
    (
        dict,
        NullObjectRef<phaseInterface>(),
        listInterfaces,
        listModels
    );

    // Transfer to a keyed table
    forAll(listInterfaces, i)
    {
        models.insert(listInterfaces[i], listModels.set(i, nullptr));
    }
}


template<class ModelType>
void Foam::phaseSystem::generateInterfacialModels
(
    HashTable
    <
        autoPtr<ModelType>,
        phaseInterfaceKey,
        phaseInterfaceKey::hash
    >& models
) const
{
    generateInterfacialModels
    (
        interfacialDict<dictionary>(modelName<ModelType>()),
        models
    );
}


template<class ValueType>
void Foam::phaseSystem::generateInterfacialValues
(
    const dictionary& dict,
    HashTable<ValueType, phaseInterfaceKey, phaseInterfaceKey::hash>& values
) const
{
    forAllConstIter(dictionary, dict, iter)
    {
        autoPtr<phaseInterface> interfacePtr =
            phaseInterface::New(*this, iter().keyword());

        const ValueType value(pTraits<ValueType>(iter().stream()));

        values.insert(interfacePtr(), value);
    }
}


template<class ValueType>
void Foam::phaseSystem::generateInterfacialValues
(
    const word& valueName,
    HashTable<ValueType, phaseInterfaceKey, phaseInterfaceKey::hash>& values
) const
{
    generateInterfacialValues(interfacialDict<ValueType>(valueName), values);
}


template<class ModelType>
const Foam::dictionary& Foam::phaseSystem::modelSubDict
(
    const dictionary& dict
)
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


template<class ModelType>
void Foam::phaseSystem::validateMassTransfer
(
    const phaseInterface& interface
) const
{
    if (interface.phase1().stationary() || interface.phase2().stationary())
    {
        FatalErrorInFunction
            << "A " << ModelType::typeName << " was specified for pair "
            << interface.name() << ", but one of these phases is stationary. "
            << "Mass transfer is not supported on stationary phases"
            << exit(FatalError);
    }
}


template<class ModelType>
bool Foam::phaseSystem::foundInterfacialModel
(
    const phaseInterface& interface
) const
{
    return
        mesh().foundObject<ModelType>
        (
            IOobject::groupName(ModelType::typeName, interface.name())
        );
}


template<class ModelType>
const ModelType& Foam::phaseSystem::lookupInterfacialModel
(
    const phaseInterface& interface
) const
{
    return
        mesh().lookupObject<ModelType>
        (
            IOobject::groupName(ModelType::typeName, interface.name())
        );
}


// ************************************************************************* //
