/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
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

#include "BlendedInterfacialModel.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class modelType>
void Foam::phaseSystem::createSubModels
(
    const dictTable& modelDicts,
    HashTable
    <
        autoPtr<modelType>,
        phasePairKey,
        phasePairKey::hash
    >& models
)
{
    forAllConstIter(dictTable, modelDicts, iter)
    {
        const phasePairKey& key = iter.key();

        models.insert
        (
            key,
            modelType::New
            (
               *iter,
                phasePairs_[key]
            )
        );
    }
}


template<class modelType>
void Foam::phaseSystem::generatePairsAndSubModels
(
    const word& modelName,
    HashTable
    <
        autoPtr<modelType>,
        phasePairKey,
        phasePairKey::hash
    >& models,
    const bool correctFixedFluxBCs
)
{
    dictTable modelDicts(lookup(modelName));

    generatePairs(modelDicts);

    createSubModels(modelDicts, models);
}


template<class modelType>
void Foam::phaseSystem::generatePairsAndSubModels
(
    const word& modelName,
    HashTable
    <
        autoPtr<BlendedInterfacialModel<modelType>>,
        phasePairKey,
        phasePairKey::hash
    >& models,
    const bool correctFixedFluxBCs
)
{
    typedef
        HashTable<autoPtr<modelType>, phasePairKey, phasePairKey::hash>
        modelTypeTable;

    modelTypeTable tempModels;
    generatePairsAndSubModels(modelName, tempModels);

    const blendingMethod& blending
    (
        blendingMethods_.found(modelName)
      ? blendingMethods_[modelName]
      : blendingMethods_.found(member(modelName))
      ? blendingMethods_[member(modelName)]
      : blendingMethods_["default"]
    );

    autoPtr<modelType> noModel(nullptr);

    forAllConstIter(typename modelTypeTable, tempModels, iter)
    {
        if (!iter().valid())
        {
            continue;
        }

        const phasePairKey key(iter.key().first(), iter.key().second());

        if (!phasePairs_.found(key))
        {
            phasePairs_.insert
            (
                key,
                autoPtr<phasePair>
                (
                    new phasePair
                    (
                        phaseModels_[key.first()],
                        phaseModels_[key.second()]
                    )
                )
            );
        }

        const phasePair& pair = phasePairs_[key];
        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();
        const phasePairKey key1In2(phase1.name(), phase2.name(), true);
        const phasePairKey key2In1(phase2.name(), phase1.name(), true);

        models.insert
        (
            key,
            autoPtr<BlendedInterfacialModel<modelType>>
            (
                new BlendedInterfacialModel<modelType>
                (
                    phase1,
                    phase2,
                    blending,
                    tempModels.found(key) ? tempModels[key] : noModel,
                    tempModels.found(key1In2) ? tempModels[key1In2] : noModel,
                    tempModels.found(key2In1) ? tempModels[key2In1] : noModel,
                    correctFixedFluxBCs
                )
            )
        );
    }
}


template<class modelType>
void Foam::phaseSystem::generatePairsAndSubModels
(
    const word& modelName,
    HashTable
    <
        Pair<autoPtr<modelType>>,
        phasePairKey,
        phasePairKey::hash
    >& models,
    const bool correctFixedFluxBCs
)
{
    typedef
        HashTable<autoPtr<modelType>, phasePairKey, phasePairKey::hash>
        modelTypeTable;

    forAll(phaseModels_, phasei)
    {
        const phaseModel& phase = phaseModels_[phasei];

        modelTypeTable tempModels;
        generatePairsAndSubModels
        (
            IOobject::groupName(modelName, phase.name()),
            tempModels,
            correctFixedFluxBCs
        );

        forAllIter(typename modelTypeTable, tempModels, tempModelIter)
        {
            const phasePairKey& key(tempModelIter.key());

            if (!models.found(key))
            {
                models.insert
                (
                    key,
                    Pair<autoPtr<modelType>>()
                );
            }

            const phasePair& pair = phasePairs_[key];

            if (!pair.contains(phase))
            {
                FatalErrorInFunction
                    << "A two-sided " << modelType::typeName << " was "
                    << "specified for the " << phase.name() << " side of the "
                    << pair << " pair, but that phase is not part of that pair."
                    << exit(FatalError);
            }

            models[key][pair.index(phase)] = tempModelIter().ptr();
        }
    }
}


template<class modelType>
void Foam::phaseSystem::validateMassTransfer(const phasePair& key) const
{
    if (key.phase1().stationary() || key.phase2().stationary())
    {
        FatalErrorInFunction
            << "A " << modelType::typeName << " was specified for pair "
            << key.name() << ", but one of these phases is stationary. "
            << "Mass transfer is not supported on stationary phases"
            << exit(FatalError);
    }
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
                    this->mesh_.time().timeName(),
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
                    this->mesh_.time().timeName(),
                    this->mesh_
                ),
                this->mesh_,
                dimensioned<Type>("zero", dims, pTraits<Type>::zero)
            )
        );
    }
}


template<class modelType>
bool Foam::phaseSystem::foundSubModel(const phasePair& key) const
{
    const word name(IOobject::groupName(modelType::typeName, key.name()));

    if (key.ordered())
    {
        if (mesh().foundObject<modelType>(name))
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        if
        (
            mesh().foundObject<modelType>(name)
         ||
            mesh().foundObject<modelType>
            (
                IOobject::groupName(modelType::typeName, key.otherName())
            )
        )
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}


template<class modelType>
const modelType& Foam::phaseSystem::lookupSubModel(const phasePair& key) const
{
    const word name(IOobject::groupName(modelType::typeName, key.name()));

    if (key.ordered() || mesh().foundObject<modelType>(name))
    {
        return mesh().lookupObject<modelType>(name);
    }
    else
    {
        return
            mesh().lookupObject<modelType>
            (
                IOobject::groupName(modelType::typeName, key.otherName())
            );
    }
}


template<class modelType>
bool Foam::phaseSystem::foundSubModel
(
    const phaseModel& dispersed,
    const phaseModel& continuous
) const
{
    return foundSubModel<modelType>(orderedPhasePair(dispersed, continuous));
}


template<class modelType>
const modelType& Foam::phaseSystem::lookupSubModel
(
    const phaseModel& dispersed,
    const phaseModel& continuous
) const
{
    return lookupSubModel<modelType>(orderedPhasePair(dispersed, continuous));
}


template<class modelType>
bool Foam::phaseSystem::foundBlendedSubModel(const phasePair& key) const
{
    if
    (
        mesh().foundObject<BlendedInterfacialModel<modelType>>
        (
            IOobject::groupName
            (
                BlendedInterfacialModel<modelType>::typeName,
                key.name()
            )
        )
     || mesh().foundObject<BlendedInterfacialModel<modelType>>
        (
            IOobject::groupName
            (
                BlendedInterfacialModel<modelType>::typeName,
                key.otherName()
            )
        )
    )
    {
        return true;
    }
    else
    {
        return false;
    }
}


template<class modelType>
const Foam::BlendedInterfacialModel<modelType>&
Foam::phaseSystem::lookupBlendedSubModel(const phasePair& key) const
{
    const word name
    (
        IOobject::groupName
        (
            BlendedInterfacialModel<modelType>::typeName,
            key.name()
        )
    );

    if (mesh().foundObject<BlendedInterfacialModel<modelType>>(name))
    {
        return mesh().lookupObject<BlendedInterfacialModel<modelType>>(name);
    }
    else
    {
        return
            mesh().lookupObject<BlendedInterfacialModel<modelType>>
            (
                IOobject::groupName
                (
                    BlendedInterfacialModel<modelType>::typeName,
                    key.otherName()
                )
            );
    }
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField, class Group>
inline void addField
(
    const Group& group,
    const word& name,
    tmp<GeoField> field,
    PtrList<GeoField>& fieldList
)
{
    if (fieldList.set(group.index()))
    {
        fieldList[group.index()] += field;
    }
    else
    {
        fieldList.set
        (
            group.index(),
            new GeoField
            (
                IOobject::groupName(name, group.name()),
                field
            )
        );
    }
}


template<class GeoField, class Group>
inline void addField
(
    const Group& group,
    const word& name,
    const GeoField& field,
    PtrList<GeoField>& fieldList
)
{
    addField(group, name, tmp<GeoField>(field), fieldList);
}


template<class GeoField, class Group>
inline void addField
(
    const Group& group,
    const word& name,
    tmp<GeoField> field,
    HashPtrTable<GeoField>& fieldTable
)
{
    if (fieldTable.found(group.name()))
    {
        *fieldTable[group.name()] += field;
    }
    else
    {
        fieldTable.set
        (
            group.name(),
            new GeoField
            (
                IOobject::groupName(name, group.name()),
                field
            )
        );
    }
}


template<class GeoField, class Group>
inline void addField
(
    const Group& group,
    const word& name,
    const GeoField& field,
    HashPtrTable<GeoField>& fieldTable
)
{
    addField(group, name, tmp<GeoField>(field), fieldTable);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
