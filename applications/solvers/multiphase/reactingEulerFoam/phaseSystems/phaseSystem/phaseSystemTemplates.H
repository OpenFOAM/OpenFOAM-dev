/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    >& models
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
        autoPtr<BlendedInterfacialModel<modelType> >,
        phasePairKey,
        phasePairKey::hash
    >& models
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
      : blendingMethods_["default"]
    );

    autoPtr<modelType> noModel(NULL);

    forAllConstIter(typename modelTypeTable, tempModels, iter)
    {
        if (!iter().valid())
        {
            continue;
        }

        const phasePairKey key(iter.key().first(), iter.key().second());
        const phasePairKey key1In2(key.first(), key.second(), true);
        const phasePairKey key2In1(key.second(), key.first(), true);

        models.insert
        (
            key,
            autoPtr<BlendedInterfacialModel<modelType> >
            (
                new BlendedInterfacialModel<modelType>
                (
                    phaseModels_[key.first()],
                    phaseModels_[key.second()],
                    blending,
                    tempModels.found(key    ) ? tempModels[key    ] : noModel,
                    tempModels.found(key1In2) ? tempModels[key1In2] : noModel,
                    tempModels.found(key2In1) ? tempModels[key2In1] : noModel
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
        HashTable<autoPtr<modelType> >,
        phasePairKey,
        phasePairKey::hash
    >& models
)
{
    typedef
        HashTable<autoPtr<modelType>, phasePairKey, phasePairKey::hash>
        modelTypeTable;

    forAll(phaseModels_, phasei)
    {
        modelTypeTable tempModels;
        generatePairsAndSubModels
        (
            IOobject::groupName(modelName, phaseModels_[phasei].name()),
            tempModels
        );

        forAllConstIter(typename modelTypeTable, tempModels, tempModelIter)
        {
            const phasePairKey key(tempModelIter.key());

            if (!models.found(key))
            {
                models.insert
                (
                    key,
                    HashTable<autoPtr<modelType> >()
                );
            }

            models[tempModelIter.key()].insert
            (
                phaseModels_[phasei].name(),
                *tempModelIter
            );
        }
    }
}

template <class modelType>
const modelType& Foam::phaseSystem::lookupSubModel(const phasePair& key) const
{
    return
        mesh().lookupObject<modelType>
        (
            IOobject::groupName(modelType::typeName, key.name())
        );
}


template <class modelType>
const modelType& Foam::phaseSystem::lookupSubModel
(
    const phaseModel& dispersed,
    const phaseModel& continuous
) const
{
    return lookupSubModel<modelType>(orderedPhasePair(dispersed, continuous));
}


// ************************************************************************* //
