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
    Functions for generating tables of blended interfacial models

\*---------------------------------------------------------------------------*/

#ifndef generateBlendedInterfacialModels_H
#define generateBlendedInterfacialModels_H

#include "generateInterfacialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class ModelType>
const dictionary& blendingDict
(
    const phaseSystem& fluid,
    const dictionary& dict
)
{
    if (dict.found("blending"))
    {
        if (dict.isDict("blending"))
        {
            return dict.subDict("blending");
        }

        const dictionary& fluidBlendingDict = fluid.subDict("blending");

        return fluidBlendingDict.subDict(dict.lookup<word>("blending"));
    }

    if (fluid.found("blending"))
    {
        const dictionary& fluidBlendingDict = fluid.subDict("blending");

        if (fluidBlendingDict.found(modelName<ModelType>()))
        {
            return fluidBlendingDict.subDict(modelName<ModelType>());
        }

        if (fluidBlendingDict.found("default"))
        {
            return fluidBlendingDict.subDict("default");
        }
    }

    return dict.subDict("blending");
}


template<class ModelType>
void checkBlendedInterfacialModelsDict
(
    const phaseSystem& fluid,
    const dictionary& dict,
    const wordHashSet& ignoreKeys = wordHashSet()
)
{
    checkInterfacialModelsDict<ModelType>
    (
        fluid,
        dict,
        ignoreKeys | wordHashSet({"blending"})
    );
}


template<class ModelType, class ... Args>
Foam::HashPtrTable
<
    ModelType,
    Foam::phaseInterfaceKey,
    Foam::phaseInterfaceKey::hash
> generateBlendedInterfacialModels
(
    const phaseSystem& fluid,
    const dictionary& dict,
    const wordHashSet& ignoreKeys,
    const bool ignoreNonModelPhaseInterfaceTypes,
    const Args& ... args
)
{
    if (dict.empty())
    {
        return
            HashPtrTable
            <
                ModelType,
                phaseInterfaceKey,
                phaseInterfaceKey::hash
            >();
    }

    return
        generateInterfacialModels<ModelType, dictionary, Args ...>
        (
            fluid,
            dict,
            ignoreKeys | wordHashSet({"blending"}),
            ignoreNonModelPhaseInterfaceTypes,
            blendingDict<ModelType>(fluid, dict),
            args ...
        );
}


template<class ModelType>
Foam::HashPtrTable
<
    ModelType,
    Foam::phaseInterfaceKey,
    Foam::phaseInterfaceKey::hash
> generateBlendedInterfacialModels
(
    const phaseSystem& fluid,
    const dictionary& dict
)
{
    return
        generateBlendedInterfacialModels<ModelType>
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
