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
    Functions for generating tables of interfacial values

\*---------------------------------------------------------------------------*/

#ifndef generateInterfacialValues_H
#define generateInterfacialValues_H

#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class ValueType>
Foam::HashTable
<
    ValueType,
    Foam::phaseInterfaceKey,
    Foam::phaseInterfaceKey::hash
> generateInterfacialValues
(
    const phaseSystem& fluid,
    const dictionary& dict,
    const wordHashSet& ignoreKeys = wordHashSet()
)
{
    HashTable<ValueType, phaseInterfaceKey, phaseInterfaceKey::hash> values;

    forAllConstIter(dictionary, dict, iter)
    {
        // Get the keyword name and skip if it is ignored
        const word& modelName = iter().keyword();
        if (ignoreKeys.found(modelName)) continue;

        // Construct the associated interface
        autoPtr<phaseInterface> interfacePtr =
            phaseInterface::New(fluid, iter().keyword());

        // Read the vale and store in the table
        const ValueType value(pTraits<ValueType>(iter().stream()));
        values.insert(interfacePtr(), value);
    }

    return values;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
