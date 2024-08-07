/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "liquidProperties.H"

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::liquidProperties> Foam::liquidProperties::New
(
    const word& name
)
{
    if (debug)
    {
        InfoInFunction << "Constructing liquidProperties" << endl;
    }

    ConstructorTable::iterator cstrIter = ConstructorTablePtr_->find(name);

    if (cstrIter == ConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown liquidProperties type "
            << name << nl << nl
            << "Valid liquidProperties types are:" << nl
            << ConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<liquidProperties>(cstrIter()());
}


Foam::autoPtr<Foam::liquidProperties> Foam::liquidProperties::New
(
    const dictionary& dict
)
{
    if (debug)
    {
        InfoInFunction << "Constructing liquidProperties" << endl;
    }

    // If the type is not specified use the name as the liquid type name
    const word& liquidPropertiesTypeName =
        dict.found("type") ? dict.lookup("type") : dict.dictName();

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(liquidPropertiesTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown liquidProperties type "
            << liquidPropertiesTypeName << nl << nl
            << "Valid liquidProperties types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<liquidProperties>(cstrIter()(dict));
}


// ************************************************************************* //
