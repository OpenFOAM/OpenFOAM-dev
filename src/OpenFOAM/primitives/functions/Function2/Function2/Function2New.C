/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#include "Constant.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::Function2<Type>> Foam::Function2<Type>::New
(
    const word& name,
    const dictionary& dict
)
{
    if (dict.isDict(name))
    {
        const dictionary& coeffsDict(dict.subDict(name));

        const word Function2Type(coeffsDict.lookup("type"));

        typename dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(Function2Type);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown Function2 type "
                << Function2Type << " for Function2 "
                << name << nl << nl
                << "Valid Function2 types are:" << nl
                << dictionaryConstructorTablePtr_->sortedToc() << nl
                << exit(FatalError);
        }

        return cstrIter()(name, coeffsDict);
    }
    else
    {
        Istream& is(dict.lookup(name, false));

        token firstToken(is);
        word Function2Type;

        if (!firstToken.isWord())
        {
            is.putBack(firstToken);
            return autoPtr<Function2<Type>>
            (
                new Function2s::Constant<Type>(name, is)
            );
        }
        else
        {
            Function2Type = firstToken.wordToken();
        }

        typename dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(Function2Type);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown Function2 type "
                << Function2Type << " for Function2 "
                << name << nl << nl
                << "Valid Function2 types are:" << nl
                << dictionaryConstructorTablePtr_->sortedToc() << nl
                << exit(FatalError);
        }

        return cstrIter()(name, dict);
    }
}


// ************************************************************************* //
