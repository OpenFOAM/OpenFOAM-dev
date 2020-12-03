/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
Foam::autoPtr<Foam::Function1<Type>> Foam::Function1<Type>::New
(
    const word& name,
    const dictionary& dict
)
{
    if (dict.isDict(name))
    {
        const dictionary& coeffsDict(dict.subDict(name));

        const word Function1Type(coeffsDict.lookup("type"));

        typename dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(Function1Type);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown Function1 type "
                << Function1Type << " for Function1 "
                << name << nl << nl
                << "Valid Function1 types are:" << nl
                << dictionaryConstructorTablePtr_->sortedToc() << nl
                << exit(FatalError);
        }

        return cstrIter()(name, coeffsDict);
    }
    else
    {
        Istream& is(dict.lookup(name, false));

        token firstToken(is);
        word Function1Type;

        if (!firstToken.isWord())
        {
            is.putBack(firstToken);
            return autoPtr<Function1<Type>>
            (
                new Function1s::Constant<Type>(name, is)
            );
        }
        else
        {
            Function1Type = firstToken.wordToken();
        }

        typename dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(Function1Type);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown Function1 type "
                << Function1Type << " for Function1 "
                << name << nl << nl
                << "Valid Function1 types are:" << nl
                << dictionaryConstructorTablePtr_->sortedToc() << nl
                << exit(FatalError);
        }

        autoPtr<Function1<Type>> funcPtr
        (
            cstrIter()
            (
                name,
                dict.found(name + "Coeffs")
              ? dict.subDict(name + "Coeffs")
              : dict
            )
        );

        if (dict.found(name + "Coeffs"))
        {
            IOWarningInFunction(dict)
                << "Using deprecated "
                << (name + "Coeffs") << " sub-dictionary."<< nl
                << "    Please use the simpler form" << endl;
            funcPtr->write(Info);
        }

        return funcPtr;
    }
}


// ************************************************************************* //
