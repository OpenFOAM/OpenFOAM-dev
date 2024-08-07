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

#include "Constant.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::Function1<Type>> Foam::Function1<Type>::New
(
    const word& name,
    const Function1s::unitConversions& units,
    const dictionary& dict
)
{
    // If the function is a dictionary (preferred) then read straightforwardly
    if (dict.isDict(name))
    {
        const dictionary& coeffDict(dict.subDict(name));

        const word Function1Type(coeffDict.lookup("type"));

        typename dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(Function1Type);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalIOErrorInFunction(dict)
                << "Unknown Function1 type "
                << Function1Type << " for Function1 "
                << name << nl << nl
                << "Valid Function1 types are:" << nl
                << dictionaryConstructorTablePtr_->sortedToc() << nl
                << exit(FatalIOError);
        }

        return cstrIter()(name, units, coeffDict);
    }

    // Find the entry
    Istream& is(dict.lookup(name));

    // Peek at the first token
    token firstToken(is);
    is.putBack(firstToken);

    // Read the type, or assume constant
    const word Function1Type =
        firstToken.isWord() ? word(is) : Function1s::Constant<Type>::typeName;

    // If the entry is not a type followed by a end statement then
    // construct the function from the stream
    if (!firstToken.isWord() || !is.eof())
    {
        return New(name, units, Function1Type, is);
    }

    // Otherwise, construct from the current or coeffs dictionary
    typename dictionaryConstructorTable::iterator dictCstrIter =
        dictionaryConstructorTablePtr_->find(Function1Type);

    if (dictCstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown Function1 type "
            << Function1Type << " for Function1 "
            << name << nl << nl
            << "Valid Function1 types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalIOError);
    }

    const bool haveCoeffsDict = dict.found(name + "Coeffs");

    autoPtr<Function1<Type>> funcPtr
    (
        dictCstrIter()
        (
            name,
            units,
            haveCoeffsDict ? dict.subDict(name + "Coeffs") : dict
        )
    );

    if (haveCoeffsDict)
    {
        IOWarningInFunction(dict)
            << "Using deprecated "
            << (name + "Coeffs") << " sub-dictionary."<< nl
            << "    Please use the simpler form" << endl;
        funcPtr->write(Info, units);
    }

    return funcPtr;
}


template<class Type>
Foam::autoPtr<Foam::Function1<Type>> Foam::Function1<Type>::New
(
    const word& name,
    const unitConversion& xUnits,
    const unitConversion& valueUnits,
    const dictionary& dict
)
{
    return New(name, {xUnits, valueUnits}, dict);
}


template<class Type>
Foam::autoPtr<Foam::Function1<Type>> Foam::Function1<Type>::New
(
    const word& name,
    const Function1s::unitConversions& units,
    const word& Function1Type,
    Istream& is
)
{
    typename dictionaryConstructorTable::iterator dictCstrIter =
        dictionaryConstructorTablePtr_->find(Function1Type);
    const bool haveDictCstrIter =
        dictCstrIter != dictionaryConstructorTablePtr_->end();

    typename IstreamConstructorTable::iterator isCstrIter =
        IstreamConstructorTablePtr_->find(Function1Type);
    const bool haveIstreamCstrIter =
        isCstrIter != IstreamConstructorTablePtr_->end();

    if (!haveDictCstrIter && !haveIstreamCstrIter)
    {
        FatalErrorInFunction
            << "Unknown Function1 type "
            << Function1Type << " for Function1 "
            << name << nl << nl
            << "Valid Function1 types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    if (!haveIstreamCstrIter)
    {
        FatalErrorInFunction
            << "Function1 type "
            << Function1Type << " for Function1 "
            << name << " cannot be specified inline" << nl << nl
            << "Make " << name << " a sub-dictionary"
            << exit(FatalError);
    }

    return isCstrIter()(name, units, is);
}


template<class Type>
Foam::autoPtr<Foam::Function1<Type>> Foam::Function1<Type>::New
(
    const word& name,
    const unitConversion& xUnits,
    const unitConversion& valueUnits,
    const word& Function1Type,
    Istream& is
)
{
    return New(name, {xUnits, valueUnits}, Function1Type, is);
}


// ************************************************************************* //
