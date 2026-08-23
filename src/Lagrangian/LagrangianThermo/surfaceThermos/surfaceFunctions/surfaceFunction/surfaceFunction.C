/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "surfaceFunction.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfaceFunction, 0);
    defineRunTimeSelectionTable(surfaceFunction, word);
    defineRunTimeSelectionTable(surfaceFunction, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceFunction::surfaceFunction(const LagrangianMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::surfaceFunction>
Foam::surfaceFunction::New
(
    const LagrangianMesh& mesh,
    const dictionary& dict
)
{
    const bool foundDefault = dict.found(typeName);

    const bool foundSpecific =
        dict.found(typeName + 's')
     && dict.subDict(typeName + 's').found(mesh.name());

    if (!foundDefault && !foundSpecific)
    {
        FatalIOErrorInFunction(dict)
            << "neither of keywords " << typeName << " or "
            << word(typeName + 's') << '/' << mesh.name()
            << " defined in dictionary " << dict.name()
            << exit(FatalIOError);
    }

    const dictionary& subDict =
        foundDefault
      ? (
            dict.isDict(typeName)
          ? dict.subDict(typeName)
          : NullObjectRef<dictionary>()
        )
      : (
            dict.subDict(typeName + 's').isDict(mesh.name())
          ? dict.subDict(typeName + 's').subDict(mesh.name())
          : NullObjectRef<dictionary>()
        );

    const word type =
        isNull(subDict)
      ? (
            foundDefault
          ? dict.lookup<word>(typeName)
          : dict.subDict(typeName + 's').lookup<word>(mesh.name())
        )
      : subDict.lookup<word>("type");

    static wordConstructorTable::iterator wordConstructorTableEnd =
        wordConstructorTable::end();
    static dictionaryConstructorTable::iterator dictionaryConstructorTableEnd =
        dictionaryConstructorTable::end();

    wordConstructorTable::iterator wordCstrIter =
        wordConstructorTablePtr_
      ? wordConstructorTablePtr_->find(type)
      : wordConstructorTableEnd;
    dictionaryConstructorTable::iterator dictCstrIter =
        dictionaryConstructorTablePtr_
      ? dictionaryConstructorTablePtr_->find(type)
      : dictionaryConstructorTableEnd;

    if
    (
        wordCstrIter == wordConstructorTableEnd
     && dictCstrIter == dictionaryConstructorTableEnd
    )
    {
        wordList types;
        if (wordConstructorTablePtr_)
        {
            types.append(wordConstructorTablePtr_->toc());
        }
        if (dictionaryConstructorTablePtr_)
        {
            types.append(dictionaryConstructorTablePtr_->sortedToc());
        }
        sort(types);

        FatalIOErrorInFunction(dict)
            << "Unknown " << typeName << " type "
            << type << endl << endl
            << "Valid " << typeName << " types are : " << endl
            << types << exit(FatalIOError);
    }

    // If we don't have a sub-dictionary, but this model is only on the
    // dictionary table, then raise a "tried to read a primitive entry as a
    // dictionary" error
    if
    (
        isNull(subDict)
     && wordCstrIter == wordConstructorTableEnd
     && dictCstrIter != dictionaryConstructorTableEnd
    )
    {
        foundDefault
      ? dict.subDict(typeName)
      : dict.subDict(typeName + 's').subDict(mesh.name());
    }

    if
    (
        isNull(subDict)
     || dictCstrIter == dictionaryConstructorTableEnd
    )
    {
        return wordCstrIter()(mesh);
    }
    else
    {
        printDictionary print(subDict);

        return dictCstrIter()(mesh, subDict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceFunction::~surfaceFunction()
{}


// ************************************************************************* //
