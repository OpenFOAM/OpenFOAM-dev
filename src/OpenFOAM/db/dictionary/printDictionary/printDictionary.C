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

#include "printDictionary.H"
#include "dictionary.H"
#include "stringOps.H"

// * * * * * * * * * * * * Private Static Data Members * * * * * * * * * * * //

namespace Foam
{
    HashTable<Tuple2<const dictionary*, label>, fileName, Hash<fileName>>
        printDictionary::dictNameToDictPtrAndCount_;

    HashTable<tmpNrc<dictionary>, const dictionary*, Hash<void*>>
        printDictionary::dictPtrToDefaults_;
}


// * * * * * * * * * * * Private Static Member Functions * * * * * * * * * * //

void Foam::printDictionary::removeDefaults
(
    const dictionary* dictPtr,
    HashSet<const dictionary*, Hash<void*>>& removeDictPtrs
)
{
    typedef
        HashTable<tmpNrc<dictionary>, const dictionary*, Hash<void*>>
        dictPtrToDefaultsType;

    if (!dictPtrToDefaults_.found(dictPtr)) return;

    const dictionary& defaults = dictPtrToDefaults_[dictPtr]();

    forAllConstIter(dictionary, defaults, iter)
    {
        if (!iter().isDict()) continue;

        const dictionary& subDefaults = iter().dict();

        forAllConstIter(dictPtrToDefaultsType, dictPtrToDefaults_, jter)
        {
            if (&jter()() == &subDefaults)
            {
                removeDefaults(jter.key(), removeDictPtrs);
                break;
            }
        }
    }

    removeDictPtrs.insert(dictPtr);
}


void Foam::printDictionary::removeDefaults(const dictionary* dictPtr)
{
    typedef
        HashSet<const dictionary*, Hash<void*>>
        removeDictPtrsType;

    removeDictPtrsType removeDictPtrs;

    removeDefaults(dictPtr, removeDictPtrs);

    forAllConstIter(removeDictPtrsType, removeDictPtrs, iter)
    {
        dictPtrToDefaults_.erase(iter.key());
    }
}


void Foam::printDictionary::setDefaults(const dictionary& dict)
{
    removeDefaults(&dict);

    dictPtrToDefaults_.set
    (
        &dict,
        tmpNrc<dictionary>(new dictionary(dict.parent(), dictionary()))
    );

    setSubDefaults(dict, printDictionary::defaults(dict));
}


void Foam::printDictionary::setSubDefaults
(
    const dictionary& dict,
    dictionary& defaults
)
{
    forAllConstIter(dictionary, dict, iter)
    {
        if (!iter().isDict()) continue;

        const dictionary& subDict = iter().dict();

        defaults.set(iter().keyword(), dictionary());

        dictionary& subDefaults = defaults.subDict(iter().keyword());

        if (!dictPtrToDefaults_.found(&subDict))
        {
            dictPtrToDefaults_.set
            (
                &subDict,
                tmpNrc<dictionary>(subDefaults)
            );
        }

        setSubDefaults(subDict, subDefaults);
    }
}


void Foam::printDictionary::print
(
    const dictionary& dict,
    const dictionary& defaults
)
{
    Info<< nl << indent << token::BEGIN_BLOCK << nl << incrIndent;

    forAllConstIter(dictionary, dict, iter)
    {
        if (!iter().isDict())
        {
            Info<< iter();
        }
        else if (defaults.isDict(iter().keyword()))
        {
            Info<< indent << iter().keyword();
            print(iter().dict(), defaults.subDict(iter().keyword()));
        }
        else
        {
            Info<< indent << iter().keyword()
                << nl << indent << token::BEGIN_BLOCK << nl << incrIndent;

            iter().dict().write(Info, false);

            OStringStream oss;
            oss << "/* #print must be specified after the " << iter().keyword()
                << " sub-dictionary in order for its defaults to be printed */";

            Info<<
                stringOps::breakIntoIndentedLines
                (
                    oss.str(),
                    80,
                    Info().indentSize()
                ).c_str() << endl
                << decrIndent << indent << token::END_BLOCK << nl;
        }
    }

    label nDefaultsEntries = 0;
    forAllConstIter(dictionary, defaults, iter)
    {
        nDefaultsEntries += !iter().isDict();
    }

    if (nDefaultsEntries) Info<< indent << "/* defaults */" << endl;

    forAllConstIter(dictionary, defaults, iter)
    {
        if (!iter().isDict())
        {
            Info<< iter();
        }
    }

    Info<< decrIndent << indent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::printDictionary::add(const dictionary& dict)
{
    if (dictNameToDictPtrAndCount_.found(dict.name()))
    {
        dicts_.append(&dict);
        dictNames_.append(fileName::null);

        dictNameToDictPtrAndCount_[dict.name()].second() ++;

        setDefaults(dict);
    }
}


void Foam::printDictionary::add(const fileName& dictName)
{
    dicts_.append(nullptr);
    dictNames_.append(dictName);

    dictNameToDictPtrAndCount_.set
    (
        dictName,
        Tuple2<const dictionary*, label>(nullptr, 1)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::printDictionary::printDictionary()
:
    dicts_(),
    dictNames_()
{
    Info<< incrIndent;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::printDictionary::~printDictionary()
{
    forAll(dicts_, i)
    {
        const word& dictName =
            dicts_.set(i) ? dicts_[i].name() : dictNames_[i];

        if (!dictNameToDictPtrAndCount_.found(dictName)) continue;

        Tuple2<const dictionary*, label>& dictPtrAndCount =
            dictNameToDictPtrAndCount_[dictName];
        const dictionary* dictPtr = dictPtrAndCount.first();
        label& count = dictPtrAndCount.second();

        count --;

        if (!dictPtrToDefaults_.found(dictPtr)) continue;

        if (dictPtr && count == 0 && dictPtrToDefaults_[dictPtr].isTmp())
        {
            Info<< indent << dictPtr->name().relativePath().c_str();

            print(*dictPtr, dictPtrToDefaults_[dictPtr]());
        }
    }

    Info<< decrIndent;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::printDictionary::set(const dictionary& dict)
{
    if (dictNameToDictPtrAndCount_.found(dict.name()))
    {
        setDefaults(dict);
    }

    const label count =
        dictNameToDictPtrAndCount_.found(dict.name())
      ? dictNameToDictPtrAndCount_[dict.name()].second()
      : 0;

    dictNameToDictPtrAndCount_.set
    (
        dict.name(),
        Tuple2<const dictionary*, label>(&dict, count)
    );
}


// ************************************************************************* //
