/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "dictionary.H"
#include "primitiveEntry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, class ... KeysAndTs>
Foam::dictionary::dictionary
(
    const keyType& k,
    const T& t,
    const KeysAndTs& ... keysAndTs
)
:
    parent_(dictionary::null)
{
    set(k, t, keysAndTs ...);
}


template<class T, class ... KeysAndTs>
Foam::dictionary::dictionary
(
    const fileName& name,
    const keyType& k,
    const T& t,
    const KeysAndTs& ... keysAndTs
)
:
    dictionaryName(name),
    parent_(dictionary::null)
{
    set(k, t, keysAndTs ...);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
T Foam::dictionary::lookup
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    const entry* entryPtr = lookupEntryPtr(keyword, recursive, patternMatch);

    if (entryPtr == nullptr)
    {
        FatalIOErrorInFunction
        (
            *this
        )   << "keyword " << keyword << " is undefined in dictionary "
            << name()
            << exit(FatalIOError);
    }

    return pTraits<T>(entryPtr->stream());
}


template<class T>
T Foam::dictionary::lookupBackwardsCompatible
(
    const wordList& keywords,
    bool recursive,
    bool patternMatch
) const
{
    const entry* entryPtr =
        lookupEntryPtrBackwardsCompatible(keywords, recursive, patternMatch);

    if (entryPtr)
    {
        return pTraits<T>(entryPtr->stream());
    }
    else
    {
        // Generate error message using the first keyword
        return lookup<T>(keywords[0], recursive, patternMatch);
    }
}


template<class T>
T Foam::dictionary::lookupOrDefault
(
    const word& keyword,
    const T& deflt,
    bool recursive,
    bool patternMatch
) const
{
    const entry* entryPtr = lookupEntryPtr(keyword, recursive, patternMatch);

    if (entryPtr)
    {
        return pTraits<T>(entryPtr->stream());
    }
    else
    {
        if (writeOptionalEntries)
        {
            IOInfoInFunction(*this)
                << "Optional entry '" << keyword << "' is not present,"
                << " returning the default value '" << deflt << "'"
                << endl;
        }

        return deflt;
    }
}


template<class T>
T Foam::dictionary::lookupOrDefaultBackwardsCompatible
(
    const wordList& keywords,
    const T& deflt,
    bool recursive,
    bool patternMatch
) const
{
    const entry* entryPtr =
        lookupEntryPtrBackwardsCompatible(keywords, recursive, patternMatch);

    if (entryPtr)
    {
        return pTraits<T>(entryPtr->stream());
    }
    else
    {
        // Generate debugging messages using the first keyword
        return lookupOrDefault<T>(keywords[0], deflt, recursive, patternMatch);
    }
}


template<class T>
T Foam::dictionary::lookupOrAddDefault
(
    const word& keyword,
    const T& deflt,
    bool recursive,
    bool patternMatch
)
{
    const entry* entryPtr = lookupEntryPtr(keyword, recursive, patternMatch);

    if (entryPtr)
    {
        return pTraits<T>(entryPtr->stream());
    }
    else
    {
        if (writeOptionalEntries)
        {
            IOInfoInFunction(*this)
                << "Optional entry '" << keyword << "' is not present,"
                << " adding and returning the default value '" << deflt << "'"
                << endl;
        }

        add(new primitiveEntry(keyword, deflt));
        return deflt;
    }
}


template<class T>
bool Foam::dictionary::readIfPresent
(
    const word& keyword,
    T& val,
    bool recursive,
    bool patternMatch
) const
{
    const entry* entryPtr = lookupEntryPtr(keyword, recursive, patternMatch);

    if (entryPtr)
    {
        entryPtr->stream() >> val;
        return true;
    }
    else
    {
        if (writeOptionalEntries)
        {
            IOInfoInFunction(*this)
                << "Optional entry '" << keyword << "' is not present,"
                << " the default value '" << val << "' will be used."
                << endl;
        }

        return false;
    }
}


template<class T>
T Foam::dictionary::lookupScoped
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    const entry* entryPtr =
        lookupScopedEntryPtr(keyword, recursive, patternMatch);

    if (entryPtr == nullptr)
    {
        FatalIOErrorInFunction
        (
            *this
        )   << "keyword " << keyword << " is undefined in dictionary "
            << name()
            << exit(FatalIOError);
    }

    return pTraits<T>(entryPtr->stream());
}


template<class T>
void Foam::dictionary::add(const keyType& k, const T& t, bool overwrite)
{
    add(new primitiveEntry(k, t), overwrite);
}


template<class T>
void Foam::dictionary::set(const keyType& k, const T& t)
{
    set(new primitiveEntry(k, t));
}


template<class T, class ... KeysAndTs>
void Foam::dictionary::set
(
    const keyType& k,
    const T& t,
    const KeysAndTs& ... keysAndTs
)
{
    set(new primitiveEntry(k, t));
    set(keysAndTs ...);
}


// * * * * * * * * * * * * * * * IOstream Functions  * * * * * * * * * * * * //

template<class EntryType>
void Foam::writeEntry
(
    Ostream& os,
    const word& entryName,
    const EntryType& value
)
{
    writeKeyword(os, entryName);
    writeEntry(os, value);
    os << token::END_STATEMENT << endl;
}


template<class EntryType>
void Foam::writeEntryIfDifferent
(
    Ostream& os,
    const word& entryName,
    const EntryType& value1,
    const EntryType& value2
)
{
    if (value1 != value2)
    {
        writeEntry(os, entryName, value2);
    }
}


// ************************************************************************* //
