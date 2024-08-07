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

#include "dictionary.H"
#include "primitiveEntry.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
T Foam::dictionary::readType
(
    const word& keyword,
    const unitConversion& defaultUnits,
    ITstream& is
) const
{
    assertNoConvertUnits(pTraits<T>::typeName, keyword, defaultUnits, is);

    return pTraits<T>(is);
}


template<class T>
T Foam::dictionary::readType(const word& keyword, ITstream& is) const
{
    return pTraits<T>(is);
}


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
        FatalIOErrorInFunction(*this)
            << "keyword " << keyword << " is undefined in dictionary "
            << name() << exit(FatalIOError);
    }

    return readType<T>(keyword, entryPtr->stream());
}


template<class T>
T Foam::dictionary::lookup
(
    const word& keyword,
    const unitConversion& defaultUnits,
    bool recursive,
    bool patternMatch
) const
{
    const entry* entryPtr = lookupEntryPtr(keyword, recursive, patternMatch);

    if (entryPtr == nullptr)
    {
        FatalIOErrorInFunction(*this)
            << "keyword " << keyword << " is undefined in dictionary "
            << name() << exit(FatalIOError);
    }

    return readType<T>(keyword, defaultUnits, entryPtr->stream());
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
        return readType<T>(entryPtr->keyword(), entryPtr->stream());
    }
    else
    {
        // Generate error message using the first keyword
        return lookup<T>(keywords[0], recursive, patternMatch);
    }
}


template<class T>
T Foam::dictionary::lookupBackwardsCompatible
(
    const wordList& keywords,
    const unitConversion& defaultUnits,
    bool recursive,
    bool patternMatch
) const
{
    const entry* entryPtr =
        lookupEntryPtrBackwardsCompatible(keywords, recursive, patternMatch);

    if (entryPtr)
    {
        return
            readType<T>(entryPtr->keyword(), defaultUnits, entryPtr->stream());
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
    const T& defaultValue,
    const bool writeDefault
) const
{
    const entry* entryPtr = lookupEntryPtr(keyword, false, false);

    if (entryPtr)
    {
        return readType<T>(keyword, entryPtr->stream());
    }
    else
    {
        if (writeDefault)
        {
            Info<< indent << "Default: " << keyword
                << " " << defaultValue
                << " in " << name().relativePath() << endl;
        }

        return defaultValue;
    }
}


template<class T>
T Foam::dictionary::lookupOrDefault
(
    const word& keyword,
    const unitConversion& defaultUnits,
    const T& defaultValue,
    const bool writeDefault
) const
{
    const entry* entryPtr = lookupEntryPtr(keyword, false, false);

    if (entryPtr)
    {
        return readType<T>(keyword, defaultUnits, entryPtr->stream());
    }
    else
    {
        if (writeDefault)
        {
            Info<< indent << "Default: " << keyword
                << " " << defaultValue
                << " in " << name().relativePath() << endl;
        }

        return defaultValue;
    }
}


template<class T>
T Foam::dictionary::lookupOrDefaultBackwardsCompatible
(
    const wordList& keywords,
    const T& defaultValue
) const
{
    const entry* entryPtr =
        lookupEntryPtrBackwardsCompatible(keywords, false, false);

    if (entryPtr)
    {
        return readType<T>(entryPtr->keyword(), entryPtr->stream());
    }
    else
    {
        // Generate debugging messages using the first keyword
        return lookupOrDefault<T>(keywords[0], defaultValue);
    }
}


template<class T>
T Foam::dictionary::lookupOrDefaultBackwardsCompatible
(
    const wordList& keywords,
    const unitConversion& defaultUnits,
    const T& defaultValue
) const
{
    const entry* entryPtr =
        lookupEntryPtrBackwardsCompatible(keywords, false, false);

    if (entryPtr)
    {
        return
            readType<T>(entryPtr->keyword(), defaultUnits, entryPtr->stream());
    }
    else
    {
        // Generate debugging messages using the first keyword
        return lookupOrDefault<T>(keywords[0], defaultValue);
    }
}


template<class T>
T Foam::dictionary::lookupOrAddDefault
(
    const word& keyword,
    const T& defaultValue
)
{
    const entry* entryPtr = lookupEntryPtr(keyword, false, false);

    if (entryPtr)
    {
        return readType<T>(keyword, entryPtr->stream());
    }
    else
    {
        if (writeOptionalEntries > 1)
        {
            Info<< indent << "Added default: " << keyword
                << " " << defaultValue
                << " to " << name().relativePath() << endl;
        }

        add(new primitiveEntry(keyword, defaultValue));
        return defaultValue;
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
        val = readType<T>(keyword, entryPtr->stream());
        return true;
    }
    else
    {
        if (writeOptionalEntries > 1)
        {
            Info<< indent << "Default: " << keyword
                << " " << val
                << " in " << name().relativePath() << endl;
        }

        return false;
    }
}


template<class T>
bool Foam::dictionary::readIfPresent
(
    const word& keyword,
    const unitConversion& defaultUnits,
    T& val,
    bool recursive,
    bool patternMatch
) const
{
    const entry* entryPtr = lookupEntryPtr(keyword, recursive, patternMatch);

    if (entryPtr)
    {
        val = readType<T>(keyword, defaultUnits, entryPtr->stream());
        return true;
    }
    else
    {
        if (writeOptionalEntries > 1)
        {
            Info<< indent << "Default: " << keyword
                << " " << val
                << " in " << name().relativePath() << endl;
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
        FatalIOErrorInFunction(*this)
            << "keyword " << keyword << " is undefined in dictionary "
            << name() << exit(FatalIOError);
    }

    return pTraits<T>(entryPtr->stream());
}


template<class T>
const T& Foam::dictionary::lookupCompoundScoped
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

    token firstToken(entryPtr->stream());

    return dynamicCast<const token::Compound<T>>(firstToken.compoundToken());
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
void Foam::writeEntry
(
    Ostream& os,
    const word& entryName,
    const unitConversion& defaultUnits,
    const EntryType& value
)
{
    writeKeyword(os, entryName);
    writeEntry(os, defaultUnits, value);
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


template<class EntryType>
void Foam::writeEntryIfDifferent
(
    Ostream& os,
    const word& entryName,
    const unitConversion& defaultUnits,
    const EntryType& value1,
    const EntryType& value2
)
{
    if (value1 != value2)
    {
        writeEntry(os, entryName, defaultUnits, value2);
    }
}


// ************************************************************************* //
