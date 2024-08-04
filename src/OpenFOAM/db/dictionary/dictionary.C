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
#include "dictionaryEntry.H"
#include "regExp.H"
#include "OSHA1stream.H"
#include "unitConversion.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(dictionary, 0);
}

const Foam::dictionary Foam::dictionary::null;

bool Foam::dictionary::writeOptionalEntries
(
    Foam::debug::infoSwitch("writeOptionalEntries", 0)
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::entry* Foam::dictionary::lookupScopedSubEntryPtr
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    // Check for the dictionary boundary marker
    const string::size_type emarkPos = keyword.find('!');

    if (emarkPos == string::npos || emarkPos == 0)
    {
        // Lookup in this dictionary

        string::size_type slashPos = keyword.find('/');

        if (slashPos == string::npos)
        {
            // Non-scoped lookup
            return lookupEntryPtr(keyword, recursive, patternMatch);
        }
        else
        {
            // Extract the first word
            word firstWord = keyword.substr(0, slashPos);
            slashPos++;

            if (firstWord == ".")
            {
                return lookupScopedSubEntryPtr
                (
                    keyword.substr(slashPos),
                    false,
                    patternMatch
                );
            }
            else if (firstWord == "..")
            {
                // Go to parent
                if (&parent_ == &dictionary::null)
                {
                    FatalIOErrorInFunction(*this)
                        << "No parent of current dictionary"
                        << " when searching for "
                        << keyword.substr(slashPos, keyword.size() - slashPos)
                        << exit(FatalIOError);
                }

                return parent_.lookupScopedSubEntryPtr
                (
                    keyword.substr(slashPos),
                    false,
                    patternMatch
                );
            }
            else
            {
                const entry* entPtr = lookupScopedSubEntryPtr
                (
                    firstWord,
                    recursive,
                    patternMatch
                );

                if (entPtr && entPtr->isDict())
                {
                    return entPtr->dict().lookupScopedSubEntryPtr
                    (
                        keyword.substr(slashPos, keyword.size() - slashPos),
                        false,
                        patternMatch
                    );
                }
                else
                {
                    return nullptr;
                }
            }
        }
    }
    else
    {
        // Lookup in the dictionary specified by file name
        // created from the part of the keyword before the '!'

        fileName fName = keyword.substr(0, emarkPos);

        if (!fName.isAbsolute())
        {
            fName = topDict().name().path()/fName;
        }

        if (fName == topDict().name())
        {
            FatalIOErrorInFunction(*this)
                << "Attempt to re-read current dictionary " << fName
                << " for keyword "
                << keyword
                << exit(FatalIOError);
        }

        const word localKeyword = keyword.substr
        (
            emarkPos + 1,
            keyword.size() - emarkPos - 1
        );

        includedDictionary dict(fName, *this);

        const Foam::entry* entryPtr = dict.lookupScopedEntryPtr
        (
            localKeyword,
            recursive,
            patternMatch
        );

        if (!entryPtr)
        {
            FatalIOErrorInFunction(dict)
                << "keyword " << localKeyword
                << " is undefined in dictionary "
                << dict.name()
                << exit(FatalIOError);
        }

        return entryPtr->clone(*this).ptr();
    }
}


bool Foam::dictionary::findInPatterns
(
    const bool patternMatch,
    const word& Keyword,
    DLList<entry*>::const_iterator& wcLink,
    DLList<autoPtr<regExp>>::const_iterator& reLink
) const
{
    if (patternEntries_.size())
    {
        while (wcLink != patternEntries_.end())
        {
            if
            (
                patternMatch
              ? reLink()->match(Keyword)
              : wcLink()->keyword() == Keyword
            )
            {
                return true;
            }

            ++reLink;
            ++wcLink;
        }
    }

    return false;
}


bool Foam::dictionary::findInPatterns
(
    const bool patternMatch,
    const word& Keyword,
    DLList<entry*>::iterator& wcLink,
    DLList<autoPtr<regExp>>::iterator& reLink
)
{
    if (patternEntries_.size())
    {
        while (wcLink != patternEntries_.end())
        {
            if
            (
                patternMatch
              ? reLink()->match(Keyword)
              : wcLink()->keyword() == Keyword
            )
            {
                return true;
            }

            ++reLink;
            ++wcLink;
        }
    }

    return false;
}


void Foam::dictionary::assertNoConvertUnits
(
    const char* typeName,
    const word& keyword,
    const unitConversion& defaultUnits,
    ITstream& is
) const
{
    if (!defaultUnits.standard())
    {
        FatalIOErrorInFunction(is)
            << "Unit conversions are not supported when reading "
            << typeName << " types" << abort(FatalError);
    }
}


template<class T>
T Foam::dictionary::readTypeAndConvertUnits
(
    const word& keyword,
    const unitConversion& defaultUnits,
    ITstream& is
) const
{
    // Read the units if they are before the value
    unitConversion units(defaultUnits);
    const bool haveUnits = units.readIfPresent(keyword, *this, is);

    // Read the value
    T value = pTraits<T>(is);

    // Read the units if they are after the value
    if (!haveUnits && !is.eof())
    {
        units.readIfPresent(keyword, *this, is);
    }

    // Modify the value by the unit conversion
    units.makeStandard(value);

    return value;
}


#define IMPLEMENT_SPECIALISED_READ_TYPE(T, nullArg)                            \
                                                                               \
    template<>                                                                 \
    Foam::T Foam::dictionary::readType                                         \
    (                                                                          \
        const word& keyword,                                                   \
        const unitConversion& defaultUnits,                                    \
        ITstream& is                                                           \
    ) const                                                                    \
    {                                                                          \
        return readTypeAndConvertUnits<T>(keyword, defaultUnits, is);          \
    }                                                                          \
                                                                               \
    template<>                                                                 \
    Foam::T Foam::dictionary::readType                                         \
    (                                                                          \
        const word& keyword,                                                   \
        ITstream& is                                                           \
    ) const                                                                    \
    {                                                                          \
        return readTypeAndConvertUnits<T>(keyword, unitAny, is);               \
    }

#define IMPLEMENT_SPECIALISED_READ_LIST_TYPE(T, nullArg)                       \
    IMPLEMENT_SPECIALISED_READ_TYPE(List<Foam::T>, nullArg)

FOR_ALL_FIELD_TYPES(IMPLEMENT_SPECIALISED_READ_TYPE)
FOR_ALL_FIELD_TYPES(IMPLEMENT_SPECIALISED_READ_LIST_TYPE)

#undef IMPLEMENT_SPECIALISED_READ_TYPE


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dictionary::dictionary()
:
    parent_(dictionary::null)
{}


Foam::dictionary::dictionary(const fileName& name)
:
    dictionaryName(name),
    parent_(dictionary::null)
{}


Foam::dictionary::dictionary
(
    const word& name,
    const dictionary& parentDict
)
:
    dictionaryName(name),
    parent_(parentDict)
{}


Foam::dictionary::dictionary
(
    const dictionary& parentDict,
    const dictionary& dict
)
:
    dictionaryName(dict.name()),
    IDLList<entry>(dict, *this),
    parent_(parentDict)
{
    forAllIter(IDLList<entry>, *this, iter)
    {
        hashedEntries_.insert(iter().keyword(), &iter());

        if (iter().keyword().isPattern())
        {
            patternEntries_.insert(&iter());
            patternRegexps_.insert
            (
                autoPtr<regExp>(new regExp(iter().keyword()))
            );
        }
    }
}


Foam::dictionary::dictionary(const dictionary& dict)
:
    dictionaryName(dict.name()),
    IDLList<entry>(dict, *this),
    parent_(dictionary::null)
{
    forAllIter(IDLList<entry>, *this, iter)
    {
        hashedEntries_.insert(iter().keyword(), &iter());

        if (iter().keyword().isPattern())
        {
            patternEntries_.insert(&iter());
            patternRegexps_.insert
            (
                autoPtr<regExp>(new regExp(iter().keyword()))
            );
        }
    }
}


Foam::dictionary::dictionary(const dictionary* dictPtr)
:
    parent_(dictionary::null)
{
    if (dictPtr)
    {
        operator=(*dictPtr);
    }
}


Foam::autoPtr<Foam::dictionary> Foam::dictionary::clone() const
{
    return autoPtr<dictionary>(new dictionary(*this));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dictionary::~dictionary()
{
    // cerr<< "~dictionary() " << name() << " " << long(this) << std::endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary& Foam::dictionary::topDict() const
{
    const dictionary& p = parent();

    if (&p != this && !p.name().empty())
    {
        return p.topDict();
    }
    else
    {
        return *this;
    }
}


Foam::word Foam::dictionary::topDictKeyword() const
{
    const dictionary& p = parent();

    if (&p != this && !p.name().empty())
    {
        const word pKeyword = p.topDictKeyword();
        const char pSeparator = '/';
        return
            pKeyword == word::null
          ? dictName()
          : word(pKeyword + pSeparator + dictName());
    }
    else
    {
        return word::null;
    }
}


Foam::label Foam::dictionary::startLineNumber() const
{
    if (size())
    {
        return first()->startLineNumber();
    }
    else
    {
        return -1;
    }
}


Foam::label Foam::dictionary::endLineNumber() const
{
    if (size())
    {
        return last()->endLineNumber();
    }
    else
    {
        return -1;
    }
}


Foam::SHA1Digest Foam::dictionary::digest() const
{
    OSHA1stream os;

    // Process entries
    forAllConstIter(IDLList<entry>, *this, iter)
    {
        os << *iter;
    }

    return os.digest();
}


Foam::tokenList Foam::dictionary::tokens() const
{
    // Serialise dictionary into a string
    OStringStream os;
    write(os, false);
    IStringStream is(os.str());

    // Parse string as tokens
    DynamicList<token> tokens;
    token t;
    while (is.read(t))
    {
        tokens.append(t);
    }

    return tokenList(move(tokens));
}


bool Foam::dictionary::found
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    if (hashedEntries_.found(keyword))
    {
        return true;
    }
    else
    {
        if (patternMatch && patternEntries_.size())
        {
            DLList<entry*>::const_iterator wcLink =
                patternEntries_.begin();
            DLList<autoPtr<regExp>>::const_iterator reLink =
                patternRegexps_.begin();

            // Find in patterns using regular expressions only
            if (findInPatterns(patternMatch, keyword, wcLink, reLink))
            {
                return true;
            }
        }

        if (recursive && &parent_ != &dictionary::null)
        {
            return parent_.found(keyword, recursive, patternMatch);
        }
        else
        {
            return false;
        }
    }
}


const Foam::entry* Foam::dictionary::lookupEntryPtr
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    HashTable<entry*>::const_iterator iter = hashedEntries_.find(keyword);

    if (iter == hashedEntries_.end())
    {
        if (patternMatch && patternEntries_.size())
        {
            DLList<entry*>::const_iterator wcLink =
                patternEntries_.begin();
            DLList<autoPtr<regExp>>::const_iterator reLink =
                patternRegexps_.begin();

            // Find in patterns using regular expressions only
            if (findInPatterns(patternMatch, keyword, wcLink, reLink))
            {
                return wcLink();
            }
        }

        if (recursive && &parent_ != &dictionary::null)
        {
            return parent_.lookupEntryPtr(keyword, recursive, patternMatch);
        }
        else
        {
            return nullptr;
        }
    }

    return iter();
}


Foam::entry* Foam::dictionary::lookupEntryPtr
(
    const word& keyword,
    bool recursive,
    bool patternMatch
)
{
    HashTable<entry*>::iterator iter = hashedEntries_.find(keyword);

    if (iter == hashedEntries_.end())
    {
        if (patternMatch && patternEntries_.size())
        {
            DLList<entry*>::iterator wcLink =
                patternEntries_.begin();
            DLList<autoPtr<regExp>>::iterator reLink =
                patternRegexps_.begin();

            // Find in patterns using regular expressions only
            if (findInPatterns(patternMatch, keyword, wcLink, reLink))
            {
                return wcLink();
            }
        }

        if (recursive && &parent_ != &dictionary::null)
        {
            return const_cast<dictionary&>(parent_).lookupEntryPtr
            (
                keyword,
                recursive,
                patternMatch
            );
        }
        else
        {
            return nullptr;
        }
    }

    return iter();
}


const Foam::entry* Foam::dictionary::lookupEntryPtrBackwardsCompatible
(
    const wordList& keywords,
    bool recursive,
    bool patternMatch
) const
{
    const entry* result = nullptr;

    forAll(keywords, keywordi)
    {
        const entry* entryPtr =
            lookupEntryPtr(keywords[keywordi], recursive, patternMatch);

        if (entryPtr)
        {
            if (result)
            {
                IOWarningInFunction((*this))
                    << "Duplicate backwards compatible keywords \""
                    << result->keyword() << "\" and \"" << entryPtr->keyword()
                    << "\" are defined in dictionary " << name() << endl
                    << "The preferred keyword for this entry is \""
                    << keywords[0] << "\"" << endl;
            }
            else
            {
                result = entryPtr;
            }
        }
    }

    return result;
}


const Foam::entry& Foam::dictionary::lookupEntry
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
            << name()
            << exit(FatalIOError);
    }

    return *entryPtr;
}


const Foam::entry& Foam::dictionary::lookupEntryBackwardsCompatible
(
    const wordList& keywords,
    bool recursive,
    bool patternMatch
) const
{
    const entry* entryPtr =
        lookupEntryPtrBackwardsCompatible(keywords, recursive, patternMatch);

    if (entryPtr == nullptr)
    {
        // Generate error message using the first keyword
        return lookupEntry(keywords[0], recursive, patternMatch);
    }
    else
    {
        return *entryPtr;
    }
}


Foam::ITstream& Foam::dictionary::lookup
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    return lookupEntry(keyword, recursive, patternMatch).stream();
}


Foam::ITstream& Foam::dictionary::lookupBackwardsCompatible
(
    const wordList& keywords,
    bool recursive,
    bool patternMatch
) const
{
    return lookupEntryBackwardsCompatible
    (
        keywords,
        recursive,
        patternMatch
    ).stream();
}


const Foam::entry* Foam::dictionary::lookupScopedEntryPtr
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    // '!' indicates the top-level directory
    if (keyword[0] == '!')
    {
        // Go up to top level
        const dictionary* dictPtr = this;
        while (&dictPtr->parent_ != &dictionary::null)
        {
            dictPtr = &dictPtr->parent_;
        }

        // At top. Recurse to find entries
        return dictPtr->lookupScopedSubEntryPtr
        (
            keyword.substr(1, keyword.size() - 1),
            false,
            patternMatch
        );
    }
    else
    {
        return lookupScopedSubEntryPtr
        (
            keyword,
            recursive,
            patternMatch
        );
    }
}


bool Foam::dictionary::substituteScopedKeyword(const word& keyword)
{
    word varName = keyword(1, keyword.size() - 1);

    // Lookup the variable name in the given dictionary
    const entry* ePtr = lookupScopedEntryPtr(varName, true, true);

    // If defined insert its entries into this dictionary
    if (ePtr != nullptr)
    {
        const dictionary& addDict = ePtr->dict();

        forAllConstIter(IDLList<entry>, addDict, iter)
        {
            add(iter());
        }

        return true;
    }

    return false;
}


bool Foam::dictionary::isDict(const word& keyword) const
{
    // Find non-recursive with patterns
    const entry* entryPtr = lookupEntryPtr(keyword, false, true);

    if (entryPtr)
    {
        return entryPtr->isDict();
    }
    else
    {
        return false;
    }
}


const Foam::dictionary* Foam::dictionary::subDictPtr(const word& keyword) const
{
    const entry* entryPtr = lookupEntryPtr(keyword, false, true);

    if (entryPtr)
    {
        return &entryPtr->dict();
    }
    else
    {
        return nullptr;
    }
}


Foam::dictionary* Foam::dictionary::subDictPtr(const word& keyword)
{
    entry* entryPtr = lookupEntryPtr(keyword, false, true);

    if (entryPtr)
    {
        return &entryPtr->dict();
    }
    else
    {
        return nullptr;
    }
}


const Foam::dictionary& Foam::dictionary::subDict(const word& keyword) const
{
    const entry* entryPtr = lookupEntryPtr(keyword, false, true);

    if (entryPtr == nullptr)
    {
        FatalIOErrorInFunction(*this)
            << "keyword " << keyword << " is undefined in dictionary "
            << name()
            << exit(FatalIOError);
    }
    return entryPtr->dict();
}


Foam::dictionary& Foam::dictionary::subDict(const word& keyword)
{
    entry* entryPtr = lookupEntryPtr(keyword, false, true);

    if (entryPtr == nullptr)
    {
        FatalIOErrorInFunction(*this)
            << "keyword " << keyword << " is undefined in dictionary "
            << name()
            << exit(FatalIOError);
    }
    return entryPtr->dict();
}


const Foam::dictionary& Foam::dictionary::subDictBackwardsCompatible
(
    const wordList& keywords
) const
{
    const entry* entryPtr =
        lookupEntryPtrBackwardsCompatible(keywords, false, true);

    if (entryPtr == nullptr)
    {
        // Generate error message using the first keyword
        return subDict(keywords[0]);
    }
    else
    {
        return entryPtr->dict();
    }
}


const Foam::dictionary& Foam::dictionary::subOrEmptyDict
(
    const word& keyword,
    const bool mustRead
) const
{
    const entry* entryPtr = lookupEntryPtr(keyword, false, true);

    if (entryPtr == nullptr)
    {
        if (mustRead)
        {
            FatalIOErrorInFunction(*this)
                << "keyword " << keyword << " is undefined in dictionary "
                << name()
                << exit(FatalIOError);
        }

        return null;
    }
    else
    {
        return entryPtr->dict();
    }
}


const Foam::dictionary& Foam::dictionary::optionalSubDict
(
    const word& keyword
) const
{
    const entry* entryPtr = lookupEntryPtr(keyword, false, true);

    if (entryPtr)
    {
        return entryPtr->dict();
    }
    else
    {
        return *this;
    }
}


const Foam::dictionary& Foam::dictionary::scopedDict(const word& keyword) const
{
    if (keyword == "")
    {
        return *this;
    }
    else
    {
        const entry* entPtr = lookupScopedEntryPtr
        (
            keyword,
            false,
            false
        );
        if (!entPtr || !entPtr->isDict())
        {
            FatalIOErrorInFunction(*this)
                << "keyword " << keyword
                << " is undefined in dictionary "
                << name() << " or is not a dictionary"
                << endl
                << "Valid keywords are " << keys()
                << exit(FatalIOError);
        }
        return entPtr->dict();
    }
}


Foam::dictionary& Foam::dictionary::scopedDict(const word& keyword)
{
    return const_cast<dictionary&>
    (
        const_cast<const dictionary*>(this)->scopedDict(keyword)
    );
}


Foam::wordList Foam::dictionary::toc() const
{
    wordList keys(size());

    label nKeys = 0;
    forAllConstIter(IDLList<entry>, *this, iter)
    {
        keys[nKeys++] = iter().keyword();
    }

    return keys;
}


Foam::wordList Foam::dictionary::sortedToc() const
{
    return hashedEntries_.sortedToc();
}


Foam::List<Foam::keyType> Foam::dictionary::keys(bool patterns) const
{
    List<keyType> keys(size());

    label nKeys = 0;
    forAllConstIter(IDLList<entry>, *this, iter)
    {
        if (iter().keyword().isPattern() ? patterns : !patterns)
        {
            keys[nKeys++] = iter().keyword();
        }
    }
    keys.setSize(nKeys);

    return keys;
}


bool Foam::dictionary::add(entry* entryPtr, bool mergeEntry)
{
    HashTable<entry*>::iterator iter = hashedEntries_.find
    (
        entryPtr->keyword()
    );

    if (mergeEntry && iter != hashedEntries_.end())
    {
        // Merge dictionary with dictionary
        if (iter()->isDict() && entryPtr->isDict())
        {
            iter()->dict().merge(entryPtr->dict());
            delete entryPtr;

            return true;
        }
        else
        {
            // Replace existing dictionary with entry or vice versa
            IDLList<entry>::replace(iter(), entryPtr);
            delete iter();
            hashedEntries_.erase(iter);

            if (hashedEntries_.insert(entryPtr->keyword(), entryPtr))
            {
                entryPtr->name() = name() + '/' + entryPtr->keyword();

                if (entryPtr->keyword().isPattern())
                {
                    patternEntries_.insert(entryPtr);
                    patternRegexps_.insert
                    (
                        autoPtr<regExp>(new regExp(entryPtr->keyword()))
                    );
                }

                return true;
            }
            else
            {
                IOWarningInFunction((*this))
                    << "problem replacing entry "<< entryPtr->keyword()
                    << " in dictionary " << name() << endl;

                IDLList<entry>::remove(entryPtr);
                delete entryPtr;
                return false;
            }
        }
    }

    if (hashedEntries_.insert(entryPtr->keyword(), entryPtr))
    {
        entryPtr->name() = name() + '/' + entryPtr->keyword();
        IDLList<entry>::append(entryPtr);

        if (entryPtr->keyword().isPattern())
        {
            patternEntries_.insert(entryPtr);
            patternRegexps_.insert
            (
                autoPtr<regExp>(new regExp(entryPtr->keyword()))
            );
        }

        return true;
    }
    else
    {
        // If function entries are disabled allow duplicate entries
        if (entry::disableFunctionEntries)
        {
            entryPtr->name() = name() + '/' + entryPtr->keyword();
            IDLList<entry>::append(entryPtr);

            return true;
        }
        else
        {
            IOWarningInFunction((*this))
                << "attempt to add entry "<< entryPtr->keyword()
                << " which already exists in dictionary " << name()
                << endl;

            delete entryPtr;
            return false;
        }
    }
}


void Foam::dictionary::add(const entry& e, bool mergeEntry)
{
    add(e.clone(*this).ptr(), mergeEntry);
}


void Foam::dictionary::add(const keyType& k, const word& w, bool overwrite)
{
    add(new primitiveEntry(k, token(w)), overwrite);
}


void Foam::dictionary::add
(
    const keyType& k,
    const Foam::string& s,
    bool overwrite
)
{
    add(new primitiveEntry(k, token(s)), overwrite);
}


void Foam::dictionary::add(const keyType& k, const label l, bool overwrite)
{
    add(new primitiveEntry(k, token(l)), overwrite);
}


void Foam::dictionary::add(const keyType& k, const scalar s, bool overwrite)
{
    add(new primitiveEntry(k, token(s)), overwrite);
}


void Foam::dictionary::add
(
    const keyType& k,
    const dictionary& d,
    bool mergeEntry
)
{
    add(new dictionaryEntry(k, *this, d), mergeEntry);
}


void Foam::dictionary::set(entry* entryPtr)
{
    entry* existingPtr = lookupEntryPtr(entryPtr->keyword(), false, true);

    // Clear dictionary so merge acts like overwrite
    if (existingPtr && existingPtr->isDict())
    {
        existingPtr->dict().clear();
    }
    add(entryPtr, true);
}


void Foam::dictionary::set(const entry& e)
{
    set(e.clone(*this).ptr());
}


void Foam::dictionary::set(const keyType& k, const dictionary& d)
{
    set(new dictionaryEntry(k, *this, d));
}


bool Foam::dictionary::remove(const word& Keyword)
{
    HashTable<entry*>::iterator iter = hashedEntries_.find(Keyword);

    if (iter != hashedEntries_.end())
    {
        // Delete from patterns first
        DLList<entry*>::iterator wcLink =
            patternEntries_.begin();
        DLList<autoPtr<regExp>>::iterator reLink =
            patternRegexps_.begin();

        // Find in pattern using exact match only
        if (findInPatterns(false, Keyword, wcLink, reLink))
        {
            patternEntries_.remove(wcLink);
            patternRegexps_.remove(reLink);
        }

        IDLList<entry>::remove(iter());
        delete iter();
        hashedEntries_.erase(iter);

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::dictionary::remove(const wordList& Keywords)
{
    forAll(Keywords, i)
    {
        remove(Keywords[i]);
    }
}


bool Foam::dictionary::changeKeyword
(
    const keyType& oldKeyword,
    const keyType& newKeyword,
    bool forceOverwrite
)
{
    // No change
    if (oldKeyword == newKeyword)
    {
        return false;
    }

    HashTable<entry*>::iterator iter = hashedEntries_.find(oldKeyword);

    // oldKeyword not found - do nothing
    if (iter == hashedEntries_.end())
    {
        return false;
    }

    if (iter()->keyword().isPattern())
    {
        FatalIOErrorInFunction(*this)
            << "Old keyword "<< oldKeyword
            << " is a pattern."
            << "Pattern replacement not yet implemented."
            << exit(FatalIOError);
    }


    HashTable<entry*>::iterator iter2 = hashedEntries_.find(newKeyword);

    // newKeyword already exists
    if (iter2 != hashedEntries_.end())
    {
        if (forceOverwrite)
        {
            if (iter2()->keyword().isPattern())
            {
                // Delete from patterns first
                DLList<entry*>::iterator wcLink =
                    patternEntries_.begin();
                DLList<autoPtr<regExp>>::iterator reLink =
                    patternRegexps_.begin();

                // Find in patterns using exact match only
                if (findInPatterns(false, iter2()->keyword(), wcLink, reLink))
                {
                    patternEntries_.remove(wcLink);
                    patternRegexps_.remove(reLink);
                }
            }

            IDLList<entry>::replace(iter2(), iter());
            delete iter2();
            hashedEntries_.erase(iter2);

        }
        else
        {
            IOWarningInFunction
            (
                *this
            )   << "cannot rename keyword "<< oldKeyword
                << " to existing keyword " << newKeyword
                << " in dictionary " << name() << endl;
            return false;
        }
    }

    // Change name and HashTable, but leave DL-List untouched
    iter()->keyword() = newKeyword;
    iter()->name() = name() + '/' + string::validate<word>(newKeyword);
    hashedEntries_.erase(oldKeyword);
    hashedEntries_.insert(newKeyword, iter());

    if (newKeyword.isPattern())
    {
        patternEntries_.insert(iter());
        patternRegexps_.insert
        (
            autoPtr<regExp>(new regExp(newKeyword))
        );
    }

    return true;
}


bool Foam::dictionary::merge(const dictionary& dict)
{
    // Check for assignment to self
    if (this == &dict)
    {
        FatalIOErrorInFunction(*this)
            << "attempted merge to self for dictionary " << name()
            << abort(FatalIOError);
    }

    bool changed = false;

    forAllConstIter(IDLList<entry>, dict, iter)
    {
        HashTable<entry*>::iterator fnd = hashedEntries_.find(iter().keyword());

        if (fnd != hashedEntries_.end())
        {
            // Recursively merge sub-dictionaries
            // TODO: merge without copying
            if (fnd()->isDict() && iter().isDict())
            {
                if (fnd()->dict().merge(iter().dict()))
                {
                    changed = true;
                }
            }
            else
            {
                add(iter().clone(*this).ptr(), true);
                changed = true;
            }
        }
        else
        {
            // Not found - just add
            add(iter().clone(*this).ptr());
            changed = true;
        }
    }

    return changed;
}


void Foam::dictionary::clear()
{
    IDLList<entry>::clear();
    hashedEntries_.clear();
    patternEntries_.clear();
    patternRegexps_.clear();
}


void Foam::dictionary::transfer(dictionary& dict)
{
    // Changing parents probably doesn't make much sense,
    // but what about the names?
    name() = dict.name();

    IDLList<entry>::transfer(dict);
    hashedEntries_.transfer(dict.hashedEntries_);
    patternEntries_.transfer(dict.patternEntries_);
    patternRegexps_.transfer(dict.patternRegexps_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::ITstream& Foam::dictionary::operator[](const word& keyword) const
{
    return lookup(keyword);
}


void Foam::dictionary::operator=(const dictionary& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalIOErrorInFunction(*this)
            << "attempted assignment to self for dictionary " << name()
            << abort(FatalIOError);
    }

    name() = rhs.name();
    clear();

    // Create clones of the entries in the given dictionary
    // resetting the parentDict to this dictionary

    forAllConstIter(IDLList<entry>, rhs, iter)
    {
        add(iter().clone(*this).ptr());
    }
}


void Foam::dictionary::operator+=(const dictionary& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalIOErrorInFunction(*this)
            << "attempted addition assignment to self for dictionary " << name()
            << abort(FatalIOError);
    }

    forAllConstIter(IDLList<entry>, rhs, iter)
    {
        add(iter().clone(*this).ptr());
    }
}


void Foam::dictionary::operator|=(const dictionary& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalIOErrorInFunction(*this)
            << "attempted assignment to self for dictionary " << name()
            << abort(FatalIOError);
    }

    forAllConstIter(IDLList<entry>, rhs, iter)
    {
        if (!found(iter().keyword()))
        {
            add(iter().clone(*this).ptr());
        }
    }
}


void Foam::dictionary::operator<<=(const dictionary& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalIOErrorInFunction(*this)
            << "attempted assignment to self for dictionary " << name()
            << abort(FatalIOError);
    }

    forAllConstIter(IDLList<entry>, rhs, iter)
    {
        set(iter().clone(*this).ptr());
    }
}


/* * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * */

Foam::dictionary Foam::operator+
(
    const dictionary& dict1,
    const dictionary& dict2
)
{
    dictionary sum(dict1);
    sum += dict2;
    return sum;
}


Foam::dictionary Foam::operator|
(
    const dictionary& dict1,
    const dictionary& dict2
)
{
    dictionary sum(dict1);
    sum |= dict2;
    return sum;
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::dictArgList
(
    const string& argString,
    word& funcName,
    wordReList& args,
    List<Tuple2<word, string>>& namedArgs
)
{
    funcName = argString;

    int argLevel = 0;
    bool namedArg = false;
    word argName;

    word::size_type start = 0;
    word::size_type i = 0;

    for
    (
        word::const_iterator iter = argString.begin();
        iter != argString.end();
        ++iter
    )
    {
        char c = *iter;

        if (c == '(')
        {
            if (argLevel == 0)
            {
                funcName = argString(start, i - start);
                start = i+1;
            }
            ++argLevel;
        }
        else if (c == ',' || c == ')')
        {
            if (argLevel == 1)
            {
                if (namedArg)
                {
                    namedArgs.append
                    (
                        Tuple2<word, string>
                        (
                            argName,
                            argString(start, i - start)
                        )
                    );
                    namedArg = false;
                }
                else
                {
                    args.append(wordRe(argString(start, i - start)));
                }
                start = i+1;
            }

            if (c == ')')
            {
                if (argLevel == 1)
                {
                    break;
                }
                --argLevel;
            }
        }
        else if (c == '=')
        {
            argName = argString(start, i - start);
            string::stripInvalid<variable>(argName);
            start = i+1;
            namedArg = true;
        }

        ++i;
    }

    // Strip whitespace from the function name
    string::stripInvalid<word>(funcName);
}


void Foam::dictArgList
(
    const string& argString,
    wordReList& args,
    List<Tuple2<word, string>>& namedArgs
)
{
    int argLevel = 0;
    bool namedArg = false;
    word argName;

    word::size_type start = 0;
    word::size_type i = 0;

    for
    (
        word::const_iterator iter = argString.begin();
        iter != argString.end();
        ++iter
    )
    {
        char c = *iter;

        if (c == '(')
        {
            ++argLevel;
        }
        else if (c == ',' || std::next(iter) == argString.end())
        {
            if (std::next(iter) == argString.end())
            {
                if (c == ')')
                {
                    --argLevel;
                }

                ++i;
            }

            if (argLevel == 0)
            {
                if (namedArg)
                {
                    namedArgs.append
                    (
                        Tuple2<word, string>
                        (
                            argName,
                            argString(start, i - start)
                        )
                    );
                    namedArg = false;
                }
                else
                {
                    args.append(wordRe(argString(start, i - start)));
                }
                start = i+1;
            }
        }
        else if (c == '=')
        {
            argName = argString(start, i - start);
            string::stripInvalid<variable>(argName);
            start = i+1;
            namedArg = true;
        }
        else if (c == ')')
        {
            --argLevel;
        }

        ++i;
    }
}


Foam::Pair<Foam::word> Foam::dictAndKeyword(const word& scopedName)
{
    string::size_type i = scopedName.find_last_of('/');

    if (i != string::npos)
    {
        return Pair<word>
        (
            scopedName.substr(0, i),
            scopedName.substr(i+1, string::npos)
        );
    }
    else
    {
        return Pair<word>("", scopedName);
    }
}


// ************************************************************************* //
