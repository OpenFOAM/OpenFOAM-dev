/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "DictionaryBase.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class IDLListType, class T>
void Foam::DictionaryBase<IDLListType, T>::addEntries()
{
    for
    (
        typename IDLListType::iterator iter = this->begin();
        iter != this->end();
        ++iter
    )
    {
        this->hashedTs_.insert((*iter).keyword(), &(*iter));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class IDLListType, class T>
Foam::DictionaryBase<IDLListType, T>::DictionaryBase(const label size)
:
    hashedTs_(size)
{}


template<class IDLListType, class T>
Foam::DictionaryBase<IDLListType, T>::DictionaryBase
(
    const DictionaryBase& dict
)
:
    IDLListType(dict)
{
    addEntries();
}


template<class IDLListType, class T>
template<class INew>
Foam::DictionaryBase<IDLListType, T>::DictionaryBase
(
    Istream& is,
    const INew& iNew
)
:
    IDLListType(is, iNew)
{
    addEntries();
}


template<class IDLListType, class T>
Foam::DictionaryBase<IDLListType, T>::DictionaryBase(Istream& is)
:
    IDLListType(is)
{
    addEntries();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Find and return T
template<class IDLListType, class T>
bool Foam::DictionaryBase<IDLListType, T>::found(const word& keyword) const
{
    return hashedTs_.found(keyword);
}


// Find and return T*, return NULL if not found
template<class IDLListType, class T>
const T* Foam::DictionaryBase<IDLListType, T>::lookupPtr
(
    const word& keyword
) const
{
    typename HashTable<T*>::const_iterator iter = hashedTs_.find(keyword);

    if (iter != hashedTs_.end())
    {
        return *iter;
    }
    else
    {
        return NULL;
    }
}


// Find and return T*, return NULL if not found
template<class IDLListType, class T>
T* Foam::DictionaryBase<IDLListType, T>::lookupPtr(const word& keyword)
{
    typename HashTable<T*>::iterator iter = hashedTs_.find(keyword);

    if (iter != hashedTs_.end())
    {
        return *iter;
    }
    else
    {
        return NULL;
    }
}


// Find and return T*, FatalError if keyword not found
template<class IDLListType, class T>
const T* Foam::DictionaryBase<IDLListType, T>::lookup(const word& keyword) const
{
    typename HashTable<T*>::const_iterator iter = hashedTs_.find(keyword);

    if (iter == hashedTs_.end())
    {
        FatalErrorIn
        (
            "DictionaryBase<IDLListType, T>::lookup(const word&) const"
        )   << keyword << " is undefined"
            << exit(FatalError);
    }

    return *iter;
}


// Find and return T*, FatalError if keyword not found
template<class IDLListType, class T>
T* Foam::DictionaryBase<IDLListType, T>::lookup(const word& keyword)
{
    typename HashTable<T*>::iterator iter = hashedTs_.find(keyword);

    if (iter == hashedTs_.end())
    {
        FatalErrorIn
        (
            "DictionaryBase<IDLListType, T>::lookup(const word&)"
        )   << keyword << " is undefined"
            << exit(FatalError);
    }

    return *iter;
}


// Return the table of contents
template<class IDLListType, class T>
Foam::wordList Foam::DictionaryBase<IDLListType, T>::toc() const
{
    wordList keywords(this->size());

    label i = 0;
    for
    (
        typename IDLListType::const_iterator iter = this->begin();
        iter != this->end();
        ++iter
    )
    {
        keywords[i++] = iter().keyword();
    }

    return keywords;
}


// Add at head of dictionary
template<class IDLListType, class T>
void Foam::DictionaryBase<IDLListType, T>::insert(const word& keyword, T* tPtr)
{
    // NOTE: we should probably check that HashTable::insert actually worked
    hashedTs_.insert(keyword, tPtr);
    IDLListType::insert(tPtr);
}


// Add at tail of dictionary
template<class IDLListType, class T>
void Foam::DictionaryBase<IDLListType, T>::append(const word& keyword, T* tPtr)
{
    // NOTE: we should probably check that HashTable::insert actually worked
    hashedTs_.insert(keyword, tPtr);
    IDLListType::append(tPtr);
}


template<class IDLListType, class T>
T* Foam::DictionaryBase<IDLListType, T>::remove(const word& keyword)
{
    typename HashTable<T*>::iterator iter = hashedTs_.find(keyword);

    if (iter != hashedTs_.end())
    {
        T* tPtr = IDLListType::remove(iter());
        hashedTs_.erase(iter);
        return tPtr;
    }
    else
    {
        return NULL;
    }
}


template<class IDLListType, class T>
void Foam::DictionaryBase<IDLListType, T>::clear()
{
    IDLListType::clear();
    hashedTs_.clear();
}


template<class IDLListType, class T>
void Foam::DictionaryBase<IDLListType, T>::transfer
(
    DictionaryBase<IDLListType, T>& dict
)
{
    IDLListType::transfer(dict);
    hashedTs_.transfer(dict.hashedTs_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class IDLListType, class T>
void Foam::DictionaryBase<IDLListType, T>::operator=
(
    const DictionaryBase<IDLListType, T>& dict
)
{
    // Check for assignment to self
    if (this == &dict)
    {
        FatalErrorIn("DictionaryBase::operator=(const DictionaryBase&)")
            << "attempted assignment to self"
            << abort(FatalError);
    }

    IDLListType::operator=(dict);
    this->hashedTs_.clear();
    this->addEntries();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DictionaryBaseIO.C"

// ************************************************************************* //
