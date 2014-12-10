/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "error.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Classes * * * * * * * * * * * * //

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline Foam::label
Foam::StaticHashTable<T, Key, Hash>::hashKeyIndex(const Key& key) const
{
    // size is power of two - this is the modulus
    return Hash()(key) & (keys_.size() - 1);
}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline Foam::label Foam::StaticHashTable<T, Key, Hash>::size() const
{
    return nElmts_;
}


template<class T, class Key, class Hash>
inline bool Foam::StaticHashTable<T, Key, Hash>::empty() const
{
    return !nElmts_;
}


template<class T, class Key, class Hash>
inline bool Foam::StaticHashTable<T, Key, Hash>::insert
(
    const Key& key,
    const T& newEntry
)
{
    return set(key, newEntry, true);
}


template<class T, class Key, class Hash>
inline bool Foam::StaticHashTable<T, Key, Hash>::set
(
    const Key& key,
    const T& newEntry
)
{
    return set(key, newEntry, false);
}


template<class T, class Key, class Hash>
inline Foam::Xfer<Foam::StaticHashTable<T, Key, Hash> >
Foam::StaticHashTable<T, Key, Hash>::xfer()
{
    return xferMove(*this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline T& Foam::StaticHashTable<T, Key, Hash>::operator[](const Key& key)
{
    iterator iter = find(key);

    if (iter == end())
    {
        FatalErrorIn("StaticHashTable<T, Key, Hash>::operator[](const Key&)")
            << key << " not found in table.  Valid entries: "
            << toc()
            << exit(FatalError);
    }

    return *iter;
}


template<class T, class Key, class Hash>
inline const T& Foam::StaticHashTable<T, Key, Hash>::operator[]
(
    const Key& key
) const
{
    const_iterator iter = find(key);

    if (iter == cend())
    {
        FatalErrorIn
        (
            "StaticHashTable<T, Key, Hash>::operator[](const Key&) const"
        )   << key << " not found in table.  Valid entries: "
            << toc()
            << exit(FatalError);
    }

    return *iter;
}


template<class T, class Key, class Hash>
inline T& Foam::StaticHashTable<T, Key, Hash>::operator()(const Key& key)
{
    iterator iter = find(key);

    if (iter == end())
    {
        insert(key, T());
        return *find(key);
    }
    else
    {
        return *iter;
    }
}


// * * * * * * * * * * * * * * * * STL iterator  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
template<class TRef, class TableRef>
inline Foam::StaticHashTable<T, Key, Hash>::Iterator<TRef, TableRef>::Iterator
(
    TableRef hashTbl,
    label hashIndex,
    label elemIndex
)
:
    hashTable_(hashTbl),
    hashIndex_(hashIndex),
    elemIndex_(elemIndex)
{}


template<class T, class Key, class Hash>
template<class TRef, class TableRef>
inline Foam::StaticHashTable<T, Key, Hash>::Iterator<TRef, TableRef>::Iterator
(
    const iterator& iter
)
:
    hashTable_(iter.hashTable_),
    hashIndex_(iter.hashIndex_),
    elemIndex_(iter.elemIndex_)
{}


template<class T, class Key, class Hash>
template<class TRef, class TableRef>
inline void
Foam::StaticHashTable<T, Key, Hash>::Iterator<TRef, TableRef>::operator=
(
    const iterator& iter
)
{
    this->hashIndex_ = iter.hashIndex_;
    this->elemIndex_ = iter.elemIndex_;
}


template<class T, class Key, class Hash>
template<class TRef, class TableRef>
inline bool
Foam::StaticHashTable<T, Key, Hash>::Iterator<TRef, TableRef>::operator==
(
    const iterator& iter
) const
{
    return hashIndex_ == iter.hashIndex_ && elemIndex_ == iter.elemIndex_;
}


template<class T, class Key, class Hash>
template<class TRef, class TableRef>
inline bool
Foam::StaticHashTable<T, Key, Hash>::Iterator<TRef, TableRef>::operator==
(
    const const_iterator& iter
) const
{
    return hashIndex_ == iter.hashIndex_ && elemIndex_ == iter.elemIndex_;
}


template<class T, class Key, class Hash>
template<class TRef, class TableRef>
inline bool
Foam::StaticHashTable<T, Key, Hash>::Iterator<TRef, TableRef>::operator!=
(
    const iterator& iter
) const
{
    return !operator==(iter);
}


template<class T, class Key, class Hash>
template<class TRef, class TableRef>
inline bool
Foam::StaticHashTable<T, Key, Hash>::Iterator<TRef, TableRef>::operator!=
(
    const const_iterator& iter
) const
{
    return !operator==(iter);
}


template<class T, class Key, class Hash>
template<class TRef, class TableRef>
inline TRef
Foam::StaticHashTable<T, Key, Hash>::Iterator<TRef, TableRef>::operator*()
{
    return hashTable_.objects_[hashIndex_][elemIndex_];
}


template<class T, class Key, class Hash>
template<class TRef, class TableRef>
inline TRef
Foam::StaticHashTable<T, Key, Hash>::Iterator<TRef, TableRef>::operator()()
{
    return operator*();
}


template<class T, class Key, class Hash>
template<class TRef, class TableRef>
inline
typename Foam::StaticHashTable<T, Key, Hash>::template Iterator
<
    TRef,
    TableRef
>&
Foam::StaticHashTable<T, Key, Hash>::Iterator
<
    TRef,
    TableRef
>::operator++()
{
    // A negative index is a special value from erase
    // (see notes in HashTable)
    if (hashIndex_ < 0)
    {
        hashIndex_ = -(hashIndex_+1) - 1;
    }
    else
    {
        // Try the next element on the local list
        elemIndex_++;

        if (elemIndex_ < hashTable_.objects_[hashIndex_].size())
        {
            return *this;
        }
    }

    // Step to the next table entry
    elemIndex_ = 0;

    while
    (
        ++hashIndex_ < hashTable_.objects_.size()
     && !hashTable_.objects_[hashIndex_].size()
    )
    {}


    if (hashIndex_ >= hashTable_.objects_.size())
    {
        // make end iterator
        hashIndex_ = hashTable_.keys_.size();
    }

    return *this;
}


template<class T, class Key, class Hash>
template<class TRef, class TableRef>
inline
typename Foam::StaticHashTable<T, Key, Hash>::template Iterator
<
    TRef,
    TableRef
>
Foam::StaticHashTable<T, Key, Hash>::Iterator
<
    TRef,
    TableRef
>::operator++
(
    int
)
{
    iterator tmp = *this;
    ++*this;
    return tmp;
}


template<class T, class Key, class Hash>
template<class TRef, class TableRef>
inline const Key&
Foam::StaticHashTable<T, Key, Hash>::Iterator<TRef, TableRef>::key() const
{
    return hashTable_.keys_[hashIndex_][elemIndex_];
}


template<class T, class Key, class Hash>
inline typename Foam::StaticHashTable<T, Key, Hash>::iterator
Foam::StaticHashTable<T, Key, Hash>::begin()
{
    // Find first non-empty entry
    forAll(keys_, hashIdx)
    {
        if (keys_[hashIdx].size())
        {
            return iterator(*this, hashIdx, 0);
        }
    }

#   ifdef FULLDEBUG
    if (debug)
    {
        Info<< "StaticHashTable is empty\n";
    }
#   endif

    return StaticHashTable<T, Key, Hash>::endIter_;
}


template<class T, class Key, class Hash>
inline const typename Foam::StaticHashTable<T, Key, Hash>::iterator&
Foam::StaticHashTable<T, Key, Hash>::end()
{
    return StaticHashTable<T, Key, Hash>::endIter_;
}


template<class T, class Key, class Hash>
inline typename Foam::StaticHashTable<T, Key, Hash>::const_iterator
Foam::StaticHashTable<T, Key, Hash>::cbegin() const
{
    // Find first non-empty entry
    forAll(keys_, hashIdx)
    {
        if (keys_[hashIdx].size())
        {
            return const_iterator(*this, hashIdx, 0);
        }
    }

#   ifdef FULLDEBUG
    if (debug)
    {
        Info<< "StaticHashTable is empty\n";
    }
#   endif

    return StaticHashTable<T, Key, Hash>::endConstIter_;
}


template<class T, class Key, class Hash>
inline const typename Foam::StaticHashTable<T, Key, Hash>::const_iterator&
Foam::StaticHashTable<T, Key, Hash>::cend() const
{
    return StaticHashTable<T, Key, Hash>::endConstIter_;
}


template<class T, class Key, class Hash>
inline typename Foam::StaticHashTable<T, Key, Hash>::const_iterator
Foam::StaticHashTable<T, Key, Hash>::begin() const
{
    return this->cbegin();
}


template<class T, class Key, class Hash>
inline const typename Foam::StaticHashTable<T, Key, Hash>::const_iterator&
Foam::StaticHashTable<T, Key, Hash>::end() const
{
    return StaticHashTable<T, Key, Hash>::endConstIter_;
}


// ************************************************************************* //
