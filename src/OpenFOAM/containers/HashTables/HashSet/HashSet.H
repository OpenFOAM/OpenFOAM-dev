/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

Class
    Foam::HashSet

Description
    A HashTable with keys but without contents.

Typedef
    Foam::wordHashSet

Description
    A HashSet with (the default) word keys.

Typedef
    Foam::labelHashSet

Description
    A HashSet with label keys.

\*---------------------------------------------------------------------------*/

#ifndef HashSet_H
#define HashSet_H

#include "HashTable.H"
#include "nil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class HashSet Declaration
\*---------------------------------------------------------------------------*/

template<class Key=word, class Hash=string::hash>
class HashSet
:
    public HashTable<nil, Key, Hash>
{

public:

    typedef typename HashTable<nil, Key, Hash>::iterator iterator;
    typedef typename HashTable<nil, Key, Hash>::const_iterator const_iterator;


    // Constructors

        //- Construct given initial size
        HashSet(const label size = 128)
        :
            HashTable<nil, Key, Hash>(size)
        {}

        //- Construct from Istream
        HashSet(Istream& is)
        :
            HashTable<nil, Key, Hash>(is)
        {}

        //- Construct from UList of Key
        HashSet(const UList<Key>&);

        //- Construct from FixedList of Key
        template<unsigned Size>
        HashSet(const FixedList<Key, Size>&);

        //- Copy constructor
        HashSet(const HashSet<Key, Hash>& hs) = default;

        //- Move constructor
        HashSet(HashSet<Key, Hash>&& hs) = default;

        //- Construct from the keys of another HashTable,
        //  the type of values held is arbitrary.
        template<class AnyType, class AnyHash>
        HashSet(const HashTable<AnyType, Key, AnyHash>&);

        //- Construct from an initialiser list
        HashSet(std::initializer_list<Key>);


    // Member Functions

        // Edit

        //- Insert a new entry
        bool insert(const Key& key)
        {
            return HashTable<nil, Key, Hash>::insert(key, nil());
        }

        //- Insert keys from a UList of Key
        //  Return the number of new elements inserted
        label insert(const UList<Key>&);

        //- Same as insert (cannot overwrite nil content)
        bool set(const Key& key)
        {
            return insert(key);
        }

        //- Same as insert (cannot overwrite nil content)
        label set(const UList<Key>& lst)
        {
            return insert(lst);
        }

        //- Unset the specified key - same as erase
        bool unset(const Key& key)
        {
            return HashTable<nil, Key, Hash>::erase(key);
        }


    // Member Operators

        //- Return true if the entry exists, same as found()
        inline bool operator[](const Key&) const;

        //- Assignment operator
        void operator=(const HashSet<Key, Hash>&);

        //- Move assignment operator
        void operator=(HashSet<Key, Hash>&&);

        //- Equality. Two hashtables are equal when their contents are equal.
        //  Independent of table size or order.
        bool operator==(const HashSet<Key, Hash>&) const;

        //- The opposite of the equality operation.
        bool operator!=(const HashSet<Key, Hash>&) const;


        //- Combine entries from HashSets
        void operator|=(const HashSet<Key, Hash>&);

        //- Only retain entries found in both HashSets
        void operator&=(const HashSet<Key, Hash>&);

        //- Only retain unique entries (xor)
        void operator^=(const HashSet<Key, Hash>&);

        //- Add entries listed in the given HashSet to this HashSet
        inline void operator+=(const HashSet<Key, Hash>& rhs)
        {
            this->operator|=(rhs);
        }

        //- Remove entries listed in the given HashSet from this HashSet
        void operator-=(const HashSet<Key, Hash>&);
};


// Global Operators

//- Combine entries from HashSets
template<class Key, class Hash>
HashSet<Key,Hash> operator|
(
    const HashSet<Key,Hash>& hash1,
    const HashSet<Key,Hash>& hash2
);


//- Create a HashSet that only contains entries found in both HashSets
template<class Key, class Hash>
HashSet<Key,Hash> operator&
(
    const HashSet<Key,Hash>& hash1,
    const HashSet<Key,Hash>& hash2
);


//- Create a HashSet that only contains unique entries (xor)
template<class Key, class Hash>
HashSet<Key,Hash> operator^
(
    const HashSet<Key,Hash>& hash1,
    const HashSet<Key,Hash>& hash2
);


//- A HashSet with word keys.
typedef HashSet<> wordHashSet;

//- A HashSet with label keys.
typedef HashSet<label, Hash<label>> labelHashSet;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HashSet.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
