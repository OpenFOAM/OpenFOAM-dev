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
    Foam::HashTable

Description
    An STL-conforming hash table.

Note
    Hashing index collisions are handled via chaining using a singly-linked
    list with the colliding entry being added to the head of the linked
    list. Thus copying the hash table (or indeed even resizing it) will
    often result in a different hash order. Use a sorted table-of-contents
    when the hash order is important.

SourceFiles
    HashTableI.H
    HashTable.C
    HashTableIO.C

\*---------------------------------------------------------------------------*/

#ifndef HashTable_H
#define HashTable_H

#include "label.H"
#include "uLabel.H"
#include "word.H"
#include "className.H"
#include <initializer_list>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class T> class List;
template<class T> class UList;
template<class T, class Key, class Hash> class HashTable;
template<class T, class Key, class Hash> class HashPtrTable;

template<class Type1, class Type2>
class Tuple2;

template<class T, class Key, class Hash>
void writeEntry(Ostream& os, const HashTable<T, Key, Hash>& ht);

template<class T, class Key, class Hash>
Istream& operator>>(Istream&, HashTable<T, Key, Hash>&);

template<class T, class Key, class Hash>
Ostream& operator<<(Ostream&, const HashTable<T, Key, Hash>&);


/*---------------------------------------------------------------------------*\
                        Class HashTableCore Declaration
\*---------------------------------------------------------------------------*/

//- Template-invariant bits for HashTable
struct HashTableCore
{
    //- Return a canonical (power-of-two) size
    static label canonicalSize(const label);

    //- Maximum allowable table size
    static const label maxTableSize;

    //- Construct null
    HashTableCore()
    {}

    //- Define template name and debug
    ClassName("HashTable");

    //- A zero-sized end iterator
    struct iteratorEnd
    {
        //- Construct null
        iteratorEnd()
        {}
    };

    //- iteratorEnd set to beyond the end of any HashTable
    inline static iteratorEnd cend()
    {
        return iteratorEnd();
    }

    //- iteratorEnd set to beyond the end of any HashTable
    inline static iteratorEnd end()
    {
        return iteratorEnd();
    }
};


/*---------------------------------------------------------------------------*\
                          Class HashTable Declaration
\*---------------------------------------------------------------------------*/

template<class T, class Key=word, class Hash=string::hash>
class HashTable
:
    public HashTableCore
{
    // Private Data type for table entries

        //- Structure to hold a hashed entry with SLList for collisions
        struct hashedEntry
        {
            //- The lookup key
            Key key_;

            //- Pointer to next hashedEntry in sub-list
            hashedEntry* next_;

            //- The data object
            T obj_;

            //- Construct from key, next pointer and object
            inline hashedEntry(const Key&, hashedEntry* next, const T&);

            //- Disallow default bitwise copy construction
            hashedEntry(const hashedEntry&) = delete;

            //- Disallow default bitwise assignment
            void operator=(const hashedEntry&) = delete;
        };


    // Private Data: size of table, the table and current number of elements

        //- The current number of elements in table
        label nElmts_;

        //- Number of primary entries allocated in table
        label tableSize_;

        //- The table of primary entries
        hashedEntry** table_;


    // Private Member Functions

        //- Return a canonical (power-of-two) size
        static label canonicalSize(const label);

        //- Return the hash index of the Key within the current table size.
        //  No checks for zero-sized tables.
        inline label hashKeyIndex(const Key&) const;

        //- Assign a new hashedEntry to a possibly already existing key
        bool set(const Key&, const T& newElmt, bool protect);


public:

    // Forward declaration of iterators

        class iteratorBase;
        class iterator;
        class const_iterator;

        //- Declare friendship with the HashPtrTable class
        template<class T2, class Key2, class Hash2>
        friend class HashPtrTable;

        //- Declare friendship with the iteratorBase
        friend class iteratorBase;

        //- Declare friendship with the iterator
        friend class iterator;

        //- Declare friendship with the const_iterator
        friend class const_iterator;


    // Constructors

        //- Construct given initial table size
        HashTable(const label size = 128);

        //- Construct from Istream
        HashTable(Istream&, const label size = 128);

        //- Copy constructor
        HashTable(const HashTable<T, Key, Hash>&);

        //- Move constructor
        HashTable(HashTable<T, Key, Hash>&&);

        //- Construct from an initialiser list
        HashTable(std::initializer_list<Tuple2<Key, T>>);


    //- Destructor
    ~HashTable();


    // Member Functions

        // Access

            //- The size of the underlying table
            inline label capacity() const;

            //- Return number of elements in table
            inline label size() const;

            //- Return true if the hash table is empty
            inline bool empty() const;

            //- Return true if hashedEntry is found in table
            bool found(const Key&) const;

            //- Find and return an iterator set at the hashedEntry
            //  If not found iterator = end()
            iterator find(const Key&);

            //- Find and return an const_iterator set at the hashedEntry
            //  If not found iterator = end()
            const_iterator find(const Key&) const;

            //- Return the table of contents
            List<Key> toc() const;

            //- Return the table of contents as a sorted list
            List<Key> sortedToc() const;

            //- Print information
            Ostream& printInfo(Ostream&) const;


        // Edit

            //- Insert a new hashedEntry
            inline bool insert(const Key&, const T& newElmt);

            //- Assign a new hashedEntry, overwriting existing entries
            inline bool set(const Key&, const T& newElmt);

            //- Erase a hashedEntry specified by given iterator
            //  This invalidates the iterator until the next operator++
            bool erase(const iterator&);

            //- Erase a hashedEntry specified by the given key
            bool erase(const Key&);

            //- Remove entries given by the listed keys from this HashTable
            //  Return the number of elements removed
            label erase(const UList<Key>&);

            //- Remove entries given by the given keys from this HashTable
            //  Return the number of elements removed.
            //  The parameter HashTable needs the same type of key, but the
            //  type of values held and the hashing function are arbitrary.
            template<class AnyType, class AnyHash>
            label erase(const HashTable<AnyType, Key, AnyHash>&);

            //- Resize the hash table for efficiency
            void resize(const label newSize);

            //- Clear all entries from table
            void clear();

            //- Clear the table entries and the table itself.
            //  Equivalent to clear() followed by resize(0)
            void clearStorage();

            //- Shrink the allocated table to approx. twice number of elements
            void shrink();

            //- Transfer the contents of the argument table into this table
            //  and annul the argument table.
            void transfer(HashTable<T, Key, Hash>&);


    // Member Operators

        //- Find and return a hashedEntry
        inline T& operator[](const Key&);

        //- Find and return a hashedEntry
        inline const T& operator[](const Key&) const;

        //- Find and return a hashedEntry, create it null if not present
        inline T& operator()(const Key&);

        //- Assignment operator
        void operator=(const HashTable<T, Key, Hash>&);

        //- Move assignment operator
        void operator=(HashTable<T, Key, Hash>&&);

        //- Assignment to an initialiser list
        void operator=(std::initializer_list<Tuple2<Key, T>>);

        //- Equality. Hash tables are equal if the keys and values are equal.
        //  Independent of table storage size and table order.
        bool operator==(const HashTable<T, Key, Hash>&) const;

        //- The opposite of the equality operation. Takes linear time.
        bool operator!=(const HashTable<T, Key, Hash>&) const;



    // STL type definitions

        //- Type of values the HashTable contains.
        typedef T value_type;

        //- Type that can be used for storing into HashTable::value_type
        //  objects.  This type is usually List::value_type&.
        typedef T& reference;

        //- Type that can be used for storing into constant
        //  HashTable::value_type objects.  This type is usually const
        //  HashTable::value_type&.
        typedef const T& const_reference;

        //- The type that can represent the size of a HashTable.
        typedef label size_type;


    // Iterators and helpers

        //- The iterator base for HashTable
        //  Note: data and functions are protected, to allow reuse by iterator
        //  and prevent most external usage.
        //  iterator and const_iterator have the same size, allowing
        //  us to reinterpret_cast between them (if desired)
        class iteratorBase
        {
            // Private Data

                //- Pointer to the HashTable for which this is an iterator
                //  This also lets us use the default bitwise copy/assignment
                HashTable<T, Key, Hash>* hashTable_;

                //- Current element
                hashedEntry* entryPtr_;

                //- Current hash index
                label hashIndex_;


        protected:

            // Constructors

                //- Construct null - equivalent to an 'end' position
                inline iteratorBase();

                //- Construct from hash table, moving to its 'begin' position
                inline explicit iteratorBase
                (
                    const HashTable<T, Key, Hash>* curHashTable
                );

                //- Construct from hash table, element and hash index
                inline iteratorBase
                (
                    const HashTable<T, Key, Hash>* curHashTable,
                    const hashedEntry* elmt,
                    const label hashIndex
                );


            // Protected Member Functions

                //- Increment to the next position
                inline void increment();

                //- Erase the HashTable element at the current position
                bool erase();

                //- Return non-const access to referenced object
                inline T& object();

                //- Return const access to referenced object
                inline const T& cobject() const;


        public:

            // Member Operators

                // Access

                //- Return the Key corresponding to the iterator
                inline const Key& key() const;

                //- Compare hashedEntry element pointers
                inline bool operator==(const iteratorBase&) const;
                inline bool operator!=(const iteratorBase&) const;

                //- Compare hashedEntry to iteratorEnd pointers
                inline bool operator==(const iteratorEnd& unused) const;
                inline bool operator!=(const iteratorEnd& unused) const;
        };


        //- An STL-conforming iterator
        class iterator
        :
            public iteratorBase
        {
            friend class HashTable;

            // Private Member Functions

                //- Construct from hash table, moving to its 'begin' position
                inline explicit iterator
                (
                    HashTable<T, Key, Hash>* curHashTable
                );

                //- Construct from hash table, element and hash index
                inline iterator
                (
                    HashTable<T, Key, Hash>* curHashTable,
                    hashedEntry* elmt,
                    const label hashIndex
                );


        public:

            // Constructors

                //- Construct null (end iterator)
                inline iterator();

                //- Construct end iterator
                inline iterator(const iteratorEnd& unused);


            // Member Operators

                //- Return referenced hash value
                inline T& operator*();
                inline T& operator()();

                //- Return referenced hash value
                inline const T& operator*() const;
                inline const T& operator()() const;

                inline iterator& operator++();
                inline iterator operator++(int);
        };

        //- Iterator set to the beginning of the HashTable
        inline iterator begin();


    // STL const_iterator

        //- An STL-conforming const_iterator
        class const_iterator
        :
            public iteratorBase
        {
            friend class HashTable;

            // Private Member Functions

                //- Construct from hash table, moving to its 'begin' position
                inline explicit const_iterator
                (
                    const HashTable<T, Key, Hash>* curHashTable
                );

                //- Construct from hash table, element and hash index
                inline const_iterator
                (
                    const HashTable<T, Key, Hash>* curHashTable,
                    const hashedEntry* elmt,
                    const label hashIndex
                );


        public:

            // Constructors

                //- Construct null (end iterator)
                inline const_iterator();

                //- Construct from iterator
                inline const_iterator(const iterator&);

                //- Construct end iterator
                inline const_iterator(const iteratorEnd& unused);


            // Member Operators

                //- Return referenced hash value
                inline const T& operator*() const;
                inline const T& operator()() const;

                inline const_iterator& operator++();
                inline const_iterator operator++(int);

        };


        //- const_iterator set to the beginning of the HashTable
        inline const_iterator cbegin() const;

        //- const_iterator set to the beginning of the HashTable
        inline const_iterator begin() const;


    // IOstream Operator

        friend Istream& operator>> <T, Key, Hash>
        (
            Istream&,
            HashTable<T, Key, Hash>&
        );

        friend Ostream& operator<< <T, Key, Hash>
        (
            Ostream&,
            const HashTable<T, Key, Hash>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "HashTableI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef NoHashTableC
#ifdef NoRepository
    #include "HashTable.C"
#endif
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
