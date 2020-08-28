/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::ListHashTable

Description
    STL conforming hash table using contiguous lists rather than linked lists.

Note
    Is slower to insert than the standard HashTable, but should be more
    memory efficient and faster to access.

SourceFiles
    ListHashTableI.H
    ListHashTable.C
    ListHashTableIO.C

\*---------------------------------------------------------------------------*/

#ifndef ListHashTable_H
#define ListHashTable_H

#include "label.H"
#include "uLabel.H"
#include "word.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class T> class List;
template<class T, class Key, class Hash> class ListHashTable;

template<class T, class Key, class Hash> Istream& operator>>
(
    Istream&,
    ListHashTable<T, Key, Hash>&
);

template<class T, class Key, class Hash> Ostream& operator<<
(
    Ostream&,
    const ListHashTable<T, Key, Hash>&
);


/*---------------------------------------------------------------------------*\
                      Class ListHashTableCore Declaration
\*---------------------------------------------------------------------------*/

//- Template-invariant bits for ListHashTable
struct ListHashTableCore
{
    //- Return a canonical (power-of-two) size
    static label canonicalSize(const label);

    //- Construct null
    ListHashTableCore()
    {}

    //- Define template name and debug
    ClassName("ListHashTable");

    //- A zero-sized end iterator
    struct iteratorEnd
    {
        //- Construct null
        iteratorEnd()
        {}
    };
};



/*---------------------------------------------------------------------------*\
                        Class ListHashTable Declaration
\*---------------------------------------------------------------------------*/

template<class T, class Key=word, class Hash=string::hash>
class ListHashTable
:
    public ListHashTableCore
{
    // Private Data type for table entries

        //- The lookup keys, ordered per hash value
        List<List<Key>> keys_;

        //- For each key the corresponding object.
        List<List<T>> objects_;

        //- The current number of elements in table
        label nElmts_;

        //- Return a canonical (power-of-two) size
        static label canonicalSize(const label);

        //- Return the hash index of the Key within the current table size.
        //  No checks for zero-sized tables.
        inline label hashKeyIndex(const Key&) const;

        //- Assign a new hashed entry to a possibly already existing key
        bool set(const Key&, const T& newElmt, bool protect);


public:


    // Forward declaration of STL iterators

        template<class TRef, class TableRef>
        class Iterator;

        typedef Iterator
        <
            T&,
            ListHashTable<T, Key, Hash>&
        > iterator;

        typedef Iterator
        <
            const T&,
            const ListHashTable<T, Key, Hash>&
        > const_iterator;


    // Declare friendship with the iterators

        friend class Iterator
        <
            T&,
            ListHashTable<T, Key, Hash>&
        >;

        friend class Iterator
        <
            const T&,
            const ListHashTable<T, Key, Hash>&
        >;


    // Constructors

        //- Construct given initial table size
        ListHashTable(const label size = 128);

        //- Construct from Istream
        ListHashTable(Istream&, const label size = 128);

        //- Copy constructor
        ListHashTable(const ListHashTable<T, Key, Hash>&);

        //- Move constructor
        ListHashTable(ListHashTable<T, Key, Hash>&&);


    //- Destructor
    ~ListHashTable();


    // Member Functions

        // Access

            //- Return number of elements in table.
            inline label size() const;

            //- Return true if the hash table is empty
            inline bool empty() const;

            //- Return true if hashed entry is found in table
            bool found(const Key& key) const;

            //- Find and return an iterator set at the hashed entry
            //  If not found iterator = end()
            iterator find(const Key& key);

            //- Find and return an const_iterator set at the hashed entry
            //  If not found iterator = end()
            const_iterator find(const Key& key) const;

            //- Return the table of contents
            List<Key> toc() const;

            //- Print information
            Ostream& printInfo(Ostream&) const;


        // Edit

            //- Insert a new hashed entry
            bool insert(const Key& key, const T& newElmt);

            //- Assign a new hashed entry, overwriting existing entries
            inline bool set(const Key&, const T& newElmt);

            //- Erase an hashed entry specified by given iterator
            bool erase(const iterator& it);

            //- Erase an hashed entry specified by given key if in table
            bool erase(const Key& key);

            //- Resize the hash table for efficiency
            void resize(const label newSize);

            //- Remove entries in the given hash table from this hash table
            //  Return the number of elements removed
            label erase(const ListHashTable<T, Key, Hash>&);

            //- Clear all entries from table
            void clear();

            //- Clear the table entries and the table itself.
            //  Equivalent to clear() followed by resize(1)
            void clearStorage();

            //- Transfer the contents of the argument table into this table
            //  and annul the argument table.
            void transfer(ListHashTable<T, Key, Hash>&);


    // Member Operators

        //- Find and return an hashed entry
        inline T& operator[](const Key&);

        //- Find and return an hashed entry
        inline const T& operator[](const Key&) const;

        //- Find and return an hashed entry, create it null if not present.
        inline T& operator()(const Key&);

        //- Assignment operator
        void operator=(const ListHashTable<T, Key, Hash>&);

        //- Move assignment operator
        void operator=(ListHashTable<T, Key, Hash>&&);

        //- Equality. Two hash tables are equal if all contents of first are
        //  also in second and vice versa.
        bool operator==(const ListHashTable<T, Key, Hash>&) const;

        //- The opposite of the equality operation.
        bool operator!=(const ListHashTable<T, Key, Hash>&) const;


    // STL type definitions

        //- Type of values the ListHashTable contains.
        typedef T value_type;

        //- Type that can be used for storing into ListHashTable::value_type
        //  objects.  This type is usually List::value_type&.
        typedef T& reference;

        //- Type that can be used for storing into constant
        //  ListHashTable::value_type objects.  This type is usually const
        //  ListHashTable::value_type&.
        typedef const T& const_reference;

        //- The type that can represent the size of a ListHashTable.
        typedef label size_type;


    // STL iterator

        //- An STL iterator
        template<class TRef, class TableRef>
        class Iterator
        {
            friend class ListHashTable;

            template<class TRef2, class TableRef2>
            friend class Iterator;

            // Private Data

                //- Reference to the ListHashTable this is an iterator for
                TableRef hashTable_;

                //- Current hash index
                label hashIndex_;

                //- Index of current element at hashIndex
                label elemIndex_;


        public:

            // Constructors

                //- Construct from hash table, hash index and element index
                inline Iterator
                (
                    TableRef,
                    label hashIndex_,
                    label elemIndex_
                );

                //- Construct from the non-const iterator
                inline Iterator(const iterator&);


            // Member Operators

                inline void operator=(const iterator&);

                inline bool operator==(const iterator&) const;
                inline bool operator==(const const_iterator&) const;

                inline bool operator!=(const iterator&) const;
                inline bool operator!=(const const_iterator&) const;

                inline TRef operator*();
                inline TRef operator()();

                inline Iterator& operator++();
                inline Iterator operator++(int);

                inline const Key& key() const;
        };


        //- Iterator set to the beginning of the ListHashTable
        inline iterator begin();

        //- Iterator set to beyond the end of the ListHashTable
        inline const iterator& end();

        //- const_iterator set to the beginning of the ListHashTable
        inline const_iterator cbegin() const;

        //- const_iterator set to beyond the end of the ListHashTable
        inline const const_iterator& cend() const;

        //- const_iterator set to the beginning of the ListHashTable
        inline const_iterator begin() const;

        //- const_iterator set to beyond the end of the ListHashTable
        inline const const_iterator& end() const;

    // IOstream Operator

        friend Istream& operator>> <T, Key, Hash>
        (
            Istream&,
            ListHashTable<T, Key, Hash>&
        );

        friend Ostream& operator<< <T, Key, Hash>
        (
            Ostream&,
            const ListHashTable<T, Key, Hash>&
        );


private:

        //- Iterator returned by end()
        iterator endIter_;

        //- const_iterator returned by end()
        const_iterator endConstIter_;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "ListHashTableI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef NoListHashTableC
#ifdef NoRepository
    #include "ListHashTable.C"
#endif
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
