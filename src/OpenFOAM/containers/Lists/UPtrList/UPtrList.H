/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::UPtrList

Description
    A templated 1D list of pointers to objects of type \<T\>, where the
    size of the array is known and used for subscript bounds checking, etc.

    The element operator [] returns a reference to the object rather than a
    pointer.  Storage is not allocated during construction or use but is
    supplied to the constructor as an argument.

SourceFiles
    UPtrList.C
    UPtrListIO.C

\*---------------------------------------------------------------------------*/

#ifndef UPtrList_H
#define UPtrList_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend classes
template<class T> class PtrList;

// Forward declaration of friend functions and operators
template<class T> class UPtrList;
template<class T> void writeEntry(Ostream& os, const UPtrList<T>&);
template<class T> Istream& operator>>(Istream&, UPtrList<T>&);
template<class T> Ostream& operator<<(Ostream&, const UPtrList<T>&);


/*---------------------------------------------------------------------------*\
                          Class UPtrList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class UPtrList
{
    // Private Data

        List<T*> ptrs_;


public:

    // Related types

        //- Declare friendship with the PtrList class
        friend class PtrList<T>;


    // Constructors

        //- Null Constructor
        UPtrList();

        //- Construct with size specified
        explicit UPtrList(const label);

        //- Construct as copy or reuse as specified
        UPtrList(UPtrList<T>&, bool reuse);

        //- Construct from an initialiser list of pointers
        UPtrList(std::initializer_list<T*>);


    // Member Functions

        // Access

            //- Return the number of elements in the UPtrList
            inline label size() const;

            //- Return true if the UPtrList is empty (ie, size() is zero)
            inline bool empty() const;

            //- Return reference to the first element of the list
            inline T& first();

            //- Return reference to first element of the list
            inline const T& first() const;

            //- Return reference to the last element of the list
            inline T& last();

            //- Return reference to the last element of the list
            inline const T& last() const;


        // Edit

            //- Reset size of UPtrList.  This can only be used to set the size
            //  of an empty UPtrList, extend a UPtrList, remove entries from
            //  the end of a UPtrList
            void setSize(const label);

            //- Reset size of UPtrList.  This can only be used to set the size
            //  of an empty UPtrList, extend a UPtrList, remove entries from
            //  the end of a UPtrList
            inline void resize(const label);

            //- Clear the UPtrList, i.e. set size to zero
            void clear();

            //- Append an element at the end of the list
            inline void append(T*);

            //- Transfer the contents of the argument UPtrList into this
            //  UPtrList and annul the argument list
            void transfer(UPtrList<T>&);

            //- Is element set
            inline bool set(const label) const;

            //- Set element. Return old element (can be nullptr).
            //  No checks on new element
            inline T* set(const label, T*);

            //- Reorders elements. Ordering does not have to be done in
            //  ascending or descending order. Reordering has to be unique.
            //  (is shuffle)
            void reorder(const labelUList& oldToNew);

            //- Reorders elements. Ordering does not have to be done in
            //  ascending or descending order. Reordering has to be unique.
            //  Note: can create unset elements
            void shuffle(const labelUList& newToOld);


        // Conversion

            //- Convert to list of different pointer type
            template<class T2>
            UPtrList<T2> convert();

            //- Convert to list of different pointer type
            template<class T2>
            UPtrList<const T2> convert() const;


    // Member Operators

        //- Return element const reference
        inline const T& operator[](const label) const;

        //- Return element reference
        inline T& operator[](const label);

        //- Return element const pointer
        inline const T* operator()(const label) const;

        //- Return element const pointer
        inline T* operator()(const label);


    // STL type definitions

        //- Type of values the UPtrList contains
        typedef T value_type;

        //- Type that can be used for storing into UPtrList::value_type objects
        typedef T& reference;

        //- Type that can be used for storing into constant UPtrList::value_type
        //  objects
        typedef const T& const_reference;


    // STL iterator
    // Random access iterator for traversing UPtrList

        class iterator;
        class const_iterator;
        friend class iterator;

        //- An STL iterator
        class iterator
        {
            T** ptr_;

        public:

            friend class const_iterator;

            //- Construct for a given UPtrList entry
            inline iterator(T**);

            // Member Operators

                inline bool operator==(const iterator&) const;
                inline bool operator!=(const iterator&) const;

                inline T& operator*();
                inline T& operator()();

                inline iterator operator++();
                inline iterator operator++(const int);

                inline iterator operator--();
                inline iterator operator--(const int);

                inline iterator operator+=(const label);
                inline iterator operator-=(const label);

                inline iterator operator+(const label) const;
                inline iterator operator-(const label) const;

                inline label operator-(const iterator&) const;

                inline T& operator[](const label);

                inline bool operator<(const iterator&) const;
                inline bool operator>(const iterator&) const;

                inline bool operator<=(const iterator&) const;
                inline bool operator>=(const iterator&) const;
        };

        //- Return an iterator to begin traversing the UPtrList
        inline iterator begin();

        //- Return an iterator to end traversing the UPtrList
        inline iterator end();


    // STL const_iterator
    // Random access iterator for traversing UPtrList

        //- An STL-conforming const_iterator
        class const_iterator
        {
            const T* const* ptr_;

        public:

            //- Construct for a given UPtrList entry
            inline const_iterator(const T* const*);

            //- Construct from an iterator
            inline const_iterator(const iterator&);


            // Member Operators

                inline bool operator==(const const_iterator&) const;
                inline bool operator!=(const const_iterator&) const;

                typedef const T& Tref;
                inline Tref operator*();
                inline Tref operator()();

                inline const_iterator operator++();
                inline const_iterator operator++(const int);

                inline const_iterator operator--();
                inline const_iterator operator--(const int);

                inline const_iterator operator+=(const label);
                inline const_iterator operator-=(const label);

                inline const_iterator operator+(const label) const;
                inline const_iterator operator-(const label) const;

                inline label operator-(const const_iterator&) const;

                inline const T& operator[](const label);

                inline bool operator<(const const_iterator&) const;
                inline bool operator>(const const_iterator&) const;

                inline bool operator<=(const const_iterator&) const;
                inline bool operator>=(const const_iterator&) const;
        };

        //- Return an const_iterator to begin traversing the UPtrList
        inline const_iterator cbegin() const;

        //- Return an const_iterator to end traversing the UPtrList
        inline const_iterator cend() const;

        //- Return an const_iterator to begin traversing the UPtrList
        inline const_iterator begin() const;

        //- Return an const_iterator to end traversing the UPtrList
        inline const_iterator end() const;


    // IOstream operator

        //- Write UPtrList to Ostream
        friend Ostream& operator<< <T>(Ostream&, const UPtrList<T>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "UPtrListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UPtrList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
