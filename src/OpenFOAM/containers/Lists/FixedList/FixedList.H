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

Class
    Foam::FixedList

Description
    A 1D vector of objects of type \<T\> with a fixed size \<Size\>.

SourceFiles
    FixedList.C
    FixedListI.H
    FixedListIO.C

\*---------------------------------------------------------------------------*/

#ifndef FixedList_H
#define FixedList_H

#include "bool.H"
#include "label.H"
#include "uLabel.H"
#include "Hash.H"
#include "autoPtr.H"
#include "StaticAssert.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class T, unsigned Size> class FixedList;

template<class T, unsigned Size>
Istream& operator>>(Istream&, FixedList<T, Size>&);

template<class T, unsigned Size>
Ostream& operator<<(Ostream&, const FixedList<T, Size>&);

template<class T> class UList;
template<class T> class SLList;


/*---------------------------------------------------------------------------*\
                           Class FixedList Declaration
\*---------------------------------------------------------------------------*/

template<class T, unsigned Size>
class FixedList
{
    //- Size must be positive (non-zero) and also fit as a signed value
    StaticAssert(Size && Size <= INT_MAX);

    // Private data

        //- Vector of values of type T of size Size.
        T v_[Size];


public:

    //- Hashing function class.
    //  Use Hasher directly for contiguous data. Otherwise hash incrementally.
    template< class HashT=Hash<T> >
    class Hash
    {
    public:
        Hash()
        {}

        inline unsigned operator()
        (
            const FixedList<T, Size>&,
            unsigned seed = 0
        ) const;
    };

    // Static Member Functions

        //- Return a null FixedList
        inline static const FixedList<T, Size>& null();


    // Constructors

        //- Null constructor.
        inline FixedList();

        //- Construct from C-array.
        explicit inline FixedList(const T v[Size]);

        //- Construct from value
        explicit inline FixedList(const T&);

        //- Construct from UList.
        explicit inline FixedList(const UList<T>&);

        //- Construct from SLList.
        explicit inline FixedList(const SLList<T>&);

        //- Copy constructor.
        inline FixedList(const FixedList<T, Size>&);

        //- Construct from Istream.
        FixedList(Istream&);

        //- Clone
        inline autoPtr< FixedList<T, Size> > clone() const;


    // Member Functions

        // Access

            //- Return the forward circular index, i.e. the next index
            //  which returns to the first at the end of the list
            inline label fcIndex(const label i) const;

            //- Return the reverse circular index, i.e. the previous index
            //  which returns to the last at the beginning of the list
            inline label rcIndex(const label i) const;


            //- Return a const pointer to the first data element,
            //  similar to the STL front() method and the string::data() method
            //  This can be used (with caution) when interfacing with C code.
            inline const T* cdata() const;

            //- Return a pointer to the first data element,
            //  similar to the STL front() method and the string::data() method
            //  This can be used (with caution) when interfacing with C code.
            inline T* data();

            //- Return the first element of the list.
            inline T& first();

            //- Return first element of the list.
            inline const T& first() const;

            //- Return the last element of the list.
            inline T& last();

            //- Return the last element of the list.
            inline const T& last() const;


        // Check

            //- Check start is within valid range (0 ... size-1).
            inline void checkStart(const label start) const;

            //- Check size is within valid range (0 ... size).
            inline void checkSize(const label size) const;

            //- Check index i is within valid range (0 ... size-1).
            inline void checkIndex(const label i) const;


        // Edit

            //- Dummy resize function
            //  needed to make FixedList consistent with List
            inline void resize(const label);

            //- Dummy setSize function
            //  needed to make FixedList consistent with List
            inline void setSize(const label);

            //- Copy (not transfer) the argument contents
            //  needed to make FixedList consistent with List
            void transfer(const FixedList<T, Size>&);

        //- Write the FixedList as a dictionary entry
        void writeEntry(Ostream&) const;

        //- Write the FixedList as a dictionary entry with keyword
        void writeEntry(const word& keyword, Ostream&) const;


    // Member operators

        //- Return element of FixedList.
        inline T& operator[](const label);

        //- Return element of constant FixedList.
        inline const T& operator[](const label) const;

        //- Assignment from array operator. Takes linear time.
        inline void operator=(const T v[Size]);

        //- Assignment from UList operator. Takes linear time.
        inline void operator=(const UList<T>&);

        //- Assignment from SLList operator. Takes linear time.
        inline void operator=(const SLList<T>&);

        //- Assignment of all entries to the given value
        inline void operator=(const T&);


    // STL type definitions

        //- Type of values the FixedList contains.
        typedef T value_type;

        //- Type that can be used for storing into
        //  FixedList::value_type objects.
        typedef T& reference;

        //- Type that can be used for storing into
        //  constant FixedList::value_type objects
        typedef const T& const_reference;

        //- The type that can represent the difference between any two
        //  FixedList iterator objects.
        typedef label difference_type;

        //- The type that can represent the size of a FixedList.
        typedef label size_type;


    // STL iterator

        //- Random access iterator for traversing FixedList.
        typedef T* iterator;

        //- Return an iterator to begin traversing the FixedList.
        inline iterator begin();

        //- Return an iterator to end traversing the FixedList.
        inline iterator end();


    // STL const_iterator

        //- Random access iterator for traversing FixedList.
        typedef const T* const_iterator;

        //- Return const_iterator to begin traversing the constant FixedList.
        inline const_iterator cbegin() const;

        //- Return const_iterator to end traversing the constant FixedList.
        inline const_iterator cend() const;

        //- Return const_iterator to begin traversing the constant FixedList.
        inline const_iterator begin() const;

        //- Return const_iterator to end traversing the constant FixedList.
        inline const_iterator end() const;


    // STL reverse_iterator

        //- Reverse iterator for reverse traversal of FixedList.
        typedef T* reverse_iterator;

        //- Return reverse_iterator to begin reverse traversing the FixedList.
        inline reverse_iterator rbegin();

        //- Return reverse_iterator to end reverse traversing the FixedList.
        inline reverse_iterator rend();


    // STL const_reverse_iterator

        //- Reverse iterator for reverse traversal of constant FixedList.
        typedef const T* const_reverse_iterator;

        //- Return const_reverse_iterator to begin reverse traversing FixedList.
        inline const_reverse_iterator crbegin() const;

        //- Return const_reverse_iterator to end reverse traversing FixedList.
        inline const_reverse_iterator crend() const;

        //- Return const_reverse_iterator to begin reverse traversing FixedList.
        inline const_reverse_iterator rbegin() const;

        //- Return const_reverse_iterator to end reverse traversing FixedList.
        inline const_reverse_iterator rend() const;


    // STL member functions

        //- Return the number of elements in the FixedList.
        inline label size() const;

        //- Return size of the largest possible FixedList.
        inline label max_size() const;

        //- Return true if the FixedList is empty (ie, size() is zero).
        inline bool empty() const;

        //- Swap two FixedLists of the same type in constant time.
        void swap(FixedList<T, Size>&);


    // STL member operators

        //- Equality operation on FixedLists of the same type.
        //  Returns true when the FixedLists are elementwise equal
        //  (using FixedList::value_type::operator==).  Takes linear time.
        bool operator==(const FixedList<T, Size>&) const;

        //- The opposite of the equality operation. Takes linear time.
        bool operator!=(const FixedList<T, Size>&) const;

        //- Compare two FixedLists lexicographically. Takes linear time.
        bool operator<(const FixedList<T, Size>&) const;

        //- Compare two FixedLists lexicographically. Takes linear time.
        bool operator>(const FixedList<T, Size>&) const;

        //- Return true if !(a > b). Takes linear time.
        bool operator<=(const FixedList<T, Size>&) const;

        //- Return true if !(a < b). Takes linear time.
        bool operator>=(const FixedList<T, Size>&) const;


    // IOstream operators

        //- Read List from Istream, discarding contents of existing List.
        friend Istream& operator>> <T, Size>
        (Istream&, FixedList<T, Size>&);

        // Write FixedList to Ostream.
        friend Ostream& operator<< <T, Size>
        (
            Ostream&,
            const FixedList<T, Size>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FixedListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "FixedList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
