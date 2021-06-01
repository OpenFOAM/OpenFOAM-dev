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
    Foam::UList

Description
    A 1D vector of objects of type \<T\>, where the size of the vector is
    known and can be used for subscript bounds checking, etc.

    Storage is not allocated during construction or use but is supplied to
    the constructor as an argument.  This type of list is particularly useful
    for lists that refer to parts of existing lists such as SubList.

SourceFiles
    UList.C
    UListI.H
    UListIO.C

\*---------------------------------------------------------------------------*/

#ifndef UList_H
#define UList_H

#include "bool.H"
#include "label.H"
#include "uLabel.H"
#include "nullObject.H"
#include "zero.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend classes
template<class T> class List;
template<class T> class SubList;

// Forward declaration of friend functions and operators
template<class T> class UList;
template<class T> void writeEntry(Ostream& os, const UList<T>&);
template<class T> Ostream& operator<<(Ostream&, const UList<T>&);
template<class T> Istream& operator>>(Istream&, UList<T>&);

typedef UList<label> labelUList;


/*---------------------------------------------------------------------------*\
                            Class UList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class UList
{
    // Private Data

        //- Number of elements in UList
        label size_;

        //- Vector of values of type T
        T* __restrict__ v_;


    // Private Member Functions

        //- Disallow default shallow-copy assignment
        //
        //  Assignment of UList<T> may need to be either shallow (copy pointer)
        //  or deep (copy elements) depending on context or the particular type
        //  of list derived from UList and it is confusing and prone to error
        //  for the default assignment to be either.  The solution is to
        //  disallow default assignment and provide separate 'shallowCopy' and
        //  'deepCopy' member functions
        void operator=(const UList<T>&) = delete;


public:

    // Related types

        //- Declare friendship with the List class
        friend class List<T>;

        //- Declare friendship with the SubList class
        friend class SubList<T>;


    // Static Member Functions

        //- Return a null UList
        inline static const UList<T>& null();


    // Public classes

        //- Less function class that can be used for sorting
        class less
        {
            const UList<T>& values_;

        public:

            less(const UList<T>& values)
            :
                values_(values)
            {}

            bool operator()(const label a, const label b)
            {
                return values_[a] < values_[b];
            }
        };

        //- Greater function class that can be used for sorting
        class greater
        {
            const UList<T>& values_;

        public:

            greater(const UList<T>& values)
            :
                values_(values)
            {}

            bool operator()(const label a, const label b)
            {
                return values_[a] > values_[b];
            }
        };


    // Constructors

        //- Null constructor
        inline UList();

        //- Construct from components
        inline UList(T* __restrict__ v, label size);


    // Member Functions

        // Access

            //- Return the forward circular index, i.e. the next index
            //  which returns to the first at the end of the list
            inline label fcIndex(const label i) const;

            //- Return the reverse circular index, i.e. the previous index
            //  which returns to the last at the beginning of the list
            inline label rcIndex(const label i) const;

            //- Return the binary size in number of characters of the UList
            //  if the element is a primitive type
            //  i.e. contiguous<T>() == true.
            //  Note that is of type streamsize since used in stream ops
            std::streamsize byteSize() const;


            //- Return a const pointer to the first data element,
            //  similar to the STL front() method and the string::data() method
            //  This can be used (with caution) when interfacing with C code
            inline const T* cdata() const;

            //- Return a pointer to the first data element,
            //  similar to the STL front() method and the string::data() method
            //  This can be used (with caution) when interfacing with C code
            inline T* data();

            //- Return the first element of the list
            inline T& first();

            //- Return first element of the list
            inline const T& first() const;

            //- Return the last element of the list
            inline T& last();

            //- Return the last element of the list
            inline const T& last() const;


        // Check

            //- Check start is within valid range (0 ... size-1)
            inline void checkStart(const label start) const;

            //- Check size is within valid range (0 ... size)
            inline void checkSize(const label size) const;

            //- Check index i is within valid range (0 ... size-1)
            inline void checkIndex(const label i) const;


        // Edit

            //- Copy the pointer held by the given UList
            inline void shallowCopy(const UList<T>&);

            //- Copy elements of the given UList
            void deepCopy(const UList<T>&);


    // Member Operators

        //- Return element of UList
        inline T& operator[](const label);

        //- Return element of constant UList
        //  Note that the bool specialisation adds lazy evaluation so reading
        //  an out-of-range element returns false without any ill-effects
        inline const T& operator[](const label) const;

        //- Allow cast to a const List<T>&
        inline operator const Foam::List<T>&() const;

        //- Assignment of all entries to the given value
        void operator=(const T&);

        //- Assignment of all entries to zero
        void operator=(const zero);


    // STL type definitions

        //- Type of values the UList contains
        typedef T value_type;

        //- Type that can be used for storing into
        //  UList::value_type objects
        typedef T& reference;

        //- Type that can be used for storing into
        //  constant UList::value_type objects
        typedef const T& const_reference;

        //- The type that can represent the difference between any two
        //  UList iterator objects
        typedef label difference_type;

        //- The type that can represent the size of a UList
        typedef label size_type;


    // STL iterator

        //- Random access iterator for traversing UList
        typedef T* iterator;

        //- Return an iterator to begin traversing the UList
        inline iterator begin();

        //- Return an iterator to end traversing the UList
        inline iterator end();


    // STL const_iterator

        //- Random access iterator for traversing UList
        typedef const T* const_iterator;

        //- Return const_iterator to begin traversing the constant UList
        inline const_iterator cbegin() const;

        //- Return const_iterator to end traversing the constant UList
        inline const_iterator cend() const;

        //- Return const_iterator to begin traversing the constant UList
        inline const_iterator begin() const;

        //- Return const_iterator to end traversing the constant UList
        inline const_iterator end() const;


    // STL reverse_iterator

        //- Reverse iterator for reverse traversal of UList
        typedef std::reverse_iterator<iterator> reverse_iterator;

        //- Return reverse_iterator to begin reverse traversing the UList
        inline reverse_iterator rbegin();

        //- Return reverse_iterator to end reverse traversing the UList
        inline reverse_iterator rend();


    // STL const_reverse_iterator

        //- Reverse iterator for reverse traversal of constant UList
        typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

        //- Return const_reverse_iterator to begin reverse traversing the UList
        inline const_reverse_iterator crbegin() const;

        //- Return const_reverse_iterator to end reverse traversing the UList
        inline const_reverse_iterator crend() const;

        //- Return const_reverse_iterator to begin reverse traversing the UList
        inline const_reverse_iterator rbegin() const;

        //- Return const_reverse_iterator to end reverse traversing the UList
        inline const_reverse_iterator rend() const;


    // STL member functions

        //- Return the number of elements in the UList
        inline label size() const;

        //- Return size of the largest possible UList
        inline label max_size() const;

        //- Return true if the UList is empty (ie, size() is zero)
        inline bool empty() const;

        //- Swap two ULists of the same type in constant time
        void swap(UList<T>&);


    // STL member operators

        //- Equality operation on ULists of the same type.
        //  Returns true when the ULists are element-wise equal
        //  (using UList::value_type::operator==).  Takes linear time
        bool operator==(const UList<T>&) const;

        //- The opposite of the equality operation. Takes linear time
        bool operator!=(const UList<T>&) const;

        //- Compare two ULists lexicographically. Takes linear time
        bool operator<(const UList<T>&) const;

        //- Compare two ULists lexicographically. Takes linear time
        bool operator>(const UList<T>&) const;

        //- Return true if !(a > b). Takes linear time
        bool operator<=(const UList<T>&) const;

        //- Return true if !(a < b). Takes linear time
        bool operator>=(const UList<T>&) const;


    // Ostream operator

        // Write UList to Ostream
        friend Ostream& operator<< <T>
        (
            Ostream&,
            const UList<T>&
        );

        //- Read UList contents from Istream. Requires size to have been set
        //  before
        friend Istream& operator>> <T>
        (
            Istream&,
            UList<T>&
        );
};

template<class T>
void sort(UList<T>&);

template<class T, class Cmp>
void sort(UList<T>&, const Cmp&);

template<class T>
void stableSort(UList<T>&);

template<class T, class Cmp>
void stableSort(UList<T>&, const Cmp&);

template<class T>
void shuffle(UList<T>&);

// Reverse the first n elements of the list
template<class T>
inline void reverse(UList<T>&, const label n);

// Reverse all the elements of the list
template<class T>
inline void reverse(UList<T>&);

template<class ListType>
void writeListEntry(Ostream& os, const ListType& l);

template<class ListType>
void writeListEntries(Ostream& os, const ListType& l);

template<class ListType>
void writeListEntries(Ostream& os, const word& keyword, const ListType& l);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "UListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Loop across all elements in \a list
// \par Usage
// \code
// forAll(anyList, i)
// {
//      statements;
// }
// \endcode
// \sa forAllReverse
#define forAll(list, i) \
    for (Foam::label i=0; i<(list).size(); i++)

//- Reverse loop across all elements in \a list
//  \par Usage
//  \code
//  forAllReverse(anyList, i)
//  {
//       statements;
//  }
//  \endcode
//  \sa forAll
#define forAllReverse(list, i) \
    for (Foam::label i=(list).size()-1; i>=0; i--)

//- Iterate across all elements in the \a container object of type
//  \a Container.
//  \par Usage
//  \code
//  forAll(ContainerType, container, iter)
//  {
//      statements;
//  }
//  \endcode
//  \sa forAllConstIter
#define forAllIter(Container,container,iter)                                   \
    for                                                                        \
    (                                                                          \
        Container::iterator iter = (container).begin();                        \
        iter != (container).end();                                             \
        ++iter                                                                 \
    )

//- Iterate across all elements in the \a container object of type
//  \a Container with const access.
//  \par Usage
//  \code
//  forAllConstIter(ContainerType, container, iter)
//  {
//      statements;
//  }
//  \endcode
//  \sa forAllIter
#define forAllConstIter(Container,container,iter)                              \
    for                                                                        \
    (                                                                          \
        Container::const_iterator iter = (container).begin();                  \
        iter != (container).end();                                             \
        ++iter                                                                 \
    )


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
