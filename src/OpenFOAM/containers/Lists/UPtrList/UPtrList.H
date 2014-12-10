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

// Forward declaration of friend functions and operators

template<class T> class UPtrList;

template<class T>
inline typename UPtrList<T>::iterator operator+
(
    const typename UPtrList<T>::iterator&,
    label
);

template<class T>
inline typename UPtrList<T>::iterator operator+
(
    label,
    const typename UPtrList<T>::iterator&
);

template<class T>
inline typename UPtrList<T>::iterator operator-
(
    const typename UPtrList<T>::iterator&,
    label
);

template<class T>
inline label operator-
(
    const typename UPtrList<T>::iterator&,
    const typename UPtrList<T>::iterator&
);

template<class T>
Istream& operator>>(Istream&, UPtrList<T>&);

template<class T>
Ostream& operator<<(Ostream&, const UPtrList<T>&);


/*---------------------------------------------------------------------------*\
                           Class UPtrList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class UPtrList
{
    // Private data

        List<T*> ptrs_;


public:

    // Constructors

        //- Null Constructor.
        UPtrList();

        //- Construct with size specified.
        explicit UPtrList(const label);

        //- Construct by transferring the parameter contents
        UPtrList(const Xfer<UPtrList<T> >&);

        //- Construct as copy or re-use as specified.
        UPtrList(UPtrList<T>&, bool reUse);


    // Member functions

        // Access

            //- Return the number of elements in the UPtrList
            inline label size() const;

            //- Return true if the UPtrList is empty (ie, size() is zero).
            inline bool empty() const;

            //- Return reference to the first element of the list.
            inline T& first();

            //- Return reference to first element of the list.
            inline const T& first() const;

            //- Return reference to the last element of the list.
            inline T& last();

            //- Return reference to the last element of the list.
            inline const T& last() const;


        // Edit

            //- Reset size of UPtrList.  This can only be used to set the size
            //  of an empty UPtrList, extend a UPtrList, remove entries from
            //  the end of a UPtrList.
            void setSize(const label);

            //- Reset size of UPtrList.  This can only be used to set the size
            //  of an empty UPtrList, extend a UPtrList, remove entries from
            //  the end of a UPtrList.
            inline void resize(const label);

            //- Clear the UPtrList, i.e. set size to zero
            void clear();

            //- Transfer the contents of the argument UPtrList into this
            //  UPtrList and annul the argument list.
            void transfer(UPtrList<T>&);

            //- Transfer contents to the Xfer container
            inline Xfer<UPtrList<T> > xfer();

            //- Is element set
            inline bool set(const label) const;

            //- Set element. Return old element (can be NULL).
            //  No checks on new element.
            inline T* set(const label, T*);

            //- Reorders elements. Ordering does not have to be done in
            //  ascending or descending order. Reordering has to be unique.
            //  (is shuffle)
            void reorder(const labelUList&);


    // Member operators

        //- Return element const reference.
        inline const T& operator[](const label) const;

        //- Return element reference.
        inline T& operator[](const label);

        //- Return element const pointer.
        inline const T* operator()(const label) const;


    // STL type definitions

        //- Type of values the UPtrList contains.
        typedef T value_type;

        //- Type that can be used for storing into UPtrList::value_type objects.
        typedef T& reference;

        //- Type that can be used for storing into constant UPtrList::value_type
        //  objects.
        typedef const T& const_reference;


    // STL iterator
    // Random access iterator for traversing UPtrList.

        class iterator;
        friend class iterator;

        //- An STL iterator
        class iterator
        {
            T** ptr_;

        public:

            //- Construct for a given UPtrList entry
            inline iterator(T**);

            // Member operators

                inline bool operator==(const iterator&) const;
                inline bool operator!=(const iterator&) const;

                inline T& operator*();
                inline T& operator()();

                inline iterator operator++();
                inline iterator operator++(int);

                inline iterator operator--();
                inline iterator operator--(int);

                inline iterator operator+=(label);

                friend iterator operator+ <T>(const iterator&, label);
                friend iterator operator+ <T>(label, const iterator&);

                inline iterator operator-=(label);

                friend iterator operator- <T>(const iterator&, label);

                friend label operator- <T>
                (
                    const iterator&,
                    const iterator&
                );

                inline T& operator[](label);

                inline bool operator<(const iterator&) const;
                inline bool operator>(const iterator&) const;

                inline bool operator<=(const iterator&) const;
                inline bool operator>=(const iterator&) const;
        };

        //- Return an iterator to begin traversing the UPtrList.
        inline iterator begin();

        //- Return an iterator to end traversing the UPtrList.
        inline iterator end();


    // IOstream operator

        // Write List to Ostream.
        friend Ostream& operator<< <T>(Ostream&, const UPtrList<T>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "UPtrListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "UPtrList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
