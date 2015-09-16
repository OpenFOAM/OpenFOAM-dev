/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::PtrList

Description
    A templated 1D list of pointers to objects of type \<T\>, where the
    size of the array is known and used for subscript bounds checking, etc.

    The element operator [] returns a reference to the object rather than a
    pointer.

SourceFiles
    PtrListI.H
    PtrList.C
    PtrListIO.C

\*---------------------------------------------------------------------------*/

#ifndef PtrList_H
#define PtrList_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class T> class PtrList;
template<class T> class SLPtrList;

template<class T>
inline typename PtrList<T>::iterator operator+
(
    const typename PtrList<T>::iterator&,
    label
);

template<class T>
inline typename PtrList<T>::iterator operator+
(
    label,
    const typename PtrList<T>::iterator&
);

template<class T>
inline typename PtrList<T>::iterator operator-
(
    const typename PtrList<T>::iterator&,
    label
);

template<class T>
inline label operator-
(
    const typename PtrList<T>::iterator&,
    const typename PtrList<T>::iterator&
);

template<class T>
inline typename PtrList<T>::const_iterator operator+
(
    const typename PtrList<T>::const_iterator&,
    label
);

template<class T>
inline typename PtrList<T>::const_iterator operator+
(
    label,
    const typename PtrList<T>::const_iterator&
);

template<class T>
inline typename PtrList<T>::const_iterator operator-
(
    const typename PtrList<T>::const_iterator&,
    label
);

template<class T>
inline label operator-
(
    const typename PtrList<T>::const_iterator&,
    const typename PtrList<T>::const_iterator&
);

template<class T>
Istream& operator>>(Istream&, PtrList<T>&);

template<class T>
Ostream& operator<<(Ostream&, const PtrList<T>&);

template<class T> class autoPtr;
template<class T> class tmp;


/*---------------------------------------------------------------------------*\
                           Class PtrList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class PtrList
{
    // Private data

        List<T*> ptrs_;


protected:

    // Protected Member Functions

        //- Read from Istream using given Istream constructor class
        template<class INew>
        void read(Istream&, const INew& inewt);


public:

    // Constructors

        //- Null Constructor.
        PtrList();

        //- Construct with size specified.
        explicit PtrList(const label);

        //- Copy constructor.
        PtrList(const PtrList<T>&);

        //- Copy constructor with additional argument for clone
        template<class CloneArg>
        PtrList(const PtrList<T>&, const CloneArg&);

        //- Construct by transferring the parameter contents
        PtrList(const Xfer<PtrList<T> >&);

        //- Construct as copy or re-use as specified.
        PtrList(PtrList<T>&, bool reUse);

        //- Construct as copy of SLPtrList<T>
        explicit PtrList(const SLPtrList<T>&);

        //- Construct from Istream using given Istream constructor class
        template<class INew>
        PtrList(Istream&, const INew&);

        //- Construct from Istream using default Istream constructor class
        PtrList(Istream&);


    //- Destructor
    ~PtrList();


    // Member functions

        // Access

            //- Return the number of elements in the PtrList
            inline label size() const;

            //- Return true if the PtrList is empty (ie, size() is zero).
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

            //- Reset size of PtrList. If extending the PtrList, new entries are
            //  set to NULL. If truncating the PtrList, removed entries are
            //  deleted.
            void setSize(const label);

            //- Alias for setSize(const label)
            inline void resize(const label);

            //- Clear the PtrList, i.e. set size to zero deleting all the
            //  allocated entries.
            void clear();

            //- Append an element at the end of the list
            inline void append(T*);
            inline void append(const autoPtr<T>&);
            inline void append(const tmp<T>&);

            //- Transfer the contents of the argument PtrList into this PtrList
            //  and annul the argument list.
            void transfer(PtrList<T>&);

            //- Transfer contents to the Xfer container
            inline Xfer<PtrList<T> > xfer();

            //- Is element set
            inline bool set(const label) const;

            //- Set element. Return old element (can be NULL).
            //  No checks on new element.
            inline autoPtr<T> set(const label, T*);
            inline autoPtr<T> set(const label, const autoPtr<T>&);
            inline autoPtr<T> set(const label, const tmp<T>&);

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


        //- Assignment.
        PtrList<T>& operator=(const PtrList<T>&);


    // STL type definitions

        //- Type of values the PtrList contains.
        typedef T value_type;

        //- Type that can be used for storing into PtrList::value_type objects.
        typedef T& reference;

        //- Type that can be used for storing into constant PtrList::value_type
        //  objects.
        typedef const T& const_reference;


    // STL iterator
    // Random access iterator for traversing PtrList.

        class iterator;
        class const_iterator;
        friend class iterator;
        friend class const_iterator;

        //- An STL-conforming iterator
        class iterator
        {
            T** ptr_;

        public:

            friend class const_iterator;

            //- Construct for a given PtrList entry
            inline iterator(T**);

            // Member operators

                inline bool operator==(const iterator&) const;
                inline bool operator!=(const iterator&) const;

                typedef T& Tref;
                inline Tref operator*();
                inline Tref operator()();

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

        //- Return an iterator to begin traversing the PtrList.
        inline iterator begin();

        //- Return an iterator to end traversing the PtrList.
        inline iterator end();


    // STL const_iterator
    // Random access iterator for traversing PtrList.

        //- An STL-conforming const_iterator
        class const_iterator
        {
            const T* const* ptr_;

        public:

            //- Construct for a given PtrList entry
            inline const_iterator(const T* const*);

            //- Construct from an iterator
            inline const_iterator(const iterator&);


            // Member operators

                inline bool operator==(const const_iterator&) const;
                inline bool operator!=(const const_iterator&) const;

                typedef const T& Tref;
                inline Tref operator*();
                inline Tref operator()();

                inline const_iterator operator++();
                inline const_iterator operator++(int);

                inline const_iterator operator--();
                inline const_iterator operator--(int);

                inline const_iterator operator+=(label);

                friend const_iterator operator+ <T>
                (
                    const const_iterator&,
                    label
                );
                friend const_iterator operator+ <T>
                (
                    label,
                    const const_iterator&
                );

                inline const_iterator operator-=(label);

                friend const_iterator operator- <T>
                (
                    const const_iterator&,
                    label
                );

                friend label operator- <T>
                (
                    const const_iterator&,
                    const const_iterator&
                );

                inline const T& operator[](label);

                inline bool operator<(const const_iterator&) const;
                inline bool operator>(const const_iterator&) const;

                inline bool operator<=(const const_iterator&) const;
                inline bool operator>=(const const_iterator&) const;
        };

        //- Return an const_iterator to begin traversing the PtrList.
        inline const_iterator cbegin() const;

        //- Return an const_iterator to end traversing the PtrList.
        inline const_iterator cend() const;

        //- Return an const_iterator to begin traversing the PtrList.
        inline const_iterator begin() const;

        //- Return an const_iterator to end traversing the PtrList.
        inline const_iterator end() const;


    // IOstream operator

        //- Read List from Istream, discarding contents of existing List.
        friend Istream& operator>> <T>(Istream&, PtrList<T>&);

        // Write List to Ostream.
        friend Ostream& operator<< <T>(Ostream&, const PtrList<T>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "PtrListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PtrList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
