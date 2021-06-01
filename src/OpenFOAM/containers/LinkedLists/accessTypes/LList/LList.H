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
    Foam::LList

Description
    Template class for non-intrusive linked lists.

SourceFiles
    LList.C
    LListIO.C

\*---------------------------------------------------------------------------*/

#ifndef LList_H
#define LList_H

#include "label.H"
#include <initializer_list>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

// Forward declaration of friend functions and operators

template<class LListBase, class T> class LList;

template<class LListBase, class T>
Istream& operator>>
(
    Istream&,
    LList<LListBase, T>&
);

template<class LListBase, class T>
Ostream& operator<<
(
    Ostream&,
    const LList<LListBase, T>&
);


/*---------------------------------------------------------------------------*\
                            Class LList Declaration
\*---------------------------------------------------------------------------*/

template<class LListBase, class T>
class LList
:
    public LListBase
{

public:

    // Forward declaration of STL iterators

        class iterator;
        friend class iterator;

        class const_iterator;
        friend class const_iterator;


    //- Link structure
    struct link
    :
        public LListBase::link
    {
        //- Stored object
        T obj_;

        //- Construct given object
        link(T a)
        :
            obj_(a)
        {}
    };


    // Constructors

        //- Null construct
        LList()
        {}

        //- Construct given initial T
        LList(T a)
        :
            LListBase(new link(a))
        {}

        //- Construct from Istream
        LList(Istream&);

        //- Copy constructor
        LList(const LList<LListBase, T>&);

        //- Move constructor
        LList(LList<LListBase, T>&&);

        //- Construct from an initialiser list
        LList(std::initializer_list<T>);


    //- Destructor
    ~LList();


    // Member Functions

        // Access

            //- Return the first entry added
            T& first()
            {
                return static_cast<link*>(LListBase::first())->obj_;
            }

            //- Return const access to the first entry added
            const T& first() const
            {
                return static_cast<const link*>(LListBase::first())->obj_;
            }

            //- Return the last entry added
            T& last()
            {
                return static_cast<link*>(LListBase::last())->obj_;
            }

            //- Return const access to the last entry added
            const T& last() const
            {
                return static_cast<const link*>(LListBase::last())->obj_;
            }


        // Edit

            //- Add at head of list
            void insert(const T& a)
            {
                LListBase::insert(new link(a));
            }

            //- Add at tail of list
            void append(const T& a)
            {
                LListBase::append(new link(a));
            }

            //- Remove and return head
            T removeHead()
            {
                link* elmtPtr = static_cast<link*>(LListBase::removeHead());
                T data = elmtPtr->obj_;
                delete elmtPtr;
                return data;
            }

            //- Remove and return element
            T remove(link* l)
            {
                link* elmtPtr = static_cast<link*>(LListBase::remove(l));
                T data = elmtPtr->obj_;
                delete elmtPtr;
                return data;
            }

            //- Remove and return element specified by iterator
            T remove(iterator& it)
            {
                link* elmtPtr = static_cast<link*>(LListBase::remove(it));
                T data = elmtPtr->obj_;
                delete elmtPtr;
                return data;
            }

            //- Delete contents of list
            void clear();

            //- Transfer the contents of the argument into this List
            //  and annul the argument list.
            void transfer(LList<LListBase, T>&);


    // Member Operators

        //- Assignment operator
        void operator=(const LList<LListBase, T>&);

        //- Move assignment operator
        void operator=(LList<LListBase, T>&&);

        //- Assignment to an initialiser list
        void operator=(std::initializer_list<T>);


    // STL type definitions

        //- Type of values the LList contains.
        typedef T value_type;

        //- Type that can be used for storing into value_type
        //  objects.
        typedef T& reference;

        //- Type that can be used for storing into constant
        //  LList::value_type objects.
        typedef const T& const_reference;

        //- The type that can represent the size of a LList.
        typedef label size_type;


    // STL iterator

        typedef typename LListBase::iterator LListBase_iterator;

        //- An STL-conforming iterator
        class iterator
        :
            public LListBase_iterator
        {

        public:

            //- Construct from base iterator
            iterator(LListBase_iterator baseIter)
            :
                LListBase_iterator(baseIter)
            {}


            // Member Operators

                T& operator*()
                {
                    return
                        static_cast<link&>
                        (LListBase_iterator::operator*()).obj_;
                }

                T& operator()()
                {
                    return operator*();
                }

                iterator& operator++()
                {
                    LListBase_iterator::operator++();
                    return *this;
                }
        };

        inline iterator begin()
        {
            return LListBase::begin();
        }

        inline const iterator& end()
        {
            return static_cast<const iterator&>(LListBase::end());
        }


    // STL const_iterator

        typedef typename LListBase::const_iterator LListBase_const_iterator;

        //- An STL-conforming const_iterator
        class const_iterator
        :
            public LListBase_const_iterator
        {

        public:

            //- Construct from base const_iterator
            const_iterator(LListBase_const_iterator baseIter)
            :
                LListBase_const_iterator(baseIter)
            {}


            //- Construct from base iterator
            const_iterator(LListBase_iterator baseIter)
            :
                LListBase_const_iterator(baseIter)
            {}


            // Member Operators

                const T& operator*()
                {
                    return
                        static_cast<const link&>
                        (LListBase_const_iterator::operator*()).obj_;
                }

                const T& operator()()
                {
                    return operator*();
                }

                const_iterator& operator++()
                {
                    LListBase_const_iterator::operator++();
                    return *this;
                }
        };

        inline const_iterator cbegin() const
        {
            return LListBase::cbegin();
        }

        inline const const_iterator& cend() const
        {
            return static_cast<const const_iterator&>(LListBase::cend());
        }

        inline const_iterator begin() const
        {
            return LListBase::begin();
        }

        inline const const_iterator& end() const
        {
            return static_cast<const const_iterator&>(LListBase::end());
        }


    // IOstream Operators

        friend Istream& operator>> <LListBase, T>
        (
            Istream&,
            LList<LListBase, T>&
        );

        friend Ostream& operator<< <LListBase, T>
        (
            Ostream&,
            const LList<LListBase, T>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
