/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::UILList

Description
    Template class for intrusive linked lists.

SourceFiles
    UILList.C
    UILListIO.C

\*---------------------------------------------------------------------------*/

#ifndef UILList_H
#define UILList_H

#include "label.H"
#include "uLabel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Ostream;

// Forward declaration of friend functions and operators

template<class LListBase, class T>
class UILList;

template<class LListBase, class T>
Ostream& operator<<
(
    Ostream&,
    const UILList<LListBase, T>&
);


/*---------------------------------------------------------------------------*\
                           Class UILList Declaration
\*---------------------------------------------------------------------------*/

template<class LListBase, class T>
class UILList
:
    public LListBase
{

public:

    // Forward declaration of STL iterators

        class iterator;
        friend class iterator;

        class const_iterator;
        friend class const_iterator;

        class const_reverse_iterator;
        friend class const_reverse_iterator;


    // Constructors

        //- Null construct
        UILList()
        {}

        //- Construct given initial T
        UILList(T* a)
        :
            LListBase(a)
        {}

        //- Construct as copy
        UILList(const UILList<LListBase, T>&);


    // Member Functions

        // Access

            //- Return the first entry
            T* first()
            {
                return static_cast<T*>(LListBase::first());
            }

            //- Return the first entry
            const T* first() const
            {
                return static_cast<const T*>(LListBase::first());
            }

            //- Return the last entry
            T* last()
            {
                return static_cast<T*>(LListBase::last());
            }

            //- Return the last entry
            const T* last() const
            {
                return static_cast<const T*>(LListBase::last());
            }


        // Edit

            //- Remove and return head
            T* removeHead()
            {
                return static_cast<T*>(LListBase::removeHead());
            }

            //- Remove and return element
            T* remove(T* p)
            {
                return static_cast<T*>(LListBase::remove(p));
            }

            //- Remove and return specified by iterator
            T* remove(iterator& it)
            {
                return static_cast<T*>(LListBase::remove(it));
            }


    // Member operators

        void operator=(const UILList<LListBase, T>&);


    // STL type definitions

        //- Type of values the DLList contains.
        typedef T value_type;

        //- Type that can be used for storing into DLList::value_type
        //  objects.
        typedef T& reference;

        //- Type that can be used for storing into constant
        //  DLList::value_type objects.
        typedef const T& const_reference;

        //- The type that can represent the size of a DLList.
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


            // Member operators

                T& operator*()
                {
                    return static_cast<T&>(LListBase_iterator::operator*());
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


            // Member operators

                const T& operator*()
                {
                    return
                        static_cast<const T&>
                        (LListBase_const_iterator::operator*());
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


    // STL const_reverse_iterator

        //- An STL-conforming const_reverse_iterator
        class const_reverse_iterator
        :
            public LListBase::const_reverse_iterator
        {

        public:

            //- Construct from base const_reverse_iterator
            const_reverse_iterator
            (
                typename LListBase::const_reverse_iterator baseIter
            )
            :
                LListBase::const_reverse_iterator(baseIter)
            {}


            // Member operators

                const T& operator*()
                {
                    return
                        static_cast<const T&>
                        (LListBase::const_reverse_iterator::operator*());
                }

                const T& operator()()
                {
                    return operator*();
                }

                const_reverse_iterator& operator++()
                {
                    LListBase::const_reverse_iterator::operator++();
                    return *this;
                }
        };

        inline const_reverse_iterator crbegin() const
        {
            return LListBase::crbegin();
        }

        inline const const_reverse_iterator& crend() const
        {
            return
                static_cast<const const_reverse_iterator&>(LListBase::crend());
        }


        inline const_reverse_iterator rbegin() const
        {
            return LListBase::rbegin();
        }

        inline const const_reverse_iterator& rend() const
        {
            return
                static_cast<const const_reverse_iterator&>(LListBase::rend());
        }


    // STL member operators

        //- Equality operation on ULists of the same type.
        //  Returns true when the ULists are element-wise equal
        //  (using UList::value_type::operator==).  Takes linear time.
        bool operator==(const UILList<LListBase, T>&) const;

        //- The opposite of the equality operation. Takes linear time.
        bool operator!=(const UILList<LListBase, T>&) const;


    // Ostream operator

        friend Ostream& operator<< <LListBase, T>
        (
            Ostream&,
            const UILList<LListBase, T>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UILList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
