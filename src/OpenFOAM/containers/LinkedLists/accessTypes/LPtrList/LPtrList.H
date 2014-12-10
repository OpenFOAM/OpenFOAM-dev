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
    Foam::LPtrList

Description
    Template class for non-intrusive linked PtrLists.

SourceFiles
    LPtrList.C
    LPtrListIO.C

\*---------------------------------------------------------------------------*/

#ifndef LPtrList_H
#define LPtrList_H

#include "LList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class LListBase, class T> class LPtrList;

template<class LListBase, class T>
Istream& operator>>
(
    Istream&,
    LPtrList<LListBase, T>&
);

template<class LListBase, class T>
Ostream& operator<<
(
    Ostream&,
    const LPtrList<LListBase, T>&
);


/*---------------------------------------------------------------------------*\
                          Class LPtrList Declaration
\*---------------------------------------------------------------------------*/

template<class LListBase, class T>
class LPtrList
:
    public LList<LListBase, T*>
{
    // Private Member Functions

        //- Read from Istream using given Istream constructor class
        template<class INew>
        void read(Istream&, const INew&);


public:

    // Forward declaration of STL iterators

        class iterator;
        friend class iterator;

        class const_iterator;
        friend class const_iterator;


    // Constructors

        //- Null construct
        LPtrList()
        {}

        //- Construct given initial T
        LPtrList(T* a)
        :
            LList<LListBase, T*>(a)
        {}

        //- Construct from Istream using given Istream constructor class
        template<class INew>
        LPtrList(Istream&, const INew&);

        //- Construct from Istream using default Istream constructor class
        LPtrList(Istream&);

        //- Construct as copy
        LPtrList(const LPtrList&);


    //- Destructor
    ~LPtrList();


    // Member Functions

        // Access

            //- Return the first entry added
            T& first()
            {
                return *LList<LListBase, T*>::first();
            }

            //- Return const access to the first entry added
            const T& first() const
            {
                return *LList<LListBase, T*>::first();
            }

            //- Return the last entry added
            T& last()
            {
                return *LList<LListBase, T*>::last();
            }

            //- Return const access to the last entry added
            const T& last() const
            {
                return *LList<LListBase, T*>::last();
            }


        // Edit

            //- Remove the head element from the list and delete the pointer
            bool eraseHead();

            //- Clear the contents of the list
            void clear();

            //- Transfer the contents of the argument into this List
            //  and annul the argument list.
            void transfer(LPtrList<LListBase, T>&);


    // Member operators

        //- Assign copy
        void operator=(const LPtrList<LListBase, T>&);


    // STL type definitions

        //- Type that can be used for storing into LPtrList::value_type
        //  objects.
        typedef T& reference;

        //- Type that can be used for storing into constant
        //  LPtrList::value_type objects.
        typedef T& const_reference;


    // STL iterator

        typedef typename LListBase::iterator LListBase_iterator;

        //- An STL-conforming iterator
        class iterator
        :
            public LList<LListBase, T*>::iterator
        {

        public:

            //- Construct from base iterator
            iterator(LListBase_iterator baseIter)
            :
                LList<LListBase, T*>::iterator(baseIter)
            {}


            // Member operators

                T& operator*()
                {
                    return *(LList<LListBase, T*>::iterator::operator*());
                }

                T& operator()()
                {
                    return operator*();
                }
        };


    // STL const_iterator

        typedef typename LListBase::const_iterator LListBase_const_iterator;

        //- An STL-conforming const_iterator
        class const_iterator
        :
            public LList<LListBase, T*>::const_iterator
        {

        public:

            //- Construct from base const_iterator
            const_iterator(LListBase_const_iterator baseIter)
            :
                LList<LListBase, T*>::const_iterator(baseIter)
            {}

            //- Construct from base iterator
            const_iterator(LListBase_iterator baseIter)
            :
                LList<LListBase, T*>::const_iterator(baseIter)
            {}


            // Member operators

                const T& operator*()
                {
                    return *(LList<LListBase, T*>::const_iterator::operator*());
                }

                const T& operator()()
                {
                    return operator*();
                }
        };


    // IOstream operators

        friend Istream& operator>> <LListBase, T>
        (
            Istream&,
            LPtrList<LListBase, T>&
        );

        friend Ostream& operator<< <LListBase, T>
        (
            Ostream&,
            const LPtrList<LListBase, T>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "LPtrList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
