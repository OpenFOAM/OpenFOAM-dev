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
    Foam::DLListBase

Description
    Base doubly-linked list.

SourceFiles
    DLListBase.C

\*---------------------------------------------------------------------------*/

#ifndef DLListBase_H
#define DLListBase_H

#include "bool.H"
#include "label.H"
#include "uLabel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class DLListBase Declaration
\*---------------------------------------------------------------------------*/

class DLListBase
{

public:

    //- Link structure
    struct link
    {
        //- Pointer to next entry in list
        link *prev_, *next_;

        //- Null construct
        inline link();

        //- Check if the link is registered with the DLListBase
        inline bool registered() const;

        //- Deregister the link after removal
        inline void deregister();
    };


private:

    // Private data

       //- first_ points to first element and last_ points to last element.
       link *first_, *last_;

       //- Number of elements in in list
       label nElmts_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        DLListBase(const DLListBase&);

        //- Disallow default bitwise assignment
        void operator=(const DLListBase&);


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
        inline DLListBase();

        //- Construct given initial entry
        inline DLListBase(link*);


    //- Destructor
    ~DLListBase();


    // Member Functions

        // Access

            //- Return number of elements in list
            inline label size() const;

            //- Return true if the list is empty
            inline bool empty() const;

            //- Return first entry
            inline link* first();

            //- Return const access to first entry
            inline const link* first() const;

            //- Return last entry
            inline link* last();

            //- Return const access to last entry
            inline const link* last() const;


        // Edit

            //- Add at head of list
            void insert(link*);

            //- Add at tail of list
            void append(link*);

            //- Swap this element with the one above unless it is at the top
            bool swapUp(link*);

            //- Swap this element with the one below unless it is at the bottom
            bool swapDown(link*);

            //- Remove and return head
            link* removeHead();

            //- Remove and return element
            link* remove(link*);

            // Remove and return element specified by iterator
            inline link* remove(iterator&);

            //- Replace oldLink with newLink and return element
            link* replace(link* oldLink, link* newLink);

            //- Replace oldIter with newLink and return element
            inline link* replace(iterator& oldIter, link* newLink);

            //- Clear the list
            inline void clear();

            //- Transfer the contents of the argument into this List
            //  and annul the argument list.
            inline void transfer(DLListBase&);

    // STL iterator

        //- An STL-conforming iterator
        class iterator
        {
            friend class DLListBase;
            friend class const_iterator;

            // Private data

                //- Reference to the list this is an iterator for
                DLListBase& curList_;

                //- Current element
                link* curElmt_;

                //- Copy of the link
                link curLink_;

            // Private Member Functions

            //- Construct for a given SLListBase with nullptr element and link.
            //  Only used to create endIter
            inline iterator(DLListBase&);

        public:

            //- Construct for a given DLListBase and link
            inline iterator(DLListBase&, link*);

            // Member operators

                inline void operator=(const iterator&);

                inline bool operator==(const iterator&) const;
                inline bool operator!=(const iterator&) const;

                inline link& operator*();

                inline iterator& operator++();
                inline iterator operator++(int);
        };

        inline iterator begin();
        inline const iterator& end();


    // STL const_iterator

        //- An STL-conforming const_iterator
        class const_iterator
        {
            // Private data

                //- Reference to the list this is an iterator for
                const DLListBase& curList_;

                //- Current element
                const link* curElmt_;

        public:

            //- Construct for a given DLListBase and link
            inline const_iterator(const DLListBase&, const link*);

            //- Construct from a non-const iterator
            inline const_iterator(const iterator&);

            // Member operators

                inline void operator=(const const_iterator&);

                inline bool operator==(const const_iterator&) const;
                inline bool operator!=(const const_iterator&) const;

                inline const link& operator*();

                inline const_iterator& operator++();
                inline const_iterator operator++(int);
        };

        inline const_iterator cbegin() const;
        inline const const_iterator& cend() const;

        inline const_iterator begin() const;
        inline const const_iterator& end() const;


    // STL const_reverse_iterator

        //- An STL-conforming const_reverse_iterator
        class const_reverse_iterator
        {
            // Private data

                //- Reference to the list this is an reverse_iterator for
                const DLListBase& curList_;

                //- Current element
                const link* curElmt_;

        public:

            //- Construct for a given DLListBase and link
            inline const_reverse_iterator(const DLListBase&, const link*);

            // Member operators

                inline void operator=(const const_reverse_iterator&);

                inline bool operator==(const const_reverse_iterator&) const;
                inline bool operator!=(const const_reverse_iterator&) const;

                inline const link& operator*();

                inline const_reverse_iterator& operator++();
                inline const_reverse_iterator operator++(int);
        };

        inline const_reverse_iterator crbegin() const;
        inline const const_reverse_iterator& crend() const;

        inline const_reverse_iterator rbegin() const;
        inline const const_reverse_iterator& rend() const;


private:

        //- Iterator returned by end()
        static iterator endIter_;

        //- const_iterator returned by end()
        static const_iterator endConstIter_;

        //- const_reverse_iterator returned by end()
        static const_reverse_iterator endConstRevIter_;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DLListBaseI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
