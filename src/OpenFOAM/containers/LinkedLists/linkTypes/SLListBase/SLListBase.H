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
    Foam::SLListBase

Description
    Base singly-linked list.

SourceFiles
    SLListBase.C

\*---------------------------------------------------------------------------*/

#ifndef SLListBase_H
#define SLListBase_H

#include "bool.H"
#include "label.H"
#include "uLabel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class SLListBase Declaration
\*---------------------------------------------------------------------------*/

class SLListBase
{

public:

    //- Link structure
    struct link
    {
        //- Pointer to next entry in list
        link* next_;

        //- Null construct
        inline link();

        //- Construct given pointer to another link
        inline link(link* p);
    };


private:

    // Private Data

       //- last_ points to last element
       //  last_->next_ points to first element, i.e. circular storage
       link* last_;

       //- Number of elements in in list
       label nElmts_;


public:

    // Forward declaration of STL iterators

        class iterator;
        friend class iterator;

        class const_iterator;
        friend class const_iterator;


    // Constructors

        //- Null construct
        inline SLListBase();

        //- Construct given initial entry
        inline SLListBase(link*);

        //- Disallow default bitwise copy construction
        SLListBase(const SLListBase&) = delete;


    //- Destructor
    ~SLListBase();


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

            //- Remove and return head
            link* removeHead();

            // Remove and return element
            link* remove(link*);

            // Remove and return element specified by iterator
            inline link* remove(iterator&);

            //- Clear the list
            inline void clear();

            //- Transfer the contents of the argument into this List
            //  and annul the argument list.
            inline void transfer(SLListBase&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const SLListBase&) = delete;


    // STL iterator

        //- An STL-conforming iterator
        class iterator
        {
            friend class SLListBase;
            friend class const_iterator;

            // Private Data

                //- Reference to the list this is an iterator for
                SLListBase& curList_;

                //- Current element
                link* curElmt_;

                //- Copy of the link
                link curLink_;

            // Private Member Functions

            //- Construct for a given SLListBase with nullptr element and link.
            //  Only used to create endIter
            inline iterator(SLListBase&);

        public:

            //- Construct for a given SLListBase and link
            inline iterator(SLListBase&, link*);

            //- Copy constructor
            inline iterator(const iterator&) = default;

            // Member Operators

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
            // Private Data

                //- Reference to the list this is an iterator for
                const SLListBase& curList_;

                //- Current element
                const link* curElmt_;

        public:

            //- Construct for a given SLListBase and link
            inline const_iterator(const SLListBase&, const link*);

            //- Construct from a non-const iterator
            inline const_iterator(const iterator&);

            //- Copy constructor
            inline const_iterator(const const_iterator&) = default;

            // Member Operators

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


private:

        //- Iterator returned by end()
        static iterator endIter_;

        //- const_iterator returned by end()
        static const_iterator endConstIter_;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SLListBaseI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
