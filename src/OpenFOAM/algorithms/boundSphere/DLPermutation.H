/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::DLPermutation

Description
    A permutation stored in the same manner as a doubly-linked list in order to
    facilitate reordering in constant time

SourceFiles
    DLPermutationI.H

\*---------------------------------------------------------------------------*/

#ifndef DLPermutation_H
#define DLPermutation_H

#include "randomGenerator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class DLPermutation Declaration
\*---------------------------------------------------------------------------*/

template<class LabelList>
class DLPermutation
{
    // Private Data

        //- The first index
        label head_;

        //- For every index, the next index
        LabelList next_;

        //- The last index
        label tail_;

        //- For every index, the previous index
        LabelList previous_;


public:

    // Public Classes

        //- A bi-directional constant iterator
        class const_iterator
        {
            // Private Data

                //- Reference to the permutation
                const DLPermutation& p_;

                //- The current index
                label i_;


        public:

            // Constructors

                //- Construct at an index within a permutation
                inline const_iterator(const DLPermutation&, const label);

                //- Use the default copy constructor
                inline const_iterator(const const_iterator&) = default;


            // Member Operators

                //- Dereference
                inline label operator*() const;

                //- Assignment
                inline const_iterator& operator=(const const_iterator& iter);

                //- Increment
                inline const_iterator& operator++();

                //- Increment
                inline const_iterator operator++(int);

                //- Decrement
                inline const_iterator& operator--();

                //- Decrement
                inline const_iterator operator--(int);

                //- Equality comparison
                inline bool operator==(const const_iterator& iter) const;

                //- Inequality comparison
                inline bool operator!=(const const_iterator& iter) const;


            // Member Functions

                //- Return the next iterator without incrementing
                inline const_iterator next() const;

                //- Return the previous iterator without incrementing
                inline const_iterator previous() const;
        };


    // Constructors

        //- Construct for a given permutation
        inline DLPermutation(const List<label>& permutation);

        //- Construct a random permutation
        inline DLPermutation(const label n, randomGenerator& rndGen);


    // Member Functions

        //- Return the permutation
        inline List<label> permutation(const bool reverse);

        //- Return the starting forward iterator
        inline const_iterator begin() const;

        //- Return the ending forward iterator
        inline const_iterator end() const;

        //- Return the starting reverse iterator
        inline const_iterator rbegin() const;

        //- Return the ending reverse iterator
        inline const_iterator rend() const;

        //- Add a new index to the end of the permutation
        void append();

        //- Add a new index to the start of the permutation
        void prepend();

        //- Raise a given index to the start of the permutation
        void raise(const label i);

        //- Sink a given index to the end of the permutation
        void sink(const label i);

        // !!! More operations are possible; e.g., swap two elements or rotate
        // a sub-section. Implement as necessary.
};


//- Define the type of a permutation based on a list
typedef DLPermutation<List<label>> dLPermutation;


//- Define the type of a permutation based on a dynamic list
typedef DLPermutation<DynamicList<label>> dynamicDLPermutation;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forAllConstIterReverse(Container,container,iter)                       \
    for                                                                        \
    (                                                                          \
        Container::const_iterator iter = (container).rbegin();                 \
        iter != (container).rend();                                            \
        -- iter                                                                \
    )

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DLPermutationI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
