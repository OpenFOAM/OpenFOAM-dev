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
    Foam::FIFOStack

Description
    A FIFO stack based on a singly-linked list.

    Operations are push(), pop(), top(), bottom() and empty().

SourceFiles
    FIFOStack.C

\*---------------------------------------------------------------------------*/

#ifndef FIFOStack_H
#define FIFOStack_H

#include "SLList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class FIFOStack Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class FIFOStack
:
    public SLList<T>
{

public:

    // Constructors

        //- Construct null
        FIFOStack()
        {}

        //- Construct given initial T
        FIFOStack(T a)
        :
            SLList<T>(a)
        {}

        //- Construct from Istream
        FIFOStack(Istream& is)
        :
            SLList<T>(is)
        {}


    // Member Functions

        // Access

            //- Return a copy of the top element
            T top() const
            {
                return this->last();
            }

            //- Return a copy of the bottom element
            T bottom() const
            {
                return this->first();
            }


        // Edit

            //- Push an element onto the stack
            void push(const T& a)
            {
                this->append(a);
            }

            //- Pop the bottom element off the stack
            T pop()
            {
                return this->removeHead();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
