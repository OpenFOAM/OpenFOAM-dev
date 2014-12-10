/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::ParSortableList

Description
    Implementation of PSRS parallel sorting routine.

    From "On the Versatility of Parallel Sorting by Regular Sampling"
    Xiaobo Li et. all.

    Construct from list of things to sort (uses SortableList, 'thing' should
    implement >, ==).

    Will contain sorted data and in
        - indices() the original indices
        - procs() the original processor id.

    Can also be constructed from size, filled at ease and then sort()'ed.

SourceFiles
    ParSortableList.C

\*---------------------------------------------------------------------------*/

#ifndef ParSortableList_H
#define ParSortableList_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                        Class ParSortableListName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(ParSortableList);


/*---------------------------------------------------------------------------*\
                           Class ParSortableList Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class ParSortableList
:
    public ParSortableListName,
    public List<Type>
{
    // Private data

        //- Original indices
        labelList indices_;

        //- Processor numbers
        labelList procs_;


    // Private class

        //- Private class for sorting. Sorts on value_.
        class taggedValue
        {
            // Private data

                //- Value
                Type value_;

                //- Index
                label index_;

                //- Proc
                label procID_;


        public:

            // Constructors

                taggedValue()
                {}

            // Member Functions

                const Type& value() const
                {
                    return value_;
                }

                Type& value()
                {
                    return value_;
                }

                label index() const
                {
                    return index_;
                }

                label& index()
                {
                    return index_;
                }

                label procID() const
                {
                    return procID_;
                }

                label& procID()
                {
                    return procID_;
                }
        };


    // Private Member Functions

        //- Write for debugging
        void write(const List<Type>& elems, Ostream& os) const;

        //- Copy into dest, setting indices and processor no.
        void copyInto
        (
            const List<Type>& values,
            const labelList& indices,
            const label fromProcNo,
            label& destI,
            List<taggedValue>& dest
        ) const;

        //- Calculate pivots.
        void getPivots(const List<Type>& elems, List<Type>& pivots) const;

        //- If destProc != myProcNo send & clear values and indices
        void checkAndSend
        (
            List<Type>& values,
            labelList& indices,
            const label bufSize,
            const label destProcI
        ) const;


public:

    // Constructors

        //- Construct from List, sorting the elements.
        ParSortableList(const UList<Type>&);

        //- Construct given size. Sort later on.
        ParSortableList(const label size);


    // Member Functions

        //- (stable) sort the list (if changed after construction time)
        void sort();

        //- Return the list of sorted point indices
        const labelList& indices() const
        {
            return indices_;
        }

        //- Return the list of processor number
        const labelList& procs() const
        {
            return procs_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ParSortableList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
