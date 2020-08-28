/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2020 OpenFOAM Foundation
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
    Foam::SortableListEFA

Description
    A list that is sorted upon construction or when explicitly requested
    with the sort() method.

SourceFiles
    SortableListEFA.C

\*---------------------------------------------------------------------------*/

#ifndef SortableListEFA_H
#define SortableListEFA_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class SortableListEFA Declaration
\*---------------------------------------------------------------------------*/

template <class Type>
class SortableListEFA
:
    public List<Type>
{
    // Private Data

        //- Original indices
        labelList indices_;


public:

    // Public classes

        //- Less function class used by the sort function
        class more
        {
            const UList<Type>& values_;

        public:

            more(const UList<Type>& values)
            :
                values_(values)
            {}

            bool operator()(const label a, const label b)
            {
                return values_[a] > values_[b];
            }
        };


    // Constructors

        //- Construct from List, sorting the elements. Starts with indices set
        //  to index in argument
        explicit SortableListEFA(const List<Type>&);

        //- Construct given size. Sort later on.
        explicit SortableListEFA(const label size);

        //- Construct given size and initial value. Sort later on.
        SortableListEFA(const label size, const Type&);

        //- Construct as copy.
        SortableListEFA(const SortableListEFA<Type>&);


    // Member Functions

        //- Return the list of sorted indices. Updated every sort.
        const labelList& indices() const
        {
            return indices_;
        }

        //- Size the list. If grow can cause undefined indices (until next sort)
        void setSize(const label);

        //- Sort the list (if changed after construction time)
        void sort();

        //- Partial sort the list (if changed after construction time)
        void partialSort(int M, int start);

        //- Sort the list (if changed after construction time)
        void stableSort();


    // Member Operators

        void operator=(const SortableListEFA<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SortableListEFA.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
