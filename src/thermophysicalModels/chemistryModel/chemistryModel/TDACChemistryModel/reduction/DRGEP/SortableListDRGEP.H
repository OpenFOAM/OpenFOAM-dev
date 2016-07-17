/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::SortableListDRGEP

Description
    A list that is sorted upon construction or when explicitly requested
    with the sort() method.

SourceFiles
    SortableListDRGEP.C

\*---------------------------------------------------------------------------*/

#ifndef SortableListDRGEP_H
#define SortableListDRGEP_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class SortableListDRGEP Declaration
\*---------------------------------------------------------------------------*/

template <class Type>
class SortableListDRGEP
:
    public List<Type>
{
    // Private data

        //- Original indices
        labelList indices_;


public:

    // Public classes

        //- Less function class used by the sort function
        class less
        {
            const UList<Type>& values_;

        public:

            less(const UList<Type>& values)
            :
                values_(values)
            {}

            bool operator()(const label a, const label b)
            {
                return values_[a] < values_[b];
            }
        };


    // Constructors

        //- Construct from List, sorting the elements. Starts with indices set
        //  to index in argument
        explicit SortableListDRGEP(const List<Type>&);

        //- Construct given size. Sort later on.
        explicit SortableListDRGEP(const label size);

        //- Construct given size and initial value. Sort later on.
        SortableListDRGEP(const label size, const Type&);

        //- Construct as copy.
        SortableListDRGEP(const SortableListDRGEP<Type>&);


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
        void partialSort(int M);

        //- Sort the list (if changed after construction time)
        void stableSort();


    // Member Operators

        void operator=(const SortableListDRGEP<Type>&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SortableListDRGEP.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
