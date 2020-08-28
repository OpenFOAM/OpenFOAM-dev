/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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
    Foam::BinSum

Description
    Sums into bins

SourceFiles
    BinSum.C

\*---------------------------------------------------------------------------*/

#ifndef BinSum_H
#define BinSum_H

#include "ops.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class BinSum Declaration
\*---------------------------------------------------------------------------*/

template
<
    class IndexType,
    class List,
    class CombineOp = plusEqOp<typename List::value_type>
>
class BinSum
:
    public List
{
    // Private Data

        const IndexType min_;

        const IndexType max_;

        const IndexType delta_;


        //- Sum < lowest bin
        typename List::value_type lowSum_;

        //- Sum of >= highest bin
        typename List::value_type highSum_;


public:

    // Constructors

        //- Construct given min, max, delta
        BinSum
        (
            const IndexType min,
            const IndexType max,
            const IndexType delta
        );

        //- Construct given min, max, delta and data
        BinSum
        (
            const IndexType min,
            const IndexType max,
            const IndexType delta,
            const UList<IndexType>& indexVals,
            const List& vals,
            const CombineOp& cop = plusEqOp<typename List::value_type>()
        );


        // Access

            //- Return the delta
            inline IndexType delta() const
            {
                return delta_;
            }

            //- Return the sum of all added elements < min
            inline const IndexType& lowSum() const
            {
                return lowSum_;
            }

            //- Return the sum of all added elements >= max
            inline const IndexType& highSum() const
            {
                return highSum_;
            }

            void add
            (
                const IndexType& indexVal,
                const typename List::const_reference val,
                const CombineOp& cop = plusEqOp<typename List::value_type>()
            );

            void add
            (
                const UList<IndexType>& indexVals,
                const List& vals,
                const CombineOp& cop = plusEqOp<typename List::value_type>()
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BinSum.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
