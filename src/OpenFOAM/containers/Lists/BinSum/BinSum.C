/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "BinSum.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class IndexType, class List, class CombineOp>
Foam::BinSum<IndexType, List, CombineOp>::BinSum
(
    const IndexType min,
    const IndexType max,
    const IndexType delta
)
:
    List(ceil((max-min)/delta), pTraits<typename List::value_type>::zero),
    min_(min),
    max_(max),
    delta_(delta),
    lowSum_(pTraits<typename List::value_type>::zero),
    highSum_(pTraits<typename List::value_type>::zero)
{}


template<class IndexType, class List, class CombineOp>
Foam::BinSum<IndexType, List, CombineOp>::BinSum
(
    const IndexType min,
    const IndexType max,
    const IndexType delta,
    const UList<IndexType>& indexVals,
    const List& vals,
    const CombineOp& cop
)
:
    List(ceil((max-min)/delta), pTraits<typename List::value_type>::zero),
    min_(min),
    max_(max),
    delta_(delta),
    lowSum_(pTraits<typename List::value_type>::zero),
    highSum_(pTraits<typename List::value_type>::zero)
{
    forAll(indexVals, i)
    {
        add(indexVals[i], vals[i], cop);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class IndexType, class List, class CombineOp>
void Foam::BinSum<IndexType, List, CombineOp>::add
(
    const IndexType& indexVal,
    const typename List::const_reference val,
    const CombineOp& cop
)
{
    if (indexVal < min_)
    {
        cop(lowSum_, val);
    }
    else if (indexVal >= max_)
    {
        cop(highSum_, val);
    }
    else
    {
        label index = (indexVal-min_)/delta_;
        cop(this->operator[](index), val);
    }
}


template<class IndexType, class List, class CombineOp>
void Foam::BinSum<IndexType, List, CombineOp>::add
(
    const UList<IndexType>& indexVals,
    const List& vals,
    const CombineOp& cop
)
{
    forAll(indexVals, i)
    {
        add(indexVals[i], vals[i], cop);
    }
}


// ************************************************************************* //
