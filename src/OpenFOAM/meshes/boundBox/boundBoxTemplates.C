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

\*---------------------------------------------------------------------------*/

#include "boundBox.H"
#include "FixedList.H"
#include "PstreamReduceOps.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<unsigned Size>
Foam::boundBox::boundBox
(
    const UList<point>& points,
    const FixedList<label, Size>& indices,
    const bool doReduce
)
:
    min_(Zero),
    max_(Zero)
{
    // a FixedList is never empty
    if (points.empty())
    {
        if (doReduce && Pstream::parRun())
        {
            // Use values that get overwritten by reduce minOp, maxOp below
            min_ = point(vGreat, vGreat, vGreat);
            max_ = point(-vGreat, -vGreat, -vGreat);
        }
    }
    else
    {
        min_ = points[indices[0]];
        max_ = points[indices[0]];

        for (unsigned i=1; i < Size; ++i)
        {
            min_ = ::Foam::min(min_, points[indices[i]]);
            max_ = ::Foam::max(max_, points[indices[i]]);
        }
    }

    // Reduce parallel information
    if (doReduce)
    {
        reduce(min_, minOp<point>());
        reduce(max_, maxOp<point>());
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<unsigned Size>
bool Foam::boundBox::contains
(
    const UList<point>& points,
    const FixedList<label, Size>& indices
) const
{
    // a FixedList is never empty
    if (points.empty())
    {
        return false;
    }

    forAll(indices, i)
    {
        if (!contains(points[indices[i]]))
        {
            return false;
        }
    }

    return true;
}


template<unsigned Size>
bool Foam::boundBox::containsAny
(
    const UList<point>& points,
    const FixedList<label, Size>& indices
) const
{
    // a FixedList is never empty
    if (points.empty())
    {
        return false;
    }

    forAll(indices, i)
    {
        if (contains(points[indices[i]]))
        {
            return true;
        }
    }

    return false;
}


// ************************************************************************* //
