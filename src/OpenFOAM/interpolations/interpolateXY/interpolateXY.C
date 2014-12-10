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

\*---------------------------------------------------------------------------*/

#include "interpolateXY.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Field<Type> interpolateXY
(
    const scalarField& xNew,
    const scalarField& xOld,
    const Field<Type>& yOld
)
{
    Field<Type> yNew(xNew.size());

    forAll(xNew, i)
    {
        yNew[i] = interpolateXY(xNew[i], xOld, yOld);
    }

    return yNew;
}


template<class Type>
Type interpolateXY
(
    const scalar x,
    const scalarField& xOld,
    const Field<Type>& yOld
)
{
    label n = xOld.size();

    label lo = 0;
    for (lo=0; lo<n && xOld[lo]>x; ++lo)
    {}

    label low = lo;
    if (low < n)
    {
        for (label i=low; i<n; ++i)
        {
            if (xOld[i] > xOld[lo] && xOld[i] <= x)
            {
                lo = i;
            }
        }
    }

    label hi = 0;
    for (hi=0; hi<n && xOld[hi]<x; ++hi)
    {}

    label high = hi;
    if (high < n)
    {
        for (label i=high; i<n; ++i)
        {
            if (xOld[i] < xOld[hi] && xOld[i] >= x)
            {
                hi = i;
            }
        }
    }


    if (lo<n && hi<n && lo != hi)
    {
        return yOld[lo]
            + ((x - xOld[lo])/(xOld[hi] - xOld[lo]))*(yOld[hi] - yOld[lo]);
    }
    else if (lo == hi)
    {
        return yOld[lo];
    }
    else if (lo == n)
    {
        return yOld[hi];
    }
    else
    {
        return yOld[lo];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
