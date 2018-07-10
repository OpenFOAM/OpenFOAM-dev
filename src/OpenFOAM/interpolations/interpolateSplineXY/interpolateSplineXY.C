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

#include "interpolateSplineXY.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::Field<Type> Foam::interpolateSplineXY
(
    const scalarField& xNew,
    const scalarField& xOld,
    const Field<Type>& yOld
)
{
    Field<Type> yNew(xNew.size());

    forAll(xNew, i)
    {
        yNew[i] = interpolateSplineXY(xNew[i], xOld, yOld);
    }

    return yNew;
}


template<class Type>
Type Foam::interpolateSplineXY
(
    const scalar x,
    const scalarField& xOld,
    const Field<Type>& yOld
)
{
    label n = xOld.size();

    // early exit if out of bounds or only one value
    if (n == 1 || x < xOld[0])
    {
        return yOld[0];
    }
    if (x > xOld[n - 1])
    {
        return yOld[n - 1];
    }

    // linear interpolation if only two values
    if (n == 2)
    {
        return (x - xOld[0])/(xOld[1] - xOld[0])*(yOld[1] - yOld[0]) + yOld[0];
    }

    // find bounding knots
    label hi = 0;
    while (hi < n && xOld[hi] < x)
    {
        hi++;
    }

    label lo = hi - 1;

    const Type& y1 = yOld[lo];
    const Type& y2 = yOld[hi];

    Type y0;
    if (lo == 0)
    {
        y0 = 2*y1 - y2;
    }
    else
    {
        y0 = yOld[lo - 1];
    }

    Type y3;
    if (hi + 1 == n)
    {
        y3 = 2*y2 - y1;
    }
    else
    {
        y3 = yOld[hi + 1];
    }

    // weighting
    scalar mu = (x - xOld[lo])/(xOld[hi] - xOld[lo]);

    // interpolate
    return
        0.5
       *(
            2*y1
          + mu
           *(
               -y0 + y2
              + mu*((2*y0 - 5*y1 + 4*y2 - y3) + mu*(-y0 + 3*y1 - 3*y2 + y3))
            )
        );
}


// ************************************************************************* //
