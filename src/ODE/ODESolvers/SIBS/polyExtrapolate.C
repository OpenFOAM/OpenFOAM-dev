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

#include "SIBS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::SIBS::polyExtrapolate
(
    const label iest,
    const scalar xest,
    const scalarField& yest,
    scalarField& yz,
    scalarField& dy,
    scalarField& x,
    scalarRectangularMatrix& d
) const
{
    label n = yz.size();

    x[iest] = xest;

    for (label j=0; j<n; j++)
    {
        dy[j] = yz[j] = yest[j];
    }

    if (iest == 0)
    {
        for (label j=0; j<n; j++)
        {
            d(j, 0) = yest[j];
        }
    }
    else
    {
        scalarField c(yest);

        for (label k1=0; k1<iest; k1++)
        {
            scalar delta = 1.0/(x[iest - k1 - 1] - xest);
            scalar f1 = xest*delta;
            scalar f2 = x[iest - k1 - 1]*delta;

            for (label j=0; j<n; j++)
            {
                scalar q = d[j][k1];
                d[j][k1] = dy[j];
                delta = c[j] - q;
                dy[j] = f1*delta;
                c[j] = f2*delta;
                yz[j] += dy[j];
            }
        }

        for (label j=0; j<n; j++)
        {
            d[j][iest] = dy[j];
        }
    }
}


// ************************************************************************* //
