/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

void Foam::SIBS::SIMPR
(
    const scalar xStart,
    const scalarField& y,
    const label li,
    const scalarField& dydx,
    const scalarField& dfdx,
    const scalarSquareMatrix& dfdy,
    const scalar deltaX,
    const label nSteps,
    scalarField& yEnd
) const
{
    scalar h = deltaX/nSteps;

    scalarSquareMatrix a(n_);
    for (label i=0; i<n_; i++)
    {
        for (label j=0; j<n_; j++)
        {
            a(i, j) = -h*dfdy(i, j);
        }
        ++a(i, i);
    }

    labelList pivotIndices(n_);
    LUDecompose(a, pivotIndices);

    for (label i=0; i<n_; i++)
    {
        yEnd[i] = h*(dydx[i] + h*dfdx[i]);
    }

    LUBacksubstitute(a, pivotIndices, yEnd);

    scalarField del(yEnd);
    scalarField ytemp(n_);

    for (label i=0; i<n_; i++)
    {
        ytemp[i] = y[i] + del[i];
    }

    scalar x = xStart + h;

    odes_.derivatives(x, ytemp, li, yEnd);

    for (label nn=2; nn<=nSteps; nn++)
    {
        for (label i=0; i<n_; i++)
        {
            yEnd[i] = h*yEnd[i] - del[i];
        }

        LUBacksubstitute(a, pivotIndices, yEnd);

        for (label i=0; i<n_; i++)
        {
            ytemp[i] += (del[i] += 2.0*yEnd[i]);
        }

        x += h;

        odes_.derivatives(x, ytemp, li, yEnd);
    }
    for (label i=0; i<n_; i++)
    {
        yEnd[i] = h*yEnd[i] - del[i];
    }

    LUBacksubstitute(a, pivotIndices, yEnd);

    for (label i=0; i<n_; i++)
    {
        yEnd[i] += ytemp[i];
    }
}


// ************************************************************************* //
