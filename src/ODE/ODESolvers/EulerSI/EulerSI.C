/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "EulerSI.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(EulerSI, 0);
    addToRunTimeSelectionTable(ODESolver, EulerSI, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::EulerSI::EulerSI(const ODESystem& ode, const dictionary& dict)
:
    ODESolver(ode, dict),
    adaptiveSolver(ode, dict),
    err_(n_),
    dydx_(n_),
    dfdx_(n_),
    dfdy_(n_, n_),
    a_(n_, n_),
    pivotIndices_(n_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::EulerSI::resize()
{
    if (ODESolver::resize())
    {
        adaptiveSolver::resize(n_);

        resizeField(err_);
        resizeField(dydx_);
        resizeField(dfdx_);
        resizeMatrix(dfdy_);
        resizeMatrix(a_);
        resizeField(pivotIndices_);

        return true;
    }
    else
    {
        return false;
    }
}


Foam::scalar Foam::EulerSI::solve
(
    const scalar x0,
    const scalarField& y0,
    const scalarField& dydx0,
    const scalar dx,
    scalarField& y
) const
{
    odes_.jacobian(x0, y0, dfdx_, dfdy_);

    for (label i=0; i<n_; i++)
    {
        for (label j=0; j<n_; j++)
        {
            a_(i, j) = -dfdy_(i, j);
        }

        a_(i, i) += 1.0/dx;
    }

    LUDecompose(a_, pivotIndices_);

    // Calculate error estimate from the change in state:
    forAll(err_, i)
    {
        err_[i] = dydx0[i] + dx*dfdx_[i];
    }

    LUBacksubstitute(a_, pivotIndices_, err_);

    forAll(y, i)
    {
        y[i] = y0[i] + err_[i];
    }

    return normalizeError(y0, y, err_);
}


void Foam::EulerSI::solve
(
    scalar& x,
    scalarField& y,
    scalar& dxTry
) const
{
    adaptiveSolver::solve(odes_, x, y, dxTry);
}


// ************************************************************************* //
