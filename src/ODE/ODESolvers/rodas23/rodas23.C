/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2021 OpenFOAM Foundation
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

#include "rodas23.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rodas23, 0);
    addToRunTimeSelectionTable(ODESolver, rodas23, dictionary);

const scalar
    rodas23::c3 = 1,
    rodas23::d1 = 1.0/2.0,
    rodas23::d2 = 3.0/2.0,
    rodas23::a31 = 2,
    rodas23::a41 = 2,
    rodas23::c21 = 4,
    rodas23::c31 = 1,
    rodas23::c32 = -1,
    rodas23::c41 = 1,
    rodas23::c42 = -1,
    rodas23::c43 = -8.0/3.0,
    rodas23::gamma = 1.0/2.0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rodas23::rodas23(const ODESystem& ode, const dictionary& dict)
:
    ODESolver(ode, dict),
    adaptiveSolver(ode, dict),
    k1_(n_),
    k2_(n_),
    k3_(n_),
    dy_(n_),
    err_(n_),
    dydx_(n_),
    dfdx_(n_),
    dfdy_(n_, n_),
    a_(n_, n_),
    pivotIndices_(n_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::rodas23::resize()
{
    if (ODESolver::resize())
    {
        adaptiveSolver::resize(n_);

        resizeField(k1_);
        resizeField(k2_);
        resizeField(k3_);
        resizeField(dy_);
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


Foam::scalar Foam::rodas23::solve
(
    const scalar x0,
    const scalarField& y0,
    const label li,
    const scalarField& dydx0,
    const scalar dx,
    scalarField& y
) const
{
    odes_.jacobian(x0, y0, li, dfdx_, dfdy_);

    for (label i=0; i<n_; i++)
    {
        for (label j=0; j<n_; j++)
        {
            a_(i, j) = -dfdy_(i, j);
        }

        a_(i, i) += 1.0/(gamma*dx);
    }

    LUDecompose(a_, pivotIndices_);

    // Calculate k1:
    forAll(k1_, i)
    {
        k1_[i] = dydx0[i] + dx*d1*dfdx_[i];
    }

    LUBacksubstitute(a_, pivotIndices_, k1_);

    // Calculate k2:
    forAll(k2_, i)
    {
        k2_[i] = dydx0[i] + dx*d2*dfdx_[i] + c21*k1_[i]/dx;
    }

    LUBacksubstitute(a_, pivotIndices_, k2_);

    // Calculate k3:
    forAll(y, i)
    {
        dy_[i] = a31*k1_[i];
        y[i] = y0[i] + dy_[i];
    }

    odes_.derivatives(x0 + dx, y, li, dydx_);

    forAll(k3_, i)
    {
        k3_[i] = dydx_[i] + (c31*k1_[i] + c32*k2_[i])/dx;
    }

    LUBacksubstitute(a_, pivotIndices_, k3_);

    // Calculate new state and error
    forAll(y, i)
    {
        dy_[i] += k3_[i];
        y[i] = y0[i] + dy_[i];
    }

    odes_.derivatives(x0 + dx, y, li, dydx_);

    forAll(err_, i)
    {
        err_[i] = dydx_[i] + (c41*k1_[i] + c42*k2_[i] + c43*k3_[i])/dx;
    }

    LUBacksubstitute(a_, pivotIndices_, err_);

    forAll(y, i)
    {
        y[i] = y0[i] + dy_[i] + err_[i];
    }

    return normaliseError(y0, y, err_);
}


void Foam::rodas23::solve
(
    scalar& x,
    scalarField& y,
    const label li,
    scalar& dxTry
) const
{
    adaptiveSolver::solve(odes_, x, y, li, dxTry);
}


// ************************************************************************* //
