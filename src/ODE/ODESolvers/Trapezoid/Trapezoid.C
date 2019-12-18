/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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

#include "Trapezoid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Trapezoid, 0);
    addToRunTimeSelectionTable(ODESolver, Trapezoid, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Trapezoid::Trapezoid(const ODESystem& ode, const dictionary& dict)
:
    ODESolver(ode, dict),
    adaptiveSolver(ode, dict),
    err_(n_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::Trapezoid::resize()
{
    if (ODESolver::resize())
    {
        adaptiveSolver::resize(n_);

        resizeField(err_);

        return true;
    }
    else
    {
        return false;
    }
}


Foam::scalar Foam::Trapezoid::solve
(
    const scalar x0,
    const scalarField& y0,
    const label li,
    const scalarField& dydx0,
    const scalar dx,
    scalarField& y
) const
{
    // Predict the state using 1st-order Trapezoid method
    forAll(y, i)
    {
        y[i] = y0[i] + dx*dydx0[i];
    }

    // Evaluate the system for the predicted state
    odes_.derivatives(x0 + dx, y, li, err_);

    // Update the state as the average between the prediction and the correction
    // and estimate the error from the difference
    forAll(y, i)
    {
        y[i] = y0[i] + 0.5*dx*(dydx0[i] + err_[i]);
        err_[i] = 0.5*dx*(err_[i] - dydx0[i]);
    }

    return normalizeError(y0, y, err_);
}


void Foam::Trapezoid::solve
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
