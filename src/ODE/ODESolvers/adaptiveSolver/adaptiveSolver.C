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

#include "adaptiveSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adaptiveSolver::adaptiveSolver
(
    const ODESystem& ode,
    const dictionary& dict
)
:
    safeScale_(dict.lookupOrDefault<scalar>("safeScale", 0.9)),
    alphaInc_(dict.lookupOrDefault<scalar>("alphaIncrease", 0.2)),
    alphaDec_(dict.lookupOrDefault<scalar>("alphaDecrease", 0.25)),
    minScale_(dict.lookupOrDefault<scalar>("minScale", 0.2)),
    maxScale_(dict.lookupOrDefault<scalar>("maxScale", 10)),
    dydx0_(ode.nEqns()),
    yTemp_(ode.nEqns())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::adaptiveSolver::resize(const label n)
{
    ODESolver::resizeField(dydx0_, n);
    ODESolver::resizeField(yTemp_, n);

    return true;
}


void Foam::adaptiveSolver::solve
(
    const ODESystem& odes,
    scalar& x,
    scalarField& y,
    scalar& dxTry
) const
{
    scalar dx = dxTry;
    scalar err = 0.0;

    odes.derivatives(x, y, dydx0_);

    // Loop over solver and adjust step-size as necessary
    // to achieve desired error
    do
    {
        // Solve step and provide error estimate
        err = solve(x, y, dydx0_, dx, yTemp_);

        // If error is large reduce dx
        if (err > 1)
        {
            scalar scale = max(safeScale_*pow(err, -alphaDec_), minScale_);
            dx *= scale;

            if (dx < vSmall)
            {
                FatalErrorInFunction
                    << "stepsize underflow"
                    << exit(FatalError);
            }
        }
    } while (err > 1);

    // Update the state
    x += dx;
    y = yTemp_;

    // If the error is small increase the step-size
    if (err > pow(maxScale_/safeScale_, -1.0/alphaInc_))
    {
        dxTry =
            min(max(safeScale_*pow(err, -alphaInc_), minScale_), maxScale_)*dx;
    }
    else
    {
        dxTry = safeScale_*maxScale_*dx;
    }
}


// ************************************************************************* //
