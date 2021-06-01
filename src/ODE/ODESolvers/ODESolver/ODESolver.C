/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "ODESolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ODESolver, 0);
    defineRunTimeSelectionTable(ODESolver, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::ODESolver::normaliseError
(
    const scalarField& y0,
    const scalarField& y,
    const scalarField& err
) const
{
    // Calculate the maximum error
    scalar maxErr = 0.0;
    forAll(err, i)
    {
        scalar tol = absTol_[i] + relTol_[i]*max(mag(y0[i]), mag(y[i]));
        maxErr = max(maxErr, mag(err[i])/tol);
    }

    return maxErr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ODESolver::ODESolver(const ODESystem& ode, const dictionary& dict)
:
    odes_(ode),
    maxN_(ode.nEqns()),
    n_(ode.nEqns()),
    absTol_(n_, dict.lookupOrDefault<scalar>("absTol", small)),
    relTol_(n_, dict.lookupOrDefault<scalar>("relTol", 1e-4)),
    maxSteps_(dict.lookupOrDefault<scalar>("maxSteps", 10000))
{}


Foam::ODESolver::ODESolver
(
    const ODESystem& ode,
    const scalarField& absTol,
    const scalarField& relTol
)
:
    odes_(ode),
    maxN_(ode.nEqns()),
    n_(ode.nEqns()),
    absTol_(absTol),
    relTol_(relTol),
    maxSteps_(10000)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::ODESolver::resize()
{
    if (odes_.nEqns() != n_)
    {
        if (odes_.nEqns() > maxN_)
        {
            FatalErrorInFunction
                << "Specified number of equations " << odes_.nEqns()
                << " greater than maximum " << maxN_
                << abort(FatalError);
        }

        n_ = odes_.nEqns();

        resizeField(absTol_);
        resizeField(relTol_);

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::ODESolver::solve
(
    scalar& x,
    scalarField& y,
    const label li,
    scalar& dxTry
) const
{
    stepState step(dxTry);
    solve(x, y, li, step);
    dxTry = step.dxTry;
}


void Foam::ODESolver::solve
(
    scalar& x,
    scalarField& y,
    const label li,
    stepState& step
) const
{
    scalar x0 = x;
    solve(x, y, li, step.dxTry);
    step.dxDid = x - x0;
}


void Foam::ODESolver::solve
(
    const scalar xStart,
    const scalar xEnd,
    scalarField& y,
    const label li,
    scalar& dxTry
) const
{
    stepState step(dxTry);
    scalar x = xStart;

    for (label nStep=0; nStep<maxSteps_; nStep++)
    {
        // Store previous iteration dxTry
        scalar dxTry0 = step.dxTry;

        step.reject = false;

        // Check if this is a truncated step and set dxTry to integrate to xEnd
        if ((x + step.dxTry - xEnd)*(x + step.dxTry - xStart) > 0)
        {
            step.last = true;
            step.dxTry = xEnd - x;
        }

        // Integrate as far as possible up to step.dxTry
        solve(x, y, li, step);

        // Check if reached xEnd
        if ((x - xEnd)*(xEnd - xStart) >= 0)
        {
            if (nStep > 0 && step.last)
            {
                step.dxTry = dxTry0;
            }

            dxTry = step.dxTry;

            return;
        }

        step.first = false;

        // If the step.dxTry was reject set step.prevReject
        if (step.reject)
        {
            step.prevReject = true;
        }
    }

    FatalErrorInFunction
        << "Integration steps greater than maximum " << maxSteps_ << nl
        << "    xStart = " << xStart << ", xEnd = " << xEnd
        << ", x = " << x << ", dxDid = " << step.dxDid << nl
        << "    y = " << y
        << exit(FatalError);
}


// ************************************************************************* //
