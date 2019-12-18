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

#include "RKF45.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RKF45, 0);
    addToRunTimeSelectionTable(ODESolver, RKF45, dictionary);

const scalar
    RKF45::c2  = 1.0/4.0,
    RKF45::c3  = 3.0/8.0,
    RKF45::c4  = 12.0/13.0,
    RKF45::c5  = 1.0,
    RKF45::c6  = 1.0/2.0,

    RKF45::a21 = 1.0/4.0,
    RKF45::a31 = 3.0/32.0,
    RKF45::a32 = 9.0/32.0,
    RKF45::a41 = 1932.0/2197.0,
    RKF45::a42 = -7200.0/2197.0,
    RKF45::a43 = 7296.0/2197.0,
    RKF45::a51 = 439.0/216.0,
    RKF45::a52 = -8.0,
    RKF45::a53 = 3680.0/513.0,
    RKF45::a54 = -845.0/4104.0,
    RKF45::a61 = -8.0/27.0,
    RKF45::a62 = 2.0,
    RKF45::a63 = -3544.0/2565.0,
    RKF45::a64 = 1859.0/4104.0,
    RKF45::a65 = -11.0/40.0,

    RKF45::b1  = 16.0/135.0,
    RKF45::b3  = 6656.0/12825.0,
    RKF45::b4  = 28561.0/56430.0,
    RKF45::b5  = -9.0/50.0,
    RKF45::b6  = 2.0/55.0,

    RKF45::e1  = 25.0/216.0 - RKF45::b1,
    RKF45::e3  = 1408.0/2565.0 - RKF45::b3,
    RKF45::e4  = 2197.0/4104.0 - RKF45::b4,
    RKF45::e5  = -1.0/5.0 - RKF45::b5,
    RKF45::e6  = -RKF45::b6;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RKF45::RKF45(const ODESystem& ode, const dictionary& dict)
:
    ODESolver(ode, dict),
    adaptiveSolver(ode, dict),
    yTemp_(n_),
    k2_(n_),
    k3_(n_),
    k4_(n_),
    k5_(n_),
    k6_(n_),
    err_(n_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::RKF45::resize()
{
    if (ODESolver::resize())
    {
        adaptiveSolver::resize(n_);

        resizeField(yTemp_);
        resizeField(k2_);
        resizeField(k3_);
        resizeField(k4_);
        resizeField(k5_);
        resizeField(k6_);
        resizeField(err_);

        return true;
    }
    else
    {
        return false;
    }
}


Foam::scalar Foam::RKF45::solve
(
    const scalar x0,
    const scalarField& y0,
    const label li,
    const scalarField& dydx0,
    const scalar dx,
    scalarField& y
) const
{
    forAll(yTemp_, i)
    {
        yTemp_[i] = y0[i] + a21*dx*dydx0[i];
    }

    odes_.derivatives(x0 + c2*dx, yTemp_, li, k2_);

    forAll(yTemp_, i)
    {
        yTemp_[i] = y0[i] + dx*(a31*dydx0[i] + a32*k2_[i]);
    }

    odes_.derivatives(x0 + c3*dx, yTemp_, li, k3_);

    forAll(yTemp_, i)
    {
        yTemp_[i] = y0[i] + dx*(a41*dydx0[i] + a42*k2_[i] + a43*k3_[i]);
    }

    odes_.derivatives(x0 + c4*dx, yTemp_, li, k4_);

    forAll(yTemp_, i)
    {
        yTemp_[i] = y0[i]
          + dx*(a51*dydx0[i] + a52*k2_[i] + a53*k3_[i] + a54*k4_[i]);
    }

    odes_.derivatives(x0 + c5*dx, yTemp_, li, k5_);

    forAll(yTemp_, i)
    {
        yTemp_[i] = y0[i]
          + dx
           *(a61*dydx0[i] + a62*k2_[i] + a63*k3_[i] + a64*k4_[i] + a65*k5_[i]);
    }

    odes_.derivatives(x0 + c6*dx, yTemp_, li, k6_);

    // Calculate the 5th-order solution
    forAll(y, i)
    {
        y[i] = y0[i]
          + dx
           *(b1*dydx0[i] + b3*k3_[i] + b4*k4_[i] + b5*k5_[i] + b6*k6_[i]);
    }

    // Calculate the error estimate from the difference between the
    // 4th-order and 5th-order solutions
    forAll(err_, i)
    {
        err_[i] =
            dx
           *(e1*dydx0[i] + e3*k3_[i] + e4*k4_[i] + e5*k5_[i] + e6*k6_[i]);
    }

    return normalizeError(y0, y, err_);
}


void Foam::RKF45::solve
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
