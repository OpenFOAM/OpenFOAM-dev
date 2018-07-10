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

#include "rodas34.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rodas34, 0);
    addToRunTimeSelectionTable(ODESolver, rodas34, dictionary);

const scalar
    rodas34::c2 = 0.386,
    rodas34::c3 = 0.21,
    rodas34::c4 = 0.63,
    rodas34::d1 =  0.25,
    rodas34::d2 = -0.1043,
    rodas34::d3 =  0.1035,
    rodas34::d4 = -0.3620000000000023e-01,
    rodas34::a21 =  0.1544e1,
    rodas34::a31 =  0.9466785280815826,
    rodas34::a32 =  0.2557011698983284,
    rodas34::a41 =  0.3314825187068521e1,
    rodas34::a42 =  0.2896124015972201e1,
    rodas34::a43 =  0.9986419139977817,
    rodas34::a51 =  0.1221224509226641e1,
    rodas34::a52 =  0.6019134481288629e1,
    rodas34::a53 =  0.1253708332932087e2,
    rodas34::a54 = -0.6878860361058950,
    rodas34::c21 = -0.56688e1,
    rodas34::c31 = -0.2430093356833875e1,
    rodas34::c32 = -0.2063599157091915,
    rodas34::c41 = -0.1073529058151375,
    rodas34::c42 = -0.9594562251023355e1,
    rodas34::c43 = -0.2047028614809616e2,
    rodas34::c51 =  0.7496443313967647e1,
    rodas34::c52 = -0.1024680431464352e2,
    rodas34::c53 = -0.3399990352819905e2,
    rodas34::c54 =  0.1170890893206160e2,
    rodas34::c61 =  0.8083246795921522e1,
    rodas34::c62 = -0.7981132988064893e1,
    rodas34::c63 = -0.3152159432874371e2,
    rodas34::c64 =  0.1631930543123136e2,
    rodas34::c65 = -0.6058818238834054e1,
    rodas34::gamma = 0.25;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rodas34::rodas34(const ODESystem& ode, const dictionary& dict)
:
    ODESolver(ode, dict),
    adaptiveSolver(ode, dict),
    k1_(n_),
    k2_(n_),
    k3_(n_),
    k4_(n_),
    k5_(n_),
    dy_(n_),
    err_(n_),
    dydx_(n_),
    dfdx_(n_),
    dfdy_(n_, n_),
    a_(n_, n_),
    pivotIndices_(n_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::rodas34::resize()
{
    if (ODESolver::resize())
    {
        adaptiveSolver::resize(n_);

        resizeField(k1_);
        resizeField(k2_);
        resizeField(k3_);
        resizeField(k4_);
        resizeField(k5_);
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


Foam::scalar Foam::rodas34::solve
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
    forAll(y, i)
    {
        y[i] = y0[i] + a21*k1_[i];
    }

    odes_.derivatives(x0 + c2*dx, y, dydx_);

    forAll(k2_, i)
    {
        k2_[i] = dydx_[i] + dx*d2*dfdx_[i] + c21*k1_[i]/dx;
    }

    LUBacksubstitute(a_, pivotIndices_, k2_);

    // Calculate k3:
    forAll(y, i)
    {
        y[i] = y0[i] + a31*k1_[i] + a32*k2_[i];
    }

    odes_.derivatives(x0 + c3*dx, y, dydx_);

    forAll(k3_, i)
    {
        k3_[i] = dydx_[i] + dx*d3*dfdx_[i] + (c31*k1_[i] + c32*k2_[i])/dx;
    }

    LUBacksubstitute(a_, pivotIndices_, k3_);

    // Calculate k4:
    forAll(y, i)
    {
        y[i] = y0[i] + a41*k1_[i] + a42*k2_[i] + a43*k3_[i];
    }

    odes_.derivatives(x0 + c4*dx, y, dydx_);

    forAll(k4_, i)
    {
        k4_[i] = dydx_[i] + dx*d4*dfdx_[i]
          + (c41*k1_[i] + c42*k2_[i] + c43*k3_[i])/dx;
    }

    LUBacksubstitute(a_, pivotIndices_, k4_);

    // Calculate k5:
    forAll(y, i)
    {
        dy_[i] = a51*k1_[i] + a52*k2_[i] + a53*k3_[i] + a54*k4_[i];
        y[i] = y0[i] + dy_[i];
    }

    odes_.derivatives(x0 + dx, y, dydx_);

    forAll(k5_, i)
    {
        k5_[i] = dydx_[i]
          + (c51*k1_[i] + c52*k2_[i] + c53*k3_[i] + c54*k4_[i])/dx;
    }

    LUBacksubstitute(a_, pivotIndices_, k5_);

    // Calculate new state and error
    forAll(y, i)
    {
        dy_[i] += k5_[i];
        y[i] = y0[i] + dy_[i];
    }

    odes_.derivatives(x0 + dx, y, dydx_);

    forAll(err_, i)
    {
        err_[i] = dydx_[i]
          + (c61*k1_[i] + c62*k2_[i] + c63*k3_[i] + c64*k4_[i] + c65*k5_[i])/dx;
    }

    LUBacksubstitute(a_, pivotIndices_, err_);

    forAll(y, i)
    {
        y[i] = y0[i] + dy_[i] + err_[i];
    }

    return normalizeError(y0, y, err_);
}


void Foam::rodas34::solve
(
    scalar& x,
    scalarField& y,
    scalar& dxTry
) const
{
    adaptiveSolver::solve(odes_, x, y, dxTry);
}


// ************************************************************************* //
