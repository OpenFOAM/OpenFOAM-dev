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
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SIBS, 0);
    addToRunTimeSelectionTable(ODESolver, SIBS, dictionary);

    const label SIBS::nSeq_[iMaxX_] = {2, 6, 10, 14, 22, 34, 50, 70};

    const scalar
        SIBS::safe1 = 0.25,
        SIBS::safe2 = 0.7,
        SIBS::redMax = 1.0e-5,
        SIBS::redMin = 0.7,
        SIBS::scaleMX = 0.1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SIBS::SIBS(const ODESystem& ode, const dictionary& dict)
:
    ODESolver(ode, dict),
    a_(iMaxX_, 0.0),
    alpha_(kMaxX_, 0.0),
    d_p_(n_, kMaxX_, 0.0),
    x_p_(kMaxX_, 0.0),
    err_(kMaxX_, 0.0),

    yTemp_(n_, 0.0),
    ySeq_(n_, 0.0),
    yErr_(n_, 0.0),
    dydx0_(n_),
    dfdx_(n_, 0.0),
    dfdy_(n_, 0.0),
    first_(1),
    epsOld_(-1.0)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::SIBS::resize()
{
    if (ODESolver::resize())
    {
        resizeField(yTemp_);
        resizeField(ySeq_);
        resizeField(yErr_);
        resizeField(dydx0_);
        resizeField(dfdx_);
        resizeMatrix(dfdy_);

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::SIBS::solve
(
    scalar& x,
    scalarField& y,
    const label li,
    scalar& dxTry
) const
{
    odes_.derivatives(x, y, li, dydx0_);

    scalar h = dxTry;
    bool exitflag = false;

    if (relTol_[0] != epsOld_)
    {
        dxTry = xNew_ = -great;
        scalar eps1 = safe1*relTol_[0];
        a_[0] = nSeq_[0] + 1;

        for (label k=0; k<kMaxX_; k++)
        {
            a_[k + 1] = a_[k] + nSeq_[k + 1];
        }

        for (label iq = 1; iq<kMaxX_; iq++)
        {
            for (label k=0; k<iq; k++)
            {
                alpha_[k][iq] =
                    pow(eps1, (a_[k + 1] - a_[iq + 1])
                   /((a_[iq + 1] - a_[0] + 1.0)*(2*k + 3)));
            }
        }

        epsOld_ = relTol_[0];
        a_[0] += n_;

        for (label k=0; k<kMaxX_; k++)
        {
            a_[k + 1] = a_[k] + nSeq_[k + 1];
        }

        for (kOpt_ = 1; kOpt_<kMaxX_ - 1; kOpt_++)
        {
            if (a_[kOpt_ + 1] > a_[kOpt_]*alpha_[kOpt_ - 1][kOpt_])
            {
                break;
            }
        }

        kMax_ = kOpt_;
    }

    label k = 0;
    yTemp_ = y;

    odes_.jacobian(x, y, li, dfdx_, dfdy_);

    if (x != xNew_ || h != dxTry)
    {
        first_ = 1;
        kOpt_ = kMax_;
    }

    label km=0;
    label reduct=0;
    scalar maxErr = small;

    for (;;)
    {
        scalar red = 1.0;

        for (k=0; k <= kMax_; k++)
        {
            xNew_ = x + h;

            if (xNew_ == x)
            {
                FatalErrorInFunction
                    << "step size underflow"
                    << exit(FatalError);
            }

            SIMPR(x, yTemp_, li, dydx0_, dfdx_, dfdy_, h, nSeq_[k], ySeq_);
            scalar xest = sqr(h/nSeq_[k]);

            polyExtrapolate(k, xest, ySeq_, y, yErr_, x_p_, d_p_);

            if (k != 0)
            {
                maxErr = small;
                for (label i=0; i<n_; i++)
                {
                    maxErr = max
                    (
                        maxErr,
                        mag(yErr_[i])/(absTol_[i] + relTol_[i]*mag(yTemp_[i]))
                    );
                }
                km = k - 1;
                err_[km] = pow(maxErr/safe1, 1.0/(2*km + 3));
            }

            if (k != 0 && (k >= kOpt_ - 1 || first_))
            {
                if (maxErr < 1.0)
                {
                    exitflag = true;
                    break;
                }

                if (k == kMax_ || k == kOpt_ + 1)
                {
                    red = safe2/err_[km];
                    break;
                }
                else if (k == kOpt_ && alpha_[kOpt_ - 1][kOpt_] < err_[km])
                {
                    red = 1.0/err_[km];
                    break;
                }
                else if (kOpt_ == kMax_ && alpha_[km][kMax_ - 1] < err_[km])
                {
                    red = alpha_[km][kMax_ - 1]*safe2/err_[km];
                    break;
                }
                else if (alpha_[km][kOpt_] < err_[km])
                {
                    red = alpha_[km][kOpt_ - 1]/err_[km];
                    break;
                }
            }
        }

        if (exitflag)
        {
            break;
        }

        red = min(red, redMin);
        red = max(red, redMax);
        h *= red;
        reduct = 1;
    }

    x = xNew_;
    first_ = 0;
    scalar wrkmin = great;
    scalar scale = 1.0;

    for (label kk=0; kk<=km; kk++)
    {
        scalar fact = max(err_[kk], scaleMX);
        scalar work = fact*a_[kk + 1];
        if (work < wrkmin)
        {
            scale = fact;
            wrkmin = work;
            kOpt_ = kk + 1;
        }
    }

    dxTry = h/scale;

    if (kOpt_ >= k && kOpt_ != kMax_ && !reduct)
    {
        scalar fact = max(scale/alpha_[kOpt_ - 1][kOpt_], scaleMX);
        if (a_[kOpt_ + 1]*fact <= wrkmin)
        {
            dxTry = h/fact;
            kOpt_++;
        }
    }
}


// ************************************************************************* //
