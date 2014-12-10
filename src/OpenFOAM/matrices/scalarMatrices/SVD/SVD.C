/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "SVD.H"
#include "scalarList.H"
#include "scalarMatrices.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * Constructor  * * * * * * * * * * * * * * //

Foam::SVD::SVD(const scalarRectangularMatrix& A, const scalar minCondition)
:
    U_(A),
    V_(A.m(), A.m()),
    S_(A.m()),
    VSinvUt_(A.m(), A.n()),
    nZeros_(0)
{
    // SVDcomp to find U_, V_ and S_ - the singular values

    const label Um = U_.m();
    const label Un = U_.n();

    scalarList rv1(Um);
    scalar g = 0;
    scalar scale = 0;
    scalar s = 0;
    scalar anorm = 0;
    label l = 0;

    for (label i = 0; i < Um; i++)
    {
        l = i+2;
        rv1[i] = scale*g;
        g = s = scale = 0;

        if (i < Un)
        {
            for (label k = i; k < Un; k++)
            {
                scale += mag(U_[k][i]);
            }

            if (scale != 0)
            {
                for (label k = i; k < Un; k++)
                {
                    U_[k][i] /= scale;
                    s += U_[k][i]*U_[k][i];
                }

                scalar f = U_[i][i];
                g = -sign(Foam::sqrt(s), f);
                scalar h = f*g - s;
                U_[i][i] = f - g;

                for (label j = l-1; j < Um; j++)
                {
                    s = 0;
                    for (label k = i; k < Un; k++)
                    {
                        s += U_[k][i]*U_[k][j];
                    }

                    f = s/h;
                    for (label k = i; k < A.n(); k++)
                    {
                        U_[k][j] += f*U_[k][i];
                    }
                }

                for (label k = i; k < Un; k++)
                {
                    U_[k][i] *= scale;
                }
            }
        }

        S_[i] = scale*g;

        g = s = scale = 0;

        if (i+1 <= Un && i != Um)
        {
            for (label k = l-1; k < Um; k++)
            {
                scale += mag(U_[i][k]);
            }

            if (scale != 0)
            {
                for (label k=l-1; k < Um; k++)
                {
                    U_[i][k] /= scale;
                    s += U_[i][k]*U_[i][k];
                }

                scalar f = U_[i][l-1];
                g = -sign(Foam::sqrt(s),f);
                scalar h = f*g - s;
                U_[i][l-1] = f - g;

                for (label k = l-1; k < Um; k++)
                {
                    rv1[k] = U_[i][k]/h;
                }

                for (label j = l-1; j < Un; j++)
                {
                    s = 0;
                    for (label k = l-1; k < Um; k++)
                    {
                        s += U_[j][k]*U_[i][k];
                    }

                    for (label k = l-1; k < Um; k++)
                    {
                        U_[j][k] += s*rv1[k];
                    }
                }
                for (label k = l-1; k < Um; k++)
                {
                    U_[i][k] *= scale;
                }
            }
        }

        anorm = max(anorm, mag(S_[i]) + mag(rv1[i]));
    }

    for (label i = Um-1; i >= 0; i--)
    {
        if (i < Um-1)
        {
            if (g != 0)
            {
                for (label j = l; j < Um; j++)
                {
                    V_[j][i] = (U_[i][j]/U_[i][l])/g;
                }

                for (label j=l; j < Um; j++)
                {
                    s = 0;
                    for (label k = l; k < Um; k++)
                    {
                        s += U_[i][k]*V_[k][j];
                    }

                    for (label k = l; k < Um; k++)
                    {
                        V_[k][j] += s*V_[k][i];
                    }
                }
            }

            for (label j = l; j < Um;j++)
            {
                V_[i][j] = V_[j][i] = 0.0;
            }
        }

        V_[i][i] = 1;
        g = rv1[i];
        l = i;
    }

    for (label i = min(Um, Un) - 1; i >= 0; i--)
    {
        l = i+1;
        g = S_[i];

        for (label j = l; j < Um; j++)
        {
            U_[i][j] = 0.0;
        }

        if (g != 0)
        {
            g = 1.0/g;

            for (label j = l; j < Um; j++)
            {
                s = 0;
                for (label k = l; k < Un; k++)
                {
                    s += U_[k][i]*U_[k][j];
                }

                scalar f = (s/U_[i][i])*g;

                for (label k = i; k < Un; k++)
                {
                    U_[k][j] += f*U_[k][i];
                }
            }

            for (label j = i; j < Un; j++)
            {
                U_[j][i] *= g;
            }
        }
        else
        {
            for (label j = i; j < Un; j++)
            {
                U_[j][i] = 0.0;
            }
        }

        ++U_[i][i];
    }

    for (label k = Um-1; k >= 0; k--)
    {
        for (label its = 0; its < 35; its++)
        {
            bool flag = true;

            label nm;
            for (l = k; l >= 0; l--)
            {
                nm = l-1;
                if (mag(rv1[l]) + anorm == anorm)
                {
                    flag = false;
                    break;
                }
                if (mag(S_[nm]) + anorm == anorm) break;
            }

            if (flag)
            {
                scalar c = 0.0;
                s = 1.0;
                for (label i = l; i < k+1; i++)
                {
                    scalar f = s*rv1[i];
                    rv1[i] = c*rv1[i];

                    if (mag(f) + anorm == anorm) break;

                    g = S_[i];
                    scalar h = sqrtSumSqr(f, g);
                    S_[i] = h;
                    h = 1.0/h;
                    c = g*h;
                    s = -f*h;

                    for (label j = 0; j < Un; j++)
                    {
                        scalar y = U_[j][nm];
                        scalar z = U_[j][i];
                        U_[j][nm] = y*c + z*s;
                        U_[j][i] = z*c - y*s;
                    }
                }
            }

            scalar z = S_[k];

            if (l == k)
            {
                if (z < 0.0)
                {
                    S_[k] = -z;

                    for (label j = 0; j < Um; j++) V_[j][k] = -V_[j][k];
                }
                break;
            }
            if (its == 34)
            {
                WarningIn
                (
                    "SVD::SVD"
                    "(scalarRectangularMatrix& A, const scalar minCondition)"
                )   << "no convergence in 35 SVD iterations"
                    << endl;
            }

            scalar x = S_[l];
            nm = k-1;
            scalar y = S_[nm];
            g = rv1[nm];
            scalar h = rv1[k];
            scalar f = ((y - z)*(y + z) + (g - h)*(g + h))/(2.0*h*y);
            g = sqrtSumSqr(f, scalar(1));
            f = ((x - z)*(x + z) + h*((y/(f + sign(g, f))) - h))/x;
            scalar c = 1.0;
            s = 1.0;

            for (label j = l; j <= nm; j++)
            {
                label i = j + 1;
                g = rv1[i];
                y = S_[i];
                h = s*g;
                g = c*g;
                scalar z = sqrtSumSqr(f, h);
                rv1[j] = z;
                c = f/z;
                s = h/z;
                f = x*c + g*s;
                g = g*c - x*s;
                h = y*s;
                y *= c;

                for (label jj = 0; jj < Um; jj++)
                {
                    x = V_[jj][j];
                    z = V_[jj][i];
                    V_[jj][j] = x*c + z*s;
                    V_[jj][i] = z*c - x*s;
                }

                z = sqrtSumSqr(f, h);
                S_[j] = z;
                if (z)
                {
                    z = 1.0/z;
                    c = f*z;
                    s = h*z;
                }
                f = c*g + s*y;
                x = c*y - s*g;

                for (label jj=0; jj < Un; jj++)
                {
                    y = U_[jj][j];
                    z = U_[jj][i];
                    U_[jj][j] = y*c + z*s;
                    U_[jj][i] = z*c - y*s;
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            S_[k] = x;
        }
    }

    // zero singular values that are less than minCondition*maxS
    const scalar minS = minCondition*S_[findMax(S_)];
    forAll(S_, i)
    {
        if (S_[i] <= minS)
        {
            //Info<< "Removing " << S_[i] << " < " << minS << endl;
            S_[i] = 0;
            nZeros_++;
        }
    }

    // now multiply out to find the pseudo inverse of A, VSinvUt_
    multiply(VSinvUt_, V_, inv(S_), U_.T());

    // test SVD
    /*scalarRectangularMatrix SVDA(A.n(), A.m());
    multiply(SVDA, U_, S_, transpose(V_));
    scalar maxDiff = 0;
    scalar diff = 0;
    for (label i = 0; i < A.n(); i++)
    {
        for (label j = 0; j < A.m(); j++)
        {
            diff = mag(A[i][j] - SVDA[i][j]);
            if (diff > maxDiff) maxDiff = diff;
        }
    }
    Info<< "Maximum discrepancy between A and svd(A) = " << maxDiff << endl;

    if (maxDiff > 4)
    {
        Info<< "singular values " << S_ << endl;
    }
    */
}


// ************************************************************************* //
