/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

Global
    Foam::incGamma

Description
    Calculates the upper and lower incomplete gamma functions as well as their
    normalized versions.

    The algorithm is described in detail in DiDonato et al. (1986).

    \verbatim
        DiDonato, A. R., & Morris Jr, A. H. (1986).
        Computation of the incomplete gamma function ratios and their inverse.
        ACM Transactions on Mathematical Software (TOMS), 12(4), 377-393.
    \endverbatim

    All equation numbers in the following code refer to the above paper.
    The algorithm in function 'incGammaRatio_Q' is described in section 3.
    The accuracy parameter IND is set to a value of 1.

\*---------------------------------------------------------------------------*/

#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Eqn. (13)
static scalar calcQE11(const scalar a, const scalar x, const int e = 30)
{
    scalar a_2n = 0;
    scalar b_2n = 1;

    scalar a_2np1 = 1;
    scalar b_2np1 = x;

    int n = 1;
    for (n = 1; (2*n) <= e; n++)
    {
        const scalar a_2nm1 = a_2np1;
        const scalar b_2nm1 = b_2np1;

        a_2n = a_2nm1 + (n - a)*a_2n;
        b_2n = b_2nm1 + (n - a)*b_2n;

        a_2np1 = x*a_2n + n*a_2nm1;
        b_2np1 = x*b_2n + n*b_2nm1;
    }

    if (2*(n - 1) < e)
    {
        return a_2np1/b_2np1;
    }
    else
    {
        return a_2n/b_2n;
    }
}


// Eqn. (15)
static scalar calcPE15(const scalar a, const scalar x, const int nmax = 20)
{
    scalar prod = 1;
    scalar sum = 0;

    for (int n = 1; n <= nmax; n++)
    {
        prod *= (a + n);
        sum += pow(x, n)/prod;
    }

    const scalar R = (exp(-x)*pow(x, a))/tgamma(a);

    return R/a*(1 + sum);
}


// Eq. (16)
static scalar calcQE16(const scalar a, const scalar x, const int N = 20)
{
    scalar an = 1;
    scalar sum = 0;

    for (int n = 1; n <= (N - 1); n++)
    {
        an *= (a - n);
        sum += an/pow(x, n);
    }

    const scalar R = (exp(-x)*pow(x, a))/tgamma(a);

    return R/x*(1 + sum);
}


// Eq. (18)
static scalar calcTE18
(
    const scalar a,
    const scalar e0,
    const scalar x,
    const scalar lambda,
    const scalar sigma,
    const scalar phi
)
{
    static const scalar D0[] =
    {
       -0.333333333333333E-00,
        0.833333333333333E-01,
       -0.148148148148148E-01,
        0.115740740740741E-02,
        0.352733686067019E-03,
       -0.178755144032922E-03,
        0.391926317852244E-04,
       -0.218544851067999E-05,
       -0.185406221071516E-05,
        0.829671134095309E-06,
       -0.176659527368261E-06,
        0.670785354340150E-08,
        0.102618097842403E-07,
       -0.438203601845335E-08
    };

    static const scalar D1[] =
    {
       -0.185185185185185E-02,
       -0.347222222222222E-02,
        0.264550264550265E-02,
       -0.990226337448560E-03,
        0.205761316872428E-03,
       -0.401877572016461E-06,
       -0.180985503344900E-04,
        0.764916091608111E-05,
       -0.161209008945634E-05,
        0.464712780280743E-08,
        0.137863344691572E-06,
       -0.575254560351770E-07,
        0.119516285997781E-07
    };

    static const scalar D2[] =
    {
        0.413359788359788E-02,
       -0.268132716049383E-02,
        0.771604938271605E-03,
        0.200938786008230E-05,
       -0.107366532263652E-03,
        0.529234488291201E-04,
       -0.127606351886187E-04,
        0.342357873409614E-07,
        0.137219573090629E-05,
       -0.629899213838006E-06,
        0.142806142060642E-06
    };

    const scalar u = 1/a;
    scalar z = sqrt(2*phi);

    if (lambda < 1)
    {
        z = -z;
    }

    if (sigma > (e0/sqrt(a)))
    {
        const scalar C0 =
            D0[6]*pow6(z) + D0[5]*pow5(z) + D0[4]*pow4(z)
          + D0[3]*pow3(z) + D0[2]*sqr(z) + D0[1]*z + D0[0];

        const scalar C1 =
            D1[4]*pow4(z) + D1[3]*pow3(z) + D1[2]*sqr(z) + D1[1]*z + D1[0];

        const scalar C2 = D2[1]*z + D2[0];

        return C2*sqr(u) + C1*u + C0;
    }
    else
    {
        const scalar C0 = D0[2]*sqr(z) + D0[1]*z + D0[0];
        const scalar C1 = D1[1]*z + D1[0];
        const scalar C2 = D2[1]*z + D2[0];

        return C2*sqr(u) + C1*u + C0;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}  // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::scalar Foam::incGammaRatio_Q(const scalar a, const scalar x)
{
    const scalar BIG = 14;
    const scalar x0 = 17;
    const scalar e0 = 0.025;

    if (a < 1)
    {
        if (a == 0.5)
        {
            // Eqn. (8)
            if (x < 0.25)
            {
                return 1 - erf(sqrt(x));
            }
            else
            {
                return erfc(sqrt(x));
            }
        }
        else if ( x < 1.1)
        {
            // Eqn. (12)
            scalar alpha = x/2.59;

            if (x < 0.5)
            {
                alpha = log(sqrt(0.765))/log(x);
            }

            scalar sum = 0;

            for (int n = 1; n <= 10; n++)
            {
                sum += pow((-x), n)/((a + n)*factorial(n));
            }

            const scalar J = -a*sum;

            if (a > alpha || a == alpha)
            {
                // Eqn. (9)
                return 1 - (pow(x, a)*(1 - J))/tgamma(a + 1);
            }
            else
            {
                // Eqn. (10)
                const scalar L = exp(a*log(x)) - 1;
                const scalar H = 1/(tgamma(a + 1)) - 1;

                return (pow(x, a)*J - L)/tgamma(a + 1) - H;
            }
        }
        else
        {
            // Eqn. (11)
            const scalar R = (exp(-x)*pow(x, a))/tgamma(a);

            return R*calcQE11(a, x);
        }
    }
    else if (a >= BIG)
    {
        const scalar sigma = fabs(1 - x/a);

        if (sigma <= e0/sqrt(a))
        {
            // Eqn. (19)
            const scalar lambda = x/a;
            const scalar phi = lambda - 1 - log(lambda);
            const scalar y = a*phi;

            const scalar E = 0.5 - (1 - y/3)*sqrt(y/pi);

            if (lambda <= 1)
            {
                return
                    1
                  - (
                        E
                      - (1 - y)/sqrt(2*pi*a)
                       *calcTE18(a, e0, x, lambda, sigma, phi)
                    );
            }
            else
            {
                return
                    E
                  + (1 - y)/sqrt(2*pi*a)
                   *calcTE18(a, e0, x, lambda, sigma, phi);
            }
        }
        else
        {
            if (sigma <= 0.4)
            {
                // Eqn. (17)
                const scalar lambda = x/a;
                const scalar phi = lambda - 1 - log(lambda);
                const scalar y = a*phi;

                if (lambda <= 1)
                {
                    return
                        1
                      - (0.5*erfc(sqrt(y))
                      - exp(-y)/sqrt(2*pi*a)
                       *calcTE18(a, e0, x, lambda, sigma, phi));
                }
                else
                {
                    return
                        0.5*erfc(sqrt(y))
                      + exp(-y)/sqrt(2*pi*a)
                       *calcTE18(a, e0, x, lambda, sigma, phi);
                }
            }
            else
            {
                if (x <= max(a, log(10.0)))
                {
                    // Eqn. (15)
                    return 1 - calcPE15(a, x);
                }
                else if (x < x0)
                {
                    // Eqn. (11)
                    const scalar R = (exp(-x)*pow(x, a))/tgamma(a);

                    return R*calcQE11(a, x);
                }
                else
                {
                    // Eqn. (16)
                    return calcQE16(a, x);
                }
            }
        }
    }
    else
    {
        if (a > x || x >= x0)
        {
            if (x <= max(a, log(10.0)))
            {
                // Eqn. (15)
                return 1 - calcPE15(a, x);
            }
            else if ( x < x0)
            {
                // Eqn. (11)
                const scalar R = (exp(-x)*pow(x, a))/tgamma(a);

                return R*calcQE11(a, x);
            }
            else
            {
                // Eqn. (16)
                return calcQE16(a, x);
            }
        }
        else
        {
            if (floor(2*a) == 2*a)
            {
                // Eqn. (14)
                if (floor(a) == a)
                {
                    scalar sum = 0;

                    for (int n = 0; n <= (a - 1); n++)
                    {
                        sum += pow(x, n)/factorial(n);
                    }

                    return exp(-x)*sum;
                }
                else
                {
                    int i = a - 0.5;
                    scalar prod = 1;
                    scalar sum = 0;

                    for (int n = 1; n <= i; n++)
                    {
                        prod *= (n - 0.5);
                        sum += pow(x, n)/prod;
                    }

                    return erfc(sqrt(x)) + exp(-x)/sqrt(pi*x)*sum;
                }
            }
            else if (x <= max(a, log(10.0)))
            {
                // Eqn. (15)
                return 1 - calcPE15(a, x);
            }
            else if ( x < x0)
            {
                // Eqn. (11)
                const scalar R = (exp(-x)*pow(x, a))/tgamma(a);

                return R*calcQE11(a, x);
            }
            else
            {
                // Eqn. (16)
                return calcQE16(a, x);
            }
        }
    }
}


Foam::scalar Foam::incGammaRatio_P(const scalar a, const scalar x)
{
    return 1 - incGammaRatio_Q(a, x);
}


Foam::scalar Foam::incGamma_Q(const scalar a, const scalar x)
{
    return incGammaRatio_Q(a, x)*tgamma(a);
}


Foam::scalar Foam::incGamma_P(const scalar a, const scalar x)
{
    return incGammaRatio_P(a, x)*tgamma(a);
}


// ************************************************************************* //
