/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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
    Foam::invIncGamma

Description
    Function to calculate the inverse of the
    normalised incomplete gamma function.

    The algorithm is described in detain in reference:
    \verbatim
        DiDonato, A. R., & Morris Jr, A. H. (1986).
        Computation of the incomplete gamma function ratios and their inverse.
        ACM Transactions on Mathematical Software (TOMS), 12(4), 377-393.
    \endverbatim

    All equation numbers in the following code refer to the above text and
    the algorithm in the function 'invIncGamma' is described in section 4.

\*---------------------------------------------------------------------------*/

#include "mathematicalConstants.H"
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

static scalar minimaxs(const scalar P)
{
    // Eqn. 32

    static const scalar a[] =
    {
        3.31125922108741,
        11.6616720288968,
        4.28342155967104,
        0.213623493715853
    };

    static const scalar b[] =
    {
        6.61053765625462,
        6.40691597760039,
        1.27364489782223,
        0.03611708101884203
    };

    const scalar t = P < 0.5 ? sqrt(-2*log(P)) : sqrt(-2*log(1 - P));

    const scalar s =
        t
      - (a[0] + t*(a[1] + t*(a[2] + t*a[3])))
       /(1 + t*(b[0] + t*(b[1] + t*(b[2] + t*b[3]))));

    return P < 0.5 ? -s : s;
}


static scalar Sn(const scalar a, const scalar x)
{
    //- Eqn. 34

    scalar Sn = 1;
    scalar Si = 1;

    for (int i=1; i<100; i++)
    {
        Si *= x/(a + i);
        Sn += Si;

        if (Si < 1e-4) break;
    }

    return Sn;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Inverse normalised incomplete gamma function
Foam::scalar Foam::invIncGamma(const scalar a, const scalar P)
{
    const scalar Q = 1 - P;

    if (a == 1)
    {
        return -log(Q);
    }
    else if (a < 1)
    {
        const scalar Ga = tgamma(a);
        const scalar B = Q*Ga;

        if (B > 0.6 || (B >= 0.45 && a >= 0.3))
        {
            // Eqn. 21
            const scalar u =
                (B*Q > 1e-8) ? pow(P*Ga*a, 1/a) : exp((-Q/a) - Eu);

            return u/(1 - (u/(a + 1)));
        }
        else if (a < 0.3 && B >= 0.35)
        {
            // Eqn. 22
            const scalar t = exp(-Eu - B);
            const scalar u = t*exp(t);

            return t*exp(u);
        }
        else if (B > 0.15 || a >= 0.3)
        {
            // Eqn. 23
            const scalar y = -log(B);
            const scalar u = y - (1 - a)*log(y);

            return y - (1 - a)*log(u) - log(1 + (1 - a)/(1 + u));
        }
        else if (B > 0.1)
        {
            // Eqn. 24
            const scalar y = -log(B);
            const scalar u = y - (1 - a)*log(y);

            return y
              - (1 - a)*log(u)
              - log
                (
                    (sqr(u) + 2*(3 - a)*u + (2 - a)*(3 - a))
                   /(sqr(u) + (5 - a)*u + 2)
                );
        }
        else
        {
            // Eqn. 25:
            const scalar y = -log(B);
            const scalar c1 = (a - 1)*log(y);
            const scalar c12 = c1*c1;
            const scalar c13 = c12*c1;
            const scalar c14 = c12*c12;
            const scalar a2 = a*a;
            const scalar a3 = a2*a;
            const scalar c2 = (a - 1)*(1 + c1);
            const scalar c3 = (a - 1)*(-(c12/2) + (a - 2)*c1 + (3*a - 5)/2);
            const scalar c4 =
                (a - 1)
               *(
                    (c13/3)
                  - (3*a - 5)*c12/2
                  + (a2 - 6*a + 7)*c1
                  + (11*a2 - 46*a + 47)/6
                );
            const scalar c5 =
                (a - 1)*(-(c14/4)
              + (11*a - 17)*c13/6
              + (-3*a2 + 13*a - 13)*c12
              + (2*a3 - 25*a2 + 72*a - 61)*c1/2
              + (25*a3 - 195*a2 + 477*a - 379)/12);
            const scalar y2 = y*y;
            const scalar y3 = y2*y;
            const scalar y4 = y2*y2;

            return y + c1 + (c2/y) + (c3/y2) + (c4/y3) + (c5/y4);
        }
    }
    else
    {
        // Eqn. 31:
        scalar s = minimaxs(P);

        const scalar s2 = sqr(s);
        const scalar s3 = s*s2;
        const scalar s4 = s2*s2;
        const scalar s5 = s*s4;
        const scalar sqrta = sqrt(a);

        const scalar w =
            a + s*sqrta + (s2 - 1)/3
          + (s3 - 7*s)/(36*sqrta)
          - (3*s4 + 7*s2 - 16)/(810*a)
          + (9*s5 + 256*s3 - 433*s)/(38880*a*sqrta);

        if (a >= 500 && mag(1 - w/a) < 1e-6)
        {
            return w;
        }
        else if (P > 0.5)
        {
            if (w < 3*a)
            {
                return w;
            }
            else
            {
                const scalar D = max(scalar(2), scalar(a*(a - 1)));
                const scalar lnGa = lgamma(a);
                const scalar lnB = log(Q) + lnGa;

                if (lnB < -2.3*D)
                {
                    // Eqn. 25:
                    const scalar y = -lnB;
                    const scalar c1 = (a - 1)*log(y);
                    const scalar c12 = c1*c1;
                    const scalar c13 = c12*c1;
                    const scalar c14 = c12*c12;
                    const scalar a2 = a*a;
                    const scalar a3 = a2*a;

                    const scalar c2 = (a - 1)*(1 + c1);
                    const scalar c3 =
                        (a - 1)
                       *(
                          - (c12/2)
                          + (a - 2)*c1
                          + (3*a - 5)/2
                        );
                    const scalar c4 =
                        (a - 1)
                       *(
                            (c13/3)
                          - (3*a - 5)*c12/2
                          + (a2 - 6*a + 7)*c1
                          + (11*a2 - 46*a + 47)/6
                        );
                    const scalar c5 =
                        (a - 1)
                       *(
                          - (c14/4)
                          + (11*a - 17)*c13/6
                          + (-3*a2 + 13*a - 13)*c12
                          + (2*a3 - 25*a2 + 72*a - 61)*c1/2
                          + (25*a3 - 195*a2 + 477*a - 379)/12
                        );

                    const scalar y2 = y*y;
                    const scalar y3 = y2*y;
                    const scalar y4 = y2*y2;

                    return y + c1 + (c2/y) + (c3/y2) + (c4/y3) + (c5/y4);
                }
                else
                {
                    // Eqn. 33
                    const scalar u =
                        -lnB + (a - 1)*log(w) - log(1 + (1 - a)/(1 + w));

                    return -lnB + (a - 1)*log(u) - log(1 + (1 - a)/(1 + u));
                }
            }
        }
        else
        {
            scalar z = w;
            const scalar ap1 = a + 1;

            if (w < 0.15*ap1)
            {
                // Eqn. 35
                const scalar ap2 = a + 2;
                const scalar v = log(P) + lgamma(ap1);
                z = exp((v + w)/a);
                s = log1p(z/ap1*(1 + z/ap2));
                z = exp((v + z - s)/a);
                s = log1p(z/ap1*(1 + z/ap2));
                z = exp((v + z - s)/a);
                s = log1p(z/ap1*(1 + z/ap2*(1 + z/(a + 3))));
                z = exp((v + z - s)/a);
            }

            if (z <= 0.01*ap1 || z > 0.7*ap1)
            {
                return z;
            }
            else
            {
                // Eqn. 36
                const scalar lnSn = log(Sn(a, z));
                const scalar v = log(P) + lgamma(ap1);
                z = exp((v + z - lnSn)/a);

                return z*(1 - (a*log(z) - z - v + lnSn)/(a - z));
            }
        }
    }
}


// ************************************************************************* //
