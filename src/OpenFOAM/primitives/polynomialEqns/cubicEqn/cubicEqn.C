/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "linearEqn.H"
#include "quadraticEqn.H"
#include "cubicEqn.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Roots<3> Foam::cubicEqn::roots() const
{
    /*

    This function solves a cubic equation of the following form:

        a*x^3 + b*x^2 + c*x + d = 0
          x^3 + B*x^2 + C*x + D = 0

    The following two substitutions are used:

        x = t - B/3
        t = w - P/3/w

    This reduces the problem to a quadratic in w^3.

        w^6 + Q*w^3 - P = 0

    Where Q and P are given in the code below.

    The properties of the cubic can be related to the properties of this
    quadratic in w^3. If it has a repeated root a zero, the cubic has a tripl
    root. If it has a repeated root not at zero, the cubic has two real roots,
    one repeated and one not. If it has two complex roots, the cubic has three
    real roots. If it has two real roots, then the cubic has one real root and
    two complex roots.

    This is solved for the most numerically accurate value of w^3. See the
    quadratic function for details on how to pick a value. This single value of
    w^3 can yield up to three cube roots for w, which relate to the three
    solutions for x.

    Only a single root, or pair of conjugate roots, is directly evaluated; the
    one, or ones with the lowest relative numerical error. Root identities are
    then used to recover the remaining roots, possibly utilising a quadratic
    and/or linear solution. This seems to be a good way of maintaining the
    accuracy of roots at very different magnitudes.

    */

    const scalar a = this->a();
    const scalar b = this->b();
    const scalar c = this->c();
    const scalar d = this->d();

    if (a == 0)
    {
        return Roots<3>(quadraticEqn(b, c, d).roots(), rootType::nan, 0);
    }

    // This is assumed not to over- or under-flow. If it does, all bets are off.
    const scalar p = c*a - b*b/3;
    const scalar q = b*b*b*scalar(2)/27 - b*c*a/3 + d*a*a;
    const scalar disc = p*p*p/27 + q*q/4;

    // How many roots of what types are available?
    const bool oneReal = disc == 0 && p == 0;
    const bool twoReal = disc == 0 && p != 0;
    const bool threeReal = disc < 0;
    // const bool oneRealTwoComplex = disc > 0;

    static const scalar sqrt3 = sqrt(3.0);

    scalar x;

    if (oneReal)
    {
        const Roots<1> r = linearEqn(a, b/3).roots();
        return Roots<3>(r.type(0), r[0]);
    }
    else if (twoReal)
    {
        if (q*b > 0)
        {
            x = - 2*cbrt(q/2) - b/3;
        }
        else
        {
            x = cbrt(q/2) - b/3;
            const Roots<1> r = linearEqn(- a, x).roots();
            return Roots<3>(Roots<2>(r, r), linearEqn(x*x, a*d).roots());
        }
    }
    else if (threeReal)
    {
        const scalar wCbRe = - q/2, wCbIm = sqrt(- disc);
        const scalar wAbs = cbrt(hypot(wCbRe, wCbIm));
        const scalar wArg = atan2(wCbIm, wCbRe)/3;
        const scalar wRe = wAbs*cos(wArg), wIm = wAbs*sin(wArg);
        if (b > 0)
        {
            x = - wRe - mag(wIm)*sqrt3 - b/3;
        }
        else
        {
            x = 2*wRe - b/3;
        }
    }
    else // if (oneRealTwoComplex)
    {
        const scalar wCb = - q/2 - sign(q)*sqrt(disc);
        const scalar w = cbrt(wCb);
        const scalar t = w - p/(3*w);
        if (p + t*b < 0)
        {
            x = t - b/3;
        }
        else
        {
            const scalar xRe = - t/2 - b/3, xIm = sqrt3/2*(w + p/3/w);
            x = - a*a*d/(xRe*xRe + xIm*xIm);

            // This form of solving for the remaining roots seems more stable
            // for this case. This has been determined by trial and error.
            return
                Roots<3>
                (
                    linearEqn(- a, x).roots(),
                    quadraticEqn(a*x, x*x + b*x, - a*d).roots()
                );
        }
    }

    return
        Roots<3>
        (
            linearEqn(- a, x).roots(),
            quadraticEqn(- x*x, c*x + a*d, d*x).roots()
        );
}


// ************************************************************************* //
