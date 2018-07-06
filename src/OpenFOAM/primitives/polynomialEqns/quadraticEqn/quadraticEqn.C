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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Roots<2> Foam::quadraticEqn::roots() const
{
    /*

    This function solves a quadraticEqn equation of the following form:

        a*x^2 + b*x + c = 0
          x^2 + B*x + C = 0

    The quadraticEqn formula is as follows:

        x = - B/2 +- sqrt(B*B - 4*C)/2

    If the sqrt generates a complex number, this provides the result. If not
    then the real root with the smallest floating point error is calculated.

        x0 = - B/2 - sign(B)*sqrt(B*B - 4*C)/2

    The other root is the obtained using an identity.

        x1 = C/x0

    */

    const scalar a = this->a();
    const scalar b = this->b();
    const scalar c = this->c();

    if (a == 0)
    {
        return Roots<2>(linearEqn(b, c).roots(), roots::nan, 0);
    }

    // This is assumed not to over- or under-flow. If it does, all bets are off.
    const scalar disc = b*b/4 - a*c;

    // How many roots of what types are available?
    const bool oneReal = disc == 0;
    const bool twoReal = disc > 0;
    // const bool twoComplex = disc < 0;

    if (oneReal)
    {
        const Roots<1> r = linearEqn(- a, b/2).roots();
        return Roots<2>(r, r);
    }
    else if (twoReal)
    {
        const scalar x = - b/2 - sign(b)*sqrt(disc);
        return Roots<2>(linearEqn(- a, x).roots(), linearEqn(- x, c).roots());
    }
    else // if (twoComplex)
    {
        return Roots<2>(roots::complex, 0);
    }
}

// ************************************************************************* //
