/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "Random.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::Random::scalarNormal()
{
    // Proper inversion of the distribution. Slow. Exactly maintains
    // the random behaviour of the generator.

    /*
    using namespace constant::mathematical;

    static const scalar sqrtTwo = sqrt(scalar(2));
    static const scalar sqrtPiByTwo = sqrt(pi)/2;
    static const scalar a = 8*(pi - 3)/(3*pi*(4 - pi));

    const scalar x = 2*scalar01() - 1;
    const scalar xPos = mag(x);

    // Initial approximation
    const scalar l = log(1 - sqr(xPos));
    const scalar ll = 2/(pi*a) + l/2;
    scalar y = sqrt(sqrt(sqr(ll) - l/a) - ll);

    // Newton improvement
    label n = 0;
    while (n < 2)
    {
        const scalar dt = (erf(y) - xPos)/exp(- y*y)*sqrtPiByTwo;
        y -= dt;
        n += mag(dt) < rootSmall;
    }

    return sign(x)*sqrtTwo*y;
    */

    // Box-Muller transform. Fast. Uses rejection and caching so the
    // random sequence is not guaranteed.

    if (scalarNormalStored_)
    {
        scalarNormalStored_ = false;

        return scalarNormalValue_;
    }
    else
    {
        scalar x1, x2, rr;

        do
        {
            x1 = 2*scalar01() - 1;
            x2 = 2*scalar01() - 1;
            rr = sqr(x1) + sqr(x2);
        }
        while (rr >= 1 || rr == 0);

        const scalar f = sqrt(- 2*log(rr)/rr);

        scalarNormalValue_ = x1*f;
        scalarNormalStored_ = true;

        return x2*f;
    }
}


Foam::scalar Foam::Random::globalScalar01()
{
    scalar value = - vGreat;

    if (Pstream::master())
    {
        value = scalar01();
    }

    Pstream::scatter(value);

    return value;
}


// ************************************************************************* //
