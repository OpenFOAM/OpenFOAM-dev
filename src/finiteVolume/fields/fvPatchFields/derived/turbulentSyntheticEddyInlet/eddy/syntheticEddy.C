/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "syntheticEddy.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::syntheticEddy::syntheticEddy()
:
    position0_(Zero),
    x_(0),
    sigma_(Zero),
    alpha_(Zero),
    Rpg_(tensor::I),
    c1_(-1),
    Uconv_(0)
{}


Foam::syntheticEddy::syntheticEddy
(
    const point& position0,
    scalar x,
    const vector& sigma,
    const vector& alpha,
    const tensor& Rpg,
    scalar c1,
    scalar Uconv
)
:
    position0_(position0),
    x_(x),
    sigma_(sigma),
    alpha_(alpha),
    Rpg_(Rpg),
    c1_(c1),
    Uconv_(Uconv)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::syntheticEddy::uPrime(const point& xp, const vector& n) const
{
    if (sigma_[2] <= 0)
    {
        return vector::zero;
    }

    vector u(vector::zero);

    const vector rp(cmptDivide(Rpg_.T() & (xp - position(n)), sigma_));
    if (mag(rp) < scalar(1))
    {
        const vector q(sigma_*(1 - magSqr(rp)));
        const vector uPrimep(cmptMultiply(q, rp^alpha_));
        u += c1_*(Rpg_ & uPrimep);
    }

    return u;
}


// ************************************************************************* //
