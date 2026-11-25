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

#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::point& Foam::syntheticEddy::position0() const noexcept
{
    return position0_;
}


inline Foam::scalar Foam::syntheticEddy::x() const noexcept
{
    return x_;
}


inline Foam::scalar Foam::syntheticEddy::Uconv() const noexcept
{
    return Uconv_;
}


inline Foam::point Foam::syntheticEddy::position(const vector& n) const
{
    return position0_ + n*x_;
}


inline Foam::scalar Foam::syntheticEddy::volume() const
{
    return Foam::constant::mathematical::pi*4.0/3.0*cmptProduct(sigma_);
}


inline void Foam::syntheticEddy::move(Foam::scalar dx)
{
    x_ += dx;
}


// ************************************************************************* //
