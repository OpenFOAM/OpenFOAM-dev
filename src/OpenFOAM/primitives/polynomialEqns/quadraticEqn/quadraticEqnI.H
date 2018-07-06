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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::quadraticEqn::quadraticEqn()
{}


inline Foam::quadraticEqn::quadraticEqn(const Foam::zero)
:
    VectorSpace<quadraticEqn, scalar, 3>(Foam::zero())
{}


inline Foam::quadraticEqn::quadraticEqn
(
    const scalar a,
    const scalar b,
    const scalar c
)
{
    this->v_[A] = a;
    this->v_[B] = b;
    this->v_[C] = c;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::quadraticEqn::a() const
{
    return this->v_[A];
}


inline Foam::scalar Foam::quadraticEqn::b() const
{
    return this->v_[B];
}


inline Foam::scalar Foam::quadraticEqn::c() const
{
    return this->v_[C];
}


inline Foam::scalar& Foam::quadraticEqn::a()
{
    return this->v_[A];
}


inline Foam::scalar& Foam::quadraticEqn::b()
{
    return this->v_[B];
}


inline Foam::scalar& Foam::quadraticEqn::c()
{
    return this->v_[C];
}


inline Foam::scalar Foam::quadraticEqn::value(const scalar x) const
{
    return x*(x*a() + b()) + c();
}


inline Foam::scalar Foam::quadraticEqn::derivative(const scalar x) const
{
    return x*2*a() + b();
}


inline Foam::scalar Foam::quadraticEqn::error(const scalar x) const
{
    return mag(small*x*(x*2*a() + b()));
}


// ************************************************************************* //
