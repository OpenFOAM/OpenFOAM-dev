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

inline Foam::cubicEqn::cubicEqn()
{}


inline Foam::cubicEqn::cubicEqn(const Foam::zero)
:
    VectorSpace<cubicEqn, scalar, 4>(Foam::zero())
{}


inline Foam::cubicEqn::cubicEqn
(
    const scalar a,
    const scalar b,
    const scalar c,
    const scalar d
)
{
    this->v_[A] = a;
    this->v_[B] = b;
    this->v_[C] = c;
    this->v_[D] = d;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::cubicEqn::a() const
{
    return this->v_[A];
}


inline Foam::scalar Foam::cubicEqn::b() const
{
    return this->v_[B];
}


inline Foam::scalar Foam::cubicEqn::c() const
{
    return this->v_[C];
}


inline Foam::scalar Foam::cubicEqn::d() const
{
    return this->v_[D];
}


inline Foam::scalar& Foam::cubicEqn::a()
{
    return this->v_[A];
}


inline Foam::scalar& Foam::cubicEqn::b()
{
    return this->v_[B];
}


inline Foam::scalar& Foam::cubicEqn::c()
{
    return this->v_[C];
}


inline Foam::scalar& Foam::cubicEqn::d()
{
    return this->v_[D];
}


inline Foam::scalar Foam::cubicEqn::value(const scalar x) const
{
    return x*(x*(x*a() + b()) + c()) + d();
}


inline Foam::scalar Foam::cubicEqn::derivative(const scalar x) const
{
    return x*(x*3*a() + 2*b()) + c();
}


inline Foam::scalar Foam::cubicEqn::error(const scalar x) const
{
    return mag(small*x*(x*(x*3*a() + 2*b()) + c()));
}


// ************************************************************************* //
