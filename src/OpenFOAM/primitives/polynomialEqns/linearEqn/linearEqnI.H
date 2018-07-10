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

inline Foam::linearEqn::linearEqn()
{}


inline Foam::linearEqn::linearEqn(const Foam::zero)
:
    VectorSpace<linearEqn, scalar, 2>(Foam::zero())
{}


inline Foam::linearEqn::linearEqn(const scalar a, const scalar b)
{
    this->v_[A] = a;
    this->v_[B] = b;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::linearEqn::a() const
{
    return this->v_[A];
}


inline Foam::scalar Foam::linearEqn::b() const
{
    return this->v_[B];
}


inline Foam::scalar& Foam::linearEqn::a()
{
    return this->v_[A];
}


inline Foam::scalar& Foam::linearEqn::b()
{
    return this->v_[B];
}


inline Foam::scalar Foam::linearEqn::value(const scalar x) const
{
    return x*a() + b();
}


inline Foam::scalar Foam::linearEqn::derivative(const scalar x) const
{
    return a();
}


inline Foam::scalar Foam::linearEqn::error(const scalar x) const
{
    return mag(small*x*a());
}


inline Foam::Roots<1> Foam::linearEqn::roots() const
{
    /*

    This function solves a linear equation of the following form:

        a*x + b = 0
          x + B = 0

    */

    const scalar a = this->a();
    const scalar b = this->b();

    if (a == 0)
    {
        return Roots<1>(roots::nan, 0);
    }

    if (mag(b/vGreat) >= mag(a))
    {
        return Roots<1>(sign(a) == sign(b) ? roots::negInf : roots::posInf, 0);
    }

    return Roots<1>(roots::real, - b/a);
}


// ************************************************************************* //
