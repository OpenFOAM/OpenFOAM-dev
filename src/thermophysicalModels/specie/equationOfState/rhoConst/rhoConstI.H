/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2023 OpenFOAM Foundation
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

#include "rhoConst.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::rhoConst<Specie>::rhoConst
(
    const Specie& sp,
    const scalar rho
)
:
    Specie(sp),
    rho_(rho)
{}


template<class Specie>
inline Foam::rhoConst<Specie>::rhoConst
(
    const word& name,
    const rhoConst<Specie>& rc
)
:
    Specie(name, rc),
    rho_(rc.rho_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::rhoConst<Specie>>
Foam::rhoConst<Specie>::clone() const
{
    return autoPtr<rhoConst<Specie>>(new rhoConst<Specie>(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::rhoConst<Specie>::rho(scalar p, scalar T) const
{
    return rho_;
}


template<class Specie>
inline Foam::scalar Foam::rhoConst<Specie>::h(scalar p, scalar T) const
{
    return p/rho_;
}


template<class Specie>
inline Foam::scalar Foam::rhoConst<Specie>::Cp(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rhoConst<Specie>::e(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rhoConst<Specie>::Cv(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rhoConst<Specie>::sp(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rhoConst<Specie>::sv(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rhoConst<Specie>::psi(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rhoConst<Specie>::Z(scalar p, scalar T) const
{
    return p/(rho(p, T)*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::rhoConst<Specie>::CpMCv(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rhoConst<Specie>::alphav(scalar p, scalar T) const
{
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::rhoConst<Specie>::operator+=(const rhoConst<Specie>& rc)
{
    const scalar Y1 = this->Y();
    Specie::operator+=(rc);

    if (mag(this->Y()) > small)
    {
        rho_ = this->Y()/(Y1/rho_ + rc.Y()/rc.rho_);
    }
}


template<class Specie>
inline void Foam::rhoConst<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::rhoConst<Specie> Foam::operator+
(
    const rhoConst<Specie>& rc1,
    const rhoConst<Specie>& rc2
)
{
    Specie sp
    (
        static_cast<const Specie&>(rc1)
      + static_cast<const Specie&>(rc2)
    );

    if (mag(sp.Y()) < small)
    {
        return rhoConst<Specie>
        (
            sp,
            rc1.rho_
        );
    }
    else
    {
        return rhoConst<Specie>
        (
            sp,
            sp.Y()/(rc1.Y()/rc1.rho_ + rc2.Y()/rc2.rho_)
        );
    }
}


template<class Specie>
inline Foam::rhoConst<Specie> Foam::operator*
(
    const scalar s,
    const rhoConst<Specie>& rc
)
{
    return rhoConst<Specie>(s*static_cast<const Specie&>(rc), rc.rho_);
}


template<class Specie>
inline Foam::rhoConst<Specie> Foam::operator==
(
    const rhoConst<Specie>& rc1,
    const rhoConst<Specie>& rc2
)
{
    return rhoConst<Specie>
    (
        static_cast<const Specie&>(rc1) == static_cast<const Specie&>(rc2),
        NaN
    );
}


// ************************************************************************* //
