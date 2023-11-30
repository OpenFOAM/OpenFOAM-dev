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

#include "perfectFluid.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::perfectFluid<Specie>::perfectFluid
(
    const Specie& sp,
    const scalar R,
    const scalar rho0
)
:
    Specie(sp),
    R_(R),
    rho0_(rho0)
{}


template<class Specie>
inline Foam::perfectFluid<Specie>::perfectFluid
(
    const word& name,
    const perfectFluid<Specie>& pf
)
:
    Specie(name, pf),
    R_(pf.R_),
    rho0_(pf.rho0_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::perfectFluid<Specie>>
Foam::perfectFluid<Specie>::clone() const
{
    return autoPtr<perfectFluid<Specie>>(new perfectFluid<Specie>(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::perfectFluid<Specie>::R() const
{
    return R_;
}


template<class Specie>
inline Foam::scalar Foam::perfectFluid<Specie>::rho(scalar p, scalar T) const
{
    return rho0_ + p/(this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::perfectFluid<Specie>::h(scalar p, scalar T) const
{
    return p/rho(p, T) - Pstd/rho(Pstd, T);
}


template<class Specie>
inline Foam::scalar Foam::perfectFluid<Specie>::Cp(scalar p, scalar T) const
{
    return
        sqr(p/(rho(p, T)*T))/this->R()
      - sqr(Pstd/(rho(Pstd, T)*T))/this->R();
}


template<class Specie>
inline Foam::scalar Foam::perfectFluid<Specie>::e(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::perfectFluid<Specie>::Cv(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::perfectFluid<Specie>::sp(scalar p, scalar T) const
{
    return -this->R()*log(p/Pstd);
}


template<class Specie>
inline Foam::scalar Foam::perfectFluid<Specie>::sv(scalar p, scalar T) const
{
    NotImplemented;
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::perfectFluid<Specie>::psi(scalar p, scalar T) const
{
    return 1.0/(this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::perfectFluid<Specie>::Z(scalar p, scalar T) const
{
    return p/(rho(p, T)*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::perfectFluid<Specie>::CpMCv(scalar p, scalar T) const
{
    const scalar R = this->R();
    const scalar rho = this->rho(p, T);

    return R*sqr(p/(rho*R*T));
}


template<class Specie>
inline Foam::scalar Foam::perfectFluid<Specie>::alphav(scalar p, scalar T) const
{
    return (1 - rho0_/this->rho(p, T))/T;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::perfectFluid<Specie>::operator+=
(
    const perfectFluid<Specie>& pf
)
{
    noCoefficientMixing(perfectFluid);
}


template<class Specie>
inline void Foam::perfectFluid<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::perfectFluid<Specie> Foam::operator+
(
    const perfectFluid<Specie>& pf1,
    const perfectFluid<Specie>& pf2
)
{
    noCoefficientMixing(perfectFluid);
    return pf1;
}


template<class Specie>
inline Foam::perfectFluid<Specie> Foam::operator*
(
    const scalar s,
    const perfectFluid<Specie>& pf
)
{
    return perfectFluid<Specie>
    (
        s*static_cast<const Specie&>(pf),
        pf.R_,
        pf.rho0_
    );
}


template<class Specie>
inline Foam::perfectFluid<Specie> Foam::operator==
(
    const perfectFluid<Specie>& pf1,
    const perfectFluid<Specie>& pf2
)
{
    noCoefficientMixing(perfectFluid);
    return pf1;
}


// ************************************************************************* //
