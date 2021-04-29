/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2021 OpenFOAM Foundation
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

#include "linear.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::linear<Specie>::linear
(
    const Specie& sp,
    const scalar psi,
    const scalar rho0
)
:
    Specie(sp),
    psi_(psi),
    rho0_(rho0)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::linear<Specie>::linear
(
    const word& name,
    const linear<Specie>& pf
)
:
    Specie(name, pf),
    psi_(pf.psi_),
    rho0_(pf.rho0_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::linear<Specie>>
Foam::linear<Specie>::clone() const
{
    return autoPtr<linear<Specie>>(new linear<Specie>(*this));
}


template<class Specie>
inline Foam::autoPtr<Foam::linear<Specie>>
Foam::linear<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<linear<Specie>>(new linear<Specie>(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::linear<Specie>::rho(scalar p, scalar T) const
{
    return rho0_ + psi_*p;
}


template<class Specie>
inline Foam::scalar Foam::linear<Specie>::H(scalar p, scalar T) const
{
    return log(rho(p, T))/psi_;
}


template<class Specie>
inline Foam::scalar Foam::linear<Specie>::Cp(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::linear<Specie>::E(scalar p, scalar T) const
{
    const scalar rho = this->rho(p, T);

    return log(rho)/psi_ - p/rho;
}


template<class Specie>
inline Foam::scalar Foam::linear<Specie>::Cv(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::linear<Specie>::Sp(scalar p, scalar T) const
{
    return -log((rho0_ + psi_*p)/(rho0_ + psi_*Pstd))/(T*psi_);
}


template<class Specie>
inline Foam::scalar Foam::linear<Specie>::Sv(scalar p, scalar T) const
{
    NotImplemented;
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::linear<Specie>::psi(scalar p, scalar T) const
{
    return psi_;
}


template<class Specie>
inline Foam::scalar Foam::linear<Specie>::Z(scalar p, scalar T) const
{
    return p/(rho(p, T)*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::linear<Specie>::CpMCv(scalar p, scalar T) const
{
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::linear<Specie>::operator+=
(
    const linear<Specie>& pf
)
{
    noCoefficientMixing(linear);
}


template<class Specie>
inline void Foam::linear<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::linear<Specie> Foam::operator+
(
    const linear<Specie>& pf1,
    const linear<Specie>& pf2
)
{
    noCoefficientMixing(linear);
    return pf1;
}


template<class Specie>
inline Foam::linear<Specie> Foam::operator*
(
    const scalar s,
    const linear<Specie>& pf
)
{
    return linear<Specie>
    (
        s*static_cast<const Specie&>(pf),
        pf.psi_,
        pf.rho0_
    );
}


template<class Specie>
inline Foam::linear<Specie> Foam::operator==
(
    const linear<Specie>& pf1,
    const linear<Specie>& pf2
)
{
    noCoefficientMixing(linear);
    return pf1;
}


// ************************************************************************* //
