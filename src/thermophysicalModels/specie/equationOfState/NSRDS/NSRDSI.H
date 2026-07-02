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

#include "NSRDS.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::NSRDS<Specie>::NSRDS(const word& name, const NSRDS<Specie>& ns)
:
    Specie(name, ns),
    a_(ns.a_),
    b_(ns.b_),
    c_(ns.c_),
    d_(ns.d_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::NSRDS<Specie>>
Foam::NSRDS<Specie>::clone() const
{
    return autoPtr<NSRDS<Specie>>
    (
        new NSRDS<Specie>(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::NSRDS<Specie>::rho
(
    const scalar p,
    const scalar T
) const
{
    return a_/pow(b_, 1 + pow(1 - T/c_, d_));
}


template<class Specie>
inline Foam::scalar Foam::NSRDS<Specie>::h
(
    const scalar p,
    const scalar T
) const
{
    return p/this->rho(p, T);
}


template<class Specie>
inline Foam::scalar Foam::NSRDS<Specie>::Cp
(
    const scalar p,
    const scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::NSRDS<Specie>::e
(
    const scalar p,
    const scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::NSRDS<Specie>::Cv
(
    const scalar p,
    const scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::NSRDS<Specie>::sp
(
    const scalar p,
    const scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::NSRDS<Specie>::sv
(
    const scalar p,
    const scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::NSRDS<Specie>::psi
(
    const scalar p,
    const scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::NSRDS<Specie>::Z
(
    const scalar p,
    const scalar T
) const
{
    return p/(rho(p, T)*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::NSRDS<Specie>::CpMCv
(
    const scalar p,
    const scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::NSRDS<Specie>::alphav
(
    const scalar p,
    const scalar T
) const
{
    return - log(b_)/c_*d_*pow(1 - T/c_, d_ - 1);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::NSRDS<Specie>::operator+=(const NSRDS<Specie>& ns)
{
    noCoefficientMixing(NSRDS);
}


template<class Specie>
inline void Foam::NSRDS<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
Foam::NSRDS<Specie> Foam::operator+
(
    const NSRDS<Specie>& ns1,
    const NSRDS<Specie>& ns2
)
{
    noCoefficientMixing(NSRDS);
    return ns1;
}


template<class Specie>
Foam::NSRDS<Specie> Foam::operator*
(
    const scalar s,
    const NSRDS<Specie>& ns
)
{
    return NSRDS<Specie>
    (
        s*static_cast<const Specie&>(ns),
        ns.rhoCoeffs_
    );
}


template<class Specie>
Foam::NSRDS<Specie> Foam::operator==
(
    const NSRDS<Specie>& ns1,
    const NSRDS<Specie>& ns2
)
{
    noCoefficientMixing(NSRDS);
    return ns1;
}


// ************************************************************************* //
