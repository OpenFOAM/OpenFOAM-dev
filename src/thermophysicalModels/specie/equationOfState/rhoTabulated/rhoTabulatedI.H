/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2023 OpenFOAM Foundation
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

#include "rhoTabulated.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::rhoTabulated<Specie>::rhoTabulated
(
    const Specie& sp,
    const table2D& rho
)
:
    Specie(sp),
    rho_(rho)
{}


template<class Specie>
inline Foam::rhoTabulated<Specie>::rhoTabulated
(
    const word& name,
    const rhoTabulated<Specie>& ip
)
:
    Specie(name, ip),
    rho_(ip.rho_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::rhoTabulated<Specie>>
Foam::rhoTabulated<Specie>::clone() const
{
    return autoPtr<rhoTabulated<Specie>>
    (
        new rhoTabulated<Specie>(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::rhoTabulated<Specie>::rho
(
    scalar p,
    scalar T
) const
{
    return rho_.value(p, T);
}


template<class Specie>
inline Foam::scalar Foam::rhoTabulated<Specie>::H
(
    scalar p,
    scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rhoTabulated<Specie>::Cp
(
    scalar p,
    scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rhoTabulated<Specie>::E
(
    scalar p,
    scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rhoTabulated<Specie>::Cv
(
    scalar p,
    scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rhoTabulated<Specie>::Sp
(
    scalar p,
    scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rhoTabulated<Specie>::Sv
(
    scalar p,
    scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rhoTabulated<Specie>::psi
(
    scalar p,
    scalar T
) const
{
    return rho_.dfdp(p, T);
}


template<class Specie>
inline Foam::scalar Foam::rhoTabulated<Specie>::Z
(
    scalar p,
    scalar T
) const
{
    return p/(rho(p, T)*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::rhoTabulated<Specie>::CpMCv
(
    scalar p,
    scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rhoTabulated<Specie>::alphav
(
    scalar p,
    scalar T
) const
{
    return -rho_.dfdT(p, T)/rho(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::rhoTabulated<Specie>::operator+=
(
    const rhoTabulated<Specie>& pf
)
{
    noCoefficientMixing(rhoTabulated);
}


template<class Specie>
inline void Foam::rhoTabulated<Specie>::operator*=(const scalar s)
{
    noCoefficientMixing(rhoTabulated);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::rhoTabulated<Specie> Foam::operator+
(
    const rhoTabulated<Specie>& pf1,
    const rhoTabulated<Specie>& pf2
)
{
    noCoefficientMixing(rhoTabulated);
    return pf1;
}


template<class Specie>
inline Foam::rhoTabulated<Specie> Foam::operator*
(
    const scalar s,
    const rhoTabulated<Specie>& pf
)
{
    noCoefficientMixing(rhoTabulated);
    return pf;
}


template<class Specie>
inline Foam::rhoTabulated<Specie> Foam::operator==
(
    const rhoTabulated<Specie>& pf1,
    const rhoTabulated<Specie>& pf2
)
{
    noCoefficientMixing(rhoTabulated);
    return pf1;
}


// ************************************************************************* //
