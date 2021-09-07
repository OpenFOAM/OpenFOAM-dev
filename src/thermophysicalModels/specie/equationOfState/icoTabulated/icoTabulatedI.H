/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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

#include "icoTabulated.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::icoTabulated<Specie>::icoTabulated
(
    const Specie& sp,
    const nonUniformTable& rho
)
:
    Specie(sp),
    rho_(rho)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::icoTabulated<Specie>::icoTabulated
(
    const word& name,
    const icoTabulated<Specie>& ip
)
:
    Specie(name, ip),
    rho_(ip.rho_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::icoTabulated<Specie>>
Foam::icoTabulated<Specie>::clone() const
{
    return autoPtr<icoTabulated<Specie>>
    (
        new icoTabulated<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::icoTabulated<Specie>>
Foam::icoTabulated<Specie>::New(const dictionary& dict)
{
    return autoPtr<icoTabulated<Specie>>
    (
        new icoTabulated<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::icoTabulated<Specie>::rho
(
    scalar p,
    scalar T
) const
{
    return rho_.value(T);
}


template<class Specie>
inline Foam::scalar Foam::icoTabulated<Specie>::H
(
    scalar p,
    scalar T
) const
{
    return p/this->rho(p, T);
}


template<class Specie>
inline Foam::scalar Foam::icoTabulated<Specie>::Cp
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::icoTabulated<Specie>::E
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::icoTabulated<Specie>::Cv
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::icoTabulated<Specie>::Sp
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::icoTabulated<Specie>::Sv
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::icoTabulated<Specie>::psi
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::icoTabulated<Specie>::Z
(
    scalar p,
    scalar T
) const
{
    return p/(rho(p, T)*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::icoTabulated<Specie>::CpMCv
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::icoTabulated<Specie>::alphav
(
    scalar p,
    scalar T
) const
{
    return -rho_.dfdT(T)/rho(p, T);
}


// ************************************************************************* //
