/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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

#include "rPolynomial.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::rPolynomial<Specie>::rPolynomial
(
    const Specie& sp,
    const coeffList& coeffs
)
:
    Specie(sp),
    C_(coeffs)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::rPolynomial<Specie>::rPolynomial
(
    const word& name,
    const rPolynomial<Specie>& rp
)
:
    Specie(name, rp),
    C_(rp.C_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::rPolynomial<Specie>>
Foam::rPolynomial<Specie>::clone() const
{
    return autoPtr<rPolynomial<Specie>>(new rPolynomial<Specie>(*this));
}


template<class Specie>
inline Foam::autoPtr<Foam::rPolynomial<Specie>>
Foam::rPolynomial<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<rPolynomial<Specie>>(new rPolynomial<Specie>(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::rPolynomial<Specie>::rho(scalar p, scalar T) const
{
    return 1/(C_[0] + (C_[1] + C_[2]*T - C_[4]*p)*T - C_[3]*p);
}


template<class Specie>
inline Foam::scalar Foam::rPolynomial<Specie>::H(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rPolynomial<Specie>::Cp(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rPolynomial<Specie>::E(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rPolynomial<Specie>::Cv(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rPolynomial<Specie>::Sp(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rPolynomial<Specie>::Sv(scalar p, scalar T) const
{
    NotImplemented
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rPolynomial<Specie>::psi(scalar p, scalar T) const
{
    return sqr(rho(p, T))*(C_[3] + C_[4]*T);
}


template<class Specie>
inline Foam::scalar Foam::rPolynomial<Specie>::Z(scalar p, scalar T) const
{
    return p/(rho(p, T)*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::rPolynomial<Specie>::CpMCv(scalar p, scalar T) const
{
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::rPolynomial<Specie>::operator+=
(
    const rPolynomial<Specie>& rp
)
{
    const scalar Y1 = this->Y();
    Specie::operator+=(rp);

    if (mag(this->Y()) > small)
    {
        C_ = (Y1*C_ + rp.Y()*rp.C_)/this->Y();
    }
}


template<class Specie>
inline void Foam::rPolynomial<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::rPolynomial<Specie> Foam::operator+
(
    const rPolynomial<Specie>& rp1,
    const rPolynomial<Specie>& rp2
)
{
    Specie sp
    (
        static_cast<const Specie&>(rp1)
      + static_cast<const Specie&>(rp2)
    );

    if (mag(sp.Y()) < small)
    {
        return rPolynomial<Specie>
        (
            sp,
            rp1.C_
        );
    }
    else
    {
        return rPolynomial<Specie>
        (
            sp,
            (rp1.Y()*rp1.C_ + rp2.Y()*rp2.C_)/sp.Y()
        );
    }

    return rp1;
}


template<class Specie>
inline Foam::rPolynomial<Specie> Foam::operator*
(
    const scalar s,
    const rPolynomial<Specie>& rp
)
{
    return rPolynomial<Specie>
    (
        s*static_cast<const Specie&>(rp),
        rp.C_
    );
}


template<class Specie>
inline Foam::rPolynomial<Specie> Foam::operator==
(
    const rPolynomial<Specie>& rp1,
    const rPolynomial<Specie>& rp2
)
{
    return rPolynomial<Specie>
    (
        static_cast<const Specie&>(rp1) == static_cast<const Specie&>(rp2),
        rPolynomial<Specie>::coeffList::uniform(NaN)
    );
}


// ************************************************************************* //
