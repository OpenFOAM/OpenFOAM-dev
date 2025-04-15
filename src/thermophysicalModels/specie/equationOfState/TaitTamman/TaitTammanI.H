/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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

#include "TaitTamman.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::TaitTamman<Specie>::TaitTamman
(
    const Specie& sp,
    const scalar T0,
    const scalar b1,
    const scalar b2,
    const scalar b3,
    const scalar b4
)
:
    Specie(sp),
    T0_(T0),
    b1_(b1),
    b2_(b2),
    b3_(b3),
    b4_(b4)
{}


template<class Specie>
inline Foam::TaitTamman<Specie>::TaitTamman
(
    const word& name,
    const TaitTamman<Specie>& pf
)
:
    Specie(name, pf),
    T0_(pf.T0_),
    b1_(pf.b1_),
    b2_(pf.b2_),
    b3_(pf.b3_),
    b4_(pf.b4_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::TaitTamman<Specie>>
Foam::TaitTamman<Specie>::clone() const
{
    return autoPtr<TaitTamman<Specie>>
    (
        new TaitTamman<Specie>(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::TaitTamman<Specie>::rho
(
    scalar p,
    scalar T
) const
{
    scalar dT = T - T0_;
    scalar v0 = b1_ + b2_*dT;
    scalar B = b3_ * Foam::exp(-b4_*dT);
    scalar v = v0*(1 - C_*Foam::log(1 + p/B));

    // scalar v = (b1_+b2_*dT)*(1-C_*log(1+p*pow(b3_,-1)*exp(b4_*dT)));
    return 1/v;
}


template<class Specie>
inline Foam::scalar Foam::TaitTamman<Specie>::h
(
    scalar p,
    scalar T
) const
{
    // FIXME: is this actually zero?
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::TaitTamman<Specie>::Cp
(
    scalar p,
    scalar T
) const
{
    // FIXME: is this actually zero?
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::TaitTamman<Specie>::e
(
    scalar p,
    scalar T
) const
{
    // FIXME: Add the energy calculation
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::TaitTamman<Specie>::Cv
(
    scalar p,
    scalar T
) const
{
    // FIXME: Add the specific heat calculation
    return 0;
}


// template<class Specie>
// inline Foam::scalar Foam::TaitTamman<Specie>::sp
// (
//     scalar p,
//     scalar T
// ) const
// {
//     scalar n = 1 - 1.0/gamma_;
//     return
//        -pow(p0_ + B_, 1.0/gamma_)*(pow((p + B_), n) - pow((Pstd + B_), n))
//        /(rho0_*T*n);
// }


// template<class Specie>
// inline Foam::scalar Foam::TaitTamman<Specie>::sv
// (
//     scalar p,
//     scalar T
// ) const
// {
//     NotImplemented;
//     return 0;
// }


template<class Specie>
inline Foam::scalar Foam::TaitTamman<Specie>::psi
(
    scalar p,
    scalar T
) const
{
    scalar dT = T - T0_;
    return C_*exp(b4_*dT)*pow(b1_+b2_*dT,-1)*pow(-1+C_*log(1+p*pow(b3_,-1)*exp(b4_*dT)),-2)*pow(b3_+p*exp(b4_*dT),-1);
}


// template<class Specie>
// inline Foam::scalar Foam::TaitTamman<Specie>::Z
// (
//     scalar p,
//     scalar T
// ) const
// {
//     return p/(rho(p, T)*this->R()*T);
// }


template<class Specie>
inline Foam::scalar Foam::TaitTamman<Specie>::CpMCv
(
    scalar p,
    scalar T
) const
{
    scalar dT = T - T0_;
    return T*pow(C_,-1)*pow(b1_+b2_*dT,-1)*pow(p+b3_*exp(b4_*dT),-1)*pow(b1_*b4_*C_*p+b2_*C_*log(1+p*pow(b3_,-1)*exp(b4_*dT))*(p+b3_*exp(b4_*dT))-b2_*(p+b4_*T0_*C_*p-b4_*C_*p*T+b3_*exp(b4_*dT)),2);
}


// template<class Specie>
// inline Foam::scalar Foam::TaitTamman<Specie>::alphav
// (
//     scalar p,
//     scalar T
// ) const
// {
//     return 0;
// }


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::TaitTamman<Specie>::operator+=
(
    const TaitTamman<Specie>& pf
)
{
    noCoefficientMixing(TaitTamman);
}


template<class Specie>
inline void Foam::TaitTamman<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::TaitTamman<Specie> Foam::operator+
(
    const TaitTamman<Specie>& pf1,
    const TaitTamman<Specie>& pf2
)
{
    noCoefficientMixing(TaitTamman);
    return pf1;
}


template<class Specie>
inline Foam::TaitTamman<Specie> Foam::operator*
(
    const scalar s,
    const TaitTamman<Specie>& pf
)
{
    return TaitTamman<Specie>
    (
        s*static_cast<const Specie&>(pf),
        pf.T0_,
        pf.b1_,
        pf.b2_,
        pf.b3_,
        pf.b4_
    );
}


template<class Specie>
inline Foam::TaitTamman<Specie> Foam::operator==
(
    const TaitTamman<Specie>& pf1,
    const TaitTamman<Specie>& pf2
)
{
    noCoefficientMixing(TaitTamman);
    return pf1;
}


// ************************************************************************* //
