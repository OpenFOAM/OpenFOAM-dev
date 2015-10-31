/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

#include "Boussinesq.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::Boussinesq<Specie>::Boussinesq
(
    const Specie& sp, const scalar rho0, const scalar T0, const scalar beta
)
:
    Specie(sp),
    rho0_(rho0),
    T0_(T0),
    beta_(beta)
{}


template<class Specie>
inline Foam::Boussinesq<Specie>::Boussinesq
(
    const Boussinesq& b
)
:
    Specie(b),
    rho0_(b.rho0_),
    T0_(b.T0_),
    beta_(b.beta_)
{}


template<class Specie>
inline Foam::Boussinesq<Specie>::Boussinesq
(
    const word& name,
    const Boussinesq<Specie>& b
)
:
    Specie(name, b),
    rho0_(b.rho0_),
    T0_(b.T0_),
    beta_(b.beta_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::Boussinesq<Specie> >
Foam::Boussinesq<Specie>::clone() const
{
    return autoPtr<Boussinesq<Specie> >
    (
        new Boussinesq<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::Boussinesq<Specie> >
Foam::Boussinesq<Specie>::New
(
    Istream& is
)
{
    return autoPtr<Boussinesq<Specie> >
    (
        new Boussinesq<Specie>(is)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::Boussinesq<Specie> >
Foam::Boussinesq<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<Boussinesq<Specie> >
    (
        new Boussinesq<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::Boussinesq<Specie>::rho
(
    scalar p,
    scalar T
) const
{
    return rho0_*(1.0 - beta_*(T - T0_));
}


template<class Specie>
inline Foam::scalar Foam::Boussinesq<Specie>::s
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::Boussinesq<Specie>::psi
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::Boussinesq<Specie>::Z
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::Boussinesq<Specie>::cpMcv
(
    scalar p,
    scalar T
) const
{
    return RR;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::Boussinesq<Specie>&
Foam::Boussinesq<Specie>::operator=
(
    const Boussinesq<Specie>& b
)
{
    Specie::operator=(b);

    rho0_ = b.rho0_;
    T0_ = b.T0_;
    beta_ = b.beta_;

    return *this;
}

template<class Specie>
inline void Foam::Boussinesq<Specie>::operator+=
(
    const Boussinesq<Specie>& b
)
{
    scalar molr1 = this->nMoles();
    Specie::operator+=(b);
    molr1 /= this->nMoles();
    scalar molr2 = b.nMoles()/this->nMoles();

    rho0_ = molr1*rho0_ + molr2*b.rho0_;
    T0_ = molr1*T0_ + molr2*b.T0_;
    beta_ = molr1*beta_ + molr2*b.beta_;
}


template<class Specie>
inline void Foam::Boussinesq<Specie>::operator-=
(
    const Boussinesq<Specie>& b
)
{
    Specie::operator-=(b);
    rho0_ = b.rho0_;
    T0_ = b.T0_;
    beta_ = b.beta_;
}


template<class Specie>
inline void Foam::Boussinesq<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::Boussinesq<Specie> Foam::operator+
(
    const Boussinesq<Specie>& b1,
    const Boussinesq<Specie>& b2
)
{
    scalar nMoles = b1.nMoles() + b2.nMoles();
    scalar molr1 = b1.nMoles()/nMoles;
    scalar molr2 = b2.nMoles()/nMoles;

    return Boussinesq<Specie>
    (
        static_cast<const Specie&>(b1)
      + static_cast<const Specie&>(b2),
        molr1*b1.rho0_ + molr2*b2.rho0_,
        molr1*b1.T0_ + molr2*b2.T0_,
        molr1*b1.beta_ + molr2*b2.beta_
    );
}


template<class Specie>
inline Foam::Boussinesq<Specie> Foam::operator-
(
    const Boussinesq<Specie>& b1,
    const Boussinesq<Specie>& b2
)
{
    return Boussinesq<Specie>
    (
        static_cast<const Specie&>(b1)
      - static_cast<const Specie&>(b2),
        b1.rho0_ - b2.rho0_,
        b1.T0_ - b2.T0_,
        b1.beta_ - b2.beta_
    );
}


template<class Specie>
inline Foam::Boussinesq<Specie> Foam::operator*
(
    const scalar s,
    const Boussinesq<Specie>& b
)
{
    return Boussinesq<Specie>
    (
        s*static_cast<const Specie&>(b),
        b.rho0_,
        b.T0_,
        b.beta_
    );
}


template<class Specie>
inline Foam::Boussinesq<Specie> Foam::operator==
(
    const Boussinesq<Specie>& pg1,
    const Boussinesq<Specie>& pg2
)
{
    return pg2 - pg1;
}


// ************************************************************************* //
