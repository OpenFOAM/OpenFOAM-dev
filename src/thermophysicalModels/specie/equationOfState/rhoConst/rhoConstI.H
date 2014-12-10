/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::rhoConst<Specie>::rhoConst
(
    const word& name,
    const rhoConst<Specie>& ico
)
:
    Specie(name, ico),
    rho_(ico.rho_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::rhoConst<Specie> >
Foam::rhoConst<Specie>::clone() const
{
    return autoPtr<rhoConst<Specie> >(new rhoConst<Specie>(*this));
}


template<class Specie>
inline Foam::autoPtr<Foam::rhoConst<Specie> >
Foam::rhoConst<Specie>::New(Istream& is)
{
    return autoPtr<rhoConst<Specie> >(new rhoConst<Specie>(is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::rhoConst<Specie>::rho(scalar p, scalar T) const
{
    return rho_;
}


template<class Specie>
inline Foam::scalar Foam::rhoConst<Specie>::psi(scalar, scalar T) const
{
    return 0.0;
}


template<class Specie>
inline Foam::scalar Foam::rhoConst<Specie>::Z(scalar, scalar) const
{
    return 0.0;
}


template<class Specie>
inline Foam::scalar Foam::rhoConst<Specie>::cpMcv(scalar, scalar) const
{
    return 0.0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::rhoConst<Specie>::operator+=(const rhoConst<Specie>& ico)
{
    scalar molr1 = this->nMoles();

    Specie::operator+=(ico);

    molr1 /= this->nMoles();
    scalar molr2 = ico.nMoles()/this->nMoles();

    rho_ = molr1*rho_ + molr2*ico.rho_;
}


template<class Specie>
inline void Foam::rhoConst<Specie>::operator-=(const rhoConst<Specie>& ico)
{
    scalar molr1 = this->nMoles();

    Specie::operator-=(ico);

    molr1 /= this->nMoles();
    scalar molr2 = ico.nMoles()/this->nMoles();

    rho_ = molr1*rho_ - molr2*ico.rho_;
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
    const rhoConst<Specie>& ico1,
    const rhoConst<Specie>& ico2
)
{
    scalar nMoles = ico1.nMoles() + ico2.nMoles();
    scalar molr1 = ico1.nMoles()/nMoles;
    scalar molr2 = ico2.nMoles()/nMoles;

    return rhoConst<Specie>
    (
        static_cast<const Specie&>(ico1)
      + static_cast<const Specie&>(ico2),
        molr1*ico1.rho_ + molr2*ico2.rho_
    );
}


template<class Specie>
inline Foam::rhoConst<Specie> Foam::operator-
(
    const rhoConst<Specie>& ico1,
    const rhoConst<Specie>& ico2
)
{
    scalar nMoles = ico1.nMoles() + ico2.nMoles();
    scalar molr1 = ico1.nMoles()/nMoles;
    scalar molr2 = ico2.nMoles()/nMoles;

    return rhoConst<Specie>
    (
        static_cast<const Specie&>(ico1)
      - static_cast<const Specie&>(ico2),
        molr1*ico1.rho_ - molr2*ico2.rho_
    );
}


template<class Specie>
inline Foam::rhoConst<Specie> Foam::operator*
(
    const scalar s,
    const rhoConst<Specie>& ico
)
{
    return rhoConst<Specie>(s*static_cast<const Specie&>(ico), ico.rho_);
}


template<class Specie>
inline Foam::rhoConst<Specie> Foam::operator==
(
    const rhoConst<Specie>& ico1,
    const rhoConst<Specie>& ico2
)
{
    return ico2 - ico1;
}


// ************************************************************************* //
