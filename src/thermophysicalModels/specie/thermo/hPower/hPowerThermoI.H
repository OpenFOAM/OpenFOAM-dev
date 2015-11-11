/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2015 OpenFOAM Foundation
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

#include "hPowerThermo.H"
#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::hPowerThermo<EquationOfState>::checkT
(
    const scalar T
) const
{
    if (T < 0)
    {
        FatalErrorInFunction
            << "attempt to evaluate hPowerThermo<EquationOfState>"
               " for negative temperature " << T
            << abort(FatalError);
    }
}


template<class EquationOfState>
inline Foam::hPowerThermo<EquationOfState>::hPowerThermo
(
    const word& name,
    const hPowerThermo& jt
)
:
    EquationOfState(name, jt),
    c0_(jt.c0_),
    n0_(jt.n0_),
    Tref_(jt.Tref_),
    Hf_(jt.Hf_)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hPowerThermo<EquationOfState>::hPowerThermo
(
    const EquationOfState& st,
    const scalar c0,
    const scalar n0,
    const scalar Tref,
    const scalar Hf
)
:
    EquationOfState(st),
    c0_(c0),
    n0_(n0),
    Tref_(Tref),
    Hf_(Hf)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hPowerThermo<EquationOfState> >
Foam::hPowerThermo<EquationOfState>::clone() const
{
    return autoPtr<hPowerThermo<EquationOfState> >
    (
        new hPowerThermo<EquationOfState>(*this)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hPowerThermo<EquationOfState> >
Foam::hPowerThermo<EquationOfState>::New(Istream& is)
{
    return autoPtr<hPowerThermo<EquationOfState> >
    (
        new hPowerThermo<EquationOfState>(is)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hPowerThermo<EquationOfState> >
Foam::hPowerThermo<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<hPowerThermo<EquationOfState> >
    (
        new hPowerThermo<EquationOfState>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::hPowerThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::hPowerThermo<EquationOfState>::cp
(
    const scalar p, const scalar T
) const
{
    return c0_*pow(T/Tref_, n0_);
}


template<class EquationOfState>
inline Foam::scalar Foam::hPowerThermo<EquationOfState>::ha
(
    const scalar p, const scalar T
) const
{
    return hs(p, T) + hc();
}


template<class EquationOfState>
inline Foam::scalar Foam::hPowerThermo<EquationOfState>::hs
(
    const scalar p, const scalar T
) const
{
    return
        c0_*(pow(T, n0_ + 1) - pow(Tstd, n0_ + 1))/(pow(Tref_, n0_)*(n0_ + 1));
}


template<class EquationOfState>
inline Foam::scalar Foam::hPowerThermo<EquationOfState>::hc() const
{
    return Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hPowerThermo<EquationOfState>::s
(
    const scalar p, const scalar T
) const
{
    return
        c0_*(pow(T, n0_) - pow(Tstd, n0_))/(pow(Tref_, n0_)*n0_)
      + EquationOfState::s(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::hPowerThermo<EquationOfState>::operator+=
(
    const hPowerThermo<EquationOfState>& ct
)
{
    scalar molr1 = this->nMoles();

    EquationOfState::operator+=(ct);
    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    Hf_ = molr1*Hf_ + molr2*ct.Hf_;
    c0_ = molr1*c0_ + molr2*ct.c0_;
    n0_ = molr1*n0_ + molr2*ct.n0_;
    Tref_ = molr1*Tref_ + molr2*ct.Tref_;
}


template<class EquationOfState>
inline void Foam::hPowerThermo<EquationOfState>::operator-=
(
    const hPowerThermo<EquationOfState>& ct
)
{
    scalar molr1 = this->nMoles();

    EquationOfState::operator-=(ct);

    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    Hf_ = molr1*Hf_ - molr2*ct.Hf_;
    c0_ = (molr1*c0_ - molr2*ct.c0_);
    n0_ = (molr1*n0_ - molr2*ct.n0_);
    Tref_ = (molr1*Tref_ - molr2*ct.Tref_);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hPowerThermo<EquationOfState> Foam::operator+
(
    const hPowerThermo<EquationOfState>& ct1,
    const hPowerThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    return hPowerThermo<EquationOfState>
    (
        eofs,
        ct1.nMoles()/eofs.nMoles()*ct1.c0_
      + ct2.nMoles()/eofs.nMoles()*ct2.c0_,
        ct1.nMoles()/eofs.nMoles()*ct1.n0_
      + ct2.nMoles()/eofs.nMoles()*ct2.n0_,
        ct1.nMoles()/eofs.nMoles()*ct1.Tref_
      + ct2.nMoles()/eofs.nMoles()*ct2.Tref_,
        ct1.nMoles()/eofs.nMoles()*ct1.Hf_
      + ct2.nMoles()/eofs.nMoles()*ct2.Hf_
    );
}


template<class EquationOfState>
inline Foam::hPowerThermo<EquationOfState> Foam::operator-
(
    const hPowerThermo<EquationOfState>& ct1,
    const hPowerThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    return hPowerThermo<EquationOfState>
    (
        eofs,
        ct1.nMoles()/eofs.nMoles()*ct1.c0_
      - ct2.nMoles()/eofs.nMoles()*ct2.c0_,
        ct1.nMoles()/eofs.nMoles()*ct1.n0_
      - ct2.nMoles()/eofs.nMoles()*ct2.n0_,
        ct1.nMoles()/eofs.nMoles()*ct1.Tref_
      - ct2.nMoles()/eofs.nMoles()*ct2.Tref_,
        ct1.nMoles()/eofs.nMoles()*ct1.Hf_
      - ct2.nMoles()/eofs.nMoles()*ct2.Hf_
    );
}


template<class EquationOfState>
inline Foam::hPowerThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const hPowerThermo<EquationOfState>& ct
)
{
    return hPowerThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ct),
        ct.c0_,
        ct.n0_,
        ct.Tref_,
        ct.Hf_
    );
}


template<class EquationOfState>
inline Foam::hPowerThermo<EquationOfState> Foam::operator==
(
    const hPowerThermo<EquationOfState>& ct1,
    const hPowerThermo<EquationOfState>& ct2
)
{
    return ct2 - ct1;
}


// ************************************************************************* //
