/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "hExponentialThermo.H"
#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class equationOfState>
inline void Foam::hExponentialThermo<equationOfState>::checkT
(
    const scalar T
) const
{
    if (T < 0.0)
    {
        FatalErrorIn
        (
            "hExponentialThermo<equationOfState>::checkT(const scalar T) const"
        )   << "attempt to use hExponentialThermo<equationOfState>"
               " out of temperature range "
            << T
            << abort(FatalError);
    }
}


template<class equationOfState>
inline Foam::scalar Foam::hExponentialThermo<equationOfState>::integrateCp
(
    const scalar T
) const
{
    return
    (
        c0_*pow(T, n0_ + 1.0)/(pow(Tref_, n0_)*(n0_ + 1.0))
    );
}


template<class equationOfState>
inline Foam::hExponentialThermo<equationOfState>::hExponentialThermo
(
    const word& name,
    const hExponentialThermo& jt
)
:
    equationOfState(name, jt),
    c0_(jt.c0_),
    n0_(jt.n0_),
    Tref_(jt.Tref_),
    Hf_(jt.Hf_)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class equationOfState>
inline Foam::hExponentialThermo<equationOfState>::hExponentialThermo
(
    const equationOfState& st,
    const scalar c0,
    const scalar n0,
    const scalar Tref,
    const scalar Hf
)
:
    equationOfState(st),
    c0_(c0),
    n0_(n0),
    Tref_(Tref),
    Hf_(Hf)
{}


template<class equationOfState>
inline Foam::autoPtr<Foam::hExponentialThermo<equationOfState> >
Foam::hExponentialThermo<equationOfState>::clone() const
{
    return autoPtr<hExponentialThermo<equationOfState> >
    (
        new hExponentialThermo<equationOfState>(*this)
    );
}


template<class equationOfState>
inline Foam::autoPtr<Foam::hExponentialThermo<equationOfState> >
Foam::hExponentialThermo<equationOfState>::New(Istream& is)
{
    return autoPtr<hExponentialThermo<equationOfState> >
    (
        new hExponentialThermo<equationOfState>(is)
    );
}


template<class equationOfState>
inline Foam::autoPtr<Foam::hExponentialThermo<equationOfState> >
Foam::hExponentialThermo<equationOfState>::New(const dictionary& dict)
{
    return autoPtr<hExponentialThermo<equationOfState> >
    (
        new hExponentialThermo<equationOfState>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class equationOfState>
inline Foam::scalar Foam::hExponentialThermo<equationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class equationOfState>
inline Foam::scalar Foam::hExponentialThermo<equationOfState>::cp
(
    const scalar p, const scalar T
) const
{
    return c0_*pow(T/Tref_, n0_);
}


template<class equationOfState>
inline Foam::scalar Foam::hExponentialThermo<equationOfState>::ha
(
    const scalar p, const scalar T
) const
{
    scalar hOffset = integrateCp(specie::Tstd);

    return
    (
        (integrateCp(T) + Hf_ - hOffset)
    );
}


template<class equationOfState>
inline Foam::scalar Foam::hExponentialThermo<equationOfState>::hs
(
    const scalar p, const scalar T
) const
{
    scalar hOffset = integrateCp(specie::Tstd);
    return (integrateCp(T) - hOffset);
}


template<class equationOfState>
inline Foam::scalar Foam::hExponentialThermo<equationOfState>::hc() const
{
    return Hf_;
}


template<class equationOfState>
inline Foam::scalar Foam::hExponentialThermo<equationOfState>::s
(
    const scalar p, const scalar T
) const
{
    notImplemented
    (
        "scalar hExponentialThermo<equationOfState>::"
        "s(const scalar p, const scalar T) const"
    );
    return T;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class equationOfState>
inline void Foam::hExponentialThermo<equationOfState>::operator+=
(
    const hExponentialThermo<equationOfState>& ct
)
{
    scalar molr1 = this->nMoles();

    equationOfState::operator+=(ct);
    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    Hf_ = molr1*Hf_ + molr2*ct.Hf_;
    c0_ = molr1*c0_ + molr2*ct.c0_;
    n0_ = molr1*n0_ + molr2*ct.n0_;
    Tref_ = molr1*Tref_ + molr2*ct.Tref_;
}


template<class equationOfState>
inline void Foam::hExponentialThermo<equationOfState>::operator-=
(
    const hExponentialThermo<equationOfState>& ct
)
{
    scalar molr1 = this->nMoles();

    equationOfState::operator-=(ct);

    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    Hf_ = molr1*Hf_ - molr2*ct.Hf_;
    c0_ = (molr1*c0_ - molr2*ct.c0_);
    n0_ = (molr1*n0_ - molr2*ct.n0_);
    Tref_ = (molr1*Tref_ - molr2*ct.Tref_);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class equationOfState>
inline Foam::hExponentialThermo<equationOfState> Foam::operator+
(
    const hExponentialThermo<equationOfState>& ct1,
    const hExponentialThermo<equationOfState>& ct2
)
{
    equationOfState eofs
    (
        static_cast<const equationOfState&>(ct1)
      + static_cast<const equationOfState&>(ct2)
    );

    return hExponentialThermo<equationOfState>
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


template<class equationOfState>
inline Foam::hExponentialThermo<equationOfState> Foam::operator-
(
    const hExponentialThermo<equationOfState>& ct1,
    const hExponentialThermo<equationOfState>& ct2
)
{
    equationOfState eofs
    (
        static_cast<const equationOfState&>(ct1)
      + static_cast<const equationOfState&>(ct2)
    );

    return hExponentialThermo<equationOfState>
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


template<class equationOfState>
inline Foam::hExponentialThermo<equationOfState> Foam::operator*
(
    const scalar s,
    const hExponentialThermo<equationOfState>& ct
)
{
    return hExponentialThermo<equationOfState>
    (
        s*static_cast<const equationOfState&>(ct),
        ct.c0_,
        ct.n0_,
        ct.Tref_,
        ct.Hf_
    );
}


template<class equationOfState>
inline Foam::hExponentialThermo<equationOfState> Foam::operator==
(
    const hExponentialThermo<equationOfState>& ct1,
    const hExponentialThermo<equationOfState>& ct2
)
{
    return ct2 - ct1;
}


// ************************************************************************* //
