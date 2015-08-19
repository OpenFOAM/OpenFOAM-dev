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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hRefConstThermo<EquationOfState>::hRefConstThermo
(
    const EquationOfState& st,
    const scalar cp,
    const scalar hf,
    const scalar tref,
    const scalar href
)
:
    EquationOfState(st),
    Cp_(cp),
    Hf_(hf),
    Tref_(tref),
    Href_(href)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hRefConstThermo<EquationOfState>::hRefConstThermo
(
    const word& name,
    const hRefConstThermo& ct
)
:
    EquationOfState(name, ct),
    Cp_(ct.Cp_),
    Hf_(ct.Hf_),
    Tref_(ct.Tref_),
    Href_(ct.Href_)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hRefConstThermo<EquationOfState> >
Foam::hRefConstThermo<EquationOfState>::clone() const
{
    return autoPtr<hRefConstThermo<EquationOfState> >
    (
        new hRefConstThermo<EquationOfState>(*this)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hRefConstThermo<EquationOfState> >
Foam::hRefConstThermo<EquationOfState>::New(Istream& is)
{
    return autoPtr<hRefConstThermo<EquationOfState> >
    (
        new hRefConstThermo<EquationOfState>(is)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hRefConstThermo<EquationOfState> >
Foam::hRefConstThermo<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<hRefConstThermo<EquationOfState> >
    (
        new hRefConstThermo<EquationOfState>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::hRefConstThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::hRefConstThermo<EquationOfState>::cp
(
    const scalar p,
    const scalar T
) const
{
    return Cp_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hRefConstThermo<EquationOfState>::ha
(
    const scalar p, const scalar T
) const
{
    return Cp_*(T-Tref_) + Href_ + Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hRefConstThermo<EquationOfState>::hs
(
    const scalar p, const scalar T
) const
{
    return Cp_*(T-Tref_) + Href_ ;
}


template<class EquationOfState>
inline Foam::scalar Foam::hRefConstThermo<EquationOfState>::hc() const
{
    return Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hRefConstThermo<EquationOfState>::s
(
    const scalar p, const scalar T
) const
{
    return Cp_*log(T/Tstd) + EquationOfState::s(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::hRefConstThermo<EquationOfState>::operator+=
(
    const hRefConstThermo<EquationOfState>& ct
)
{
    scalar molr1 = this->nMoles();

    EquationOfState::operator+=(ct);

    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    Cp_ = molr1*Cp_ + molr2*ct.Cp_;
    Hf_ = molr1*Hf_ + molr2*ct.Hf_;
}


template<class EquationOfState>
inline void Foam::hRefConstThermo<EquationOfState>::operator-=
(
    const hRefConstThermo<EquationOfState>& ct
)
{
    scalar molr1 = this->nMoles();

    EquationOfState::operator-=(ct);

    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    Cp_ = molr1*Cp_ - molr2*ct.Cp_;
    Hf_ = molr1*Hf_ - molr2*ct.Hf_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hRefConstThermo<EquationOfState> Foam::operator+
(
    const hRefConstThermo<EquationOfState>& ct1,
    const hRefConstThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    return hRefConstThermo<EquationOfState>
    (
        eofs,
        ct1.nMoles()/eofs.nMoles()*ct1.Cp_
      + ct2.nMoles()/eofs.nMoles()*ct2.Cp_,
        ct1.nMoles()/eofs.nMoles()*ct1.Hf_
      + ct2.nMoles()/eofs.nMoles()*ct2.Hf_,
        ct1.nMoles()/eofs.nMoles()*ct1.Tref_
      + ct2.nMoles()/eofs.nMoles()*ct2.Tref_,
        ct1.nMoles()/eofs.nMoles()*ct1.Href_
      + ct2.nMoles()/eofs.nMoles()*ct2.Href_
    );
}


template<class EquationOfState>
inline Foam::hRefConstThermo<EquationOfState> Foam::operator-
(
    const hRefConstThermo<EquationOfState>& ct1,
    const hRefConstThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      - static_cast<const EquationOfState&>(ct2)
    );

    return hRefConstThermo<EquationOfState>
    (
        eofs,
        ct1.nMoles()/eofs.nMoles()*ct1.Cp_
      - ct2.nMoles()/eofs.nMoles()*ct2.Cp_,
        ct1.nMoles()/eofs.nMoles()*ct1.Hf_
      - ct2.nMoles()/eofs.nMoles()*ct2.Hf_
    );
}


template<class EquationOfState>
inline Foam::hRefConstThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const hRefConstThermo<EquationOfState>& ct
)
{
    return hRefConstThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ct),
        ct.Cp_,
        ct.Hf_,
        ct.Tref_,
        ct.Href_
    );
}


template<class EquationOfState>
inline Foam::hRefConstThermo<EquationOfState> Foam::operator==
(
    const hRefConstThermo<EquationOfState>& ct1,
    const hRefConstThermo<EquationOfState>& ct2
)
{
    return ct2 - ct1;
}


// ************************************************************************* //
