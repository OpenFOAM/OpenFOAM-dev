/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
inline Foam::eConstThermo<EquationOfState>::eConstThermo
(
    const EquationOfState& st,
    const scalar cv,
    const scalar hf
)
:
    EquationOfState(st),
    Cv_(cv),
    Hf_(hf)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::eConstThermo<EquationOfState>::eConstThermo
(
    const word& name,
    const eConstThermo<EquationOfState>& ct
)
:
    EquationOfState(name, ct),
    Cv_(ct.Cv_),
    Hf_(ct.Hf_)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::eConstThermo<EquationOfState> >
Foam::eConstThermo<EquationOfState>::clone() const
{
    return autoPtr<eConstThermo<EquationOfState> >
    (
        new eConstThermo<EquationOfState>(*this)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::eConstThermo<EquationOfState> >
Foam::eConstThermo<EquationOfState>::New(Istream& is)
{
    return autoPtr<eConstThermo<EquationOfState> >
    (
        new eConstThermo<EquationOfState>(is)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::eConstThermo<EquationOfState> >
Foam::eConstThermo<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<eConstThermo<EquationOfState> >
    (
        new eConstThermo<EquationOfState>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::cp
(
    const scalar p,
    const scalar T
) const
{
    return Cv_ + this->cpMcv(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::ha
(
    const scalar p,
    const scalar T
) const
{
    return cp(p, T)*T + Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::hs
(
    const scalar p,
    const scalar T
) const
{
    return cp(p, T)*T;
}


template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::hc() const
{
    return Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::eConstThermo<EquationOfState>::s
(
    const scalar p,
    const scalar T
) const
{
    notImplemented
    (
        "scalar eConstThermo<EquationOfState>::"
        "s(const scalar p, const scalar) const"
    );
    return T;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::eConstThermo<EquationOfState>::operator+=
(
    const eConstThermo<EquationOfState>& ct
)
{
    scalar molr1 = this->nMoles();

    EquationOfState::operator+=(ct);

    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    Cv_ = molr1*Cv_ + molr2*ct.Cv_;
    Hf_ = molr1*Hf_ + molr2*ct.Hf_;
}


template<class EquationOfState>
inline void Foam::eConstThermo<EquationOfState>::operator-=
(
    const eConstThermo<EquationOfState>& ct
)
{
    scalar molr1 = this->nMoles();

    EquationOfState::operator-=(ct);

    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    Cv_ = molr1*Cv_ - molr2*ct.Cv_;
    Hf_ = molr1*Hf_ - molr2*ct.Hf_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::eConstThermo<EquationOfState> Foam::operator+
(
    const eConstThermo<EquationOfState>& ct1,
    const eConstThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    return eConstThermo<EquationOfState>
    (
        eofs,
        ct1.nMoles()/eofs.nMoles()*ct1.Cv_
      + ct2.nMoles()/eofs.nMoles()*ct2.Cv_,
        ct1.nMoles()/eofs.nMoles()*ct1.Hf_
      + ct2.nMoles()/eofs.nMoles()*ct2.Hf_
    );
}


template<class EquationOfState>
inline Foam::eConstThermo<EquationOfState> Foam::operator-
(
    const eConstThermo<EquationOfState>& ct1,
    const eConstThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      - static_cast<const EquationOfState&>(ct2)
    );

    return eConstThermo<EquationOfState>
    (
        eofs,
        ct1.nMoles()/eofs.nMoles()*ct1.Cv_
      - ct2.nMoles()/eofs.nMoles()*ct2.Cv_,
        ct1.nMoles()/eofs.nMoles()*ct1.Hf_
      - ct2.nMoles()/eofs.nMoles()*ct2.Hf_
    );
}


template<class EquationOfState>
inline Foam::eConstThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const eConstThermo<EquationOfState>& ct
)
{
    return eConstThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ct),
        ct.Cv_,
        ct.Hf_
    );
}


template<class EquationOfState>
inline Foam::eConstThermo<EquationOfState> Foam::operator==
(
    const eConstThermo<EquationOfState>& ct1,
    const eConstThermo<EquationOfState>& ct2
)
{
    return ct2 - ct1;
}


// ************************************************************************* //
