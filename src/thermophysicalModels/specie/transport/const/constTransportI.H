/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::constTransport<Thermo>::constTransport
(
    const Thermo& t,
    const scalar mu,
    const scalar Pr
)
:
    Thermo(t),
    mu_(mu),
    rPr_(1.0/Pr)
{}


template<class Thermo>
inline Foam::constTransport<Thermo>::constTransport
(
    const word& name,
    const constTransport& ct
)
:
    Thermo(name, ct),
    mu_(ct.mu_),
    rPr_(ct.rPr_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::constTransport<Thermo> >
Foam::constTransport<Thermo>::clone() const
{
    return autoPtr<constTransport<Thermo> >
    (
        new constTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::constTransport<Thermo> >
Foam::constTransport<Thermo>::New
(
    Istream& is
)
{
    return autoPtr<constTransport<Thermo> >
    (
        new constTransport<Thermo>(is)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::constTransport<Thermo> >
Foam::constTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<constTransport<Thermo> >
    (
        new constTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::constTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    return mu_;
}


template<class Thermo>
inline Foam::scalar Foam::constTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    return this->Cpv(p, T)*mu(p, T)*rPr_;
}


template<class Thermo>
inline Foam::scalar Foam::constTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{
    return mu(p, T)*rPr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::constTransport<Thermo>& Foam::constTransport<Thermo>::operator=
(
    const constTransport<Thermo>& ct
)
{
    Thermo::operator=(ct);

    mu_ = ct.mu_;
    rPr_ = ct.rPr_;

    return *this;
}


template<class Thermo>
inline void Foam::constTransport<Thermo>::operator+=
(
    const constTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator+=(st);

    molr1 /= this->nMoles();
    scalar molr2 = st.nMoles()/this->nMoles();

    mu_ = molr1*mu_ + molr2*st.mu_;
    rPr_ = 1.0/(molr1/rPr_ + molr2/st.rPr_);
}


template<class Thermo>
inline void Foam::constTransport<Thermo>::operator-=
(
    const constTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator-=(st);

    molr1 /= this->nMoles();
    scalar molr2 = st.nMoles()/this->nMoles();

    mu_ = molr1*mu_ - molr2*st.mu_;
    rPr_ = 1.0/(molr1/rPr_ - molr2/st.rPr_);
}


template<class Thermo>
inline void Foam::constTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::constTransport<Thermo> Foam::operator+
(
    const constTransport<Thermo>& ct1,
    const constTransport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) + static_cast<const Thermo&>(ct2)
    );

    scalar molr1 = ct1.nMoles()/t.nMoles();
    scalar molr2 = ct2.nMoles()/t.nMoles();

    return constTransport<Thermo>
    (
        t,
        molr1*ct1.mu_ + molr2*ct2.mu_,
        1.0/(molr1/ct1.rPr_ + molr2/ct2.rPr_)
    );
}


template<class Thermo>
inline Foam::constTransport<Thermo> Foam::operator-
(
    const constTransport<Thermo>& ct1,
    const constTransport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) - static_cast<const Thermo&>(ct2)
    );

    scalar molr1 = ct1.nMoles()/t.nMoles();
    scalar molr2 = ct2.nMoles()/t.nMoles();

    return constTransport<Thermo>
    (
        t,
        molr1*ct1.mu_ - molr2*ct2.mu_,
        1.0/(molr1/ct1.rPr_ - molr2/ct2.rPr_)
    );
}


template<class Thermo>
inline Foam::constTransport<Thermo> Foam::operator*
(
    const scalar s,
    const constTransport<Thermo>& ct
)
{
    return constTransport<Thermo>
    (
        s*static_cast<const Thermo&>(ct),
        ct.mu_,
        1.0/ct.rPr_
    );
}


template<class Thermo>
inline Foam::constTransport<Thermo> Foam::operator==
(
    const constTransport<Thermo>& ct1,
    const constTransport<Thermo>& ct2
)
{
    return ct2 - ct1;
}


// ************************************************************************* //
