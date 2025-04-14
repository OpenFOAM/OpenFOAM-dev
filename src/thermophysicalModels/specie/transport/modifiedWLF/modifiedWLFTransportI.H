/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2023 OpenFOAM Foundation
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

#include "specie.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::modifiedWLFTransport<Thermo>::modifiedWLFTransport
(
    const Thermo& t,
    const scalar mu0,
    const scalar Tr0,
    const scalar Trp,
    const scalar C1,
    const scalar C2,
    const scalar muMax,
    const scalar muMin,
    const scalar Pr
)
:
    Thermo(t),
    mu0_(mu0),
    Tr0_(Tr0),
    Trp_(Trp),
    C1_(C1),
    C2_(C2),
    muMax_(muMax),
    muMin_(muMin),
    rPr_(scalar(1)/Pr)
{}


template<class Thermo>
inline Foam::modifiedWLFTransport<Thermo>::modifiedWLFTransport
(
    const word& name,
    const modifiedWLFTransport& wlft
)
:
    Thermo(name, wlft),
    mu0_(wlft.mu0_),
    Tr0_(wlft.Tr0_),
    Trp_(wlft.Trp_),
    C1_(wlft.C1_),
    C2_(wlft.C2_),
    muMax_(wlft.muMax_),
    muMin_(wlft.muMin_),
    rPr_(wlft.rPr_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::modifiedWLFTransport<Thermo>>
Foam::modifiedWLFTransport<Thermo>::clone() const
{
    return autoPtr<modifiedWLFTransport<Thermo>>
    (
        new modifiedWLFTransport<Thermo>(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::modifiedWLFTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    Foam::scalar Tr = Tr0_ + Trp_ * p;
    // Info << "WLF Tr: " << Tr << endl;
    // Info << "WLF mu0_: " << mu0_ << endl;
    return max(
        muMin_,
        min(muMax_, mu0_*exp(-C1_*(T - Tr)/(C2_ + T - Tr)))
    );
}


template<class Thermo>
inline Foam::scalar Foam::modifiedWLFTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    return this->Cp(p, T)*mu(p, T)*rPr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::modifiedWLFTransport<Thermo>::operator+=
(
    const modifiedWLFTransport<Thermo>& wlft
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(wlft);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        scalar Y2 = wlft.Y()/this->Y();

        mu0_ = Y1*mu0_ + Y2*wlft.mu0_;
        Tr0_ = Y1*Tr0_ + Y2*wlft.Tr0_;
        Trp_ = Y1*Trp_ + Y2*wlft.Trp_;
        C1_ = Y1*C1_ + Y2*wlft.C1_;
        C2_ = Y1*C2_ + Y2*wlft.C2_;
        muMax_ = Y1*muMax_ + Y2*wlft.muMax_;
        muMin_ = Y1*muMin_ + Y2*wlft.muMin_;
        rPr_ = 1.0/(Y1/rPr_ + Y2/wlft.rPr_);
    }
}


template<class Thermo>
inline void Foam::modifiedWLFTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::modifiedWLFTransport<Thermo> Foam::operator+
(
    const modifiedWLFTransport<Thermo>& wlft1,
    const modifiedWLFTransport<Thermo>& wlft2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(wlft1) + static_cast<const Thermo&>(wlft2)
    );

    if (mag(t.Y()) < small)
    {
        return modifiedWLFTransport<Thermo>
        (
            t,
            0,
            wlft1.mu0_,
            wlft1.Tr0_,
            wlft1.Trp_,
            wlft1.C1_,
            wlft1.C2_,
            wlft1.muMax_,
            wlft1.muMin_,
            wlft1.rPr_
        );
    }
    else
    {
        scalar Y1 = wlft1.Y()/t.Y();
        scalar Y2 = wlft2.Y()/t.Y();

        return modifiedWLFTransport<Thermo>
        (
            t,
            Y1*wlft1.mu0_ + Y2*wlft2.mu0_,
            Y1*wlft1.Tr0_ + Y2*wlft2.Tr0_,
            Y1*wlft1.Trp_ + Y2*wlft2.Trp_,
            Y1*wlft1.C1_ + Y2*wlft2.C1_,
            Y1*wlft1.C2_ + Y2*wlft2.C2_,
            Y1*wlft1.muMax_ + Y2*wlft2.muMax_,
            Y1*wlft1.muMin_ + Y2*wlft2.muMin_,
            1.0/(Y1/wlft1.rPr_ + Y2/wlft2.rPr_)
        );
    }
}


template<class Thermo>
inline Foam::modifiedWLFTransport<Thermo> Foam::operator*
(
    const scalar s,
    const modifiedWLFTransport<Thermo>& wlft
)
{
    return modifiedWLFTransport<Thermo>
    (
        s*static_cast<const Thermo&>(wlft),
        wlft.mu0_,
        wlft.Tr0_,
        wlft.Trp_,
        wlft.C1_,
        wlft.C2_,
        wlft.muMax_,
        wlft.muMin_,
        1.0/wlft.rPr_
    );
}


// ************************************************************************* //
