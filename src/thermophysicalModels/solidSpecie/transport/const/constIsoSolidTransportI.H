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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class thermo>
inline Foam::constIsoSolidTransport<thermo>::constIsoSolidTransport
(
    const thermo& t,
    const scalar kappa
)
:
    thermo(t),
    kappa_(kappa)
{}


template<class thermo>
inline Foam::constIsoSolidTransport<thermo>::constIsoSolidTransport
(
    const word& name,
    const constIsoSolidTransport& ct
)
:
    thermo(name, ct),
    kappa_(ct.kappa_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::constIsoSolidTransport<Thermo> >
Foam::constIsoSolidTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<constIsoSolidTransport<Thermo> >
    (
        new constIsoSolidTransport<Thermo>(dict)
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class thermo>
inline Foam::scalar Foam::constIsoSolidTransport<thermo>::
kappa(const scalar p, const scalar T) const
{
    return kappa_;
}

template<class thermo>
inline Foam::vector Foam::constIsoSolidTransport<thermo>::
Kappa(const scalar p, const scalar T) const
{
    return vector(kappa_, kappa_, kappa_);
}


template<class thermo>
inline Foam::scalar Foam::constIsoSolidTransport<thermo>::
mu(const scalar p, const scalar T) const
{
    notImplemented
    (
        "Foam::scalar Foam::constIsoSolidTransport<thermo>mu::"
        "("
        "    const scalar p, const scalar T"
        ") const"
    );
    return scalar(0);
}


template<class thermo>
inline Foam::scalar Foam::constIsoSolidTransport<thermo>::
alphah(const scalar p, const scalar T) const
{
    return kappa_/this->Cpv(p, T);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class thermo>
inline Foam::constIsoSolidTransport<thermo>&
Foam::constIsoSolidTransport<thermo>::operator=
(
    const constIsoSolidTransport<thermo>& ct
)
{
    thermo::operator=(ct);
    kappa_ = ct.kappa_;
    return *this;
}


template<class thermo>
inline void Foam::constIsoSolidTransport<thermo>::operator+=
(
    const constIsoSolidTransport<thermo>& ct
)
{
    scalar molr1 = this->nMoles();
    thermo::operator+=(ct);

    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    kappa_ = molr1*kappa_ + molr2*ct.kappa_;
}


template<class thermo>
inline void Foam::constIsoSolidTransport<thermo>::operator-=
(
    const constIsoSolidTransport<thermo>& ct
)
{
    scalar molr1 = this->nMoles();

    thermo::operator-=(ct);

    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    kappa_ = molr1*kappa_ - molr2*ct.kappa_;
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class thermo>
inline Foam::constIsoSolidTransport<thermo> Foam::operator*
(
    const scalar s,
    const constIsoSolidTransport<thermo>& ct
)
{
    return constIsoSolidTransport<thermo>
    (
        s*static_cast<const thermo&>(ct),
        ct.kappa_
    );
}


// ************************************************************************* //
