/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
inline Foam::exponentialSolidTransport<Thermo>::exponentialSolidTransport
(
    const Thermo& t,
    const scalar kappa0,
    const scalar n0,
    const scalar Tref
)
:
    Thermo(t),
    kappa0_(kappa0),
    n0_(n0),
    Tref_(Tref)
{}


template<class Thermo>
inline Foam::exponentialSolidTransport<Thermo>::exponentialSolidTransport
(
    const word& name,
    const exponentialSolidTransport& ct
)
:
    Thermo(name, ct),
    kappa0_(ct.kappa0_),
    n0_(ct.n0_),
    Tref_(ct.Tref_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::exponentialSolidTransport<Thermo>>
Foam::exponentialSolidTransport<Thermo>::clone() const
{
    return autoPtr<exponentialSolidTransport<Thermo>>
    (
        new exponentialSolidTransport<Thermo>(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::exponentialSolidTransport<Thermo>::kappa
(
    const scalar p, const scalar T
) const
{
    return (kappa0_*pow(T/Tref_, n0_));
}


template<class Thermo>
inline Foam::vector Foam::exponentialSolidTransport<Thermo>::Kappa
(
    const scalar p, const scalar T
) const
{
    NotImplemented;
    const scalar kappa(kappa0_*pow(T/Tref_, n0_));
    return vector(kappa, kappa, kappa);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::exponentialSolidTransport<Thermo>::operator+=
(
    const exponentialSolidTransport<Thermo>& ct
)
{
    scalar Y1 = this->Y();

    Y1 /= this->Y();
    scalar Y2 = ct.Y()/this->Y();

    kappa0_ = Y1*kappa0_ + Y2*ct.kappa0_;
    n0_ = (Y1*n0_ + Y2*ct.n0_);
    Tref_ = (Y1*Tref_ + Y2*ct.Tref_);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::exponentialSolidTransport<Thermo> Foam::operator*
(
    const scalar s,
    const exponentialSolidTransport<Thermo>& ct
)
{
    return exponentialSolidTransport<Thermo>
    (
        s*static_cast<const Thermo&>(ct),
        ct.kappa0_,
        ct.n0_,
        ct.Tref_
    );
}

// ************************************************************************* //
