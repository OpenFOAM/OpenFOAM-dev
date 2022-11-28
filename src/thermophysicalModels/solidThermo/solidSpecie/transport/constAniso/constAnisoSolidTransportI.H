/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
inline Foam::constAnisoSolidTransport<Thermo>::constAnisoSolidTransport
(
    const Thermo& t,
    const vector kappa
)
:
    Thermo(t),
    kappa_(kappa)
{}


template<class Thermo>
inline Foam::constAnisoSolidTransport<Thermo>::constAnisoSolidTransport
(
    const word& name,
    const constAnisoSolidTransport& ct
)
:
    Thermo(name, ct),
    kappa_(ct.kappa_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::constAnisoSolidTransport<Thermo>>
Foam::constAnisoSolidTransport<Thermo>::clone() const
{
    return autoPtr<constAnisoSolidTransport<Thermo>>
    (
        new constAnisoSolidTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::constAnisoSolidTransport<Thermo>>
Foam::constAnisoSolidTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<constAnisoSolidTransport<Thermo>>
    (
        new constAnisoSolidTransport<Thermo>(dict)
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::constAnisoSolidTransport<Thermo>::
kappa(const scalar p, const scalar T) const
{
    NotImplemented;
    return mag(kappa_);
}


template<class Thermo>
inline Foam::vector Foam::constAnisoSolidTransport<Thermo>::
Kappa(const scalar p, const scalar T) const
{
    return kappa_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::constAnisoSolidTransport<Thermo>::operator+=
(
    const constAnisoSolidTransport<Thermo>& ct
)
{
    scalar Y1 = this->Y();

    Y1 /= this->Y();
    scalar Y2 = ct.Y()/this->Y();

    kappa_ = Y1*kappa_ + Y2*ct.kappa_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class Thermo>
inline Foam::constAnisoSolidTransport<Thermo> Foam::operator*
(
    const scalar s,
    const constAnisoSolidTransport<Thermo>& ct
)
{
    return constAnisoSolidTransport<Thermo>
    (
        s*static_cast<const Thermo&>(ct),
        ct.kappa_
    );
}

// ************************************************************************* //
