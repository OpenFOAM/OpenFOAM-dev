/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "NSRDSTransport.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::NSRDSTransport<Thermo>::NSRDSTransport
(
    const word& name,
    const NSRDSTransport& ns
)
:
    Thermo(name, ns),
    muA_(ns.muA_),
    muB_(ns.muB_),
    muC_(ns.muC_),
    muD_(ns.muD_),
    muE_(ns.muE_),
    kappaCoeffs_(ns.kappaCoeffs_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::NSRDSTransport<Thermo>>
Foam::NSRDSTransport<Thermo>::clone() const
{
    return autoPtr<NSRDSTransport<Thermo>>
    (
        new NSRDSTransport<Thermo>(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::NSRDSTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    return exp(muA_ + muB_/T + muC_*log(T) + muD_*pow(T, muE_));
}


template<class Thermo>
inline Foam::scalar Foam::NSRDSTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    return kappaCoeffs_.value(T);
}


// ************************************************************************* //
