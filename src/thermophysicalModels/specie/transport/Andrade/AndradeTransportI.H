/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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
inline Foam::AndradeTransport<Thermo>::AndradeTransport
(
    const Thermo& t,
    const coeffList& muCoeffs,
    const coeffList& kappaCoeffs
)
:
    Thermo(t),
    muCoeffs_(muCoeffs),
    kappaCoeffs_(kappaCoeffs)
{}


template<class Thermo>
inline Foam::AndradeTransport<Thermo>::AndradeTransport
(
    const word& name,
    const AndradeTransport& pt
)
:
    Thermo(name, pt),
    muCoeffs_(pt.muCoeffs_),
    kappaCoeffs_(pt.kappaCoeffs_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::AndradeTransport<Thermo>>
Foam::AndradeTransport<Thermo>::clone() const
{
    return autoPtr<AndradeTransport<Thermo>>
    (
        new AndradeTransport<Thermo>(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::AndradeTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    return exp
    (
        muCoeffs_[0]
      + muCoeffs_[1]*T
      + muCoeffs_[2]*sqr(T)
      + muCoeffs_[3]/(muCoeffs_[4] + T)
    );
}


template<class Thermo>
inline Foam::scalar Foam::AndradeTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    return exp
    (
        kappaCoeffs_[0]
      + kappaCoeffs_[1]*T
      + kappaCoeffs_[2]*sqr(T)
      + kappaCoeffs_[3]/(kappaCoeffs_[4] + T)
    );
}


// ************************************************************************* //
