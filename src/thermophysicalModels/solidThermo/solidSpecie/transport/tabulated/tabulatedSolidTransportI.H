/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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
inline Foam::tabulatedSolidTransport<Thermo>::tabulatedSolidTransport
(
    const Thermo& t,
    const nonUniformTable& kappa
)
:
    Thermo(t),
    kappa_(kappa)
{}


template<class Thermo>
inline Foam::tabulatedSolidTransport<Thermo>::tabulatedSolidTransport
(
    const word& name,
    const tabulatedSolidTransport& tst
)
:
    Thermo(name, tst),
    kappa_(tst.kappa_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::tabulatedSolidTransport<Thermo>>
Foam::tabulatedSolidTransport<Thermo>::clone() const
{
    return autoPtr<tabulatedSolidTransport<Thermo>>
    (
        new tabulatedSolidTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::tabulatedSolidTransport<Thermo>>
Foam::tabulatedSolidTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<tabulatedSolidTransport<Thermo>>
    (
        new tabulatedSolidTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::tabulatedSolidTransport<Thermo>::kappa
(
    const scalar p, const scalar T
) const
{
    return kappa_.value(T);
}


template<class Thermo>
inline Foam::vector Foam::tabulatedSolidTransport<Thermo>::Kappa
(
    const scalar p, const scalar T
) const
{
    const scalar kappa = kappa_.value(T);
    return vector(kappa, kappa, kappa);
}


template<class Thermo>
inline Foam::scalar Foam::tabulatedSolidTransport<Thermo>::
mu(const scalar p, const scalar T) const
{
    NotImplemented;
    return scalar(0);
}


// ************************************************************************* //
