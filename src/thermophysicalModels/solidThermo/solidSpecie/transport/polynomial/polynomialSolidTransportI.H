/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2022 OpenFOAM Foundation
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

template<class Thermo, int PolySize>
inline Foam::polynomialSolidTransport<Thermo, PolySize>::
polynomialSolidTransport
(
    const Thermo& t,
    const Polynomial<PolySize>& kappaCoeffs
)
:
    Thermo(t),
    kappaCoeffs_(kappaCoeffs)
{}


template<class Thermo, int PolySize>
inline Foam::polynomialSolidTransport<Thermo, PolySize>::
polynomialSolidTransport
(
    const word& name,
    const polynomialSolidTransport& pt
)
:
    Thermo(name, pt),
    kappaCoeffs_(pt.kappaCoeffs_)
{}


template<class Thermo, int PolySize>
inline Foam::autoPtr<Foam::polynomialSolidTransport<Thermo, PolySize>>
Foam::polynomialSolidTransport<Thermo, PolySize>::clone() const
{
    return autoPtr<polynomialSolidTransport<Thermo, PolySize>>
    (
        new polynomialSolidTransport<Thermo, PolySize>(*this)
    );
}


template<class Thermo, int PolySize>
inline Foam::autoPtr<Foam::polynomialSolidTransport<Thermo, PolySize>>
Foam::polynomialSolidTransport<Thermo, PolySize>::New(const dictionary& dict)
{
    return autoPtr<polynomialSolidTransport<Thermo, PolySize>>
    (
        new polynomialSolidTransport<Thermo, PolySize>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline Foam::scalar Foam::polynomialSolidTransport<Thermo, PolySize>::kappa
(
    const scalar p,
    const scalar T
) const
{
    return kappaCoeffs_.value(T);
}


template<class Thermo, int PolySize>
inline Foam::vector Foam::polynomialSolidTransport<Thermo, PolySize>::Kappa
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    const scalar kappa(kappaCoeffs_.value(T));
    return vector(kappa, kappa, kappa);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline void Foam::polynomialSolidTransport<Thermo, PolySize>::operator+=
(
    const polynomialSolidTransport<Thermo, PolySize>& pt
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(pt);

    Y1 /= this->Y();
    scalar Y2 = pt.Y()/this->Y();

    kappaCoeffs_ = Y1*kappaCoeffs_ + Y2*pt.kappaCoeffs_;
}


template<class Thermo, int PolySize>
inline void Foam::polynomialSolidTransport<Thermo, PolySize>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline Foam::polynomialSolidTransport<Thermo, PolySize> Foam::operator+
(
    const polynomialSolidTransport<Thermo, PolySize>& pt1,
    const polynomialSolidTransport<Thermo, PolySize>& pt2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(pt1) + static_cast<const Thermo&>(pt2)
    );

    scalar Y1 = pt1.Y()/t.Y();
    scalar Y2 = pt2.Y()/t.Y();

    return polynomialSolidTransport<Thermo, PolySize>
    (
        t,
        Y1*pt1.kappaCoeffs_ + Y2*pt2.kappaCoeffs_
    );
}


template<class Thermo, int PolySize>
inline Foam::polynomialSolidTransport<Thermo, PolySize> Foam::operator*
(
    const scalar s,
    const polynomialSolidTransport<Thermo, PolySize>& pt
)
{
    return polynomialSolidTransport<Thermo, PolySize>
    (
        s*static_cast<const Thermo&>(pt),
        pt.kappaCoeffs_
    );
}


// ************************************************************************* //
