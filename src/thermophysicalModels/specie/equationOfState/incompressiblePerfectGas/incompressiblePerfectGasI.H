/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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

#include "incompressiblePerfectGas.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::incompressiblePerfectGas<Specie>::incompressiblePerfectGas
(
    const Specie& sp,  const scalar pRef
)
:
    Specie(sp),
    pRef_(pRef)
{}


template<class Specie>
inline Foam::incompressiblePerfectGas<Specie>::incompressiblePerfectGas
(
    const word& name,
    const incompressiblePerfectGas<Specie>& ipg
)
:
    Specie(name, ipg),
    pRef_(ipg.pRef_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::incompressiblePerfectGas<Specie>>
Foam::incompressiblePerfectGas<Specie>::clone() const
{
    return autoPtr<incompressiblePerfectGas<Specie>>
    (
        new incompressiblePerfectGas<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::incompressiblePerfectGas<Specie>>
Foam::incompressiblePerfectGas<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<incompressiblePerfectGas<Specie>>
    (
        new incompressiblePerfectGas<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::incompressiblePerfectGas<Specie>::rho
(
    scalar p,
    scalar T
) const
{
    return pRef_/(this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::incompressiblePerfectGas<Specie>::H
(
    scalar p,
    scalar T
) const
{
    return p/this->rho(p, T);
}


template<class Specie>
inline Foam::scalar Foam::incompressiblePerfectGas<Specie>::Cp
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::incompressiblePerfectGas<Specie>::E
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::incompressiblePerfectGas<Specie>::Cv
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::incompressiblePerfectGas<Specie>::Sp
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::incompressiblePerfectGas<Specie>::Sv
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::incompressiblePerfectGas<Specie>::psi
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::incompressiblePerfectGas<Specie>::Z
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::incompressiblePerfectGas<Specie>::CpMCv
(
    scalar p,
    scalar T
) const
{
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::incompressiblePerfectGas<Specie>::operator+=
(
    const incompressiblePerfectGas<Specie>& ipg
)
{
    if (notEqual(pRef_, ipg.pRef_))
    {
        FatalErrorInFunction
            << "pRef " << pRef_ << " for "
            << (this->name().size() ? this->name() : "others")
            << " != " << ipg.pRef_ << " for "
            << (ipg.name().size() ? ipg.name() : "others")
            << exit(FatalError);
    }

    Specie::operator+=(ipg);
}


template<class Specie>
inline void Foam::incompressiblePerfectGas<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::incompressiblePerfectGas<Specie> Foam::operator+
(
    const incompressiblePerfectGas<Specie>& ipg1,
    const incompressiblePerfectGas<Specie>& ipg2
)
{
    if (notEqual(ipg1.pRef_, ipg2.pRef_))
    {
        FatalErrorInFunction
            << "pRef " << ipg1.pRef_ << " for "
            << (ipg1.name().size() ? ipg1.name() : "others")
            << " != " << ipg2.pRef_ << " for "
            << (ipg2.name().size() ? ipg2.name() : "others")
            << exit(FatalError);
    }

    return incompressiblePerfectGas<Specie>
    (
        static_cast<const Specie&>(ipg1) + static_cast<const Specie&>(ipg2),
        ipg1.pRef_
    );
}


template<class Specie>
inline Foam::incompressiblePerfectGas<Specie> Foam::operator*
(
    const scalar s,
    const incompressiblePerfectGas<Specie>& ipg
)
{
    return incompressiblePerfectGas<Specie>
    (
        s*static_cast<const Specie&>(ipg),
        ipg.pRef_
    );
}


template<class Specie>
inline Foam::incompressiblePerfectGas<Specie> Foam::operator==
(
    const incompressiblePerfectGas<Specie>& ipg1,
    const incompressiblePerfectGas<Specie>& ipg2
)
{
    return incompressiblePerfectGas<Specie>
    (
        static_cast<const Specie&>(ipg1) == static_cast<const Specie&>(ipg2),
        NaN
    );
}


// ************************************************************************* //
