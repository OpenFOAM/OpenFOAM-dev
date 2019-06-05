/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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
inline Foam::scalar Foam::incompressiblePerfectGas<Specie>::S
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
    scalar Y1 = this->Y();
    Specie::operator+=(ipg);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        const scalar Y2 = ipg.Y()/this->Y();

        pRef_ = Y1*pRef_ + Y2*ipg.pRef_;
    }
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
    Specie sp
    (
        static_cast<const Specie&>(ipg1)
      + static_cast<const Specie&>(ipg2)
    );

    if (mag(sp.Y()) < small)
    {
        return incompressiblePerfectGas<Specie>
        (
            sp,
            ipg1.pRef_
        );
    }
    else
    {
        const scalar Y1 = ipg1.Y()/sp.Y();
        const scalar Y2 = ipg2.Y()/sp.Y();

        return incompressiblePerfectGas<Specie>
        (
            sp,
            Y1*ipg1.pRef_ + Y2*ipg2.pRef_
        );
    }
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
    Specie sp
    (
        static_cast<const Specie&>(ipg1)
     == static_cast<const Specie&>(ipg2)
    );

    const scalar Y1 = ipg1.Y()/sp.Y();
    const scalar Y2 = ipg2.Y()/sp.Y();

    return incompressiblePerfectGas<Specie>
    (
        sp,
        Y2*ipg2.pRef_ - Y1*ipg1.pRef_
    );
}


// ************************************************************************* //
