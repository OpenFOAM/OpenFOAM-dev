/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "Reaction.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
inline void Foam::Reaction<ThermoType>::C
(
    const state& s,
    scalar& Cf,
    scalar& Cr
) const
{
    Cf = pow(s.c[lhs()[0].index], lhs()[0].exponent);
    for(label i=1; i<lhs().size(); i++)
    {
        Cf *= pow(s.c[lhs()[i].index], lhs()[i].exponent);
    }

    Cr = pow(s.c[rhs()[0].index], rhs()[0].exponent);
    for(label i=1; i<rhs().size(); i++)
    {
        Cr *= pow(s.c[rhs()[i].index], rhs()[i].exponent);
    }
}


template<class ThermoType>
inline Foam::scalar Foam::Reaction<ThermoType>::omega
(
    const state& s,
    scalar& omegaf,
    scalar& omegar
) const
{
    C(s, omegaf, omegar);

    const scalar kf = this->kf(s);
    omegaf = kf*omegaf;
    omegar = kr(kf, s)*omegar;

    return omegaf - omegar;
}


template<class ThermoType>
inline Foam::scalar Foam::Reaction<ThermoType>::omega(const state& s) const
{
    scalar omegaf, omegar;
    return omega(s, omegaf, omegar);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class ThermoType>
inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const Reaction<ThermoType>& r
)
{
    r.write(os);
    return os;
}


// ************************************************************************* //
