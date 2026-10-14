/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#include "reaction.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word& Foam::reaction::name() const
{
    return name_;
}


inline const Foam::List<Foam::specieCoeffs>& Foam::reaction::lhs() const
{
    return lhs_;
}


inline const Foam::List<Foam::specieCoeffs>& Foam::reaction::rhs() const
{
    return rhs_;
}


inline const Foam::speciesTable& Foam::reaction::species() const
{
    return species_;
}


inline void Foam::reaction::state::reset
(
    const scalar p_,
    const scalar T_,
    const label li_
)
{
    p = p_;
    T = T_;
    li = li_;

    lnp = log(p);
    rT = 1/T;
    lnT = log(T);
    C = sum(c);
}


inline void Foam::reaction::state::reset
(
    const scalar p_,
    const scalar T_,
    const scalarField& c_,
    const label li_
)
{
    c = c_;
    reset(p_, T_, li_);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<(Foam::Ostream& os, const reaction& r)
{
    r.write(os);
    return os;
}


// ************************************************************************* //
