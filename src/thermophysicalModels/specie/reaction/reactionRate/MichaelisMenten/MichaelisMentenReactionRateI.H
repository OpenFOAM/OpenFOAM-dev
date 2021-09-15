/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2021 OpenFOAM Foundation
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

inline Foam::MichaelisMentenReactionRate::MichaelisMentenReactionRate
(
    const speciesTable& st,
    const dictionary& dict
)
:
    species_(st),
    Vmax_(dict.lookup<scalar>("Vmax")),
    Km_(dict.lookup<scalar>("Km")),
    s_(st[dict.lookup("S")])
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::MichaelisMentenReactionRate::preEvaluate() const
{}


inline void Foam::MichaelisMentenReactionRate::postEvaluate() const
{}


inline Foam::scalar Foam::MichaelisMentenReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return Vmax_/(Km_ + c[s_]);
}


inline Foam::scalar Foam::MichaelisMentenReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return 0;
}


inline bool Foam::MichaelisMentenReactionRate::hasDdc() const
{
    return true;
}


inline void Foam::MichaelisMentenReactionRate::ddc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& ddc
) const
{
    ddc = 0;
    ddc[s_] = - Vmax_/sqr(Km_ + c[s_]);
}


inline void Foam::MichaelisMentenReactionRate::write(Ostream& os) const
{
    writeEntry(os, "Vmax", Vmax_);
    writeEntry(os, "Km", Km_);
    writeEntry(os, "S", species_[s_]);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const MichaelisMentenReactionRate& mmrr
)
{
    mmrr.write(os);
    return os;
}


// ************************************************************************* //
