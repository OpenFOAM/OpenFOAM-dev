/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2026 OpenFOAM Foundation
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

#include "MichaelisMenten_reactionRate.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::reactionRates::MichaelisMenten::MichaelisMenten
(
    const speciesTable& st,
    const dimensionSet& dims,
    const dictionary& dict
)
:
    species_(st),
    Vmax_
    (
        dict.lookup<scalar>("Vmax", dims*dimensions::moles/dimensions::volume)
    ),
    Km_(dict.lookup<scalar>("Km", dimensions::moles/dimensions::volume)),
    s_(st[dict.lookup("S")])
{}


inline Foam::reactionRates::MichaelisMenten::MichaelisMenten
(
    const speciesTable& species,
    const objectRegistry& ob,
    const dimensionSet& dims,
    const dictionary& dict
)
:
    MichaelisMenten(species, dims, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::reactionRates::MichaelisMenten::preEvaluate() const
{}


inline void Foam::reactionRates::MichaelisMenten::postEvaluate() const
{}


inline Foam::scalar Foam::reactionRates::MichaelisMenten::k
(
    const reaction::state& s
) const
{
    return Vmax_/(Km_ + s.c[s_]);
}


inline Foam::scalar Foam::reactionRates::MichaelisMenten::ddT
(
    const reaction::state& s
) const
{
    return 0;
}


inline void Foam::reactionRates::MichaelisMenten::ddc
(
    const reaction::state& s,
    scalarField& ddc
) const
{
    ddc = 0;
    ddc[s_] = -Vmax_/sqr(Km_ + s.c[s_]);
}


inline void Foam::reactionRates::MichaelisMenten::write(Ostream& os) const
{
    writeEntry(os, "Vmax", Vmax_);
    writeEntry(os, "Km", Km_);
    writeEntry(os, "S", species_[s_]);
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const reactionRates::MichaelisMenten& rr
)
{
    rr.write(os);
    return os;
}


// ************************************************************************* //
