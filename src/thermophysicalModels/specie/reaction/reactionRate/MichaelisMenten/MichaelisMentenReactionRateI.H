/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    Vmax_(readScalar(dict.lookup("Vmax"))),
    Km_(readScalar(dict.lookup("Km"))),
    s_(st[dict.lookup("S")])
{
    beta_.append(Tuple2<label, scalar>(s_, 1.0));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::MichaelisMentenReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return Vmax_/(Km_ + c[s_]);
}


inline Foam::scalar Foam::MichaelisMentenReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return 0;
}


inline const Foam::List<Foam::Tuple2<Foam::label, Foam::scalar>>&
Foam::MichaelisMentenReactionRate::beta() const
{
    return beta_;
}


inline void Foam::MichaelisMentenReactionRate::dcidc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    scalarField& dcidc
) const
{
    dcidc[0] = -1.0/(Km_ + c[s_]);
}


inline Foam::scalar Foam::MichaelisMentenReactionRate::dcidT
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return 0;
}


inline void Foam::MichaelisMentenReactionRate::write(Ostream& os) const
{
    os.writeKeyword("Vmax") << Vmax_ << token::END_STATEMENT << nl;
    os.writeKeyword("Km") << Km_ << token::END_STATEMENT << nl;
    os.writeKeyword("S") << species_[s_] << token::END_STATEMENT << nl;
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
