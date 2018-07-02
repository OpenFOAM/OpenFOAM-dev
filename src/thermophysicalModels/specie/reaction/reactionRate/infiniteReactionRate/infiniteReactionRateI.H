/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

inline Foam::infiniteReactionRate::infiniteReactionRate()
{}


inline Foam::infiniteReactionRate::infiniteReactionRate
(
    const speciesTable&,
    const dictionary&
)
{}


inline void Foam::infiniteReactionRate::write(Ostream& os) const
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::infiniteReactionRate::operator()
(
    const scalar p,
    const scalar,
    const scalarField&
) const
{
    return (1);
}

inline Foam::scalar Foam::infiniteReactionRate::ddT
(
    const scalar p,
    const scalar,
    const scalarField&
) const
{
    return (0);
}


inline const Foam::List<Foam::Tuple2<Foam::label, Foam::scalar>>&
Foam::infiniteReactionRate::beta() const
{
    return NullObjectRef<List<Tuple2<label, scalar>>>();
}


inline void Foam::infiniteReactionRate::dcidc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    scalarField& dcidc
) const
{}


inline Foam::scalar Foam::infiniteReactionRate::dcidT
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return 0;
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const infiniteReactionRate& rr
)
{
    rr.write(os);
    return os;
}


// ************************************************************************* //
