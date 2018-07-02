/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "ReactionProxy.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::ReactionProxy<ReactionThermo>::ReactionProxy
(
    const speciesTable& species,
    const List<typename Reaction<ReactionThermo>::specieCoeffs>& lhs,
    const List<typename Reaction<ReactionThermo>::specieCoeffs>& rhs,
    const HashPtrTable<ReactionThermo>& thermoDatabase
)
:
    Reaction<ReactionThermo>
    (
        species,
        lhs,
        rhs,
        thermoDatabase
    )
{}


template<class ReactionThermo>
Foam::ReactionProxy<ReactionThermo>::ReactionProxy
(
    const Reaction<ReactionThermo>& r,
    const speciesTable& species
)
:
    Reaction<ReactionThermo>
    (
        r,
        species
    )
{}


template<class ReactionThermo>
Foam::ReactionProxy<ReactionThermo>::ReactionProxy
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    const dictionary& dict
)
:
    Reaction<ReactionThermo>
    (
        species,
        thermoDatabase,
        dict
    )
{}


template<class ReactionThermo>
Foam::autoPtr<Foam::Reaction<ReactionThermo>>
Foam::ReactionProxy<ReactionThermo>::clone() const
{
    NotImplemented;
    return autoPtr<Foam::Reaction<ReactionThermo>>();
}


template<class ReactionThermo>
Foam::autoPtr<Foam::Reaction<ReactionThermo>>
Foam::ReactionProxy<ReactionThermo>::clone
(
    const speciesTable& species
) const
{
    NotImplemented;
    return autoPtr<Reaction<ReactionThermo>>();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::scalar Foam::ReactionProxy<ReactionThermo>::kf
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    NotImplemented;
    return 0;
}


template<class ReactionThermo>
Foam::scalar Foam::ReactionProxy<ReactionThermo>::kr
(
    const scalar kfwd,
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    NotImplemented;
    return 0;
}


template<class ReactionThermo>
Foam::scalar Foam::ReactionProxy<ReactionThermo>::kr
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    NotImplemented;
    return 0;
}


template<class ReactionThermo>
Foam::scalar Foam::ReactionProxy<ReactionThermo>::dkfdT
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    NotImplemented;
    return 0;
}


template<class ReactionThermo>
Foam::scalar Foam::ReactionProxy<ReactionThermo>::dkrdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const scalar dkfdT,
    const scalar kr
) const
{
    NotImplemented;
    return 0;
}


template<class ReactionThermo>
const Foam::List<Foam::Tuple2<Foam::label, Foam::scalar>>&
Foam::ReactionProxy<ReactionThermo>::beta() const
{
    NotImplemented;
    return NullObjectRef<List<Tuple2<label, scalar>>>();
}


template<class ReactionThermo>
void Foam::ReactionProxy<ReactionThermo>::dcidc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    scalarField& dcidc
) const
{
    NotImplemented;
}


template<class ReactionThermo>
Foam::scalar Foam::ReactionProxy<ReactionThermo>::dcidT
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    NotImplemented;
    return 0;
}


// ************************************************************************* //
