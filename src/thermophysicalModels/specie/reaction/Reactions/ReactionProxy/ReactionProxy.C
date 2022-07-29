/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2022 OpenFOAM Foundation
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

template<class MulticomponentThermo>
Foam::ReactionProxy<MulticomponentThermo>::ReactionProxy
(
    const speciesTable& species,
    const List<specieCoeffs>& lhs,
    const List<specieCoeffs>& rhs,
    const HashPtrTable<MulticomponentThermo>& thermoDatabase
)
:
    Reaction<MulticomponentThermo>
    (
        species,
        lhs,
        rhs,
        thermoDatabase
    )
{}


template<class MulticomponentThermo>
Foam::ReactionProxy<MulticomponentThermo>::ReactionProxy
(
    const Reaction<MulticomponentThermo>& r,
    const speciesTable& species
)
:
    Reaction<MulticomponentThermo>
    (
        r,
        species
    )
{}


template<class MulticomponentThermo>
Foam::ReactionProxy<MulticomponentThermo>::ReactionProxy
(
    const speciesTable& species,
    const HashPtrTable<MulticomponentThermo>& thermoDatabase,
    const dictionary& dict
)
:
    Reaction<MulticomponentThermo>
    (
        species,
        thermoDatabase,
        dict
    )
{}


template<class MulticomponentThermo>
Foam::autoPtr<Foam::Reaction<MulticomponentThermo>>
Foam::ReactionProxy<MulticomponentThermo>::clone() const
{
    NotImplemented;
    return autoPtr<Foam::Reaction<MulticomponentThermo>>();
}


template<class MulticomponentThermo>
Foam::autoPtr<Foam::Reaction<MulticomponentThermo>>
Foam::ReactionProxy<MulticomponentThermo>::clone
(
    const speciesTable& species
) const
{
    NotImplemented;
    return autoPtr<Reaction<MulticomponentThermo>>();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MulticomponentThermo>
void Foam::ReactionProxy<MulticomponentThermo>::preEvaluate() const
{}


template<class MulticomponentThermo>
void Foam::ReactionProxy<MulticomponentThermo>::postEvaluate() const
{}


template<class MulticomponentThermo>
Foam::scalar Foam::ReactionProxy<MulticomponentThermo>::kf
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    NotImplemented;
    return 0;
}


template<class MulticomponentThermo>
Foam::scalar Foam::ReactionProxy<MulticomponentThermo>::kr
(
    const scalar kfwd,
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    NotImplemented;
    return 0;
}


template<class MulticomponentThermo>
Foam::scalar Foam::ReactionProxy<MulticomponentThermo>::kr
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    NotImplemented;
    return 0;
}


template<class MulticomponentThermo>
Foam::scalar Foam::ReactionProxy<MulticomponentThermo>::dkfdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    NotImplemented;
    return 0;
}


template<class MulticomponentThermo>
Foam::scalar Foam::ReactionProxy<MulticomponentThermo>::dkrdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    const scalar dkfdT,
    const scalar kr
) const
{
    NotImplemented;
    return 0;
}


template<class MulticomponentThermo>
bool Foam::ReactionProxy<MulticomponentThermo>::hasDkdc() const
{
    NotImplemented;
    return false;
}


template<class MulticomponentThermo>
void Foam::ReactionProxy<MulticomponentThermo>::dkfdc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dkfdc
) const
{
    NotImplemented;
}


template<class MulticomponentThermo>
void Foam::ReactionProxy<MulticomponentThermo>::dkrdc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    const scalarField& dkfdc,
    const scalar kr,
    scalarField& dkrdc
) const
{
    NotImplemented;
}


// ************************************************************************* //
