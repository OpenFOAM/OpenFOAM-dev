/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2023 OpenFOAM Foundation
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

template<class ThermoType>
Foam::ReactionProxy<ThermoType>::ReactionProxy
(
    const speciesTable& species,
    const PtrList<ThermoType>& speciesThermo,
    const List<specieCoeffs>& lhs,
    const List<specieCoeffs>& rhs
)
:
    Reaction<ThermoType>(species, speciesThermo, lhs, rhs)
{}


template<class ThermoType>
Foam::autoPtr<Foam::Reaction<ThermoType>>
Foam::ReactionProxy<ThermoType>::clone() const
{
    NotImplemented;
    return autoPtr<Foam::Reaction<ThermoType>>();
}


template<class ThermoType>
Foam::autoPtr<Foam::Reaction<ThermoType>>
Foam::ReactionProxy<ThermoType>::clone
(
    const speciesTable& species
) const
{
    NotImplemented;
    return autoPtr<Reaction<ThermoType>>();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::ReactionProxy<ThermoType>::preEvaluate() const
{}


template<class ThermoType>
void Foam::ReactionProxy<ThermoType>::postEvaluate() const
{}


template<class ThermoType>
Foam::scalar Foam::ReactionProxy<ThermoType>::kf
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


template<class ThermoType>
Foam::scalar Foam::ReactionProxy<ThermoType>::kr
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


template<class ThermoType>
Foam::scalar Foam::ReactionProxy<ThermoType>::kr
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


template<class ThermoType>
Foam::scalar Foam::ReactionProxy<ThermoType>::dkfdT
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


template<class ThermoType>
Foam::scalar Foam::ReactionProxy<ThermoType>::dkrdT
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


template<class ThermoType>
bool Foam::ReactionProxy<ThermoType>::hasDkdc() const
{
    NotImplemented;
    return false;
}


template<class ThermoType>
void Foam::ReactionProxy<ThermoType>::dkfdc
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


template<class ThermoType>
void Foam::ReactionProxy<ThermoType>::dkrdc
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
