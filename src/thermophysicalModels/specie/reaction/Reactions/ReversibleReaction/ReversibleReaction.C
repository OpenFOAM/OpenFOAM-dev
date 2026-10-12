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

#include "ReversibleReaction.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType, class ReactionRate>
Foam::ReversibleReaction<ThermoType, ReactionRate>::ReversibleReaction
(
    const Reaction<ThermoType>& reaction,
    const ReactionRate& k
)
:
    Reaction<ThermoType>(reaction),
    k_(k)
{}


template<class ThermoType, class ReactionRate>
Foam::ReversibleReaction<ThermoType, ReactionRate>::ReversibleReaction
(
    const speciesTable& species,
    const PtrList<ThermoType>& speciesThermo,
    const dictionary& dict
)
:
    Reaction<ThermoType>(species, speciesThermo, dict),
    k_(species, this->kfDims(), dict)
{}


template<class ThermoType, class ReactionRate>
Foam::ReversibleReaction<ThermoType, ReactionRate>::ReversibleReaction
(
    const speciesTable& species,
    const PtrList<ThermoType>& speciesThermo,
    const objectRegistry& ob,
    const dictionary& dict
)
:
    Reaction<ThermoType>(species, speciesThermo, dict),
    k_(species, ob, this->kfDims(), dict)
{}


template<class ThermoType, class ReactionRate>
Foam::ReversibleReaction<ThermoType, ReactionRate>::ReversibleReaction
(
    const ReversibleReaction<ThermoType, ReactionRate>& rr,
    const speciesTable& species
)
:
    Reaction<ThermoType>(rr, species),
    k_(rr.k_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType, class ReactionRate>
void Foam::ReversibleReaction<ThermoType, ReactionRate>::preEvaluate() const
{
    k_.preEvaluate();
}


template<class ThermoType, class ReactionRate>
void Foam::ReversibleReaction<ThermoType, ReactionRate>::postEvaluate() const
{
    k_.postEvaluate();
}


template<class ThermoType, class ReactionRate>
Foam::scalar Foam::ReversibleReaction<ThermoType, ReactionRate>::kf
(
    const state& s
) const
{
    return k_(s);
}


template<class ThermoType, class ReactionRate>
Foam::scalar Foam::ReversibleReaction<ThermoType, ReactionRate>::kr
(
    const scalar kfwd,
    const state& s
) const
{
    return kfwd/this->Kc(s.p, s.T);
}


template<class ThermoType, class ReactionRate>
Foam::scalar Foam::ReversibleReaction<ThermoType, ReactionRate>::kr
(
    const state& s
) const
{
    return kr(kf(s), s);
}


template<class ThermoType, class ReactionRate>
Foam::scalar Foam::ReversibleReaction<ThermoType, ReactionRate>::dkfdT
(
    const state& s
) const
{
    return k_.ddT(s);
}


template<class ThermoType, class ReactionRate>
Foam::scalar Foam::ReversibleReaction<ThermoType, ReactionRate>::dkrdT
(
    const state& s,
    const scalar dkfdT,
    const scalar kr
) const
{
    return dkfdT/this->Kc(s.p, s.T) - kr*this->dKcdTbyKc(s.p, s.T);
}


template<class ThermoType, class ReactionRate>
bool
Foam::ReversibleReaction<ThermoType, ReactionRate>::hasDkdc() const
{
    return k_.hasDdc();
}


template<class ThermoType, class ReactionRate>
void Foam::ReversibleReaction<ThermoType, ReactionRate>::dkfdc
(
    const state& s,
    scalarField& dkfdc
) const
{
    k_.ddc(s, dkfdc);
}


template<class ThermoType, class ReactionRate>
void Foam::ReversibleReaction<ThermoType, ReactionRate>::dkrdc
(
    const state& s,
    const scalarField& dkfdc,
    const scalar kr,
    scalarField& dkrdc
) const
{
    dkrdc = dkfdc/this->Kc(s.p, s.T);
}


template<class ThermoType, class ReactionRate>
void Foam::ReversibleReaction<ThermoType, ReactionRate>::write
(
    Ostream& os
) const
{
    Reaction<ThermoType>::write(os);
    k_.write(os);
}


// ************************************************************************* //
