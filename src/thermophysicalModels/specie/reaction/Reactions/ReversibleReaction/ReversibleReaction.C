/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

template<class ReactionThermo, class ReactionRate>
Foam::ReversibleReaction<ReactionThermo, ReactionRate>::
ReversibleReaction
(
    const Reaction<ReactionThermo>& reaction,
    const ReactionRate& k
)
:
    Reaction<ReactionThermo>(reaction),
    k_(k)
{}


template<class ReactionThermo, class ReactionRate>
Foam::ReversibleReaction<ReactionThermo, ReactionRate>::
ReversibleReaction
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    const dictionary& dict
)
:
    Reaction<ReactionThermo>(species, thermoDatabase, dict),
    k_(species, dict)
{}


template<class ReactionThermo, class ReactionRate>
Foam::ReversibleReaction<ReactionThermo, ReactionRate>::
ReversibleReaction
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    const objectRegistry& ob,
    const dictionary& dict
)
:
    Reaction<ReactionThermo>(species, thermoDatabase, dict),
    k_(species, ob, dict)
{}


template<class ReactionThermo, class ReactionRate>
Foam::ReversibleReaction<ReactionThermo, ReactionRate>::
ReversibleReaction
(
    const ReversibleReaction<ReactionThermo, ReactionRate>& rr,
    const speciesTable& species
)
:
    Reaction<ReactionThermo>(rr, species),
    k_(rr.k_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class ReactionRate>
void
Foam::ReversibleReaction<ReactionThermo, ReactionRate>::preEvaluate() const
{
    k_.preEvaluate();
}


template<class ReactionThermo, class ReactionRate>
void
Foam::ReversibleReaction<ReactionThermo, ReactionRate>::postEvaluate() const
{
    k_.postEvaluate();
}


template<class ReactionThermo, class ReactionRate>
Foam::scalar Foam::ReversibleReaction<ReactionThermo, ReactionRate>::kf
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return k_(p, T, c, li);
}


template<class ReactionThermo, class ReactionRate>
Foam::scalar Foam::ReversibleReaction<ReactionThermo, ReactionRate>::kr
(
    const scalar kfwd,
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return kfwd/max(this->Kc(p, T), rootSmall);
}


template<class ReactionThermo, class ReactionRate>
Foam::scalar Foam::ReversibleReaction<ReactionThermo, ReactionRate>::kr
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return kr(kf(p, T, c, li), p, T, c, li);
}


template<class ReactionThermo, class ReactionRate>
Foam::scalar Foam::ReversibleReaction<ReactionThermo, ReactionRate>::dkfdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return k_.ddT(p, T, c, li);
}


template<class ReactionThermo, class ReactionRate>
Foam::scalar Foam::ReversibleReaction<ReactionThermo, ReactionRate>::dkrdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    const scalar dkfdT,
    const scalar kr
) const
{
    scalar Kc = max(this->Kc(p, T), rootSmall);

    return dkfdT/Kc - kr*this->dKcdTbyKc(p, T);
}


template<class ReactionThermo, class ReactionRate>
const Foam::List<Foam::Tuple2<Foam::label, Foam::scalar>>&
Foam::ReversibleReaction<ReactionThermo, ReactionRate>::beta() const
{
    return k_.beta();
}


template<class ReactionThermo, class ReactionRate>
void Foam::ReversibleReaction<ReactionThermo, ReactionRate>::dcidc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dcidc
) const
{
    k_.dcidc(p, T, c, li, dcidc);
}


template<class ReactionThermo, class ReactionRate>
Foam::scalar Foam::ReversibleReaction<ReactionThermo, ReactionRate>::dcidT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return k_.dcidT(p, T, c, li);
}


template<class ReactionThermo, class ReactionRate>
void Foam::ReversibleReaction<ReactionThermo, ReactionRate>::write
(
    Ostream& os
) const
{
    Reaction<ReactionThermo>::write(os);
    k_.write(os);
}


// ************************************************************************* //
