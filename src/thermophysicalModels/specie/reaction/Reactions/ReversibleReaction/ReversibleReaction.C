/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::ReversibleReaction<ReactionType, ReactionThermo, ReactionRate>::
ReversibleReaction
(
    const ReactionType<ReactionThermo>& reaction,
    const ReactionRate& k
)
:
    ReactionType<ReactionThermo>(reaction),
    k_(k)
{}


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::ReversibleReaction<ReactionType, ReactionThermo, ReactionRate>::
ReversibleReaction
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    const dictionary& dict
)
:
    ReactionType<ReactionThermo>(species, thermoDatabase, dict),
    k_(species, dict)
{}


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::ReversibleReaction<ReactionType, ReactionThermo, ReactionRate>::
ReversibleReaction
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    const objectRegistry& ob,
    const dictionary& dict
)
:
    ReactionType<ReactionThermo>(species, thermoDatabase, dict),
    k_(species, ob, dict)
{}


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::ReversibleReaction<ReactionType, ReactionThermo, ReactionRate>::
ReversibleReaction
(
    const ReversibleReaction<ReactionType, ReactionThermo, ReactionRate>& rr,
    const speciesTable& species
)
:
    ReactionType<ReactionThermo>(rr, species),
    k_(rr.k_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::scalar Foam::ReversibleReaction
<
    ReactionType,
    ReactionThermo,
    ReactionRate
>::kf
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return k_(p, T, c, li);
}


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::scalar Foam::ReversibleReaction
<
    ReactionType,
    ReactionThermo,
    ReactionRate
>::kr
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


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::scalar Foam::ReversibleReaction
<
    ReactionType,
    ReactionThermo,
    ReactionRate
>::kr
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return kr(kf(p, T, c, li), p, T, c, li);
}


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::scalar Foam::ReversibleReaction
<
    ReactionType,
    ReactionThermo,
    ReactionRate
>::dkfdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return k_.ddT(p, T, c, li);
}


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::scalar Foam::ReversibleReaction
<
    ReactionType,
    ReactionThermo,
    ReactionRate
>::dkrdT
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


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
const Foam::List<Foam::Tuple2<Foam::label, Foam::scalar>>&
Foam::ReversibleReaction
<
    ReactionType,
    ReactionThermo,
    ReactionRate
>::beta() const
{
    return k_.beta();
}


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
void Foam::ReversibleReaction
<
    ReactionType,
    ReactionThermo,
    ReactionRate
>::dcidc
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


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::scalar Foam::ReversibleReaction
<
    ReactionType,
    ReactionThermo,
    ReactionRate
>::dcidT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return k_.dcidT(p, T, c, li);
}


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
void Foam::ReversibleReaction
<
    ReactionType,
    ReactionThermo,
    ReactionRate
>::write
(
    Ostream& os
) const
{
    Reaction<ReactionThermo>::write(os);
    k_.write(os);
}


// ************************************************************************* //
