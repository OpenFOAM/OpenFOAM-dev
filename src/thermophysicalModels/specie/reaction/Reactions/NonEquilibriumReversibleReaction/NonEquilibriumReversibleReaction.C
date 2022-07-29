/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "NonEquilibriumReversibleReaction.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MulticomponentThermo, class ReactionRate>
Foam::NonEquilibriumReversibleReaction<MulticomponentThermo, ReactionRate>::
NonEquilibriumReversibleReaction
(
    const Reaction<MulticomponentThermo>& reaction,
    const ReactionRate& forwardReactionRate,
    const ReactionRate& reverseReactionRate
)
:
    Reaction<MulticomponentThermo>(reaction),
    fk_(forwardReactionRate),
    rk_(reverseReactionRate)
{}


template<class MulticomponentThermo, class ReactionRate>
Foam::NonEquilibriumReversibleReaction<MulticomponentThermo, ReactionRate>::
NonEquilibriumReversibleReaction
(
    const speciesTable& species,
    const HashPtrTable<MulticomponentThermo>& thermoDatabase,
    const dictionary& dict
)
:
    Reaction<MulticomponentThermo>(species, thermoDatabase, dict),
    fk_(species, dict.subDict("forward")),
    rk_(species, dict.subDict("reverse"))
{}


template<class MulticomponentThermo, class ReactionRate>
Foam::NonEquilibriumReversibleReaction<MulticomponentThermo, ReactionRate>::
NonEquilibriumReversibleReaction
(
    const speciesTable& species,
    const HashPtrTable<MulticomponentThermo>& thermoDatabase,
    const objectRegistry& ob,
    const dictionary& dict
)
:
    Reaction<MulticomponentThermo>(species, thermoDatabase, dict),
    fk_(species, ob, dict.subDict("forward")),
    rk_(species, ob, dict.subDict("reverse"))
{}


template<class MulticomponentThermo, class ReactionRate>
Foam::NonEquilibriumReversibleReaction<MulticomponentThermo, ReactionRate>::
NonEquilibriumReversibleReaction
(
    const NonEquilibriumReversibleReaction<MulticomponentThermo, ReactionRate>&
        nerr,
    const speciesTable& species
)
:
    Reaction<MulticomponentThermo>(nerr, species),
    fk_(nerr.fk_),
    rk_(nerr.rk_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MulticomponentThermo, class ReactionRate>
void
Foam::NonEquilibriumReversibleReaction<MulticomponentThermo, ReactionRate>::
preEvaluate() const
{
    fk_.preEvaluate();
    rk_.preEvaluate();
}


template<class MulticomponentThermo, class ReactionRate>
void
Foam::NonEquilibriumReversibleReaction<MulticomponentThermo, ReactionRate>::
postEvaluate() const
{
    fk_.postEvaluate();
    rk_.postEvaluate();
}


template<class MulticomponentThermo, class ReactionRate>
Foam::scalar
Foam::NonEquilibriumReversibleReaction<MulticomponentThermo, ReactionRate>::kf
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return fk_(p, T, c, li);
}


template<class MulticomponentThermo, class ReactionRate>
Foam::scalar
Foam::NonEquilibriumReversibleReaction<MulticomponentThermo, ReactionRate>::kr
(
    const scalar,
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return rk_(p, T, c, li);
}


template<class MulticomponentThermo, class ReactionRate>
Foam::scalar
Foam::NonEquilibriumReversibleReaction<MulticomponentThermo, ReactionRate>::kr
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return rk_(p, T, c, li);
}


template<class MulticomponentThermo, class ReactionRate>
Foam::scalar
Foam::NonEquilibriumReversibleReaction<MulticomponentThermo, ReactionRate>::
dkfdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return fk_.ddT(p, T, c, li);
}


template<class MulticomponentThermo, class ReactionRate>
Foam::scalar
Foam::NonEquilibriumReversibleReaction<MulticomponentThermo, ReactionRate>::
dkrdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    const scalar dkfdT,
    const scalar kr
) const
{
    return rk_.ddT(p, T, c, li);
}


template<class MulticomponentThermo, class ReactionRate>
bool
Foam::NonEquilibriumReversibleReaction<MulticomponentThermo, ReactionRate>::
hasDkdc() const
{
    return fk_.hasDdc() || rk_.hasDdc();
}


template<class MulticomponentThermo, class ReactionRate>
void
Foam::NonEquilibriumReversibleReaction<MulticomponentThermo, ReactionRate>::
dkfdc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dkfdc
) const
{
    fk_.ddc(p, T, c, li, dkfdc);
}


template<class MulticomponentThermo, class ReactionRate>
void
Foam::NonEquilibriumReversibleReaction<MulticomponentThermo, ReactionRate>::
dkrdc
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
    rk_.ddc(p, T, c, li, dkrdc);
}


template<class MulticomponentThermo, class ReactionRate>
void
Foam::NonEquilibriumReversibleReaction<MulticomponentThermo, ReactionRate>::
write
(
    Ostream& os
) const
{
    Reaction<MulticomponentThermo>::write(os);

    os  << indent << "forward" << nl;
    os  << indent << token::BEGIN_BLOCK << nl;
    os  << incrIndent;
    fk_.write(os);
    os  << decrIndent;
    os  << indent << token::END_BLOCK << nl;

    os  << indent << "reverse" << nl;
    os  << indent << token::BEGIN_BLOCK << nl;
    os  << incrIndent;
    rk_.write(os);
    os  << decrIndent;
    os  << indent << token::END_BLOCK << nl;
}


// ************************************************************************* //
