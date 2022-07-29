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

#include "IrreversibleReaction.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MulticomponentThermo, class ReactionRate>
Foam::IrreversibleReaction<MulticomponentThermo, ReactionRate>::
IrreversibleReaction
(
    const Reaction<MulticomponentThermo>& reaction,
    const ReactionRate& k
)
:
    Reaction<MulticomponentThermo>(reaction),
    k_(k)
{}


template<class MulticomponentThermo, class ReactionRate>
Foam::IrreversibleReaction<MulticomponentThermo, ReactionRate>::
IrreversibleReaction
(
    const speciesTable& species,
    const HashPtrTable<MulticomponentThermo>& thermoDatabase,
    const dictionary& dict
)
:
    Reaction<MulticomponentThermo>(species, thermoDatabase, dict),
    k_(species, dict)
{}


template<class MulticomponentThermo, class ReactionRate>
Foam::IrreversibleReaction<MulticomponentThermo, ReactionRate>::
IrreversibleReaction
(
    const speciesTable& species,
    const HashPtrTable<MulticomponentThermo>& thermoDatabase,
    const objectRegistry& ob,
    const dictionary& dict
)
:
    Reaction<MulticomponentThermo>(species, thermoDatabase, dict),
    k_(species, ob, dict)
{}


template<class MulticomponentThermo, class ReactionRate>
Foam::IrreversibleReaction<MulticomponentThermo, ReactionRate>::
IrreversibleReaction
(
    const IrreversibleReaction<MulticomponentThermo,ReactionRate>& irr,
    const speciesTable& species
)
:
    Reaction<MulticomponentThermo>(irr, species),
    k_(irr.k_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MulticomponentThermo, class ReactionRate>
void Foam::IrreversibleReaction<MulticomponentThermo, ReactionRate>::
preEvaluate() const
{
    k_.preEvaluate();
}


template<class MulticomponentThermo, class ReactionRate>
void Foam::IrreversibleReaction<MulticomponentThermo, ReactionRate>::
postEvaluate() const
{
    k_.postEvaluate();
}


template<class MulticomponentThermo, class ReactionRate>
Foam::scalar Foam::IrreversibleReaction<MulticomponentThermo, ReactionRate>::kf
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return k_(p, T, c, li);
}


template<class MulticomponentThermo, class ReactionRate>
Foam::scalar Foam::IrreversibleReaction<MulticomponentThermo, ReactionRate>::kr
(
    const scalar kfwd,
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return 0;
}


template<class MulticomponentThermo, class ReactionRate>
Foam::scalar Foam::IrreversibleReaction<MulticomponentThermo, ReactionRate>::kr
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return 0;
}


template<class MulticomponentThermo, class ReactionRate>
Foam::scalar
Foam::IrreversibleReaction<MulticomponentThermo, ReactionRate>::dkfdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return k_.ddT(p, T, c, li);
}


template<class MulticomponentThermo, class ReactionRate>
Foam::scalar
Foam::IrreversibleReaction<MulticomponentThermo, ReactionRate>::dkrdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    const scalar dkfdT,
    const scalar kr
) const
{
    return 0;
}


template<class MulticomponentThermo, class ReactionRate>
bool
Foam::IrreversibleReaction<MulticomponentThermo, ReactionRate>::hasDkdc() const
{
    return k_.hasDdc();
}


template<class MulticomponentThermo, class ReactionRate>
void Foam::IrreversibleReaction<MulticomponentThermo, ReactionRate>::dkfdc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dkfdc
) const
{
    k_.ddc(p, T, c, li, dkfdc);
}


template<class MulticomponentThermo, class ReactionRate>
void Foam::IrreversibleReaction<MulticomponentThermo, ReactionRate>::dkrdc
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
    dkrdc = 0;
}


template<class MulticomponentThermo, class ReactionRate>
void Foam::IrreversibleReaction<MulticomponentThermo, ReactionRate>::write
(
    Ostream& os
) const
{
    Reaction<MulticomponentThermo>::write(os);
    k_.write(os);
}


// ************************************************************************* //
