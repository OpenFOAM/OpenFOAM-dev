/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

template<class ThermoType, class ReactionRate>
Foam::NonEquilibriumReversibleReaction<ThermoType, ReactionRate>::
NonEquilibriumReversibleReaction
(
    const Reaction<ThermoType>& reaction,
    const ReactionRate& forwardReactionRate,
    const ReactionRate& reverseReactionRate
)
:
    Reaction<ThermoType>(reaction),
    kf_(forwardReactionRate),
    kr_(reverseReactionRate)
{}


template<class ThermoType, class ReactionRate>
Foam::NonEquilibriumReversibleReaction<ThermoType, ReactionRate>::
NonEquilibriumReversibleReaction
(
    const speciesTable& species,
    const PtrList<ThermoType>& speciesThermo,
    const dictionary& dict
)
:
    Reaction<ThermoType>(species, speciesThermo, dict),
    kf_(species, this->kfDims(), dict.subDict("forward")),
    kr_(species, this->krDims(), dict.subDict("reverse"))
{}


template<class ThermoType, class ReactionRate>
Foam::NonEquilibriumReversibleReaction<ThermoType, ReactionRate>::
NonEquilibriumReversibleReaction
(
    const speciesTable& species,
    const PtrList<ThermoType>& speciesThermo,
    const objectRegistry& ob,
    const dictionary& dict
)
:
    Reaction<ThermoType>(species, speciesThermo, dict),
    kf_(species, ob, this->kfDims(), dict.subDict("forward")),
    kr_(species, ob, this->krDims(), dict.subDict("reverse"))
{}


template<class ThermoType, class ReactionRate>
Foam::NonEquilibriumReversibleReaction<ThermoType, ReactionRate>::
NonEquilibriumReversibleReaction
(
    const NonEquilibriumReversibleReaction<ThermoType, ReactionRate>&
        nerr,
    const speciesTable& species
)
:
    Reaction<ThermoType>(nerr, species),
    kf_(nerr.kf_),
    kr_(nerr.kr_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType, class ReactionRate>
void Foam::NonEquilibriumReversibleReaction<ThermoType, ReactionRate>::
preEvaluate() const
{
    kf_.preEvaluate();
    kr_.preEvaluate();
}


template<class ThermoType, class ReactionRate>
void Foam::NonEquilibriumReversibleReaction<ThermoType, ReactionRate>::
postEvaluate() const
{
    kf_.postEvaluate();
    kr_.postEvaluate();
}


template<class ThermoType, class ReactionRate>
Foam::scalar
Foam::NonEquilibriumReversibleReaction<ThermoType, ReactionRate>::kf
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return kf_(p, T, c, li);
}


template<class ThermoType, class ReactionRate>
Foam::scalar
Foam::NonEquilibriumReversibleReaction<ThermoType, ReactionRate>::kr
(
    const scalar,
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return kr_(p, T, c, li);
}


template<class ThermoType, class ReactionRate>
Foam::scalar
Foam::NonEquilibriumReversibleReaction<ThermoType, ReactionRate>::kr
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return kr_(p, T, c, li);
}


template<class ThermoType, class ReactionRate>
Foam::scalar
Foam::NonEquilibriumReversibleReaction<ThermoType, ReactionRate>::dkfdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return kf_.ddT(p, T, c, li);
}


template<class ThermoType, class ReactionRate>
Foam::scalar
Foam::NonEquilibriumReversibleReaction<ThermoType, ReactionRate>::dkrdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    const scalar dkfdT,
    const scalar kr
) const
{
    return kr_.ddT(p, T, c, li);
}


template<class ThermoType, class ReactionRate>
bool Foam::NonEquilibriumReversibleReaction<ThermoType, ReactionRate>::
hasDkdc() const
{
    return kf_.hasDdc() || kr_.hasDdc();
}


template<class ThermoType, class ReactionRate>
void Foam::NonEquilibriumReversibleReaction<ThermoType, ReactionRate>::dkfdc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dkfdc
) const
{
    kf_.ddc(p, T, c, li, dkfdc);
}


template<class ThermoType, class ReactionRate>
void Foam::NonEquilibriumReversibleReaction<ThermoType, ReactionRate>::dkrdc
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
    kr_.ddc(p, T, c, li, dkrdc);
}


template<class ThermoType, class ReactionRate>
void Foam::NonEquilibriumReversibleReaction<ThermoType, ReactionRate>::write
(
    Ostream& os
) const
{
    Reaction<ThermoType>::write(os);

    os  << indent << "forward" << nl;
    os  << indent << token::BEGIN_BLOCK << nl;
    os  << incrIndent;
    kf_.write(os);
    os  << decrIndent;
    os  << indent << token::END_BLOCK << nl;

    os  << indent << "reverse" << nl;
    os  << indent << token::BEGIN_BLOCK << nl;
    os  << incrIndent;
    kr_.write(os);
    os  << decrIndent;
    os  << indent << token::END_BLOCK << nl;
}


// ************************************************************************* //
