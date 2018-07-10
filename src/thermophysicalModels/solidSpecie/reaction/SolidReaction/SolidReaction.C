/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

#include "SolidReaction.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::SolidReaction<ReactionThermo>::SolidReaction
(
    const Reaction<ReactionThermo>& reaction,
    const speciesTable& pyrolisisGases,
    const List<specieCoeffs>& glhs,
    const List<specieCoeffs>& grhs
)
:
    Reaction<ReactionThermo>(reaction),
    pyrolisisGases_(pyrolisisGases),
    glhs_(glhs),
    grhs_(grhs)
{}


template<class ReactionThermo>
Foam::SolidReaction<ReactionThermo>::SolidReaction
(
    const SolidReaction<ReactionThermo>& r,
    const speciesTable& pyrolisisGases
)
:
    Reaction<ReactionThermo>(r),
    pyrolisisGases_(pyrolisisGases),
    glhs_(r.glhs_),
    grhs_(r.grhs_)
{}


template<class ReactionThermo>
Foam::SolidReaction<ReactionThermo>::SolidReaction
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    const dictionary& dict
)
:
    Reaction<ReactionThermo>(species, thermoDatabase, dict),
    pyrolisisGases_(dict.parent().parent().lookup("gaseousSpecies")),
    glhs_(),
    grhs_()
{
    this->setLRhs
    (
        IStringStream(dict.lookup("reaction"))(),
        pyrolisisGases_,
        glhs_,
        grhs_
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo>
const Foam::List<typename Foam::SolidReaction<ReactionThermo>::specieCoeffs>&
Foam::SolidReaction<ReactionThermo>::glhs() const
{
    return glhs_;
}


template<class ReactionThermo>
const Foam::List<typename Foam::Reaction<ReactionThermo>::specieCoeffs>&
Foam::SolidReaction<ReactionThermo>::grhs() const
{
    return grhs_;
}


template<class ReactionThermo>
const Foam::speciesTable& Foam::SolidReaction<ReactionThermo>::
gasSpecies() const
{
    return pyrolisisGases_;
}


template<class ReactionThermo>
void Foam::SolidReaction<ReactionThermo>::write(Ostream& os) const
{
    OStringStream reaction;
    os.writeKeyword("reaction") << solidReactionStr(reaction)
        << token::END_STATEMENT << nl;
}


template<class ReactionThermo>
Foam::string Foam::SolidReaction<ReactionThermo>::solidReactionStr
(
    OStringStream& reaction
) const
{
    this->reactionStrLeft(reaction);
    if (glhs().size() > 0)
    {
        reaction << " + ";
        solidReactionStrLeft(reaction);
    }
    reaction << " = ";
    this->reactionStrRight(reaction);
    if (grhs().size() > 0)
    {
        reaction << " + ";
        solidReactionStrRight(reaction);
    }
    return reaction.str();

}


template<class ReactionThermo>
void Foam::SolidReaction<ReactionThermo>::solidReactionStrLeft
(
    OStringStream& reaction
) const
{
    for (label i = 0; i < glhs().size(); ++i)
    {
        if (i > 0)
        {
            reaction << " + ";
        }
        if (mag(glhs()[i].stoichCoeff - 1) > small)
        {
            reaction << glhs()[i].stoichCoeff;
        }
        reaction << gasSpecies()[glhs()[i].index];
        if (mag(glhs()[i].exponent - glhs()[i].stoichCoeff) > small)
        {
            reaction << "^" << glhs()[i].exponent;
        }
    }
}


template<class ReactionThermo>
void Foam::SolidReaction<ReactionThermo>::solidReactionStrRight
(
    OStringStream& reaction
) const
{

    for (label i = 0; i < grhs().size(); ++i)
    {
        if (i > 0)
        {
            reaction << " + ";
        }
        if (mag(grhs()[i].stoichCoeff - 1) > small)
        {
            reaction << grhs()[i].stoichCoeff;
        }
        reaction << gasSpecies()[grhs()[i].index];
        if (mag(grhs()[i].exponent - grhs()[i].stoichCoeff) > small)
        {
            reaction << "^" << grhs()[i].exponent;
        }
    }
}

// ************************************************************************* //
