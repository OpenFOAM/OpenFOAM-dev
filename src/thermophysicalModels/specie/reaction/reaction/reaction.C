/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#include "reaction.H"
#include "dictionary.H"


// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

Foam::label Foam::reaction::nUnNamedReactions(0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::reaction::getNewReactionID()
{
    return nUnNamedReactions++;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reaction::reaction
(
    const speciesTable& species,
    const List<specieCoeffs>& lhs,
    const List<specieCoeffs>& rhs
)
:
    name_("un-named-reaction-" + Foam::name(getNewReactionID())),
    species_(species),
    lhs_(lhs),
    rhs_(rhs)
{}


Foam::reaction::reaction
(
    const reaction& r,
    const speciesTable& species
)
:
    name_(r.name() + "Copy"),
    species_(species),
    lhs_(r.lhs_),
    rhs_(r.rhs_)
{}


Foam::reaction::reaction
(
    const speciesTable& species,
    const dictionary& dict
)
:
    name_(dict.dictName()),
    species_(species)
{
    specieCoeffs::setLRhs
    (
        IStringStream(dict.lookup("reaction"))(),
        species_,
        lhs_,
        rhs_
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::reaction::write(Ostream& os) const
{
    OStringStream reaction;
    writeEntry
    (
        os,
        "reaction",
        specieCoeffs::reactionStr(reaction, species_, lhs_, rhs_)
    );
}


// ************************************************************************* //
