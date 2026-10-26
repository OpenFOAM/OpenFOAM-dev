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

#include "reactionList.H"
#include "HashPtrTable.H"
#include "delimitDictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactionList::reactionList()
:
    Tlow_(0),
    Thigh_(great)
{}


Foam::reactionList::reactionList
(
    const speciesTable& species,
    const specieThermosRef& specieThermos,
    const dictionary& dict
)
:
    Tlow_(dict.lookupOrDefault<scalar>("Tlow", 0)),
    Thigh_(dict.lookupOrDefault<scalar>("Thigh", great))
{
    const dictionary& reactions(dict.subDict("reactions"));

    this->setSize(reactions.size());

    label i = 0;

    forAllConstIter(dictionary, reactions, iter)
    {
        this->set
        (
            i++,
            reaction::New
            (
                species,
                specieThermos,
                reactions.subDict(iter().keyword())
            ).ptr()
        );
    }
}


Foam::reactionList::reactionList
(
    const speciesTable& species,
    const specieThermosRef& specieThermos,
    const objectRegistry& ob,
    const dictionary& dict
)
:
    Tlow_(dict.lookupOrDefault<scalar>("Tlow", 0)),
    Thigh_(dict.lookupOrDefault<scalar>("Thigh", great))
{
    const dictionary& reactions(dict.subDict("reactions"));

    this->setSize(reactions.size());

    label i = 0;

    forAllConstIter(dictionary, reactions, iter)
    {
        this->set
        (
            i++,
            reaction::New
            (
                species,
                specieThermos,
                ob,
                reactions.subDict(iter().keyword())
            ).ptr()
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::reactionList::write(Ostream& os) const
{
    {
        const delimitDictionary delimitReactions(os, "reactions");

        forAll(*this, i)
        {
            const reaction& r = this->operator[](i);
            const delimitDictionary delimitReaction(os, r.name());
            writeEntry(os, "type", r.selectionType());
            r.write(os);
        }
    }

    writeEntry(os, "Tlow", Tlow_);
    writeEntry(os, "Thigh", Thigh_);
}


// ************************************************************************* //
