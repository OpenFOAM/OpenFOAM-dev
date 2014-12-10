/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "ReactionList.H"
#include "IFstream.H"
#include "SLPtrList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::ReactionList<ThermoType>::ReactionList
(
    const speciesTable& species,
    const HashPtrTable<ThermoType>& thermoDb
)
:
    SLPtrList<Reaction<ThermoType> >(),
    species_(species),
    thermoDb_(thermoDb),
    dict_(dictionary::null)
{}


template<class ThermoType>
Foam::ReactionList<ThermoType>::ReactionList
(
    const speciesTable& species,
    const HashPtrTable<ThermoType>& thermoDb,
    const dictionary& dict
)
:
    SLPtrList<Reaction<ThermoType> >(),
    species_(species),
    thermoDb_(thermoDb),
    dict_(dict)
{
    readReactionDict();
}


template<class ThermoType>
Foam::ReactionList<ThermoType>::ReactionList
(
    const speciesTable& species,
    const HashPtrTable<ThermoType>& thermoDb,
    const fileName& fName
)
:
    SLPtrList<Reaction<ThermoType> >
    (
        dictionary(IFstream(fName)()).lookup("reactions"),
        Reaction<ThermoType>::iNew(species, thermoDb)
    ),
    species_(species),
    thermoDb_(thermoDb),
    dict_(dictionary::null)
{}


template<class ThermoType>
Foam::ReactionList<ThermoType>::ReactionList(const ReactionList& reactions)
:
    SLPtrList<Reaction<ThermoType> >(reactions),
    species_(reactions.species_),
    thermoDb_(reactions.thermoDb_),
    dict_(reactions.dict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::ReactionList<ThermoType>::~ReactionList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
bool Foam::ReactionList<ThermoType>::readReactionDict()
{
    const dictionary& reactions(dict_.subDict("reactions"));

    forAllConstIter(dictionary, reactions, iter)
    {
        const word reactionName = iter().keyword();

        this->append
        (
            Reaction<ThermoType>::New
            (
                species_,
                thermoDb_,
                reactions.subDict(reactionName)
            ).ptr()
        );
    }

    return true;
}


template<class ThermoType>
void Foam::ReactionList<ThermoType>::write(Ostream& os) const
{
    os  << "reactions" << nl;
    os  << token::BEGIN_BLOCK << incrIndent << nl;

    forAllConstIter(typename SLPtrList<Reaction<ThermoType> >, *this, iter)
    {
        const Reaction<ThermoType>& r = iter();
        os  << indent << r.name() << nl
            << indent << token::BEGIN_BLOCK << incrIndent << nl;
        os.writeKeyword("type") << r.type() << token::END_STATEMENT << nl;
        r.write(os);
        os  << decrIndent << indent << token::END_BLOCK << nl;
    }

    os << decrIndent << token::END_BLOCK << nl;
}


// ************************************************************************* //
