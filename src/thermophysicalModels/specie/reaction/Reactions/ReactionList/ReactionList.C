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

#include "ReactionList.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::ReactionList<ThermoType>::ReactionList
(
    const speciesTable& species,
    const PtrList<ThermoType>& speciesThermo,
    const dictionary& dict
)
{
    // Set general temperature limits from the dictionary
    Reaction<ThermoType>::TlowDefault =
        dict.lookupOrDefault<scalar>("Tlow", 0);

    Reaction<ThermoType>::ThighDefault =
        dict.lookupOrDefault<scalar>("Thigh", great);

    const dictionary& reactions(dict.subDict("reactions"));

    HashPtrTable<ThermoType> thermoDatabase;
    forAll(speciesThermo, i)
    {
        thermoDatabase.insert
        (
            speciesThermo[i].name(),
            speciesThermo[i].clone().ptr()
        );
    }

    this->setSize(reactions.size());
    label i = 0;

    forAllConstIter(dictionary, reactions, iter)
    {
        this->set
        (
            i++,
            Reaction<ThermoType>::New
            (
                species,
                thermoDatabase,
                reactions.subDict(iter().keyword())
            ).ptr()
        );
    }
}


template<class ThermoType>
Foam::ReactionList<ThermoType>::ReactionList
(
    const speciesTable& species,
    const PtrList<ThermoType>& speciesThermo,
    const objectRegistry& ob,
    const dictionary& dict
)
{
    // Set general temperature limits from the dictionary
    Reaction<ThermoType>::TlowDefault =
        dict.lookupOrDefault<scalar>("Tlow", 0);

    Reaction<ThermoType>::ThighDefault =
        dict.lookupOrDefault<scalar>("Thigh", great);

    const dictionary& reactions(dict.subDict("reactions"));

    HashPtrTable<ThermoType> thermoDatabase;
    forAll(speciesThermo, i)
    {
        thermoDatabase.insert
        (
            speciesThermo[i].name(),
            speciesThermo[i].clone().ptr()
        );
    }

    this->setSize(reactions.size());
    label i = 0;

    forAllConstIter(dictionary, reactions, iter)
    {
        this->set
        (
            i++,
            Reaction<ThermoType>::New
            (
                species,
                thermoDatabase,
                ob,
                reactions.subDict(iter().keyword())
            ).ptr()
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::ReactionList<ThermoType>::write(Ostream& os) const
{
    os  << "reactions" << nl;
    os  << token::BEGIN_BLOCK << incrIndent << nl;

    forAll(*this, i)
    {
        const Reaction<ThermoType>& r = this->operator[](i);

        os  << indent << r.name() << nl
            << indent << token::BEGIN_BLOCK << incrIndent << nl;

        writeEntry(os, "type", r.type());

        r.write(os);

        os  << decrIndent << indent << token::END_BLOCK << nl;
    }

    os << decrIndent << token::END_BLOCK << nl;

    writeEntry(os, "Tlow", Reaction<ThermoType>::TlowDefault);
    writeEntry(os, "Thigh", Reaction<ThermoType>::ThighDefault);

    os << nl;
}


// ************************************************************************* //
