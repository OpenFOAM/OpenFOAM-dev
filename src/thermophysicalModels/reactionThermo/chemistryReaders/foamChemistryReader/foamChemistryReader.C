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

#include "foamChemistryReader.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
Foam::speciesTable& Foam::foamChemistryReader<ThermoType>::setSpecies
(
    const dictionary& dict,
    speciesTable& species
)
{
    wordList s(dict.lookup("species"));
    species.transfer(s);
    return species;
}


template<class ThermoType>
void Foam::foamChemistryReader<ThermoType>::readSpeciesComposition()
{
    if (!chemDict_.found("elements"))
    {
        Info<< "    elements not defined in " << chemDict_.name() << endl;
        return;
    }

    wordList e(chemDict_.lookup("elements"));
    label currentElementIndex(0);

    DynamicList<word> elementNames_;
    HashTable<label> elementIndices_;

    forAll(e, ei)
    {
        if (!elementIndices_.found(e[ei]))
        {
            elementIndices_.insert(e[ei], currentElementIndex++);
            elementNames_.append(e[ei]);
        }
        else
        {
            IOWarningInFunction(chemDict_)
                << "element " << e[ei] << " already in table." << endl;
        }
    }

    // Loop through all species in thermoDict to retrieve
    // the species composition
    forAll(speciesTable_, si)
    {
        if (thermoDict_.subDict(speciesTable_[si]).isDict("elements"))
        {
            dictionary currentElements
            (
                thermoDict_.subDict(speciesTable_[si]).subDict("elements")
            );

            wordList currentElementsName(currentElements.toc());
            List<specieElement> currentComposition(currentElementsName.size());

            forAll(currentElementsName, eni)
            {
                currentComposition[eni].name() = currentElementsName[eni];

                currentComposition[eni].nAtoms() =
                    currentElements.lookupOrDefault
                    (
                        currentElementsName[eni],
                        0
                    );
            }

            // Add current specie composition to the hash table
            speciesCompositionTable::iterator specieCompositionIter
            (
                speciesComposition_.find(speciesTable_[si])
            );

            if (specieCompositionIter != speciesComposition_.end())
            {
                speciesComposition_.erase(specieCompositionIter);
            }

            speciesComposition_.insert(speciesTable_[si], currentComposition);
        }
        else
        {
            FatalIOErrorInFunction(thermoDict_)
                << "Specie " << speciesTable_[si]
                << " does not contain element description."
                << exit(FatalIOError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::foamChemistryReader<ThermoType>::foamChemistryReader
(
    const fileName& reactionsFileName,
    speciesTable& species,
    const fileName& thermoFileName
)
:
    chemistryReader<ThermoType>(),
    chemDict_
    (
        IFstream
        (
            fileName(reactionsFileName).expand()
        )()
    ),
    thermoDict_
    (
        IFstream
        (
            fileName(thermoFileName).expand()
        )()
    ),
    speciesTable_(setSpecies(chemDict_, species)),
    speciesThermo_(thermoDict_),
    reactions_(speciesTable_, speciesThermo_, chemDict_)
{
    readSpeciesComposition();
}


template<class ThermoType>
Foam::foamChemistryReader<ThermoType>::foamChemistryReader
(
    const dictionary& thermoDict,
    speciesTable& species
)
:
    chemistryReader<ThermoType>(),
    chemDict_
    (
        IFstream
        (
            fileName(thermoDict.lookup("foamChemistryFile")).expand()
        )()
    ),
    thermoDict_
    (
        IFstream
        (
            fileName(thermoDict.lookup("foamChemistryThermoFile")).expand()
        )()
    ),
    speciesTable_(setSpecies(chemDict_, species)),
    speciesThermo_(thermoDict_),
    reactions_(speciesTable_, speciesThermo_, chemDict_)
{
    readSpeciesComposition();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
