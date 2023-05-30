/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "multicomponentMixture.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::multicomponentMixture<ThermoType>::multicomponentMixture
(
    const dictionary& dict
)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
Foam::wordList Foam::multicomponentMixture<ThermoType>::specieNames() const
{
    wordList result(specieThermos_.size());

    forAll(specieThermos_, speciei)
    {
        result[speciei] = specieThermos_[speciei].name();
    }

    return result;
}


template<class ThermoType>
void Foam::multicomponentMixture<ThermoType>::read
(
    const dictionary& dict
)
{
    const wordList specieNames(dict.lookup("species"));

    specieThermos_.setSize(specieNames.size());
    specieCompositions_.setSize(specieNames.size());
    specieDictLocations_.setSize(specieNames.size());

    forAll(specieNames, speciei)
    {
        const dictionary& specieDict = dict.subDict(specieNames[speciei]);

        specieThermos_.set
        (
            speciei,
            new ThermoType(specieNames[speciei], specieDict)
        );

        if (specieDict.isDict("elements"))
        {
            const dictionary& specieElementsDict =
                specieDict.subDict("elements");

            const wordList elementsNames(specieElementsDict.toc());

            specieCompositions_[speciei].resize(elementsNames.size());

            forAll(elementsNames, eni)
            {
                specieCompositions_[speciei][eni].name() = elementsNames[eni];
                specieCompositions_[speciei][eni].nAtoms() =
                    specieElementsDict.lookupOrDefault(elementsNames[eni], 0);
            }
        }

        specieDictLocations_[speciei] = IOerrorLocation(specieDict);
    }
}


template<class ThermoType>
const Foam::List<Foam::specieElement>&
Foam::multicomponentMixture<ThermoType>::specieComposition
(
    const label speciei
) const
{
    if (specieCompositions_[speciei].empty())
    {
        FatalIOErrorInFunction(specieDictLocations_[speciei])
            << "Elemental composition not specified for specie "
            << specieThermos_[speciei].name()
            << exit(FatalIOError);
    }

    return specieCompositions_[speciei];
}


// ************************************************************************* //
