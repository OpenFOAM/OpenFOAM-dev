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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
Foam::PtrList<ThermoType>
Foam::multicomponentMixture<ThermoType>::readSpeciesData
(
    const dictionary& thermoDict
) const
{
    PtrList<ThermoType> specieThermos(species_.size());

    forAll(species_, i)
    {
        specieThermos.set
        (
            i,
            new ThermoType(species_[i], thermoDict.subDict(species_[i]))
        );
    }

    return specieThermos;
}


template<class ThermoType>
Foam::List<Foam::List<Foam::specieElement>>
Foam::multicomponentMixture<ThermoType>::readSpeciesComposition
(
    const dictionary& thermoDict
) const
{
    List<List<specieElement>> specieCompositions(species_.size());

    // Loop through all species in thermoDict to retrieve
    // the species composition
    forAll(species_, i)
    {
        if (thermoDict.subDict(species_[i]).isDict("elements"))
        {
            const dictionary& elements =
                thermoDict.subDict(species_[i]).subDict("elements");

            const wordList elementsNames(elements.toc());

            specieCompositions[i].resize(elementsNames.size());

            forAll(elementsNames, eni)
            {
                specieCompositions[i][eni].name() = elementsNames[eni];
                specieCompositions[i][eni].nAtoms() =
                    elements.lookupOrDefault(elementsNames[eni], 0);
            }
        }
    }

    return specieCompositions;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::multicomponentMixture<ThermoType>::multicomponentMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicSpecieMixture
    (
        thermoDict,
        thermoDict.lookup("species"),
        mesh,
        phaseName
    ),
    specieThermos_(readSpeciesData(thermoDict)),
    specieCompositions_(readSpeciesComposition(thermoDict))
{
    correctMassFractions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::multicomponentMixture<ThermoType>::read
(
    const dictionary& thermoDict
)
{
    specieThermos_ = readSpeciesData(thermoDict);
    specieCompositions_ = readSpeciesComposition(thermoDict);
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
        // Spit an error associated with the lookup of this specie's elements
        refCast<const dictionary>(*this)
            .subDict(species_[speciei])
            .subDict("elements");
    }

    return specieCompositions_[speciei];
}


// ************************************************************************* //
