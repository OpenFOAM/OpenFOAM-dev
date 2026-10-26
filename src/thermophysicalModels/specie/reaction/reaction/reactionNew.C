/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2026 OpenFOAM Foundation
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
#include "thermoName.H"
#include "compileTemplate.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::HashTable<Foam::word> Foam::reaction::backwardTypeToType()
{
    static const List<Pair<word>> reactions
    ({
        {"irreversible", "Irreversible"},
        {"reversible", "Reversible"},
        {"nonEquilibriumReversible", "NonEquilibriumReversible"}
    });

    static const List<Pair<word>> rates
    ({
        {"Arrhenius", "Arrhenius"},
        {
            "ArrheniusLindemannChemicallyActivated",
            "ChemicallyActivated<Arrhenius,Lindemann>"
        },
        {
            "ArrheniusSRIChemicallyActivated",
            "ChemicallyActivated<Arrhenius,SRI>"
        },
        {
            "ArrheniusTroeChemicallyActivated",
            "ChemicallyActivated<Arrhenius,Troe>"
        },
        {"ArrheniusLindemannFallOff", "FallOff<Arrhenius,Lindemann>"},
        {"ArrheniusSRIFallOff", "FallOff<Arrhenius,SRI>"},
        {"ArrheniusTroeFallOff", "FallOff<Arrhenius,Troe>"},
        {
            "FluxLimitedLangmuirHinshelwood",
            "fluxLimitedLangmuirHinshelwood"
        },
        {"Janev", "Janev"},
        {"LandauTeller", "LandauTeller"},
        {"LangmuirHinshelwood", "LangmuirHinshelwood"},
        {"MichaelisMenten", "MichaelisMenten"},
        {"PhaseSurfaceArrhenius", "PhaseSurface<Arrhenius>"},
        {"PLog", "pLog"},
        {"PowerSeries", "powerSeries"},
        {"SurfaceArrhenius", "Surface<Arrhenius>"},
        {"ThirdBodyArrhenius", "ThirdBody<Arrhenius>"}
    });

    HashTable<word> result;

    forAll(reactions, reactioni)
    {
        forAll(rates, ratei)
        {
            result.insert
            (
                reactions[reactioni].first() + rates[ratei].first(),
                reactions[reactioni].second()
              + '<' + rates[ratei].second() + '>'
            );
        }
    }

    return result;
}


Foam::word Foam::reaction::backwardTypeToType(const word& type)
{
    static HashTable<word> backwardTypeToType
    (
        Foam::reaction::backwardTypeToType()
    );

    return
        backwardTypeToType.found(type)
      ? backwardTypeToType[type]
      : type;
}


bool Foam::reaction::compile
(
    const specieThermosRef& specieThermos,
    ITstream& typeIts,
    const word& type,
    const word& typeThermo,
    const word& table
)
{
    if (!dynamicCode::allowSystemOperations) return false;

    const wordList typeSplit(Foam::templateName::split(type));

    List<Pair<word>> substitutions =
        thermoName::components(specieThermos.thermoName());
    substitutions.append({"table", table});
    substitutions.append({"reaction", typeSplit[0]});
    word templateName = word::null;
    switch (typeSplit.size())
    {
        case 2:
            substitutions.append({"rate", typeSplit[1]});
            templateName = "rateReaction";
            break;
        case 3:
            substitutions.append({"modifier", typeSplit[1]});
            substitutions.append({"rate", typeSplit[2]});
            templateName = "modifiedRateReaction";
            break;
        case 4:
            substitutions.append({"modifier", typeSplit[1]});
            substitutions.append({"rate", typeSplit[2]});
            substitutions.append({"function", typeSplit[3]});
            templateName = "functionModifiedRateReaction";
            break;
        default:
            break;
    }

    if
    (
        templateName.empty()
     || dynamicCode::resolveTemplate(templateName).empty()
    )
    {
        return false;
    }

    compileTemplate chemistryModel
    (
        templateName,
        typeThermo,
        substitutions
    );

    if
    (
        !dictionaryConstructorTablePtr_->found(typeThermo)
     && !objectRegistryConstructorTablePtr_->found(typeThermo)
    )
    {
        FatalIOErrorInFunction(typeIts)
            << "Compilation and linkage of " << reaction::typeName
            << " type " << type << " for thermo of type "
            << specieThermos.thermoName() << " failed." << nl << nl
            << exit(FatalIOError);
    }

    return true;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::reaction> Foam::reaction::New
(
    const speciesTable& species,
    const specieThermosRef& specieThermos,
    const dictionary& dict
)
{
    ITstream& typeIts = dict.lookup("type");
    const word type = backwardTypeToType(typeIts);
    const word typeThermo =
        type[type.size() - 1] == '>'
      ? type(type.size() - 1) + ',' + specieThermos.thermoName() + '>'
      : type + '<' + specieThermos.thermoName() + '>';

    // Try the dictionary constructor
    {
        typename dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(type);
        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            cstrIter = dictionaryConstructorTablePtr_->find(typeThermo);
        }
        if (cstrIter != dictionaryConstructorTablePtr_->end())
        {
            return cstrIter()(species, specieThermos, dict);
        }
    }

    // Try dynamic compilation
    if (compile(specieThermos, typeIts, type, typeThermo, "dictionary"))
    {
        typename dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(typeThermo);
        if (cstrIter != dictionaryConstructorTablePtr_->end())
        {
            return cstrIter()(species, specieThermos, dict);
        }
    }

    FatalIOErrorInFunction(typeIts)
        << "Unknown reaction type " << type << " (or " << typeThermo << ')'
        << nl << nl << "Valid reaction types are :" << nl
        << dictionaryConstructorTablePtr_->sortedToc()
        << exit(FatalIOError);

    return autoPtr<reaction>(nullptr);
}


Foam::autoPtr<Foam::reaction> Foam::reaction::New
(
    const speciesTable& species,
    const specieThermosRef& specieThermos,
    const objectRegistry& ob,
    const dictionary& dict
)
{
    // If the objectRegistry constructor table is empty
    // use the dictionary constructor table only
    if (!objectRegistryConstructorTablePtr_)
    {
        return New(species, specieThermos, dict);
    }

    ITstream& typeIts = dict.lookup("type");
    const word type = backwardTypeToType(typeIts);
    const word typeThermo =
        type[type.size() - 1] == '>'
      ? type(type.size() - 1) + ',' + specieThermos.thermoName() + '>'
      : type + '<' + specieThermos.thermoName() + '>';

    // Try the objectRegistry constructor
    {
        typename objectRegistryConstructorTable::iterator cstrIter =
            objectRegistryConstructorTablePtr_->find(type);
        if (cstrIter == objectRegistryConstructorTablePtr_->end())
        {
            cstrIter = objectRegistryConstructorTablePtr_->find(typeThermo);
        }
        if (cstrIter != objectRegistryConstructorTablePtr_->end())
        {
            return cstrIter()(species, specieThermos, ob, dict);
        }
    }

    // Try the dictionary constructor
    {
        typename dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(type);
        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            cstrIter = dictionaryConstructorTablePtr_->find(typeThermo);
        }
        if (cstrIter != dictionaryConstructorTablePtr_->end())
        {
            return cstrIter()(species, specieThermos, dict);
        }
    }

    // Try dynamic compilation
    if (compile(specieThermos, typeIts, type, typeThermo, "objectRegistry"))
    {
        typename objectRegistryConstructorTable::iterator cstrIter =
            objectRegistryConstructorTablePtr_->find(typeThermo);
        if (cstrIter != objectRegistryConstructorTablePtr_->end())
        {
            return cstrIter()(species, specieThermos, ob, dict);
        }
    }

    wordHashSet toc;
    toc.insert(dictionaryConstructorTablePtr_->toc());
    toc.insert(objectRegistryConstructorTablePtr_->toc());
    FatalIOErrorInFunction(typeIts)
        << "Unknown reaction type " << type << " (or " << typeThermo << ')'
        << nl << nl << "Valid reaction types are :" << nl << toc.sortedToc()
        << exit(FatalIOError);

    return autoPtr<reaction>(nullptr);
}


// ************************************************************************* //
