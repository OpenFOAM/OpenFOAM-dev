/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "zoneGenerator.H"
#include "dlLibraryTable.H"
#include "lookup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zoneGenerator, 0);
    defineRunTimeSelectionTable(zoneGenerator, dictionary);
}

template<>
const char* Foam::NamedEnum<Foam::zoneGenerator::zoneTypes, 3>::names[] =
{
    "point",
    "cell",
    "face"
};

const Foam::NamedEnum<Foam::zoneGenerator::zoneTypes, 3>
    Foam::zoneGenerator::zoneTypesNames;


template<>
const char* Foam::NamedEnum<Foam::zoneGenerator::zoneTypesAll, 4>::names[] =
{
    "point",
    "cell",
    "face",
    "all"
};

const Foam::NamedEnum<Foam::zoneGenerator::zoneTypesAll, 4>
    Foam::zoneGenerator::zoneTypesAllNames;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::zoneGenerator::indices(const boolList& selected)
{
    label nSelected = 0;
    forAll(selected, i)
    {
        if (selected[i])
        {
            nSelected++;
        }
    }

    labelList selectedIndices(nSelected);

    label ui = 0;
    forAll(selected, i)
    {
        if (selected[i])
        {
            selectedIndices[ui++] = i;
        }
    }

    return selectedIndices;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerator::zoneGenerator
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    dict_(dict),
    zoneName_(dict.lookupOrDefault("name", name)),
    mesh_(mesh),
    moveUpdate_(dict.lookupOrDefault("moveUpdate", false))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerator::~zoneGenerator()
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::zoneGenerator>
Foam::zoneGenerator::New
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
{
    const word type(dict.lookup("type"));

    if (debug)
    {
        InfoInFunction
            << "Constructing " << typeName
            << " " << name << " of type " << type << endl;
    }

    if
    (
        !dictionaryConstructorTablePtr_
     || dictionaryConstructorTablePtr_->find(type)
        == dictionaryConstructorTablePtr_->end()
    )
    {
        if
        (
           !libs.open
            (
                dict,
                "libs",
                dictionaryConstructorTablePtr_
            )
        )
        {
            libs.open("lib" + type.remove(':') + ".so", false);
        }

        if (!dictionaryConstructorTablePtr_)
        {
            FatalErrorInFunction
                << "Unknown " << typeName << " type "
                << type << nl << nl
                << "Table of " << typeName << " is empty"
                << exit(FatalError);
        }
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Unknown " << typeName << " type "
            << type << nl << nl
            << "Valid " << typeName << " types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<zoneGenerator>(cstrIter()(name, mesh, dict));
}


Foam::autoPtr<Foam::zoneGenerator>
Foam::zoneGenerator::New
(
    const word& name,
    const zoneTypes& zoneType,
    const polyMesh& mesh,
    const dictionary& dict
)
{
    // Copy the dictionary and add the zoneType entry
    dictionary zoneDict(dict);
    zoneDict.add("zoneType", zoneTypesNames[zoneType]);

    return New(name, mesh, zoneDict);
}


Foam::autoPtr<Foam::zoneGenerator>
Foam::zoneGenerator::New
(
    const polyMesh& mesh,
    const dictionary& dict
)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing " << typeName << endl;
    }

    forAllConstIter(dictionary, dict, iter)
    {
        const word& name = iter().keyword();

        if (iter().isDict())
        {
            if (iter().dict().found("type"))
            {
                return zoneGenerator::New(name, mesh, iter().dict());
            }
            else
            {
                // If an empty keyword is present assume it is a zone name
                // and add a zone lookup
                return autoPtr<zoneGenerator>
                (
                    new zoneGenerators::lookup(name, mesh, iter().dict())
                );
            }
        }
        else if (!iter().stream().size())
        {
            // If an empty keyword is present assume it is a zone name
            // and add a zone lookup
            dictionary zoneDict(name, dict);
            zoneDict.add
            (
                primitiveEntry
                (
                    "type",
                    zoneGenerators::lookup::typeName,
                    iter().startLineNumber(),
                    iter().startLineNumber()
                )
            );

            return autoPtr<zoneGenerator>
            (
                new zoneGenerators::lookup(name, mesh, zoneDict)
            );
        }
    }

    FatalIOErrorInFunction(dict)
        << "Cannot find valid " << typeName
        << " entry in dictionary " << dict << exit(FatalIOError);

    return autoPtr<zoneGenerator>(nullptr);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerator::movePoints() const
{
    if (moveUpdate_)
    {
        return generate();
    }
    else
    {
        return zoneSet();
    }
}


// ************************************************************************* //
