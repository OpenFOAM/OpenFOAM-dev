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

#include "LagrangianSchemes.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LagrangianSchemes, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LagrangianSchemes::clear()
{
    ddtSchemes_.clear();
    defaultDdtScheme_.clear();
    SpSchemes_.clear();
    defaultSpScheme_.clear();
    averagingSchemes_.clear();
    defaultAveragingScheme_.clear();
    interpolationSchemes_.clear();
    defaultInterpolationScheme_.clear();
    accumulationSchemes_.clear();
    defaultAccumulationScheme_.clear();
}


void Foam::LagrangianSchemes::read
(
    const dictionary& dict,
    const word& type,
    dictionary& typeSchemes,
    ITstream& defaultTypeScheme
)
{
    typeSchemes = dict.subDict(type + "Schemes");

    if
    (
        typeSchemes.found("default")
     && word(typeSchemes.lookup("default")) != "none"
    )
    {
        defaultTypeScheme = typeSchemes.lookup("default");
    }
}


void Foam::LagrangianSchemes::read(const dictionary& dict)
{
    read(dict, "ddt", ddtSchemes_, defaultDdtScheme_);
    read(dict, "Sp", SpSchemes_, defaultSpScheme_);
    read
    (
        dict,
        "averaging",
        averagingSchemes_,
        defaultAveragingScheme_
    );
    read
    (
        dict,
        "interpolation",
        interpolationSchemes_,
        defaultInterpolationScheme_
    );
    read
    (
        dict,
        "accumulation",
        accumulationSchemes_,
        defaultAccumulationScheme_
    );
}


Foam::ITstream& Foam::LagrangianSchemes::lookup
(
    const word& name,
    const dictionary& typeSchemes,
    const ITstream& defaultTypeScheme
)
{
    if (debug)
    {
        Info<< "Lookup scheme for " << name << " in dictionary "
            << typeSchemes.name().caseName() << endl;
    }

    // Find all the indices of the sub-strings that might need to be removed
    DynamicList<Pair<word::size_type>> nameSubIndices;
    {
        bool isSub = false;

        for
        (
            word::size_type nameChari = 0;
            nameChari < name.size();
            ++ nameChari
        )
        {
            if (isSub)
            {
                nameSubIndices.last().second() ++;
            }

            if (name[nameChari] == '.' || name[nameChari] == ':')
            {
                nameSubIndices.append({nameChari, 0});
                isSub = true;
            }
            else if (!isalnum(name[nameChari]))
            {
                isSub = false;
            }
        }
    }

    // Try looking up the name with all combinations of sub-strings removed.
    // Preferentially take the most specific entries with the smallest number
    // of removals. Fail if we get multiple successful lookups with the same
    // number of removed sub-strings.
    for (label nRemoves = 0; nRemoves < nameSubIndices.size() + 1; ++ nRemoves)
    {
        const entry* schemePtr = nullptr;

        for
        (
            label removeMask = 0;
            removeMask < (1 << nameSubIndices.size());
            ++ removeMask
        )
        {
            // Only proceed if the mask has the desired number of removes
            label thisNRemoves = 0;
            for (label r = removeMask; r > 0; r >>= 1) thisNRemoves += r % 2;
            if (thisNRemoves != nRemoves) continue;

            // Build the filtered name by appending parts of the name,
            // including or excluding the identified sub-strings as appropriate
            word filteredName;
            word::size_type nameChari = 0;
            forAll(nameSubIndices, nameSubi)
            {
                filteredName.append
                (
                    name
                    (
                        nameChari,
                        nameSubIndices[nameSubi].first() - nameChari
                    )
                );

                if (!((1 << nameSubi) & removeMask))
                {
                    filteredName.append
                    (
                        name
                        (
                            nameSubIndices[nameSubi].first(),
                            nameSubIndices[nameSubi].second()
                        )
                    );
                }

                nameChari =
                    nameSubIndices[nameSubi].first()
                  + nameSubIndices[nameSubi].second();
            }
            filteredName.append(name(nameChari, name.size() - nameChari));

            // Look up the scheme
            const entry* thisSchemePtr =
                typeSchemes.lookupEntryPtr(filteredName, false, true);

            // Fail if a scheme was found and if we already have a scheme with
            // this precedence
            if (schemePtr && thisSchemePtr)
            {
                FatalIOErrorInFunction(typeSchemes)
                    << "keyword " << name << " ambiguously matches multiple "
                    << "schemes in dictionary " << typeSchemes.name()
                    << exit(FatalIOError);
            }

            schemePtr = thisSchemePtr;
        }

        if (schemePtr)
        {
            return schemePtr->stream();
        }
    }

    // An entry was not found. Generate a lookup error if there is no default.
    if (defaultTypeScheme.empty())
    {
        return typeSchemes.lookup(name);
    }

    // Return the default scheme
    const_cast<ITstream&>(defaultTypeScheme).rewind();
    return const_cast<ITstream&>(defaultTypeScheme);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LagrangianSchemes::LagrangianSchemes(const objectRegistry& db)
:
    IOdictionary
    (
        IOobject
        (
            "LagrangianSchemes",
            db.time().system(),
            db,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    ddtSchemes_(ITstream(objectPath() + ".ddtSchemes", tokenList())()),
    defaultDdtScheme_(ddtSchemes_.name() + ".default", tokenList()),
    SpSchemes_(ITstream(objectPath() + ".SpSchemes", tokenList())()),
    defaultSpScheme_(SpSchemes_.name() + ".default", tokenList()),
    averagingSchemes_
    (
        ITstream
        (
            objectPath() + ".averagingSchemes",
            tokenList()
        )()
    ),
    defaultAveragingScheme_
    (
        averagingSchemes_.name() + ".default",
        tokenList()
    ),
    interpolationSchemes_
    (
        ITstream
        (
            objectPath() + ".interpolationSchemes",
            tokenList()
        )()
    ),
    defaultInterpolationScheme_
    (
        interpolationSchemes_.name() + ".default",
        tokenList()
    ),
    accumulationSchemes_
    (
        ITstream
        (
            objectPath() + ".accumulationSchemes",
            tokenList()
        )()
    ),
    defaultAccumulationScheme_
    (
        accumulationSchemes_.name() + ".default",
        tokenList()
    )
{
    read(schemesDict());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LagrangianSchemes::~LagrangianSchemes()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::LagrangianSchemes::read()
{
    if (regIOobject::read())
    {
        clear();

        read(schemesDict());

        return true;
    }
    else
    {
        return false;
    }
}


const Foam::dictionary& Foam::LagrangianSchemes::schemesDict() const
{
    if (found("select"))
    {
        return subDict(word(IOdictionary::lookup("select")));
    }
    else
    {
        return *this;
    }
}


Foam::ITstream& Foam::LagrangianSchemes::ddt(const word& name) const
{
    return lookup(name, ddtSchemes_, defaultDdtScheme_);
}


Foam::ITstream& Foam::LagrangianSchemes::Sp(const word& name) const
{
    return lookup(name, SpSchemes_, defaultSpScheme_);
}


Foam::ITstream& Foam::LagrangianSchemes::averaging(const word& name) const
{
    return lookup(name, averagingSchemes_, defaultAveragingScheme_);
}


Foam::ITstream& Foam::LagrangianSchemes::interpolation(const word& name) const
{
    return lookup(name, interpolationSchemes_, defaultInterpolationScheme_);
}


Foam::ITstream& Foam::LagrangianSchemes::accumulation(const word& name) const
{
    return lookup(name, accumulationSchemes_, defaultAccumulationScheme_);
}


// ************************************************************************* //
