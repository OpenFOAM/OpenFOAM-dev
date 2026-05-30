/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "coded_zoneGenerator.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(coded, 0);

        addToRunTimeSelectionTable
        (
            zoneGenerator,
            coded,
            dictionary
        );
    }
}


const Foam::wordList Foam::zoneGenerators::coded::codeKeys
{
    "code",
    "codeInclude"
};

const Foam::wordList Foam::zoneGenerators::coded::codeDictVars
{
    word::null,
    word::null
};

const Foam::word Foam::zoneGenerators::coded::codeOptions
(
    "codedZoneGeneratorOptions"
);

const Foam::wordList Foam::zoneGenerators::coded::compileFiles
{
    "codedZoneGeneratorTemplate.C"
};

const Foam::wordList Foam::zoneGenerators::coded::copyFiles
{
    "codedZoneGeneratorTemplate.H"
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::coded::coded
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    zoneGenerator(name, mesh, dict),
    codedBase
    (
        name,
        dict,
        codeKeys,
        codeDictVars,
        codeOptions,
        compileFiles,
        copyFiles
    )
{
    // Set verbose if debugging
    varSubstitutions().set("verbose", Foam::name(bool(debug)));

    updateLibrary(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::coded::~coded()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::coded::generate() const
{
    if (!redirectZoneGeneratorPtr_.valid())
    {
        dictionary redirectDict;
        redirectDict.set("type", codeName());

        redirectZoneGeneratorPtr_ = zoneGenerator::New
        (
            codeName(),
            mesh_,
            redirectDict
        );
    }

    return redirectZoneGeneratorPtr_->generate();
}


// ************************************************************************* //
