/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "stopAtFile.H"
#include "dictionary.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(stopAtFile, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        stopAtFile,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::stopAtFile::removeFile() const
{
    bool fileExists = isFile(stopAtFileFile_);
    reduce(fileExists, orOp<bool>());

    if (fileExists && Pstream::master())
    {
        // Cleanup ABORT file (on master only)
        rm(stopAtFileFile_);
    }
}


bool Foam::functionObjects::stopAtFile::condition() const
{
    bool fileExists = isFile(stopAtFileFile_);
    reduce(fileExists, orOp<bool>());
    return fileExists;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::stopAtFile::stopAtFile
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    stopAt(name, runTime, dict),
    stopAtFileFile_("$FOAM_CASE/" + name)
{
    stopAtFileFile_.expand();
    read(dict);

    // Remove any old files from previous runs
    removeFile();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::stopAtFile::~stopAtFile()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::stopAtFile::read(const dictionary& dict)
{
    if (dict.readIfPresent("file", stopAtFileFile_))
    {
        stopAtFileFile_.expand();
    }

    return true;
}


bool Foam::functionObjects::stopAtFile::end()
{
    removeFile();
    return true;
}


// ************************************************************************* //
