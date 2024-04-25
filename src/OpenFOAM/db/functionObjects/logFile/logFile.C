/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "logFile.H"
#include "Time.H"
#include "OSspecific.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::fileName Foam::functionObjects::logFile::filePathName() const
{
    const word timeName = fileObr_.time().name();
    const fileName outputDir(baseFileDir()/prefix_/timeName);
    mkDir(outputDir);
    return outputDir/(name_ + ".dat");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::logFile::logFile
(
    const objectRegistry& obr,
    const word& prefix,
    const word& name
)
:
    writeFile(obr, prefix),
    name_(name),
    file_(filePathName())
{
    initStream(file_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::logFile::~logFile()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::logFile::writeTimeColumnHeaders
(
    const wordList& titles
)
{
    writeCommented(file_, "Time");
    forAll(titles, i)
    {
        writeTabbed(file_, titles[i]);
    }
    file_ << endl;
}


// ************************************************************************* //
