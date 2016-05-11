/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "timeActivatedFileUpdate.H"
#include "objectRegistry.H"
#include "Time.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(timeActivatedFileUpdate, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::timeActivatedFileUpdate::updateFile()
{
    label i = lastIndex_;
    while
    (
        i < timeVsFile_.size()-1
     && timeVsFile_[i+1].first() < obr_.time().value()
    )
    {
        i++;
    }

    if (i > lastIndex_)
    {
        Info<< nl << type() << ": copying file" << nl << timeVsFile_[i].second()
            << nl << "to:" << nl << fileToUpdate_ << nl << endl;

        cp(timeVsFile_[i].second(), fileToUpdate_);
        lastIndex_ = i;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::timeActivatedFileUpdate::timeActivatedFileUpdate
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    fileToUpdate_(dict.lookup("fileToUpdate")),
    timeVsFile_(),
    lastIndex_(-1)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::timeActivatedFileUpdate::~timeActivatedFileUpdate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::timeActivatedFileUpdate::read
(
    const dictionary& dict
)
{
    dict.lookup("fileToUpdate") >> fileToUpdate_;
    dict.lookup("timeVsFile") >> timeVsFile_;

    lastIndex_ = -1;
    fileToUpdate_.expand();

    Info<< type() << ": time vs file list:" << nl;
    forAll(timeVsFile_, i)
    {
        timeVsFile_[i].second() = timeVsFile_[i].second().expand();
        if (!isFile(timeVsFile_[i].second()))
        {
            FatalErrorInFunction
                << "File: " << timeVsFile_[i].second() << " not found"
                << nl << exit(FatalError);
        }

        Info<< "    " << timeVsFile_[i].first() << tab
            << timeVsFile_[i].second() << endl;
    }
    Info<< endl;

    updateFile();
}


void Foam::functionObjects::timeActivatedFileUpdate::execute()
{
    updateFile();
}


void Foam::functionObjects::timeActivatedFileUpdate::end()
{
    execute();
}


void Foam::functionObjects::timeActivatedFileUpdate::timeSet()
{}


void Foam::functionObjects::timeActivatedFileUpdate::write()
{}


// ************************************************************************* //
