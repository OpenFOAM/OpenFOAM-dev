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

#include "timeActivatedFileUpdate.H"
#include "Time.H"
#include "OSspecific.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(timeActivatedFileUpdate, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        timeActivatedFileUpdate,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::timeActivatedFileUpdate::updateFile()
{
    label i = lastIndex_;
    while
    (
        i < fileVsTime_.size()-1
     && fileVsTime_[i+1].first() <= time_.userTimeValue()
    )
    {
        i++;
    }

    if (i > lastIndex_)
    {
        Info<< nl << type() << ": copying file" << nl << fileVsTime_[i].second()
            << nl << "to:" << nl << fileToUpdate_ << nl << endl;

        fileName destFile(fileToUpdate_ + Foam::name(pid()));
        cp(fileVsTime_[i].second(), destFile);
        mv(destFile, fileToUpdate_);
        lastIndex_ = i;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::timeActivatedFileUpdate::timeActivatedFileUpdate
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name, runTime),
    fileToUpdate_(dict.lookup("fileToUpdate")),
    fileVsTime_(),
    lastIndex_(-1)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::timeActivatedFileUpdate::~timeActivatedFileUpdate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::timeActivatedFileUpdate::read
(
    const dictionary& dict
)
{
    dict.lookup("fileToUpdate") >> fileToUpdate_;
    dict.lookupBackwardsCompatible({"fileVsTime", "timeVsFile"}) >> fileVsTime_;

    lastIndex_ = -1;
    fileToUpdate_.expand();

    Info<< type() << ": file vs time list:" << nl;
    forAll(fileVsTime_, i)
    {
        fileVsTime_[i].second() = fileVsTime_[i].second().expand();
        if (!isFile(fileVsTime_[i].second()))
        {
            FatalErrorInFunction
                << "File: " << fileVsTime_[i].second() << " not found"
                << nl << exit(FatalError);
        }

        Info<< "    " << fileVsTime_[i].first() << tab
            << fileVsTime_[i].second() << endl;
    }
    Info<< endl;

    updateFile();

    return true;
}


bool Foam::functionObjects::timeActivatedFileUpdate::execute()
{
    updateFile();

    return true;
}


bool Foam::functionObjects::timeActivatedFileUpdate::write()
{
    return true;
}


// ************************************************************************* //
