/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "abort.H"
#include "dictionary.H"
#include "error.H"
#include "Time.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(abort, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        abort,
        dictionary
    );
}
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::abort::actionType,
    3
>::names[] =
{
    "noWriteNow",
    "writeNow",
    "nextWrite"
};

const Foam::NamedEnum
<
    Foam::functionObjects::abort::actionType,
    3
> Foam::functionObjects::abort::actionTypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::abort::removeFile() const
{
    bool hasAbort = isFile(abortFile_);
    reduce(hasAbort, orOp<bool>());

    if (hasAbort && Pstream::master())
    {
        // Cleanup ABORT file (on master only)
        rm(abortFile_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::abort::abort
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    time_(runTime),
    abortFile_("$FOAM_CASE/" + name),
    action_(nextWrite)
{
    abortFile_.expand();
    read(dict);

    // Remove any old files from previous runs
    removeFile();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::abort::~abort()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::abort::read(const dictionary& dict)
{
    if (dict.found("action"))
    {
        action_ = actionTypeNames_.read(dict.lookup("action"));
    }
    else
    {
        action_ = nextWrite;
    }

    if (dict.readIfPresent("file", abortFile_))
    {
        abortFile_.expand();
    }

    return true;
}


bool Foam::functionObjects::abort::execute()
{
    bool hasAbort = isFile(abortFile_);
    reduce(hasAbort, orOp<bool>());

    if (hasAbort)
    {
        switch (action_)
        {
            case noWriteNow :
            {
                if (time_.stopAt(Time::saNoWriteNow))
                {
                    Info<< "USER REQUESTED ABORT (timeIndex="
                        << time_.timeIndex()
                        << "): stop without writing data"
                        << endl;
                }
                break;
            }

            case writeNow :
            {
                if (time_.stopAt(Time::saWriteNow))
                {
                    Info<< "USER REQUESTED ABORT (timeIndex="
                        << time_.timeIndex()
                        << "): stop+write data"
                        << endl;
                }
                break;
            }

            case nextWrite :
            {
                if (time_.stopAt(Time::saNextWrite))
                {
                    Info<< "USER REQUESTED ABORT (timeIndex="
                        << time_.timeIndex()
                        << "): stop after next data write"
                        << endl;
                }
                break;
            }
        }
    }

    return true;
}


bool Foam::functionObjects::abort::write()
{
    return true;
}


bool Foam::functionObjects::abort::end()
{
    removeFile();
    return true;
}


// ************************************************************************* //
