/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "abortCalculation.H"
#include "dictionary.H"
#include "error.H"
#include "Time.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(abortCalculation, 0);
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::abortCalculation::actionType,
        3
    >::names[] =
    {
        "noWriteNow",
        "writeNow",
        "nextWrite"
    };
}


const Foam::NamedEnum<Foam::abortCalculation::actionType, 3>
    Foam::abortCalculation::actionTypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::abortCalculation::removeFile() const
{
    bool hasAbort = isFile(abortFile_);
    reduce(hasAbort, orOp<bool>());

    if (hasAbort && Pstream::master())
    {
        // cleanup ABORT file (on master only)
        rm(abortFile_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::abortCalculation::abortCalculation
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    abortFile_("$FOAM_CASE/" + name),
    action_(nextWrite)
{
    abortFile_.expand();
    read(dict);

    // remove any old files from previous runs
    removeFile();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::abortCalculation::~abortCalculation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::abortCalculation::read(const dictionary& dict)
{
    if (dict.found("action"))
    {
        action_ = actionTypeNames_.read(dict.lookup("action"));
    }
    else
    {
        action_ = nextWrite;
    }

    if (dict.readIfPresent("fileName", abortFile_))
    {
        abortFile_.expand();
    }
}


void Foam::abortCalculation::execute()
{
    bool hasAbort = isFile(abortFile_);
    reduce(hasAbort, orOp<bool>());

    if (hasAbort)
    {
        switch (action_)
        {
            case noWriteNow :
            {
                if (obr_.time().stopAt(Time::saNoWriteNow))
                {
                    Info<< "USER REQUESTED ABORT (timeIndex="
                        << obr_.time().timeIndex()
                        << "): stop without writing data"
                        << endl;
                }
                break;
            }

            case writeNow :
            {
                if (obr_.time().stopAt(Time::saWriteNow))
                {
                    Info<< "USER REQUESTED ABORT (timeIndex="
                        << obr_.time().timeIndex()
                        << "): stop+write data"
                        << endl;
                }
                break;
            }

            case nextWrite :
            {
                if (obr_.time().stopAt(Time::saNextWrite))
                {
                    Info<< "USER REQUESTED ABORT (timeIndex="
                        << obr_.time().timeIndex()
                        << "): stop after next data write"
                        << endl;
                }
                break;
            }
        }
    }
}


void Foam::abortCalculation::end()
{
    removeFile();
}


void Foam::abortCalculation::timeSet()
{
    // Do nothing - only valid on execute
}


void Foam::abortCalculation::write()
{
    // Do nothing - only valid on execute
}


// ************************************************************************* //
