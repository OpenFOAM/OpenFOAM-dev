/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#include "stopAt.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(stopAt, 0);
}
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::stopAt::actionType,
    3
>::names[] =
{
    "noWriteNow",
    "writeNow",
    "nextWrite"
};

const Foam::NamedEnum
<
    Foam::functionObjects::stopAt::actionType,
    3
> Foam::functionObjects::stopAt::actionTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::stopAt::stopAt
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    time_(runTime),
    action_(actionType::nextWrite),
    stopped_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::stopAt::~stopAt()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::stopAt::read(const dictionary& dict)
{
    if (dict.found("action"))
    {
        action_ = actionTypeNames_.read(dict.lookup("action"));
    }
    else
    {
        action_ = actionType::nextWrite;
    }

    return true;
}


bool Foam::functionObjects::stopAt::execute()
{
    if (!stopped_)
    {
        bool stopCondition = condition();
        reduce(stopCondition, orOp<bool>());

        if (stopCondition)
        {
            switch (action_)
            {
                case actionType::noWriteNow :
                {
                    if (time_.stopAt(Time::stopAtControl::noWriteNow))
                    {
                        Info<< type() << "(timeIndex="
                            << time_.timeIndex()
                            << "): stopping now without writing"
                            << endl;
                    }
                    break;
                }

                case actionType::writeNow :
                {
                    if (time_.stopAt(Time::stopAtControl::writeNow))
                    {
                        Info<< type() << "(timeIndex="
                            << time_.timeIndex()
                            << "): stopping now after writing"
                            << endl;
                    }
                    break;
                }

                case actionType::nextWrite :
                {
                    if (time_.stopAt(Time::stopAtControl::nextWrite))
                    {
                        Info<< type() << "(timeIndex="
                            << time_.timeIndex()
                            << "): stopping after next write"
                            << endl;
                    }
                    break;
                }
            }

            stopped_ = true;
        }
    }

    return true;
}


bool Foam::functionObjects::stopAt::write()
{
    return true;
}


bool Foam::functionObjects::stopAt::end()
{
    return true;
}


// ************************************************************************* //
