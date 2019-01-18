/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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

#include "timeControlFunctionObject.H"
#include "polyMesh.H"
#include "mapPolyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(timeControl, 0);
}
}


// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

void Foam::functionObjects::timeControl::readControls()
{
    if (!dict_.readIfPresent("startTime", startTime_))
    {
        dict_.readIfPresent("timeStart", startTime_);
    }

    if (!dict_.readIfPresent("endTime", endTime_))
    {
        dict_.readIfPresent("timeEnd", endTime_);
    }

    dict_.readIfPresent("nStepsToStartTimeChange", nStepsToStartTimeChange_);
}


bool Foam::functionObjects::timeControl::active() const
{
    return
        time_.value() >= startTime_
     && time_.value() <= endTime_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::timeControl::timeControl
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    time_(t),
    dict_(dict),
    startTime_(-vGreat),
    endTime_(vGreat),
    nStepsToStartTimeChange_
    (
        dict.lookupOrDefault("nStepsToStartTimeChange", 3)
    ),
    executeControl_(t, dict, "execute"),
    writeControl_(t, dict, "write"),
    foPtr_(functionObject::New(name, t, dict_))
{
    readControls();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::timeControl::execute()
{
    if (active() && (postProcess || executeControl_.execute()))
    {
        foPtr_->execute();
    }

    return true;
}


bool Foam::functionObjects::timeControl::write()
{
    if (active() && (postProcess || writeControl_.execute()))
    {
        foPtr_->write();
    }

    return true;
}


bool Foam::functionObjects::timeControl::end()
{
    if (active() && (executeControl_.execute() || writeControl_.execute()))
    {
        foPtr_->end();
    }

    return true;
}



Foam::scalar Foam::functionObjects::timeControl::timeToNextWrite()
{
    if
    (
        active()
     && writeControl_.control() ==
        Foam::timeControl::timeControls::adjustableRunTime
    )
    {
        const label  writeTimeIndex = writeControl_.executionIndex();
        const scalar writeInterval = writeControl_.interval();

        return
            max
            (
                0.0,
                (writeTimeIndex + 1)*writeInterval
              - (time_.value() - time_.startTime().value())
            );
    }

    return vGreat;
}


bool Foam::functionObjects::timeControl::read
(
    const dictionary& dict
)
{
    if (dict != dict_)
    {
        dict_ = dict;

        writeControl_.read(dict);
        executeControl_.read(dict);
        readControls();

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::functionObjects::timeControl::updateMesh
(
    const mapPolyMesh& mpm
)
{
    if (active())
    {
        foPtr_->updateMesh(mpm);
    }
}


void Foam::functionObjects::timeControl::movePoints
(
    const polyMesh& mesh
)
{
    if (active())
    {
        foPtr_->movePoints(mesh);
    }
}


// ************************************************************************* //
