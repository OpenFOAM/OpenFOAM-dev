/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2024 OpenFOAM Foundation
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
#include "polyTopoChangeMap.H"
#include "polyMeshMap.H"
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

bool Foam::functionObjects::timeControl::active() const
{
    return executeControl_.active() || writeControl_.active();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::timeControl::timeControl
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name, t),
    time_(t),
    executeControl_(t, dict, "execute"),
    writeControl_(t, dict, "write"),
    foPtr_(functionObject::New(name, t, dict))
{
    writeControl_.read(dict);
    executeControl_.read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::timeControl::fields() const
{
    return foPtr_->fields();
}


bool Foam::functionObjects::timeControl::executeAtStart() const
{
    return foPtr_->executeAtStart();
}


bool Foam::functionObjects::timeControl::execute()
{
    if
    (
        executeControl_.active()
     && (
            postProcess
         || executeControl_.execute()
         || (executeAtStart() && time_.timeIndex() == time_.startTimeIndex())
        )
    )
    {
        foPtr_->execute();
    }

    return true;
}


bool Foam::functionObjects::timeControl::write()
{
    if
    (
        writeControl_.active()
     && (
            postProcess
         || writeControl_.execute()
         || (executeAtStart() && time_.timeIndex() == time_.startTimeIndex())
        )
    )
    {
        foPtr_->write();
    }

    return true;
}


bool Foam::functionObjects::timeControl::end()
{
    if (executeControl_.execute() || writeControl_.execute())
    {
        foPtr_->end();
    }

    return true;
}


Foam::scalar Foam::functionObjects::timeControl::timeToNextAction()
{
    return min
    (
        executeControl_.timeToNextAction(),
        writeControl_.timeToNextAction()
    );
}


bool Foam::functionObjects::timeControl::read(const dictionary& dict)
{
    writeControl_.read(dict);
    executeControl_.read(dict);

    if (active())
    {
        foPtr_->read(dict);
    }

    return true;
}


void Foam::functionObjects::timeControl::movePoints(const polyMesh& mesh)
{
    if (active())
    {
        foPtr_->movePoints(mesh);
    }
}


void Foam::functionObjects::timeControl::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (active())
    {
        foPtr_->topoChange(map);
    }
}


void Foam::functionObjects::timeControl::mapMesh
(
    const polyMeshMap& map
)
{
    if (active())
    {
        foPtr_->mapMesh(map);
    }
}


void Foam::functionObjects::timeControl::distribute
(
    const polyDistributionMap& map
)
{
    if (active())
    {
        foPtr_->distribute(map);
    }
}


// ************************************************************************* //
