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

#include "userTimeStepFunctionObject.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(userTimeStep, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        userTimeStep,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::userTimeStep::userTimeStep
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    regionFunctionObject(name, runTime, dict),
    logFiles(obr_, name)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::userTimeStep::~userTimeStep()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::userTimeStep::read(const dictionary& dict)
{
    functionObject::read(dict);

    resetName(typeName);

    return true;
}


void Foam::functionObjects::userTimeStep::writeFileHeader(const label i)
{
    if (Pstream::master())
    {
        writeHeader(file(), "User time step");
        writeCommented(file(), "Time");

        file() << tab << "deltaT";

        file() << endl;
    }
}


bool Foam::functionObjects::userTimeStep::execute()
{
    return true;
}


bool Foam::functionObjects::userTimeStep::write()
{
    logFiles::write();

    if (Pstream::master())
    {
        writeTime(file());

        file() << tab << time_.userDeltaTValue();

        file() << endl;
    }

    return true;
}


// ************************************************************************* //
