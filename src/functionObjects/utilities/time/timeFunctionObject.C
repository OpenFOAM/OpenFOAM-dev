/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2021 OpenFOAM Foundation
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

#include "timeFunctionObject.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(time, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        time,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::time::time
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    regionFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    perTimeStep_(false),
    cpuTime0_(time_.elapsedCpuTime()),
    clockTime0_(time_.elapsedClockTime())
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::time::~time()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::time::read(const dictionary& dict)
{
    functionObject::read(dict);

    dict.readIfPresent("perTimeStep", perTimeStep_);

    resetName(typeName);

    return true;
}


void Foam::functionObjects::time::writeFileHeader(const label i)
{
    if (Pstream::master())
    {
        writeHeader(file(), "time");
        writeCommented(file(), "Time");
        writeTabbed(file(), "cpu");
        writeTabbed(file(), "clock");

        if (perTimeStep_)
        {
            writeTabbed(file(), "cpu/step");
            writeTabbed(file(), "clock/step");
        }

        file() << endl;
    }
}


bool Foam::functionObjects::time::execute()
{
    return true;
}


bool Foam::functionObjects::time::write()
{
    logFiles::write();

    if (Pstream::master())
    {
        writeTime(file());

        const scalar cpuTime(time_.elapsedCpuTime());
        const scalar clockTime(time_.elapsedClockTime());

        file() << tab << time_.elapsedCpuTime();
        file() << tab << time_.elapsedClockTime();

        if (perTimeStep_)
        {
            file() << tab << time_.elapsedCpuTime() - cpuTime0_;
            file() << tab << time_.elapsedClockTime() - clockTime0_;
        }

        file() << endl;

        cpuTime0_ = cpuTime;
        clockTime0_ = clockTime;
    }

    return true;
}


// ************************************************************************* //
