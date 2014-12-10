/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "setTimeStepFunctionObject.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(setTimeStepFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        setTimeStepFunctionObject,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::setTimeStepFunctionObject::setTimeStepFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    time_(runTime),
    enabled_(true)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::setTimeStepFunctionObject::on()
{
    enabled_ = true;
}


void Foam::setTimeStepFunctionObject::off()
{
    enabled_ = false;
}


bool Foam::setTimeStepFunctionObject::start()
{
    return true;
}


bool Foam::setTimeStepFunctionObject::execute(const bool forceWrite)
{
    return true;
}


bool Foam::setTimeStepFunctionObject::end()
{
    return true;
}


bool Foam::setTimeStepFunctionObject::timeSet()
{
    return true;
}


bool Foam::setTimeStepFunctionObject::adjustTimeStep()
{
    if (enabled())
    {
        // Wanted timestep
        scalar newDeltaT = timeStepPtr_().value(time_.timeOutputValue());

        const_cast<Time&>(time()).setDeltaT(newDeltaT, false);

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::setTimeStepFunctionObject::read(const dictionary& dict)
{
    enabled_ = dict.lookupOrDefault("enabled", true);

    if (enabled_)
    {
        timeStepPtr_ = DataEntry<scalar>::New("deltaT", dict);

        // Check that time has adjustTimeStep
        const dictionary& controlDict = time_.controlDict();

        Switch adjust;
        if
        (
            !controlDict.readIfPresent<Switch>("adjustTimeStep", adjust)
         || !adjust
        )
        {
            FatalIOErrorIn("setTimeStep::read(const dictionary&)", dict)
                << "Need to have 'adjustTimeStep' true to enable external"
                << " timestep control" << exit(FatalIOError);
        }
    }
    return true;
}


void Foam::setTimeStepFunctionObject::updateMesh(const mapPolyMesh& mpm)
{}


void Foam::setTimeStepFunctionObject::movePoints(const polyMesh& mesh)
{}


// ************************************************************************* //
