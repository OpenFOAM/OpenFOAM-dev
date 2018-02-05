/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
namespace functionObjects
{
    defineTypeNameAndDebug(setTimeStepFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        setTimeStepFunctionObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::setTimeStepFunctionObject::setTimeStepFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    time_(runTime)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::setTimeStepFunctionObject::~setTimeStepFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::Time&
Foam::functionObjects::setTimeStepFunctionObject::time() const
{
    return time_;
}


bool Foam::functionObjects::setTimeStepFunctionObject::setTimeStep()
{
    const_cast<Time&>(time()).setDeltaTNoAdjust
    (
        timeStepPtr_().value(time_.timeOutputValue())
    );

    return true;
}


bool Foam::functionObjects::setTimeStepFunctionObject::read
(
    const dictionary& dict
)
{
    timeStepPtr_ = Function1<scalar>::New("deltaT", dict);

    return true;
}


bool Foam::functionObjects::setTimeStepFunctionObject::execute()
{
    bool adjustTimeStep =
        time().controlDict().lookupOrDefault("adjustTimeStep", false);

    if (!adjustTimeStep)
    {
        return setTimeStep();
    }

    return true;
}


bool Foam::functionObjects::setTimeStepFunctionObject::write()
{
    return true;
}


// ************************************************************************* //
