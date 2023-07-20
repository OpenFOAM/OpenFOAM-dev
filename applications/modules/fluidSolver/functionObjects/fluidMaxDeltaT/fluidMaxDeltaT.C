/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "fluidMaxDeltaT.H"
#include "fluidSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fluidMaxDeltaT, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        fluidMaxDeltaT,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fluidMaxDeltaT::fluidMaxDeltaT
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fluidMaxDeltaT::~fluidMaxDeltaT()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fluidMaxDeltaT::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    maxCoPtr_ = Function1<scalar>::New("maxCo", dict);
    maxDeltaTPtr_ = Function1<scalar>::New("maxDeltaT", dict);

    return true;
}


bool Foam::functionObjects::fluidMaxDeltaT::execute()
{
    return true;
}


bool Foam::functionObjects::fluidMaxDeltaT::write()
{
    return true;
}


Foam::scalar Foam::functionObjects::fluidMaxDeltaT::maxDeltaT() const
{
    scalar deltaT =
        time_.userTimeToTime(maxDeltaTPtr_().value(time_.userTimeValue()));

    const scalar CoNum =
        mesh_.lookupObject<solvers::fluidSolver>(solver::typeName).CoNum;

    if (CoNum > small)
    {
        const scalar maxCo = maxCoPtr_().value(time_.userTimeValue());

        deltaT = min(deltaT, maxCo/CoNum*time_.deltaTValue());
    }

    return deltaT;
}


// ************************************************************************* //
