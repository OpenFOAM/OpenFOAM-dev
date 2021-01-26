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

#include "sixDoFRigidBodyControl.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(sixDoFRigidBodyControl, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        sixDoFRigidBodyControl,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::sixDoFRigidBodyControl::sixDoFRigidBodyControl
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    sixDoFRigidBodyState(name, runTime, dict),
    time_(runTime),
    meanVelocity_(Zero),
    meanAngularVelocity_(Zero)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::sixDoFRigidBodyControl::~sixDoFRigidBodyControl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::sixDoFRigidBodyControl::read(const dictionary& dict)
{
    sixDoFRigidBodyState::read(dict);

    dict.lookup("window") >> w_;
    dict.lookup("convergedVelocity") >> convergedVelocity_;
    dict.lookup("convergedAngularVelocity") >> convergedAngularVelocity_;

    resetName(typeName);

    return true;
}


bool Foam::functionObjects::sixDoFRigidBodyControl::execute()
{
    if (time_.timeIndex() <= time_.startTimeIndex() + 1)
    {
        meanVelocity_ = cmptMag(velocity());
        meanAngularVelocity_ = cmptMag(angularVelocity());
    }
    else
    {
        const scalar dt = time_.deltaTValue();
        const scalar beta = min(dt/w_, 1);

        meanVelocity_ = (1 - beta)*meanVelocity_ + beta*cmptMag(velocity());

        meanAngularVelocity_ =
            (1 - beta)*meanAngularVelocity_ + beta*cmptMag(angularVelocity());
    }

    if
    (
        time_.value() - time_.startTime().value() > w_
     && meanVelocity_ < convergedVelocity_
     && meanAngularVelocity_ < convergedAngularVelocity_
    )
    {
        time_.stopAt(Time::stopAtControl::writeNow);
    }

    return true;
}


// ************************************************************************* //
