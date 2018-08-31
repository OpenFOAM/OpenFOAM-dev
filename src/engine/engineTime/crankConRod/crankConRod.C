/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "crankConRod.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(crankConRod, 0);
    addToRunTimeSelectionTable(engineTime, crankConRod, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::crankConRod::timeAdjustment()
{
    deltaT_  = degToTime(deltaT_);
    endTime_ = degToTime(endTime_);

    if
    (
        writeControl_ == writeControl::runTime
        || writeControl_ == writeControl::adjustableRunTime
    )
    {
        writeInterval_ = degToTime(writeInterval_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::crankConRod::crankConRod
(
    const word& name,
    const fileName& rootPath,
    const fileName& caseName,
    const fileName& systemName,
    const fileName& constantName,
    const fileName& dictName
)
:
    engineTime
    (
        name,
        rootPath,
        caseName,
        systemName,
        constantName
    ),
    rpm_(dict_.lookup("rpm")),
    conRodLength_(dimensionedScalar("conRodLength", dimLength, 0)),
    bore_(dimensionedScalar("bore", dimLength, 0)),
    stroke_(dimensionedScalar("stroke", dimLength, 0)),
    clearance_(dimensionedScalar("clearance", dimLength, 0))
{
    // geometric parameters are not strictly required for Time
    dict_.readIfPresent("conRodLength", conRodLength_);
    dict_.readIfPresent("bore", bore_);
    dict_.readIfPresent("stroke", stroke_);
    dict_.readIfPresent("clearance", clearance_);

    timeAdjustment();

    startTime_  = degToTime(startTime_);
    value()     = degToTime(value());

    deltaTSave_ = deltaT_;
    deltaT0_    = deltaT_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::crankConRod::readDict()
{
    Time::readDict();
    timeAdjustment();
}


bool Foam::crankConRod::read()
{
    if (Time::read())
    {
        timeAdjustment();
        return true;
    }
    else
    {
        return false;
    }
}


Foam::scalar Foam::crankConRod::degToTime(const scalar theta) const
{
    // 6 * rpm => deg/s
    return theta/(6.0*rpm_.value());
}


Foam::scalar Foam::crankConRod::timeToDeg(const scalar t) const
{
    // 6 * rpm => deg/s
    return t*(6.0*rpm_.value());
}


Foam::scalar Foam::crankConRod::theta() const
{
    return timeToDeg(value());
}


Foam::word Foam::crankConRod::unit() const
{
    return " CAD";
}


Foam::scalar Foam::crankConRod::thetaRevolution() const
{
    scalar t = theta();

    while (t > 180.0)
    {
        t -= 360.0;
    }

    while (t < -180.0)
    {
        t += 360.0;
    }

    return t;
}


Foam::scalar Foam::crankConRod::deltaTheta() const
{
    return timeToDeg(deltaTValue());
}


Foam::scalar Foam::crankConRod::pistonPosition(const scalar theta) const
{
    return
    (
        conRodLength_.value()
      + stroke_.value()/2.0
      + clearance_.value()
    )
  - (
        stroke_.value()*::cos(degToRad(theta))/2.0
      + ::sqrt
        (
            sqr(conRodLength_.value())
            - sqr(stroke_.value()*::sin(degToRad(theta))/2.0)
        )
    );
}



Foam::scalar Foam::crankConRod::userTimeToTime(const scalar theta) const
{
    return degToTime(theta);
}


Foam::scalar Foam::crankConRod::timeToUserTime(const scalar t) const
{
    return timeToDeg(t);
}


// ************************************************************************* //
