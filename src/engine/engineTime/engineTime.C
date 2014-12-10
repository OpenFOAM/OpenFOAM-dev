/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "engineTime.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::engineTime::timeAdjustment()
{
    deltaT_  = degToTime(deltaT_);
    endTime_ = degToTime(endTime_);

    if
    (
        writeControl_ == wcRunTime
     || writeControl_ == wcAdjustableRunTime
    )
    {
        writeInterval_ = degToTime(writeInterval_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from objectRegistry arguments
Foam::engineTime::engineTime
(
    const word& name,
    const fileName& rootPath,
    const fileName& caseName,
    const fileName& systemName,
    const fileName& constantName,
    const fileName& dictName
)
:
    Time
    (
        name,
        rootPath,
        caseName,
        systemName,
        constantName
    ),
    dict_
    (
        IOobject
        (
            "engineGeometry",
            constant(),
            *this,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
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

// Read the controlDict and set all the parameters
void Foam::engineTime::readDict()
{
    Time::readDict();
    timeAdjustment();
}


// Read the controlDict and set all the parameters
bool Foam::engineTime::read()
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


Foam::scalar Foam::engineTime::degToTime(const scalar theta) const
{
    // 6 * rpm => deg/s
    return theta/(6.0*rpm_.value());
}


Foam::scalar Foam::engineTime::timeToDeg(const scalar t) const
{
    // 6 * rpm => deg/s
    return t*(6.0*rpm_.value());
}


Foam::scalar Foam::engineTime::theta() const
{
    return timeToDeg(value());
}


// Return current crank-angle translated to a single revolution
// (value between -180 and 180 with 0 = top dead centre)
Foam::scalar Foam::engineTime::thetaRevolution() const
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


Foam::scalar Foam::engineTime::deltaTheta() const
{
    return timeToDeg(deltaTValue());
}


Foam::scalar Foam::engineTime::pistonPosition(const scalar theta) const
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


Foam::dimensionedScalar Foam::engineTime::pistonPosition() const
{
    return dimensionedScalar
    (
        "pistonPosition",
        dimLength,
        pistonPosition(theta())
    );
}


Foam::dimensionedScalar Foam::engineTime::pistonDisplacement() const
{
    return dimensionedScalar
    (
        "pistonDisplacement",
        dimLength,
        pistonPosition(theta() - deltaTheta()) - pistonPosition().value()
    );
}


Foam::dimensionedScalar Foam::engineTime::pistonSpeed() const
{
    return dimensionedScalar
    (
        "pistonSpeed",
        dimVelocity,
        pistonDisplacement().value()/(deltaTValue() + VSMALL)
    );
}


Foam::scalar Foam::engineTime::userTimeToTime(const scalar theta) const
{
    return degToTime(theta);
}


Foam::scalar Foam::engineTime::timeToUserTime(const scalar t) const
{
    return timeToDeg(t);
}


// ************************************************************************* //
