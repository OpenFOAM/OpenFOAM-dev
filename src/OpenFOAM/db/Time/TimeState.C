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

#include "TimeState.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TimeState::TimeState()
:
    dimensionedScalar(Time::timeName(0), dimTime, 0),
    timeIndex_(0),
    deltaT_(0),
    deltaTSave_(0),
    deltaT0_(0),
    deltaTchanged_(false),
    outputTimeIndex_(0),
    primaryOutputTime_(false),
    secondaryOutputTimeIndex_(0),
    secondaryOutputTime_(false),
    outputTime_(false)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * //

Foam::TimeState::~TimeState()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::TimeState::userTimeToTime(const scalar theta) const
{
    return theta;
}


Foam::scalar Foam::TimeState::timeToUserTime(const scalar t) const
{
    return t;
}


Foam::scalar Foam::TimeState::timeOutputValue() const
{
    return timeToUserTime(value());
}


Foam::label Foam::TimeState::timeIndex() const
{
    return timeIndex_;
}


Foam::dimensionedScalar Foam::TimeState::deltaT() const
{
    return dimensionedScalar("deltaT", dimTime, deltaT_);
}


Foam::dimensionedScalar Foam::TimeState::deltaT0() const
{
    return dimensionedScalar("deltaT0", dimTime, deltaT0_);
}


bool Foam::TimeState::outputTime() const
{
    return outputTime_;
}


// ************************************************************************* //
