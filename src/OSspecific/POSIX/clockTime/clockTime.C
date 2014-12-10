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

#include "clockTime.H"
#include <sys/time.h>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::clockTime::getTime(timeType& t)
{
    gettimeofday(&t, 0);
}


double Foam::clockTime::timeDifference(const timeType& beg, const timeType& end)
{
    return end.tv_sec - beg.tv_sec + 1e-6*(end.tv_usec - beg.tv_usec);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clockTime::clockTime()
{
    getTime(startTime_);
    lastTime_ = startTime_;
    newTime_ = startTime_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

double Foam::clockTime::elapsedTime() const
{
    getTime(newTime_);
    return timeDifference(startTime_, newTime_);
}


double Foam::clockTime::timeIncrement() const
{
    lastTime_ = newTime_;
    getTime(newTime_);
    return timeDifference(lastTime_, newTime_);
}


// ************************************************************************* //
