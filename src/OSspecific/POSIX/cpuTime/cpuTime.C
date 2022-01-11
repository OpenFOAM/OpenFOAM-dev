/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "cpuTime.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline double Foam::cpuTime::timeDifference
(
    const clock_t prev,
    const clock_t cur
) const
{
    return double(cur - prev)/CLOCKS_PER_SEC;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cpuTime::cpuTime()
:
    startTime_(clock()),
    prevTime_(startTime_),
    curTime_(startTime_)

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

double Foam::cpuTime::elapsedCpuTime() const
{
    curTime_ = clock();
    return timeDifference(startTime_, curTime_);
}


double Foam::cpuTime::cpuTimeIncrement() const
{
    prevTime_ = curTime_;
    curTime_ = clock();
    return timeDifference(prevTime_, curTime_);
}


// ************************************************************************* //
