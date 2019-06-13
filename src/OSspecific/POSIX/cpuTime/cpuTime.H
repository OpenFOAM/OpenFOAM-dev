/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

Class
    Foam::cpuTime

Description
    Starts timing CPU usage and return elapsed time from start.

See also
    clockTime

SourceFiles
    cpuTime.C

\*---------------------------------------------------------------------------*/

#ifndef cpuTime_H
#define cpuTime_H

#include <time.h>
#include <sys/times.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cpuTime Declaration
\*---------------------------------------------------------------------------*/

class cpuTime
{
    // Private Data

        //- Time structure used
        typedef struct tms timeType;

        //- Clock-ticks per second
        static const long Hz_;

        //- The start time
        timeType startTime_;
        mutable timeType lastTime_;
        mutable timeType newTime_;

    // Private Member Functions

        //- Retrieve the current time values from the system
        static void getTime(timeType&);

        //- Difference between two times
        static double timeDifference(const timeType& beg, const timeType& end);


public:

    // Constructors

        //- Construct with the current clock time
        cpuTime();


    // Member Functions

        //- Return CPU time (in seconds) from the start
        double elapsedCpuTime() const;

        //- Return CPU time (in seconds) since last call to cpuTimeIncrement()
        double cpuTimeIncrement() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
