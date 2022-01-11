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

Class
    Foam::cpuTime

Description
    Starts timing CPU usage and return elapsed time from start.

    Uses the POSIX clock() function which return the processor time consumed
    in clock ticks @ CLOCKS_PER_SEC clock ticks per second.

See also
    clockTime

SourceFiles
    cpuTime.C

\*---------------------------------------------------------------------------*/

#ifndef cpuTime_H
#define cpuTime_H

#include <time.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cpuTime Declaration
\*---------------------------------------------------------------------------*/

class cpuTime
{
    // Private Data

        //- Start CPU clock tick count
        clock_t startTime_;

        //- Previous CPU clock tick count
        mutable clock_t prevTime_;

        //- Current CPU clock tick count
        mutable clock_t curTime_;


    // Private Member Functions

        //- Return the difference between two times clock tick counts
        //  as elapsed CPU time (in seconds)
        inline double timeDifference
        (
            const clock_t prev,
            const clock_t cur
        ) const;


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
