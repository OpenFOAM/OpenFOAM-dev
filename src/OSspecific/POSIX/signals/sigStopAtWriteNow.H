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

Class
    Foam::sigStopAtWriteNow

Description
    Signal handler for interupt defined by
        OptimisationSwitches::stopAtWriteNowSignal

    Write and stop the job.

SourceFiles
    sigStopAtWriteNow.C

\*---------------------------------------------------------------------------*/

#ifndef sigStopAtWriteNow_H
#define sigStopAtWriteNow_H

#include <signal.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Time;

/*---------------------------------------------------------------------------*\
                           Class sigStopAtWriteNow Declaration
\*---------------------------------------------------------------------------*/

class sigStopAtWriteNow
{
    // Private data

        //- number of signal to use
        static int signal_;

        //- Saved old signal trapping setting
        static struct sigaction oldAction_;

    // Private Member Functions

        static void sigHandler(int);


public:

    //- wip. Have setter have access to signal_
    friend class addstopAtWriteNowSignalToOpt;

    // Constructors

        //- Construct null
        sigStopAtWriteNow();

        //- Construct from components
        sigStopAtWriteNow(const bool verbose, const Time& runTime);


    //- Destructor
    ~sigStopAtWriteNow();


    // Member functions

        //- (re)set signal catcher
        static void set(const bool verbose);

        //- Is active?
        bool active() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
