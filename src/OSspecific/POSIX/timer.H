/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::timer

Description
    Implements a timeout mechanism via sigalarm.

    Example usage:
    \code
        timer myTimer(5);     // 5 sec
        ..
        if (timedOut(myTimer))
        {
            // timed out
        }
        else
        {
            // do something possible blocking
        }
    \endcode

    Constructor set signal handler on sigalarm and alarm(). Destructor
    clears these.

    timedOut is macro because setjmp can't be in member function of timer.
    ?something to do with stack frames.

Warning
    The setjmp restores complete register state so including local vars
    held in regs. So if in blocking part something gets calced in a stack
    based variable make sure it is declared 'volatile'.

SourceFiles
    timer.C

\*---------------------------------------------------------------------------*/

#ifndef timer_H
#define timer_H

#include "className.H"

#include <signal.h>
#include <setjmp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Check it a timeout has occurred
//  keep setjmp in same stack frame so no function calls
#define timedOut(x) \
    (((x).newTimeOut_ > 0) ? setjmp(Foam::timer::envAlarm) : false)

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class timer Declaration
\*---------------------------------------------------------------------------*/

class timer
{
    // Private Data

        //- Old signal masks
        static struct sigaction oldAction_;

        //- Old alarm() value
        static unsigned int oldTimeOut_;


    // Private Member Functions

        //- Alarm handler
        static void signalHandler(int);


public:

    // Public data

        //- Declare name of the class and its debug switch
        ClassName("timer");

        //- Current time out value. Needed by macro timedOut
        unsigned int newTimeOut_;

        //- State for setjmp. Needed by macro timedOut
        static jmp_buf envAlarm;


    // Constructors

        //- Construct from components.
        //  newTimeOut=0 makes it do nothing.
        timer(const unsigned int newTimeOut);


    //- Destructor
    ~timer();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
