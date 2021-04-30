/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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
    Foam::functionObjects::stopAtClockTime

Description
    Stops the run when the specified clock time in seconds has been reached
    and optionally write results before stopping.

    The following actions are supported:
    - noWriteNow
    - writeNow
    - nextWrite (default)

    Examples of function object specification:
    \verbatim
    stop
    {
        type        stopAtClockTime;
        libs        ("libutilityFunctionObjects.so");
        stopTime    10;
        action      writeNow;
    }
    \endverbatim
    will write the fields and stop the run once the clock time has exceeded 10
    seconds.

Usage
    \table
        Property | Description              | Required | Default value
        type     | type name: stopAtClockTime | yes    |
        stopTime | Maximum elapsed time [s] | yes      |
        action   | Action executed          | no       | nextWrite
    \endtable

SourceFiles
    stopAtClockTime.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_stopAtClockTime_H
#define functionObjects_stopAtClockTime_H

#include "stopAt.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class stopAtClockTime Declaration
\*---------------------------------------------------------------------------*/

class stopAtClockTime
:
    public stopAt
{
    // Private Data

        //-
        time_t stopTime_;


    // Private Member Functions

        //- Return true when the stop condition is achieved
        virtual bool condition() const;


public:

    //- Runtime type information
    TypeName("stopAtClockTime");


    // Constructors

        //- Construct from Time and dictionary
        stopAtClockTime
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        stopAtClockTime(const stopAtClockTime&) = delete;


    //- Destructor
    virtual ~stopAtClockTime();


    // Member Functions

        //- Read the dictionary settings
        virtual bool read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const stopAtClockTime&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
