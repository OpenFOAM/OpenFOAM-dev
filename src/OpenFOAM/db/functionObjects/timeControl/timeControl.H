/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::timeControl

Description
    General time dependent execution controller.
    The default to execute every time-step.

SourceFiles
    timeControl.C

\*---------------------------------------------------------------------------*/

#ifndef timeControl_H
#define timeControl_H

#include "Time.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class timeControl Declaration
\*---------------------------------------------------------------------------*/

class timeControl
{
public:

    //- The time control options
    enum class timeControls
    {
        timeStep,           //!< Execute/write every 'Interval' time-steps
        writeTime,          //!< Execute/write every 'Interval' write times
        outputTime,         //!< Execute/write every 'Interval' write times
        adjustableRunTime,  //!< Execute/write every 'Interval' run time period
        runTime,            //!< Execute/write every 'Interval' run time period
        runTimes,           //!< Execute/write at specified list of run times
        clockTime,          //!< Execute/write every 'Interval' CPU time period
        cpuTime,            //!< Execute/write every 'Interval' CPU time period
        none                //!< No execution
    };


private:

    // Private Data

        //- Reference to time
        const Time& time_;

        //- Prefix
        const word prefix_;

        //- String representation of timeControls enums
        static const NamedEnum<timeControls, 9> timeControlNames_;

        //- Type of time control
        timeControls timeControl_;

        //- Start time for execution - defaults to Time::beginTime
        scalar startTime_;

        //- End time for execution - defaults to vGreat
        scalar endTime_;

        //- Begin time for execution - defaults to Time::beginTime
        scalar beginTime_;

        //- Execution interval steps for timeStep mode
        //  a value <= 1 means execute at every time step
        label intervalSteps_;

        //- Execution interval
        scalar interval_;

        //- List of times (for runTimes option)
        scalarList times_;

        //- Time delta used for time -> index (for runTimes option)
        scalar timeDelta_;

        //- Hash set of mesh mapping time indices (for runTimes option)
        HashSet<int64_t, Hash<int64_t>> timeIndices_;

        //- Index of previous execution
        mutable label executionIndex_;


    // Private Member Functions

        //- Round the given scalar down to the nearest label. As a static_cast
        //  to label, except that negative numbers are rounded down rather than
        //  towards zero.
        static label roundDown(const scalar);


public:

    // Constructors

        //- Construct from Time object and dictionary
        timeControl
        (
            const Time&,
            const dictionary&,
            const word& prefix
        );

        //- Disallow default bitwise copy construction
        timeControl(const timeControl&) = delete;


    //- Destructor
    ~timeControl();


    // Member Functions

        //- Read from dictionary
        void read(const dictionary&);

        //- Return Time
        inline const Time& time() const;

        //- Returns true the time is within the start to end period
        bool active() const;

        //- Flag to indicate whether to execute
        bool execute();

        //- Return the time to the next write
        scalar timeToNextAction();

        //- Return control
        inline timeControls control() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const timeControl&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "timeControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
