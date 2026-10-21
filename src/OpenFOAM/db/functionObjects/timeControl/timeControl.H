/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
        static const NamedEnum<timeControls, 8> timeControlNames_;

        //- Type of time control
        timeControls timeControl_;

        //- Start time for execution - defaults to Time::beginTime
        scalar startTime_;

        //- End time for execution - defaults to vGreat
        scalar endTime_;

        //- Begin time for execution - defaults to Time::beginTime
        scalar beginTime_;

        //- Counter for the number of write times
        label writeTimeIndex_;

        //- Times at which intervals apply
        scalarList intervalTimes_;

        //- Integer interval sizes. Counts of the numbers of steps between
        //  executions. For timeStep, writeTime, outputTime.
        List<label> intervalCounts_;

        //- Scalar interval sizes. Physical or compute times between
        //  executions. For adjustableRunTime, runTime, runTime, cpuTime.
        List<scalar> intervalPeriods_;

        //- When we get to the end of the list of intervals, should we start
        //  over, or should we use the last interval forever
        bool repeat_;

        //- Index of the current interval
        label intervalIndex_;

        //- Current execution index. Incremented on execution to mark the time
        //  as having been executed and prevent re-execution.
        label executionIndex_;

        //- Offset for counts
        label offsetCount_;

        //- Offset for periods
        scalar offsetPeriod_;


    // Private Member Functions

        //- Read the intervals
        template<class Type, class Units>
        void readIntervals
        (
            const word& intervalName,
            const dictionary& dict,
            const Units& units,
            List<Type>& intervalSizes
        );

        //- Read the intervals with a default single interval
        template<class Type, class Units>
        void readIntervals
        (
            const word& intervalName,
            const dictionary& dict,
            const Units& units,
            const bool allowDefault,
            const Type& defaultInterval,
            List<Type>& intervalSizes
        );

        //- Return the amount of the given time spent repeating
        scalar repeating(const scalar t) const;

        //- Return a time just after the given time
        scalar after(const scalar t) const;

        //- Update the interval
        template<class Type>
        void update
        (
            const List<Type>& intervalSizes,
            const scalar t,
            const Type& x0,
            const Type& x,
            label& intervalIndex,
            label& executionIndex,
            Type& offset
        ) const;

        //- Return whether or not to execute
        template<class Type>
        bool execute
        (
            const List<Type>& intervalSizes,
            const scalar t,
            const Type& x0,
            const Type& x,
            const Type& dx,
            Type& offset
        );


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
        scalar timeToNextAction() const;

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
