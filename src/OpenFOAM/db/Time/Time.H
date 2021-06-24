/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::Time

Description
    Class to control time during OpenFOAM simulations that is also the
    top-level objectRegistry.

SourceFiles
    Time.C
    TimeIO.C
    findInstance.C

\*---------------------------------------------------------------------------*/

#ifndef Time_H
#define Time_H

#include "TimePaths.H"
#include "objectRegistry.H"
#include "IOdictionary.H"
#include "FIFOStack.H"
#include "clock.H"
#include "cpuTime.H"
#include "TimeState.H"
#include "Switch.H"
#include "instantList.H"
#include "NamedEnum.H"
#include "typeInfo.H"
#include "dlLibraryTable.H"
#include "functionObjectList.H"
#include "sigWriteNow.H"
#include "sigStopAtWriteNow.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class argList;

/*---------------------------------------------------------------------------*\
                            Class Time Declaration
\*---------------------------------------------------------------------------*/

class Time
:
    public clock,
    public cpuTime,
    public TimePaths,
    public objectRegistry,
    public TimeState
{
    // Private Data

        //- Is runtime modification of dictionaries allowed?
        Switch runTimeModifiable_;

        //- The controlDict
        IOdictionary controlDict_;


public:

        //- Write control options
        enum class writeControl
        {
            timeStep,
            runTime,
            adjustableRunTime,
            clockTime,
            cpuTime
        };

        //- Stop-run control options
        enum class stopAtControl
        {
            endTime,    //!< stop when Time reaches the prescribed endTime
            noWriteNow, //!< set endTime to stop immediately w/o writing
            writeNow,   //!< set endTime to stop immediately w/ writing
            nextWrite   //!< stop the next time data are written
        };

        //- Supported time directory name formats
        enum class format
        {
            general    = 0,
            fixed      = ios_base::fixed,
            scientific = ios_base::scientific
        };


protected:

    // Protected data

        label  startTimeIndex_;
        scalar startTime_;
        mutable scalar endTime_;
        scalar beginTime_;

        static const NamedEnum<stopAtControl, 4> stopAtControlNames_;
        mutable stopAtControl stopAt_;

        static const NamedEnum<writeControl, 5> writeControlNames_;
        writeControl writeControl_;

        scalar writeInterval_;

        label  purgeWrite_;
        mutable FIFOStack<word> previousWriteTimes_;

        // One-shot writing
        bool writeOnce_;

        //- Is the time currently being sub-cycled?
        bool subCycling_;

        //- If time is being sub-cycled this is the previous TimeState
        autoPtr<TimeState> prevTimeState_;


        // Signal handlers for writing

            //- Enable one-shot writing upon signal
            sigWriteNow sigWriteNow_;

            //- Enable write and clean exit upon signal
            sigStopAtWriteNow sigStopAtWriteNow_;


        //- Time directory name format
        static format format_;

        //- Time directory name precision
        static int precision_;

        //- Maximum time directory name precision
        static const int maxPrecision_;

        //- Adjust the time step so that writing occurs at the specified time
        void adjustDeltaT();

        //- Set the controls from the current controlDict
        void setControls();

        //- Read the control dictionary and set the write controls etc.
        virtual void readDict();


private:

        //- Default write option
        IOstream::streamFormat writeFormat_;

        //- Default output file format version number
        IOstream::versionNumber writeVersion_;

        //- Default output compression
        IOstream::compressionType writeCompression_;

        //- Default graph format
        word graphFormat_;

        //- Is temporary object cache enabled
        mutable bool cacheTemporaryObjects_;

        //- Function objects executed at start and on ++, +=
        mutable functionObjectList functionObjects_;


public:

    TypeName("time");

    //- The default control dictionary name (normally "controlDict")
    static word controlDictName;


    // Constructors

        //- Construct given name of dictionary to read and argument list
        Time
        (
            const word& name,
            const argList& args,
            const word& systemName = "system",
            const word& constantName = "constant"
        );

        //- Construct given name of dictionary to read, rootPath and casePath
        Time
        (
            const word& name,
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant",
            const bool enableFunctionObjects = true
        );

        //- Construct given dictionary, rootPath and casePath
        Time
        (
            const dictionary& dict,
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant",
            const bool enableFunctionObjects = true
        );

        //- Construct given endTime, rootPath and casePath
        Time
        (
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant",
            const bool enableFunctionObjects = true
        );


    //- Destructor
    virtual ~Time();


    // Member Functions

        // Database functions

            //- Explicitly inherit rootPath from TimePaths to disambiguate from
            //  the corresponding method in objectRegistry
            using TimePaths::rootPath;

            //- Explicitly inherit caseName from TimePaths to disambiguate from
            //  the corresponding method in objectRegistry
            using TimePaths::caseName;

            //- Explicitly inherit path from TimePaths to disambiguate from
            //  the corresponding method in objectRegistry
            using TimePaths::path;

            //- Return the control dict
            const dictionary& controlDict() const
            {
                return controlDict_;
            }

            //- Local directory path of this objectRegistry relative to time
            virtual const fileName& dbDir() const
            {
                return fileName::null;
            }

            //- Return current time path
            fileName timePath() const
            {
                return path()/timeName();
            }

            //- Default write format
            IOstream::streamFormat writeFormat() const
            {
                return writeFormat_;
            }

            //- Default write version number
            IOstream::versionNumber writeVersion() const
            {
                return writeVersion_;
            }

            //- Default write compression
            IOstream::compressionType writeCompression() const
            {
                return writeCompression_;
            }

            //- Default graph format
            const word& graphFormat() const
            {
                return graphFormat_;
            }

            //- Supports re-reading
            const Switch& runTimeModifiable() const
            {
                return runTimeModifiable_;
            }

            //- Read control dictionary, update controls and time
            virtual bool read();

            //- Read the objects that have been modified
            void readModifiedObjects();

            //- Return the location of "dir" containing the file "name".
            //  (eg, used in reading mesh data)
            //  If name is null, search for the directory "dir" only.
            //  Does not search beyond stopInstance (if set) or constant.
            word findInstance
            (
                const fileName& dir,
                const word& name = word::null,
                const IOobject::readOption rOpt = IOobject::MUST_READ,
                const word& stopInstance = word::null
            ) const;

            //- Search the case for valid time directories
            instantList times() const;

            //- Search the case for the time directory path
            //  corresponding to the given instance
            word findInstancePath(const fileName& path, const instant&) const;

            //- Search the case for the time directory path
            //  corresponding to the given instance
            word findInstancePath(const instant&) const;

            //- Search the case for the time closest to the given time
            instant findClosestTime(const scalar) const;

            //- Search instantList for the time index closest to the given time
            static label findClosestTimeIndex
            (
                const instantList&,
                const scalar,
                const word& constantName = "constant"
            );

            //- Write time dictionary to the \<time\>/uniform directory
            virtual bool writeTimeDict() const;

            //- Write using given format, version and compression
            virtual bool writeObject
            (
                IOstream::streamFormat,
                IOstream::versionNumber,
                IOstream::compressionType,
                const bool write
            ) const;

            //- Write the objects now (not at end of iteration) and continue
            //  the run
            bool writeNow();

            //- Write the objects now (not at end of iteration) and end the run
            bool writeAndEnd();

            //- Write the objects once (one shot) and continue the run
            void writeOnce();


        // Access

            //- Return time name of given scalar time
            //  formatted with given precision
            static word timeName
            (
                const scalar,
                const int precision = precision_
            );

            //- Return current time name
            virtual word timeName() const;

            //- Search a given directory for valid time directories
            static instantList findTimes
            (
                const fileName&,
                const word& constantName = "constant"
            );

            //- Return start time index
            virtual label startTimeIndex() const;

            //- Return begin time (initial start time)
            virtual dimensionedScalar beginTime() const;

            //- Return start time
            virtual dimensionedScalar startTime() const;

            //- Return end time
            virtual dimensionedScalar endTime() const;

            //- Return the list of function objects
            const functionObjectList& functionObjects() const
            {
                return functionObjects_;
            }

            //- Return true if the run is a restart, i.e. startTime != beginTime
            bool restart() const
            {
                return startTime_ != beginTime_;
            }

            //- Return true if time currently being sub-cycled, otherwise false
            bool subCycling() const
            {
                return subCycling_;
            }

            //- Return previous TimeState if time is being sub-cycled
            const TimeState& prevTimeState() const
            {
                return prevTimeState_();
            }


        // Check

            //- Return true if run should continue without any side effects
            virtual bool running() const;

            //- Return true if run should continue,
            //  also invokes the functionObjectList::end() method
            //  when the time goes out of range
            //  \note
            //  For correct behaviour, the following style of time-loop
            //  is recommended:
            //  \code
            //      while (runTime.run())
            //      {
            //          runTime++;
            //          solve;
            //          runTime.write();
            //      }
            //  \endcode
            virtual bool run() const;

            //- Return true if run should continue and if so increment time
            //  also invokes the functionObjectList::end() method
            //  when the time goes out of range
            //  \note
            //  For correct behaviour, the following style of time-loop
            //  is recommended:
            //  \code
            //      while (runTime.loop())
            //      {
            //          solve;
            //          runTime.write();
            //      }
            //  \endcode
            virtual bool loop();

            //- Return true if end of run,
            //  does not invoke any functionObject methods
            //  \note
            //      The rounding heuristics near endTime mean that
            //      \code run() \endcode and \code !end() \endcode may
            //      not yield the same result
            virtual bool end() const;


        // Edit

            //- Adjust the current stopAtControl. Note that this value
            //  only persists until the next time the dictionary is read.
            //  Return true if the stopAtControl changed.
            virtual bool stopAt(const stopAtControl) const;

            //- Reset the time and time-index to those of the given time
            virtual void setTime(const Time&);

            //- Reset the time and time-index
            virtual void setTime(const instant&, const label newIndex);

            //- Reset the time and time-index
            virtual void setTime
            (
                const dimensionedScalar&,
                const label newIndex
            );

            //- Reset the time and time-index
            virtual void setTime(const scalar, const label newIndex);

            //- Reset end time
            virtual void setEndTime(const dimensionedScalar&);

            //- Reset end time
            virtual void setEndTime(const scalar);

            //- Reset time step
            virtual void setDeltaT(const dimensionedScalar&);

            //- Reset time step
            virtual void setDeltaT(const scalar);

            //- Reset time step without additional adjustment or modification
            //  by function objects
            virtual void setDeltaTNoAdjust(const scalar);

            //- Reset the write interval
            virtual void setWriteInterval(const scalar writeInterval);

            //- Set time to sub-cycle for the given number of steps
            virtual TimeState subCycle(const label nSubCycles);

            //- Reset time after sub-cycling back to previous TimeState
            virtual void endSubCycle();

            //- Return non-const access to the list of function objects
            functionObjectList& functionObjects()
            {
                return functionObjects_;
            }


    // Member Operators

        //- Set deltaT to that specified and increment time via operator++()
        virtual Time& operator+=(const dimensionedScalar&);

        //- Set deltaT to that specified and increment time via operator++()
        virtual Time& operator+=(const scalar);

        //- Prefix increment,
        //  also invokes the functionObjectList::start() or
        //  functionObjectList::execute() method, depending on the time-index
        virtual Time& operator++();

        //- Postfix increment, this is identical to the prefix increment
        virtual Time& operator++(int);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
