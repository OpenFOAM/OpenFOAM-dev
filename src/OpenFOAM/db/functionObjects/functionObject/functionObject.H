/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

Namespace
    Foam::functionObjects

Description
    Namespace for functionObjects.

    OpenFOAM includes a collection of functionObjects selected by the user at
    run-time to manipulate the simulation and provide mechanisms to extract
    field and derived quantities.  Alternatively, the same actions can be
    executed after the simulation using the \c -postProcess command-line option.

    functionObjects are selected by entries in the $FOAM_CASE/system/functions
    dictionary e.g. to select the \c functionObjectType functionObject the
    following entry would be specified:

    \verbatim
    <functionObjectName>
    {
        type                functionObjectType;
        libs                ("libMyFunctionObjectlib.so");
        region              defaultRegion;
        enabled             yes;
        startTime           0;
        endTime             10;
        writeControl        writeTime;
        writeInterval       1;
        ...
    }
    \endverbatim

    Where:
    \table
        Property          | Description              | Required | Default value
        type              | Type of functionObject         | yes |
        libs              | Shared libraries               | no  |
        region            | Name of region                 | no  |
        enabled           | On/off switch                  | no  | yes
        log               | Print data to log              | no  | no
        startTime         | Start time                     | no  |
        endTime           | End time                       | no  |
        executeAtStart    | Execute at start time switch   | no  | yes
        executeControl    | See time controls below        | no  | timeStep
        executeInterval   | Steps between each execution   | no  |
        executeTimes      | List of execution times        | no  |
        writeControl      | See time controls below        | no  | timeStep
        writeInterval     | Steps between each write       | no  |
        writeTimes        | List of write times            | no  |
        writeFrequencies  | List of write frequencies      | no  |
    \endtable

    Time controls:
    \table
        Option            | Description
        timeStep          | Execute/write every 'Interval' time-steps
        writeTime         | Execute/write every 'Interval' write times
        adjustableRunTime | Execute/write every 'Interval' run time period
        runTime           | Execute/write every 'Interval' run time period
        runTimes          | Execute/write at specified list of run times
        clockTime         | Execute/write every 'Interval' clock time period
        cpuTime           | Execute/write every 'Interval' CPU time period
        none              | No execution
    \endtable

    The sub-dictionary name \c \<functionObjectName\> is chosen by the user, and
    is typically used as the name of the output directory for any data written
    by the functionObject.  The \c type entry defines the type of function
    object properties that follow.  FunctionObjects are packaged into separate
    libraries and the \c libs entry is used to specify which library should be
    loaded.

    Each functionObject has two separate run phases:

      - The \c execute phase is meant to be used for updating calculations
        or for management tasks.
      - The \c write phase is meant for writing the calculated data to disk.

    For each phase the respective time controls are provided, as listed above.

Class
    Foam::functionObject

Description
    Abstract base-class for Time/database functionObjects.

See also
    Foam::functionObjectList
    Foam::functionObjects::timeControl

SourceFiles
    functionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObject_H
#define functionObject_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;
class polyMesh;
class polyTopoChangeMap;
class polyMeshMap;
class polyDistributionMap;

/*---------------------------------------------------------------------------*\
                       Class functionObject Declaration
\*---------------------------------------------------------------------------*/

class functionObject
{
    // Private Data

        //- Name
        const word name_;


protected:

        //- Reference to time
        const Time& time_;


public:

    ClassName("functionObject");

    //- Runtime type information
    virtual const word& type() const = 0;

    //- Global post-processing mode switch
    static bool postProcess;

    //- Switch write log to Info
    Switch log;

    //- Switch write log to Info
    Switch executeAtStart_;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            functionObject,
            dictionary,
            (const word& name, const Time& runTime, const dictionary& dict),
            (name, runTime, dict)
        );


    // Constructors

        //- Construct from components
        functionObject(const word& name, const Time& runTime);

        //- Return clone
        autoPtr<functionObject> clone() const
        {
            NotImplemented;
            return autoPtr<functionObject>(nullptr);
        }

        //- Disallow default bitwise copy construction
        functionObject(const functionObject&) = delete;


    // Selectors

        //- Select from dictionary, based on its "type" entry
        static autoPtr<functionObject> New
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    //- Destructor
    virtual ~functionObject();


    // Member Functions

        //- Return the name of this functionObject
        const word& name() const;

        //- Read and set the functionObject if its data have changed
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const = 0;

        //- Return true if the functionObject should be executed at the start
        virtual bool executeAtStart() const;

        //- Called at each ++ or += of the time-loop.
        //  foamPostProcess overrides the usual executeControl behaviour and
        //  forces execution (used in post-processing mode)
        virtual bool execute() = 0;

        //- Called at each ++ or += of the time-loop.
        //  foamPostProcess overrides the usual writeControl behaviour and
        //  forces writing always (used in post-processing mode)
        virtual bool write() = 0;

        //- Called when Time::run() determines that the time-loop exits.
        virtual bool end();

        //- Called by Time::adjustTimeStep(). Allows the functionObject to
        //  insert a write time earlier than that already in use by the run
        //  time. Returns the write time, or vGreat.
        virtual scalar timeToNextAction();

        //- Return the maximum time-step for stable operation
        virtual scalar maxDeltaT() const;

        //- Update topology using the given map
        virtual void movePoints(const polyMesh& mesh);

        //- Update topology using the given map
        virtual void topoChange(const polyTopoChangeMap& map);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const functionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
