/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

    \linebreak
    \subsection secFunctionObjects Using function objects

    FunctionObjects are selected by additional entries in the
    $FOAM_CASE/system/controlDict dictionary.  Each object is listed in the \c
    functions sub-dictionary, e.g. to select the \c functionObjectType
    functionObject the following entry would be specified:

    \verbatim
        functions
        {
            myFunctionObject    // Name of functionObject entry
            {
                type                functionObjectType;
                functionObjectLibs  ("libMyFunctionObjectlib.so");
                region              defaultRegion;
                enabled             yes;
                timeStart           0;
                timeEnd             10;
                writeControl        outputTime;
                writeInterval       1;
                ...
            }
        }
    \endverbatim

    Where:
    \table
        Property | Description                       | Required | Default value
        type     | Type of function object               | yes      |
        functionObjectLibs | Libraries containing implementation | yes |
        region   | Name of region for multi-region cases | no       |
        enabled  | On/off switch                         | no       | yes
        timeStart| Start time                            | no       |
        timeEnd  | End time                              | no       |
        evaluateControl | See time controls below        | no       | timeStep
        evaluateInterval   | Steps between output           | no       |
        writeControl   | See time controls below        | no       | timeStep
        writeInterval   | Steps between output           | no       |
    \endtable

    Time controls:
    \table
        Option            | Description
        timeStep          | Execute/write every 'writeInterval' time-steps
        outputTime        | Execute/write every 'writeInterval' output times
        adjustableRunTime | Execute/write every 'writeInterval' run time period
        runTime           | Execute/write every 'writeInterval' run time period
        clockTime       | Execute/write every 'writeInterval' clock time period
        cpuTime           | Execute/write every 'writeInterval' CPU time period
        none              | Execute/write every time-step
    \endtable

    The sub-dictionary name \c myFunctionObject is chosen by the user, and is
    typically used as the name of the output directory for any data written by
    the functionObject.  The \c type entry defines the type of function object
    properties that follow.  FunctionObjects are packaged into separate
    libraries and the \c libs entry is used to specify which library
    should be loaded.

Class
    Foam::functionObject

Description
    Abstract base-class for Time/database function objects.

See Also
    Foam::OutputFilterFunctionObject

SourceFiles
    functionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObject_H
#define functionObject_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                       Class functionObject Declaration
\*---------------------------------------------------------------------------*/

class functionObject
{
    // Private data

        //- Name
        const word name_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        functionObject(const functionObject&);

        //- Disallow default bitwise assignment
        void operator=(const functionObject&);


public:

    //- Runtime type information
    virtual const word& type() const = 0;

    static int debug;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            functionObject,
            dictionary,
            (const word& name, const Time& t, const dictionary& dict),
            (name, t, dict)
        );


    // Constructors

        //- Construct from components
        functionObject(const word& name);

        //- Return clone
        autoPtr<functionObject> clone() const
        {
            NotImplemented;
            return autoPtr<functionObject>(NULL);
        }


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

        //- Name
        virtual const word& name() const;

        //- Called at each ++ or += of the time-loop. forceWrite overrides
        //  the usual writeControl behaviour and forces writing always
        //  (used in post-processing mode)
        virtual bool execute(const bool forceWrite) = 0;

        //- Called when Time::run() determines that the time-loop exits.
        //  By default it simply calls execute().
        virtual bool end();

        //- Called when time was set at the end of the Time::operator++
        virtual bool timeSet();

        //- Called at the end of Time::adjustDeltaT() if adjustTime is true
        virtual bool adjustTimeStep();

        //- Read and set the function object if its data have changed
        virtual bool read(const dictionary&) = 0;

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh& mpm) = 0;

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh& mesh) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
