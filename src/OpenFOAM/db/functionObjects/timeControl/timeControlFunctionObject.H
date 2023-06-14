/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2023 OpenFOAM Foundation
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
    Foam::functionObjects::timeControl

Description
    General time control for functionObjects

SourceFiles
    timeControlFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef timeControlFunctionObject_H
#define timeControlFunctionObject_H

#include "functionObject.H"
#include "dictionary.H"
#include "timeControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class timeControl Declaration
\*---------------------------------------------------------------------------*/

class timeControl
:
    public functionObject
{
    // Private Data

        //- Reference to the time database
        const Time& time_;


        // Optional user inputs

            //- Activation time - defaults to -vGreat
            scalar startTime_;

            //- De-activation time - defaults to vGreat
            scalar endTime_;

            //- Number of steps before the dump-time during which deltaT
            //  may be changed (valid for adjustableRunTime)
            label nStepsToStartTimeChange_;


        //- Execute controls
        Foam::timeControl executeControl_;

        //- Write controls
        Foam::timeControl writeControl_;

        //- The functionObject to execute
        autoPtr<functionObject> foPtr_;


    // Private Member Functions

        //- Read relevant dictionary entries
        void readControls(const dictionary& dict);

        //- Returns true if within time bounds
        bool active() const;


public:

    //- Runtime type information
    TypeName("timeControl");


    // Constructors

        //- Construct from components
        timeControl
        (
            const word& name,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        timeControl(const timeControl&) = delete;


    // Member Functions

        // Access

            //- Return time database
            inline const Time& time() const;

            //- Return the execute control object
            inline const Foam::timeControl& executeControl() const;

            //- Return the write control object
            inline const Foam::timeControl& writeControl() const;

            //- Return the functionObject filter
            inline const functionObject& filter() const;


        // Function object control

            //- Return the list of fields required
            virtual wordList fields() const;

            //- Return true if the functionObject should be executed
            //  at the start
            virtual bool executeAtStart() const;

            //- Called at each ++ or += of the time-loop.
            //  foamPostProcess overrides the usual executeControl behaviour and
            //  forces execution (used in post-processing mode)
            virtual bool execute();

            //- Called at each ++ or += of the time-loop.
            //  foamPostProcess overrides the usual writeControl behaviour and
            //  forces writing (used in post-processing mode)
            virtual bool write();

            //- Called when Time::run() determines that the time-loop exits
            virtual bool end();

            //- Return the time to the next execute or write
            //  if time requires adjustment
            virtual scalar timeToNextAction();

            //- Read and set the function object if its data have changed
            virtual bool read(const dictionary&);

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
        void operator=(const timeControl&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "timeControlFunctionObjectI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
