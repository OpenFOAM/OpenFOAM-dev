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

Class
    Foam::OutputFilterFunctionObject

Description
    A functionObject wrapper around OutputFilter to allow them to be
    created via the functions entry within controlDict.

Note
    Since the timeIndex is used directly from Foam::Time, it is unaffected
    by user-time conversions. For example, Foam::engineTime might cause \a
    writeInterval to be degrees crank angle, but the functionObject
    execution \a interval would still be in timestep.

SourceFiles
    OutputFilterFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef OutputFilterFunctionObject_H
#define OutputFilterFunctionObject_H

#include "functionObject.H"
#include "dictionary.H"
#include "outputFilterControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class OutputFilterFunctionObject Declaration
\*---------------------------------------------------------------------------*/

template<class OutputFilter>
class OutputFilterFunctionObject
:
    public functionObject
{
    // Private data

        //- Reference to the time database
        const Time& time_;

        //- Input dictionary
        dictionary dict_;


        // Optional user inputs

            //- Activation time - defaults to -VGREAT
            scalar timeStart_;

            //- De-activation time - defaults to VGREAT
            scalar timeEnd_;

            //- Number of steps before the dump-time during which deltaT
            //  may be changed (valid for adjustableRunTime)
            label nStepsToStartTimeChange_;


        //- Output controls
        outputFilterControl writeControl_;

        //- Evaluate controls
        outputFilterControl evaluateControl_;

        //- The output filter
        OutputFilter filter_;


    // Private Member Functions

        //- Read relevant dictionary entries
        void readControls();

        //- Returns true if within time bounds
        bool active() const;

        //- Disallow default bitwise copy construct
        OutputFilterFunctionObject(const OutputFilterFunctionObject&);

        //- Disallow default bitwise assignment
        void operator=(const OutputFilterFunctionObject&);


public:

    //- Runtime type information
    TypeName(OutputFilter::typeName_());


    // Constructors

        //- Construct from components
        OutputFilterFunctionObject
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        // Access

            //- Return time database
            inline const Time& time() const;

            //- Return the input dictionary
            inline const dictionary& dict() const;

            //- Return the region name
            inline const word& regionName() const;

            //- Return the output control object
            inline const outputFilterControl& writeControl() const;

            //- Return the output filter
            inline const OutputFilter& outputFilter() const;


        // Function object control

            //- Called at each ++ or += of the time-loop.
            //  forceWrite overrides the usual writeControl behaviour and
            //  forces writing (used in post-processing mode)
            virtual bool execute(const bool forceWrite = false);

            //- Called when Time::run() determines that the time-loop exits
            virtual bool end();

            //- Called when time was set at the end of the Time::operator++
            virtual bool timeSet();

            //- Called at the end of Time::adjustDeltaT() if adjustTime is true
            virtual bool adjustTimeStep();

            //- Read and set the function object if its data have changed
            virtual bool read(const dictionary&);

            //- Update for changes of mesh
            virtual void updateMesh(const mapPolyMesh& mpm);

            //- Update for changes of mesh
            virtual void movePoints(const polyMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "OutputFilterFunctionObjectI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "OutputFilterFunctionObject.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
