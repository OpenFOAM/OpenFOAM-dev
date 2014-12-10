/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
#include "outputFilterOutputControl.H"

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

            //- Name of region - defaults to name of polyMesh::defaultRegion
            word regionName_;

            //- Dictionary name to supply required inputs
            word dictName_;

            //- Switch for the execution - defaults to 'yes/on'
            bool enabled_;

            //- Switch to store filter in between writes or use on-the-fly
            //  construction - defaults to true
            bool storeFilter_;

            //- Activation time - defaults to -VGREAT
            scalar timeStart_;

            //- De-activation time - defaults to VGREAT
            scalar timeEnd_;

            //- Number of steps before the dumping time in which the deltaT
            // will start to change (valid for ocAdjustableTime)
            label nStepsToStartTimeChange_;


        //- Output controls
        outputFilterOutputControl outputControl_;

        //- Evaluate controls
        outputFilterOutputControl evaluateControl_;

        //- Pointer to the output filter
        autoPtr<OutputFilter> ptr_;


    // Private Member Functions

        //- Read relevant dictionary entries
        void readDict();

        //- Creates most of the data associated with this object.
        void allocateFilter();

        //- Destroys most of the data associated with this object.
        void destroyFilter();

        //- Returns true if active (enabled and within time bounds)
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
            virtual const Time& time() const
            {
                return time_;
            }

            //- Return the input dictionary
            virtual const dictionary& dict() const
            {
                return dict_;
            }

            //- Return the region name
            virtual const word& regionName() const
            {
                return regionName_;
            }

            //- Return the optional dictionary name
            virtual const word& dictName() const
            {
                return dictName_;
            }

            //- Return the enabled flag
            virtual bool enabled() const
            {
                return enabled_;
            }

            //- Return the output control object
            virtual const outputFilterOutputControl& outputControl() const
            {
                return outputControl_;
            }

            //- Return the output filter
            virtual const OutputFilter& outputFilter() const
            {
                return ptr_();
            }


        // Function object control

            //- Switch the function object on
            virtual void on();

            //- Switch the function object off
            virtual void off();


            //- Called at the start of the time-loop
            virtual bool start();

            //- Called at each ++ or += of the time-loop
            virtual bool execute(const bool forceWrite);

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

#ifdef NoRepository
#   include "OutputFilterFunctionObject.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
