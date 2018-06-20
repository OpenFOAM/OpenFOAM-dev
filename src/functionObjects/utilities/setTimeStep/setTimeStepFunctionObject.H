/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Foam::functionObjects::setTimeStepFunctionObject

Description
    Overrides the timeStep. Can only be used with
    solvers with adjustTimeStep control (e.g. pimpleFoam). Makes no attempt
    to cooperate with other timeStep 'controllers' (maxCo, other
    functionObjects). Supports 'enabled' flag but none of the other ones
    'timeStart', 'timeEnd', 'writeControl' etc.

SourceFiles
    setTimeStepFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_setTimeStepFunctionObject_H
#define functionObjects_setTimeStepFunctionObject_H

#include "functionObject.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class setTimeStepFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class setTimeStepFunctionObject
:
    public functionObject
{
    // Private data

        //- Reference to the time database
        const Time& time_;

        //- Time step function/table
        autoPtr<Function1<scalar>> timeStepPtr_;


    // Private member functions

        //- Disallow default bitwise copy construct
        setTimeStepFunctionObject(const setTimeStepFunctionObject&);

        //- Disallow default bitwise assignment
        void operator=(const setTimeStepFunctionObject&);


public:

    //- Runtime type information
    TypeName("setTimeStep");


    // Constructors

        //- Construct from components
        setTimeStepFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    // Destructor
    virtual ~setTimeStepFunctionObject();


    // Member Functions

        //- Return time database
        const Time& time() const;

        //- Override the time-step value
        virtual bool setTimeStep();

        //- Read and set the function object if its data have changed
        virtual bool read(const dictionary&);

        //- Called at each ++ or += of the time-loop.
        //  postProcess overrides the usual executeControl behaviour and
        //  forces execution (used in post-processing mode)
        virtual bool execute();

        //- Called at each ++ or += of the time-loop.
        //  postProcess overrides the usual writeControl behaviour and
        //  forces writing always (used in post-processing mode)
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
