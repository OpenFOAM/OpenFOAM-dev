/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2021 OpenFOAM Foundation
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
    Updates the time step as a Function1 of time.

    Makes no attempt to cooperate with other timeStep controllers (i.e., solver
    courant number control or other functionObjects). In general any value set
    here will be overwritten by solver time-step adjustment if enabled, so for
    this to work 'adjustTime' should be switched off. This is compatible with
    'adjustableWriteTime', in which case the time-step values set will not be
    exactly as specified, but write intervals will be matched exactly.

    This function supports the 'enabled' flag but none of the other run
    controls; 'startTime', 'endTime', 'writeControl', etc...

SourceFiles
    setTimeStepFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_setTimeStepFunctionObject_H
#define functionObjects_setTimeStepFunctionObject_H

#include "functionObject.H"
#include "Function1.H"
#include "Time.H"

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
    // Private Data

        //- Reference to the time database
        const Time& time_;

        //- Time step function/table
        autoPtr<Function1<scalar>> timeStepPtr_;


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

        //- Disallow default bitwise copy construction
        setTimeStepFunctionObject(const setTimeStepFunctionObject&) = delete;


    // Destructor
    virtual ~setTimeStepFunctionObject();


    // Member Functions

        //-  Read and reset the timeStep Function1
        virtual bool read(const dictionary&);

        //- Reset the timeStep from the Function1 of time
        virtual bool execute();

        //- Do nothing
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const setTimeStepFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
