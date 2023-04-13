/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::functionObjects::stopAtTimeStep

Description
    Stops the run if the time-step drops below the specified value in seconds
    and optionally write results before stopping.

    The following actions are supported:
    - noWriteNow
    - writeNow (default)
    - nextWrite

    Examples of function object specification:
    \verbatim
    stop
    {
        type            stopAtTimeStep;

        libs            ("libutilityFunctionObjects.so");

        minDeltaT       1e-8;
        action          writeNow;
    }
    \endverbatim
    will write the fields and stop if the time-step drops below 1e-8s.

Usage
    \table
        Property  | Description               | Required | Default value
        type      | type name: stopAtTimeStep | yes      |
        minDeltaT | Minimum time-step [s]     | yes      |
        action    | Action executed           | no       | writeNow
    \endtable

SourceFiles
    stopAtTimeStep.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_stopAtTimeStep_H
#define functionObjects_stopAtTimeStep_H

#include "stopAt.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class stopAtTimeStep Declaration
\*---------------------------------------------------------------------------*/

class stopAtTimeStep
:
    public stopAt
{
    // Private Data

        //- Minimum time-step, below which the run is stopped
        scalar minDeltaT_;


    // Private Member Functions

        //- Return true when the stop condition is achieved
        virtual bool condition() const;


public:

    //- Runtime type information
    TypeName("stopAtTimeStep");


    // Constructors

        //- Construct from Time and dictionary
        stopAtTimeStep
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        stopAtTimeStep(const stopAtTimeStep&) = delete;


    //- Destructor
    virtual ~stopAtTimeStep();


    // Member Functions

        //- Read the dictionary settings
        virtual bool read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const stopAtTimeStep&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
