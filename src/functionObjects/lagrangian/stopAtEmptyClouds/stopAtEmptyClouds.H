/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::functionObjects::stopAtEmptyClouds

Description
    Stops the run when all parcel clouds are empty.

    The following actions are supported:
    - noWriteNow
    - writeNow
    - nextWrite (default)

    Example of function object specification:
    \verbatim
    stop
    {
        type            stopAtEmptyClouds;
        libs            ("liblagrangianFunctionObjects.so");

        executeControl  timeStep; // <-- Likely only to be needed if injection
        startTime       0.500001; // starts after time zero. In which case the
                                  // startTime should be slightly after the
                                  // injection start time.

        action          nextWrite;
    }
    \endverbatim

Usage
    \table
        Property | Description                  | Required | Default value
        type     | type name: stopAtEmptyClouds | yes      |
        action   | Action executed              | no       | nextWrite
    \endtable

SourceFiles
    stopAtEmptyClouds.C

\*---------------------------------------------------------------------------*/

#ifndef stopAtEmptyClouds_H
#define stopAtEmptyClouds_H

#include "stopAt.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class stopAtEmptyClouds Declaration
\*---------------------------------------------------------------------------*/

class stopAtEmptyClouds
:
    public stopAt
{
    // Private Member Functions

        //- Return true when the stop condition is achieved
        virtual bool condition() const;


public:

    //- Runtime type information
    TypeName("stopAtEmptyClouds");


    // Constructors

        //- Construct from Time and dictionary
        stopAtEmptyClouds
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~stopAtEmptyClouds();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
