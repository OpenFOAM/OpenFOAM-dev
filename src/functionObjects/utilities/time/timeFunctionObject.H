/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::functionObjects::time

Description
    Writes run time, CPU time and clock time
    and optionally the CPU and clock times per time step.

    Example of function object specification:
    \verbatim
    time
    {
        type            time;

        libs            ("libutilityFunctionObjects.so");

        writeControl    timeStep;
        writeInterval   1;

        perTimeStep     no;
    }
    \endverbatim

See also
    Foam::functionObject
    Foam::regionFunctionObject
    Foam::functionObjects::logFiles

SourceFiles
    time.C

\*---------------------------------------------------------------------------*/

#ifndef timeFunctionObject_H
#define timeFunctionObject_H

#include "regionFunctionObject.H"
#include "logFiles.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class time Declaration
\*---------------------------------------------------------------------------*/

class time
:
    public regionFunctionObject,
    public logFiles
{
    // Private member data

        //- Switch to write CPU and clock times per time-step
        Switch perTimeStep_;

        //- Previous time-step CPU time
        scalar cpuTime0_;

        //- Previous time-step clock time
        scalar clockTime0_;


    // Private member functions

        //- Disallow default bitwise copy construct
        time(const time&);

        //- Disallow default bitwise assignment
        void operator=(const time&);


protected:

    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(const label i);


public:

    //- Runtime type information
    TypeName("time");


    // Constructors

        //- Construct from Time and dictionary
        time
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~time();


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the time
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
