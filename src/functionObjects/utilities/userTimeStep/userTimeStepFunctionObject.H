/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::functionObjects::userTimeStep

Description
    Writes the time step

    Example of function object specification:
    \verbatim
    userTimeStep
    {
        type            userTimeStep;

        libs            ("libutilityFunctionObjects.so");

        writeControl    userTimeStep;
        writeInterval   1;
    }
    \endverbatim

    Output data is written to the dir postProcessing/userTimeStep/\<timeDir\>/

See also
    Foam::functionObject
    Foam::functionObjects::regionFunctionObject
    Foam::functionObjects::logFiles

SourceFiles
    userTimeStep.C

\*---------------------------------------------------------------------------*/

#ifndef userTimeStepFunctionObject_H
#define userTimeStepFunctionObject_H

#include "regionFunctionObject.H"
#include "logFiles.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class userTimeStep Declaration
\*---------------------------------------------------------------------------*/

class userTimeStep
:
    public regionFunctionObject,
    public logFiles
{

protected:

    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(const label i);


public:

    //- Runtime type information
    TypeName("userTimeStep");


    // Constructors

        //- Construct from Time and dictionary
        userTimeStep
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        userTimeStep(const userTimeStep&) = delete;


    //- Destructor
    virtual ~userTimeStep();


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the time step value
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const userTimeStep&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
