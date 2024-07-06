/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::functionObjects::stopAtFile

Description
    Stops the run when the specified file is created in the case directory
    and optionally write results before stopping.

    The default name of the trigger file is \c $FOAM_CASE/\<name\> where \c
    \<name\> is the name of the functionObject entry.

    Currently the following action types are supported:
    - noWriteNow
    - writeNow
    - nextWrite (default)

    Examples of function object specification:
    \verbatim
    stop
    {
        type            stopAtFile;

        libs            ("libutilityFunctionObjects.so");
    }
    \endverbatim
    will stop the run at the next write after the file "stop" is created in the
    case directory.

    \verbatim
    stop
    {
        type    stopAtFile;
        libs    ("libutilityFunctionObjects.so");
        file    "$FOAM_CASE/stop";
        action  writeNow;
    }
    \endverbatim
    will write the fields and stop the run when the file "stop" is created in
    the case directory.

Usage
    \table
        Property | Description            | Required | Default value
        type     | type name: stopAtFile  | yes      |
        file     | Trigger file path name | no       | $FOAM_CASE/\<name\>
        action   | Action executed        | no       | nextWrite
    \endtable

SourceFiles
    stopAtFile.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_stopAtFile_H
#define functionObjects_stopAtFile_H

#include "stopAt.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class stopAtFile Declaration
\*---------------------------------------------------------------------------*/

class stopAtFile
:
    public stopAt
{
    // Private Data

        //- The fully-qualified name of the stopAtFile file
        fileName stopAtFileFile_;


    // Private Member Functions

        //- Remove stopAtFile file.
        void removeFile() const;

        //- Return true when the stop condition is achieved
        virtual bool condition() const;


public:

    //- Runtime type information
    TypeName("stopAtFile");


    // Constructors

        //- Construct from Time and dictionary
        stopAtFile
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        stopAtFile(const stopAtFile&) = delete;


    //- Destructor
    virtual ~stopAtFile();


    // Member Functions

        //- Read the dictionary settings
        virtual bool read(const dictionary&);

        //- Execute at the final time-loop, used for cleanup
        virtual bool end();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const stopAtFile&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
