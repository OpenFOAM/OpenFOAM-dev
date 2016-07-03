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
    Foam::functionObjects::writeObjects

Group
    grpUtilitiesFunctionObjects

Description
    Allows specification of different writing frequency of objects registered to
    the database.

    It has similar functionality as the main time database through the
    writeControl setting:
      - timeStep
      - writeTime
      - adjustableRunTime
      - runTime
      - clockTime
      - cpuTime

    Example of function object specification:
    \verbatim
    writeObjects1
    {
        type        writeObjects;
        libs        ("libutilityFunctionObjects.so");
        exclusiveWriting     true;
        ...
        objects     (obj1 obj2);
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: writeObjects | yes |
        objects      | objects to write        | yes         |
        exclusiveWriting    | Takes over object writing | no | yes
    \endtable

    \c exclusiveWriting disables automatic writing (i.e through database) of the
    objects to avoid duplicate writing.

See also
    Foam::functionObject
    Foam::functionObjects::timeControl

SourceFiles
    writeObjects.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_writeObjects_H
#define functionObjects_writeObjects_H

#include "functionObject.H"
#include "wordReList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class writeObjects Declaration
\*---------------------------------------------------------------------------*/

class writeObjects
:
    public functionObject
{
    // Private data

        //- Refererence to Db
        const objectRegistry& obr_;

        //- Takes over the writing from Db
        bool exclusiveWriting_;

        //- Names of objects to control
        wordReList objectNames_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        writeObjects(const writeObjects&);

        //- Disallow default bitwise assignment
        void operator=(const writeObjects&);


public:

    //- Runtime type information
    TypeName("writeObjects");


    // Constructors

        //- Construct from Time and dictionary
        writeObjects
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~writeObjects();


    // Member Functions

        //- Read the writeObjects data
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute();

        //- Write the registered objects
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
