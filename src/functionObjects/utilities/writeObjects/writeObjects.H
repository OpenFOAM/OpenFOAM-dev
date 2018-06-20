/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Description
    Allows specification of different writing frequency of objects registered
    to the database.

    It has similar functionality as the main time database through the
    \c writeControl setting:
      - timeStep
      - writeTime
      - adjustableRunTime
      - runTime
      - clockTime
      - cpuTime

    It also has the ability to write the selected objects that were defined
    with the respective write mode for the requested \c writeOption, namely:
    \vartable
        autoWrite | objects set to write at output time
        noWrite   | objects set to not write by default
        anyWrite  | any option of the previous two
    \endvartable

    Example of function object specification:
    \verbatim
    writeObjects1
    {
        type        writeObjects;
        libs        ("libutilityFunctionObjects.so");
        ...
        objects     (obj1 obj2);
        writeOption anyWrite;
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required     | Default value
        type         | type name: writeObjects | yes          |
        objects      | objects to write        | yes          |
        writeOption  | only those with this write option | no | anyWrite
    \endtable

    Note: Regular expressions can also be used in \c objects.

See also
    Foam::functionObject
    Foam::functionObjects::writeObjectsBase
    Foam::functionObjects::timeControl

SourceFiles
    writeObjects.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_writeObjects_H
#define functionObjects_writeObjects_H

#include "functionObject.H"
#include "writeObjectsBase.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class writeObjects Declaration
\*---------------------------------------------------------------------------*/

class writeObjects
:
    public functionObject,
    public writeObjectsBase
{
public:

    // Public data types

        //- Re-enumeration defining the write options, based on the original
        //  ones at IOobject::writeOption
        enum writeOption
        {
            AUTO_WRITE,
            NO_WRITE,
            ANY_WRITE
        };

        static const NamedEnum<writeOption, 3> writeOptionNames_;

private:

    // Private data

        //- To only write objects of defined writeOption
        writeOption writeOption_;


    // Private Member Functions

        //- Write the requested registered IO object
        virtual void writeObject(const regIOobject& obj);

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
