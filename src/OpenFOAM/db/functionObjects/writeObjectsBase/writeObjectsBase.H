/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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
    Foam::functionObjects::writeObjectsBase

Description
    FunctionObject base class for writing a list of objects registered to the
    database, on behalf of the inheriting function object, on when those should
    be written to disk.

    FunctionObjects that inherit this class will receive the additional
    dictionary option \c objects which allows selecting which fields of the
    inherited function should be written to disk when \c write() is called.

    Example of function object specification:
    \verbatim
    <functionObjectName>
    {
        ...
        objects     (obj1 obj2);
        ...
    }
    \endverbatim

Usage
    \table
        Property | Description                       | Required | Default value
        objects  | List of objects to be written     | yes      |
    \endtable

    Note: Regular expressions can also be used in \c objects.

See also
    Foam::functionObject
    Foam::functionObjects::writeObjects
    Foam::functionObjects::writeLocalObjects

SourceFiles
    writeObjectsBase.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_writeObjectsBase_H
#define functionObjects_writeObjectsBase_H

#include "wordList.H"
#include "wordReList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class regIOobject;
class Switch;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class writeObjectsBase Declaration
\*---------------------------------------------------------------------------*/

class writeObjectsBase
{

protected:

    // Protected data

        //- Reference to the region objectRegistry
        const objectRegistry& writeObr_;

        //- Reference to the inheriting function object's log variable
        const Switch& log_;

        //- Object names requested by the user to be written
        wordReList writeObjectNames_;


    // Protected Member Functions

        //- Reset the list of object names to be written to a single regular
        //  expression
        void resetWriteObjectName(const wordRe& name);

        //- Reset the list of object names to be written
        void resetWriteObjectNames(const wordReList& names);

        //- Get the list of field names to be written
        virtual wordList objectNames();

        //- Write the requested registered IO object
        virtual void writeObject(const regIOobject& obj);


public:

    // Constructors

        //- Construct from objectRegistry and inheriting function object
        writeObjectsBase
        (
            const objectRegistry& obr,
            const Switch& logRef
        );

        //- Disallow default bitwise copy construction
        writeObjectsBase(const writeObjectsBase&) = delete;


    //- Destructor
    virtual ~writeObjectsBase();


    // Member Functions

        //- Return const access to the object names requested to be written
        const wordReList& writeObjectNames() const;

        //- Read the list of objects to be written
        virtual bool read(const dictionary&);

        //- Write function
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const writeObjectsBase&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
