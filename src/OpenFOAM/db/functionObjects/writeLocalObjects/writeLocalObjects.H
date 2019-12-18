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
    Foam::functionObjects::writeLocalObjects

Description
    FunctionObject base class for managing a list of objects on behalf of the
    inheriting function object, on when those should be written to disk.

    FunctionObjects that inherit this class will receive the additional
    dictionary option \c objects which allows selecting which fields of the
    inherited function should be written to disk when \c write() is called.
    When \c objects is omitted, it will write all objects and when that list is
    empty, it will not write any of the inheriting function object's managed
    objects.

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
        objects  | List of objects to be written     | no       | ".*"
    \endtable

    Note: Regular expressions can also be used in \c objects.

See also
    Foam::functionObject
    Foam::functionObjects::writeObjectsBase

SourceFiles
    writeLocalObjects.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_writeLocalObjects_H
#define functionObjects_writeLocalObjects_H

#include "wordList.H"
#include "writeObjectsBase.H"

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
                     Class writeLocalObjects Declaration
\*---------------------------------------------------------------------------*/

class writeLocalObjects
:
    public writeObjectsBase
{

protected:

    // Protected data

        //- Object names that are handled on behalf of the inheritor
        wordList localObjectNames_;


    // Protected Member Functions

        //- Reset the list of local object names from a single word
        void resetLocalObjectName(const word& name);

        //- Reset the list of local object names from a wordList
        void resetLocalObjectNames(const wordList& names);

        //- Get the list of field names to be written
        virtual wordList objectNames();


public:

    // Constructors

        //- Construct from objectRegistry and inheriting function object
        writeLocalObjects
        (
            const objectRegistry& obr,
            const Switch& logRef
        );

        //- Disallow default bitwise copy construction
        writeLocalObjects(const writeLocalObjects&) = delete;


    //- Destructor
    virtual ~writeLocalObjects();


    // Member Functions

        //- Return const access to the local object names
        const wordList& localObjectNames() const;

        //- Read the list of objects to be written
        virtual bool read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const writeLocalObjects&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
