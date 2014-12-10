/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::writeRegisteredObject

Group
    grpIOFunctionObjects

Description
    This function object allows specification of different writing frequency
    of objects registered to the database. It has similar functionality
    as the main time database through the outputControl setting:
        timeStep
        outputTime
        adjustableTime
        runTime
        clockTime
        cpuTime

    Example of function object specification:
    \verbatim
    writeRegisteredObject1
    {
        type        writeRegisteredObject;
        functionObjectLibs ("libIOFunctionObjects.so");
        exclusiveWriting     true;
        ...
        objectNames (obj1 obj2);
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: writeRegisteredObject | yes |
        objectNames  | objects to write        | yes         |
        exclusiveWriting    | Takes over object writing | no | yes
    \endtable

    exclusiveWriting disables automatic writing (i.e through database) of the
    objects to avoid duplicate writing.

SeeAlso
    Foam::functionObject
    Foam::OutputFilterFunctionObject

SourceFiles
    writeRegisteredObject.C
    IOwriteRegisteredObject.H

\*---------------------------------------------------------------------------*/

#ifndef writeRegisteredObject_H
#define writeRegisteredObject_H

#include "wordReList.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                   Class writeRegisteredObject Declaration
\*---------------------------------------------------------------------------*/

class writeRegisteredObject
{
protected:

    // Private data

        //- Name of this set of writeRegisteredObject
        word name_;

        //- Takes over the writing from Db
        bool exclusiveWriting_;

        //- Refererence to Db
        const objectRegistry& obr_;

        // Read from dictionary

            //- Names of objects to control
            wordReList objectNames_;


    // Private Member Functions


        //- Disallow default bitwise copy construct
        writeRegisteredObject(const writeRegisteredObject&);

        //- Disallow default bitwise assignment
        void operator=(const writeRegisteredObject&);


public:

    //- Runtime type information
    TypeName("writeRegisteredObject");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        writeRegisteredObject
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~writeRegisteredObject();


    // Member Functions

        //- Return name of the writeRegisteredObject
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the writeRegisteredObject data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Write the writeRegisteredObject
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
