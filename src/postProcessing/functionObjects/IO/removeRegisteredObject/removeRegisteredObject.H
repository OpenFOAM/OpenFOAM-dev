/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::removeRegisteredObject

Group
    grpIOFunctionObjects

Description
    This function object removes registered objects if present in the database

    Example of function object specification:
    \verbatim
    removeRegisteredObject1
    {
        type        removeRegisteredObject;
        functionObjectLibs ("libIOFunctionObjects.so");
        ...
        objectNames (obj1 obj2);
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: removeRegisteredObject | yes |
        objectNames  | objects to remove       | yes         |
    \endtable

SeeAlso
    Foam::functionObject
    Foam::OutputFilterFunctionObject

SourceFiles
    removeRegisteredObject.C
    IOremoveRegisteredObject.H

\*---------------------------------------------------------------------------*/

#ifndef removeRegisteredObject_H
#define removeRegisteredObject_H

#include "wordList.H"
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
                   Class removeRegisteredObject Declaration
\*---------------------------------------------------------------------------*/

class removeRegisteredObject
{
protected:

    // Private data

        //- Name of this set of removeRegisteredObject
        word name_;

        const objectRegistry& obr_;

        // Read from dictionary

            //- Names of objects to control
            wordList objectNames_;


    // Private Member Functions


        //- Disallow default bitwise copy construct
        removeRegisteredObject(const removeRegisteredObject&);

        //- Disallow default bitwise assignment
        void operator=(const removeRegisteredObject&);


public:

    //- Runtime type information
    TypeName("removeRegisteredObject");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        removeRegisteredObject
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~removeRegisteredObject();


    // Member Functions

        //- Return name of the removeRegisteredObject
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the removeRegisteredObject data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Write the removeRegisteredObject
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
