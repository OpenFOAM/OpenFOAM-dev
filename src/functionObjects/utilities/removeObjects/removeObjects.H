/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2024 OpenFOAM Foundation
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
    Foam::functionObjects::removeObjects

Description
    Removes registered objects if present in the database.

    Example of function object specification:
    \verbatim
    removeObjects1
    {
        type        removeObjects;

        libs        ("libutilityFunctionObjects.so");

        objects     (obj1 obj2);
    }
    \endverbatim

Usage
    \table
        Property     | Description              | Required    | Default value
        type         | type name: removeObjects | yes         |
        objects      | objects to remove        | yes         |
    \endtable

See also
    Foam::functionObject

SourceFiles
    removeObjects.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_removeObjects_H
#define functionObjects_removeObjects_H

#include "functionObject.H"
#include "wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class removeObjects Declaration
\*---------------------------------------------------------------------------*/

class removeObjects
:
    public functionObject
{
    // Private Data

        //- Reference to the objectRegistry
        const objectRegistry& obr_;

        //- Names of objects to control
        wordList objectNames_;


public:

    //- Runtime type information
    TypeName("removeObjects");


    // Constructors

        //- Construct from Time and dictionary
        removeObjects
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        removeObjects(const removeObjects&) = delete;


    //- Destructor
    virtual ~removeObjects();


    // Member Functions

        //- Read the removeObjects data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Remove the registered objects
        virtual bool execute();

        //- Do nothing
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const removeObjects&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
