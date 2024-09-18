/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2024 OpenFOAM Foundation
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
    Foam::functionObjects::objectRegistryFunctionObject

Description
    Specialisation of Foam::functionObject which holds a reference to an object
    registry. Provides a number of functions for cacheing and accessing objects
    from the registry.

See also
    Foam::functionObject

SourceFiles
    objectRegistryFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_objectRegistryFunctionObject_H
#define functionObjects_objectRegistryFunctionObject_H

#include "functionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                  Class objectRegistryFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class objectRegistryFunctionObject
:
    public functionObject
{

protected:

    // Protected member data

        //- Reference to the objectRegistry
        const objectRegistry& obr_;


    // Protected member functions

        //- Find field in the objectRegistry
        template<class ObjectType>
        bool foundObject(const word& fieldName) const;

        //- Prints a warning message that fieldName cannot be found
        template<class ObjectType>
        void cannotFindObject(const word& fieldName);

        //- Prints a warning message that fieldName cannot be found
        void cannotFindObject(const word& fieldName);

        //- Prints a warning message that fieldNames cannot be found
        void cannotFindObjects(const wordList& fieldNames);

        //- Lookup object from the objectRegistry
        template<class ObjectType>
        const ObjectType& lookupObject(const word& fieldName) const;

        //- Lookup non-const object reference from the objectRegistry
        template<class ObjectType>
        ObjectType& lookupObjectRef(const word& fieldName);

        //- Store the given field in the objectRegistry
        template<class ObjectType>
        bool store(const tmp<ObjectType>& tfield);

        //- Store the given field in the objectRegistry under the given name
        template<class ObjectType>
        bool store
        (
            const word& fieldName,
            const tmp<ObjectType>& tfield,
            bool cacheable = false
        );

        //- Write field if present in objectRegistry
        bool writeObject(const word& fieldName);

        //- Clear field from the objectRegistry if present
        bool clearObject(const word& fieldName);


public:

    //- Runtime type information
    TypeName("objectRegistryFunctionObject");


    // Constructors

        //- Construct from an objectRegistry and dictionary
        objectRegistryFunctionObject
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );

        //- Construct from an objectRegistry
        objectRegistryFunctionObject
        (
            const word& name,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construction
        objectRegistryFunctionObject
        (
            const objectRegistryFunctionObject&
        ) = delete;


    //- Destructor
    virtual ~objectRegistryFunctionObject();


    // Member Functions

        //- Read optional controls
        virtual bool read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const objectRegistryFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "objectRegistryFunctionObjectTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
