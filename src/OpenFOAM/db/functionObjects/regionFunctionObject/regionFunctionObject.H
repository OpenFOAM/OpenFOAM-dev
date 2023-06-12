/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2023 OpenFOAM Foundation
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
    Foam::functionObjects::regionFunctionObject

Description
    Specialisation of Foam::functionObject for a region and providing a
    reference to the region Foam::objectRegistry.

See also
    Foam::functionObject

SourceFiles
    regionFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_regionFunctionObject_H
#define functionObjects_regionFunctionObject_H

#include "functionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class regionFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class regionFunctionObject
:
    public functionObject
{

protected:

    // Protected member data

        //- Reference to the region objectRegistry
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
    TypeName("regionFunctionObject");


    // Constructors

        //- Construct from Time and dictionary.
        //  The region objectRegistry is looked-up runTime with the name
        //  looked-up from the dictionary (defaults to polyMesh::defaultRegion)
        regionFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from the region objectRegistry and dictionary
        regionFunctionObject
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        regionFunctionObject(const regionFunctionObject&) = delete;


    //- Destructor
    virtual ~regionFunctionObject();


    // Member Functions

        //- Read optional controls
        virtual bool read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const regionFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "regionFunctionObjectTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
