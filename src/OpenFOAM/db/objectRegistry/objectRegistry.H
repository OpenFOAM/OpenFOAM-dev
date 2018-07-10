/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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
    Foam::objectRegistry

Description
    Registry of regIOobjects

SourceFiles
    objectRegistry.C

\*---------------------------------------------------------------------------*/

#ifndef objectRegistry_H
#define objectRegistry_H

#include "HashTable.H"
#include "regIOobject.H"
#include "wordReList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class objectRegistry Declaration
\*---------------------------------------------------------------------------*/

class objectRegistry
:
    public regIOobject,
    public HashTable<regIOobject*>
{
    // Private Data

        //- Master time objectRegistry
        const Time& time_;

        //- Parent objectRegistry
        const objectRegistry& parent_;

        //- Local directory path of this objectRegistry relative to time
        fileName dbDir_;

        //- Current event
        mutable label event_;


    // Private Member Functions

        //- Is the objectRegistry parent_ different from time_
        //  Used to terminate searching within the ancestors
        bool parentNotTime() const;

        //- Disallow Copy constructor
        objectRegistry(const objectRegistry&);

        //- Disallow default bitwise copy construct and assignment
        void operator=(const objectRegistry&);


public:

    //- Declare type name for this IOobject
    TypeName("objectRegistry");


    // Constructors

        //- Construct the time objectRegistry given an initial estimate
        //  for the number of entries
        explicit objectRegistry
        (
            const Time& db,
            const label nIoObjects = 128
        );

        //- Construct a sub-registry given an IObject to describe the registry
        //  and an initial estimate for the number of entries
        explicit objectRegistry
        (
            const IOobject& io,
            const label nIoObjects = 128
        );


    //- Destructor
    virtual ~objectRegistry();


    // Member functions

        // Access

            //- Return time
            const Time& time() const
            {
                return time_;
            }

            //- Return the parent objectRegistry
            const objectRegistry& parent() const
            {
                return parent_;
            }

            //- Local directory path of this objectRegistry relative to the time
            virtual const fileName& dbDir() const
            {
                return dbDir_;
            }

            //- Return the list of names of the IOobjects
            wordList names() const;

            //- Return the sorted list of names of the IOobjects
            wordList sortedNames() const;

            //- Return the list of names of IOobjects of given class name
            wordList names(const word& className) const;

            //- Return the sorted list of names of IOobjects of given class name
            wordList sortedNames(const word& className) const;

            //- Return the list of names of the IOobjects of given type
            template<class Type>
            wordList names() const;

            //- Return the list of objects whose name matches the input regExp
            template<class Type>
            wordList names(const wordRe& name) const;

            //- Return the list of objects whose name matches the input regExp
            template<class Type>
            wordList names(const wordReList& name) const;

            //- Lookup and return a const sub-objectRegistry. Optionally create
            //  it if it does not exist.
            const objectRegistry& subRegistry
            (
                const word& name,
                const bool forceCreate = false
            ) const;

            //- Lookup and return all objects of the given Type
            template<class Type>
            HashTable<const Type*> lookupClass(const bool strict = false) const;

            //- Lookup and return all objects of the given Type
            template<class Type>
            HashTable<Type*> lookupClass(const bool strict = false);

            //- Is the named Type found?
            template<class Type>
            bool foundObject(const word& name) const;

            //- Lookup and return the object of the given Type
            template<class Type>
            const Type& lookupObject(const word& name) const;

            //- Lookup and return the object reference of the given Type
            template<class Type>
            Type& lookupObjectRef(const word& name) const;

            //- Return new event number.
            label getEvent() const;


        // Edit

            //- Rename
            virtual void rename(const word& newName);

            //- Add an regIOobject to registry
            bool checkIn(regIOobject&) const;

            //- Remove an regIOobject from registry
            bool checkOut(regIOobject&) const;

            //- Remove all regIOobject owned by the registry
            void clear();


        // Reading

            //- Return true if any of the object's files have been modified
            virtual bool modified() const;

            //- Read the objects that have been modified
            void readModifiedObjects();

            //- Read object if modified
            virtual bool readIfModified();


        // Writing

            //- writeData function required by regIOobject but not used
            //  for this class, write is used instead
            virtual bool writeData(Ostream&) const
            {
                NotImplemented;

                return false;
            }

            //- Write the objects
            virtual bool writeObject
            (
                IOstream::streamFormat fmt,
                IOstream::versionNumber ver,
                IOstream::compressionType cmp,
                const bool valid
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "objectRegistryTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
