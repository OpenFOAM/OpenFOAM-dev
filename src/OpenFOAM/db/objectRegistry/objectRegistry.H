/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
#include "HashSet.H"
#include "Pair.H"

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

        //- Table of temporary object names with current state
        mutable HashTable<Pair<bool>> cacheTemporaryObjects_;

        //- State of cacheTemporaryObjects_, set true after reading
        mutable bool cacheTemporaryObjectsSet_;

        //- Accumulated list of temporary objects available to cache
        //  Used to provide diagnostics in case the requested object is not
        //  available
        mutable HashSet<word> temporaryObjects_;


    // Private Member Functions

        //- Is the objectRegistry parent_ different from time_
        //  Used to terminate searching within the ancestors
        bool parentNotTime() const;

        //- Read the cacheTemporaryObjects list from controlDict
        void readCacheTemporaryObjects() const;

        //- Delete the current cached object before caching a new object
        void deleteCachedObject(regIOobject& cachedOb) const;


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

        //- Disallow default bitwise copy construction
        objectRegistry(const objectRegistry&);


    //- Destructor
    virtual ~objectRegistry();


    // Member Functions

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

            //- Return the object registry
            const objectRegistry& thisDb() const
            {
                return *this;
            }


        // Edit

            //- Rename
            virtual void rename(const word& newName);

            //- Add an regIOobject to registry
            bool checkIn(regIOobject&) const;

            //- Remove an regIOobject from registry
            bool checkOut(regIOobject&) const;

            //- Remove all regIOobject owned by the registry
            void clear();

            //- Add the given name to the set of temporary objects to cache
            void addTemporaryObject(const word& name) const;

            //- Return true if given name is in the cacheTemporaryObjects set
            bool cacheTemporaryObject(const word& name) const;

            //- Cache the given object
            template<class Object>
            bool cacheTemporaryObject(Object& ob) const;

            //- Reset the cache state of the given object
            //  in the cacheTemporaryObjects set
            void resetCacheTemporaryObject(const regIOobject& ob) const;

            //- Check that all objects in the cacheTemporaryObjects set
            //  were cached
            bool checkCacheTemporaryObjects() const;


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
                const bool write
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const objectRegistry&) = delete;
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
