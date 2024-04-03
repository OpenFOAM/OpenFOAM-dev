/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "objectRegistry.H"
#include "stringListOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::wordList Foam::objectRegistry::toc() const
{
    wordList objectNames(size());

    label count=0;
    forAllConstIter(HashTable<regIOobject*>, *this, iter)
    {
        if (isA<Type>(*iter()))
        {
            objectNames[count++] = iter()->name();
        }
    }

    objectNames.setSize(count);

    return objectNames;
}


template<class Type>
Foam::wordList Foam::objectRegistry::toc(const wordRe& name) const
{
    wordList objectNames(size());

    label count = 0;
    forAllConstIter(HashTable<regIOobject*>, *this, iter)
    {
        if (isA<Type>(*iter()))
        {
            const word& objectName = iter()->name();

            if (name.match(objectName))
            {
                objectNames[count++] = objectName;
            }
        }
    }

    objectNames.setSize(count);

    return objectNames;
}


template<class Type>
Foam::wordList Foam::objectRegistry::toc(const wordReList& patterns) const
{
    wordList names(this->toc<Type>());

    return wordList(names, findStrings(patterns, names));
}


template<class Type>
Foam::HashTable<const Type*> Foam::objectRegistry::lookupClass
(
    const bool strict
) const
{
    HashTable<const Type*> objectsOfClass(size());

    forAllConstIter(HashTable<regIOobject*>, *this, iter)
    {
        if
        (
            (strict && isType<Type>(*iter()))
         || (!strict && isA<Type>(*iter()))
        )
        {
            objectsOfClass.insert
            (
                iter()->name(),
                dynamic_cast<const Type*>(iter())
            );
        }
    }

    return objectsOfClass;
}


template<class Type>
Foam::HashTable<Type*> Foam::objectRegistry::lookupClass
(
    const bool strict
)
{
    HashTable<Type*> objectsOfClass(size());

    forAllIter(HashTable<regIOobject*>, *this, iter)
    {
        if
        (
            (strict && isType<Type>(*iter()))
         || (!strict && isA<Type>(*iter()))
        )
        {
            objectsOfClass.insert
            (
                iter()->name(),
                dynamic_cast<Type*>(iter())
            );
        }
    }

    return objectsOfClass;
}


template<class Type>
bool Foam::objectRegistry::foundObject(const word& name) const
{
    const_iterator iter = find(name);

    if (iter != end())
    {
        const Type* vpsiPtr_ = dynamic_cast<const Type*>(iter());

        if (vpsiPtr_)
        {
            return true;
        }
    }
    else if (this->parentNotTime())
    {
        return parent_.foundObject<Type>(name);
    }

    return false;
}


template<class Type>
const Type& Foam::objectRegistry::lookupObject(const word& name) const
{
    const_iterator iter = find(name);

    if (iter != end())
    {
        const Type* vpsiPtr_ = dynamic_cast<const Type*>(iter());

        if (vpsiPtr_)
        {
            return *vpsiPtr_;
        }

        FatalErrorInFunction
            << nl
            << "    lookup of " << name << " from objectRegistry "
            << this->name()
            << " successful\n    but it is not a " << Type::typeName
            << ", it is a " << iter()->type()
            << abort(FatalError);
    }
    else
    {
        if (this->parentNotTime())
        {
            return parent_.lookupObject<Type>(name);
        }

        FatalErrorInFunction
            << nl
            << "    request for " << Type::typeName
            << " " << name << " from objectRegistry " << this->name()
            << " failed\n    available objects of type " << Type::typeName
            << " are" << nl
            << toc<Type>();

        if (cacheTemporaryObject(name))
        {
            FatalErrorInFunction
                << nl
                << "    request for " << name << " from objectRegistry "
                << this->name() << " to be cached failed" << nl
                << "    available temporary objects are" << nl
                << temporaryObjects_;
        }

        FatalErrorInFunction
            << abort(FatalError);
    }

    return NullObjectRef<Type>();
}


template<class Type>
Type& Foam::objectRegistry::lookupObjectRef(const word& name) const
{
    return const_cast<Type&>(lookupObject<Type>(name));
}


template<class Type>
bool Foam::objectRegistry::foundType(const word& group) const
{
    return foundObject<Type>(IOobject::groupName(Type::typeName, group));
}


template<class Type>
const Type& Foam::objectRegistry::lookupType(const word& group) const
{
    return lookupObject<Type>(IOobject::groupName(Type::typeName, group));
}


template<class Object>
bool Foam::objectRegistry::cacheTemporaryObject(Object& ob) const
{
    readCacheTemporaryObjects();

    if (cacheTemporaryObjects_.size())
    {
        temporaryObjects_.insert(ob.name());

        HashTable<Pair<bool>>::iterator iter
        (
            cacheTemporaryObjects_.find(ob.name())
        );

        // Cache object ob if is in the cacheTemporaryObjects list
        // and hasn't been cached yet
        if (iter != cacheTemporaryObjects_.end() && iter().first() == false)
        {
            iter().first() = true;
            iter().second() = true;

            if (ob.db().template foundObject<Object>(ob.name()))
            {
                Object& cachedOb =
                    ob.db().template lookupObjectRef<Object>(ob.name());

                // If the object is already cached in the database delete it
                if (&cachedOb != &ob && cachedOb.ownedByRegistry())
                {
                    deleteCachedObject(cachedOb);
                }
            }

            if (debug)
            {
                Info<< "Caching " << ob.name()
                    << " of type " << ob.type() << endl;
            }

            ob.release();
            ob.checkOut();
            store(new Object(move(ob)));

            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
