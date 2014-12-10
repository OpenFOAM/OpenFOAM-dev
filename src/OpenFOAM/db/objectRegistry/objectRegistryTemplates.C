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

\*---------------------------------------------------------------------------*/

#include "objectRegistry.H"
#include "stringListOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::wordList Foam::objectRegistry::names() const
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
Foam::wordList Foam::objectRegistry::names(const wordRe& name) const
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
Foam::wordList Foam::objectRegistry::names(const wordReList& patterns) const
{
    wordList names(this->names<Type>());

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

        FatalErrorIn
        (
            "objectRegistry::lookupObject<Type>(const word&) const"
        )   << nl
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

        FatalErrorIn
        (
            "objectRegistry::lookupObject<Type>(const word&) const"
        )   << nl
            << "    request for " << Type::typeName
            << " " << name << " from objectRegistry " << this->name()
            << " failed\n    available objects of type " << Type::typeName
            << " are" << nl
            << names<Type>()
            << abort(FatalError);
    }

    return *reinterpret_cast< const Type* >(0);
}


// ************************************************************************* //
