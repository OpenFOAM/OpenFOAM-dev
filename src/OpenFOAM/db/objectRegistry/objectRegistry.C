/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(objectRegistry, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::objectRegistry::parentNotTime() const
{
    return (&parent_ != dynamic_cast<const objectRegistry*>(&time_));
}


void Foam::objectRegistry::readCacheTemporaryObjects() const
{
    if
    (
        !cacheTemporaryObjectsSet_
     && time_.controlDict().found("cacheTemporaryObjects")
    )
    {
        cacheTemporaryObjectsSet_ = true;

        const dictionary& controlDict = time_.controlDict();

        wordList cacheTemporaryObjects;

        if (controlDict.isDict("cacheTemporaryObjects"))
        {
            if(controlDict.subDict("cacheTemporaryObjects").found(name()))
            {
                controlDict.subDict("cacheTemporaryObjects").lookup(name())
                    >> cacheTemporaryObjects;
            }
        }
        else
        {
            controlDict.lookup("cacheTemporaryObjects")
                >> cacheTemporaryObjects;
        }

        forAll(cacheTemporaryObjects, i)
        {
            cacheTemporaryObjects_.insert
            (
                cacheTemporaryObjects[i],
                {false, false}
            );
        }
    }
}


void Foam::objectRegistry::deleteCachedObject(regIOobject& cachedOb) const
{
    cachedOb.release();
    cachedOb.checkOut();
    cachedOb.rename(cachedOb.name() + "Cached");
    delete &cachedOb;
}


// * * * * * * * * * * * * * * * * Constructors *  * * * * * * * * * * * * * //

Foam::objectRegistry::objectRegistry
(
    const Time& t,
    const label nIoObjects
)
:
    regIOobject
    (
        IOobject
        (
            string::validate<word>(t.caseName()),
            t.path(),
            t,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        ),
        true    // to flag that this is the top-level regIOobject
    ),
    HashTable<regIOobject*>(nIoObjects),
    time_(t),
    parent_(t),
    dbDir_(name()),
    event_(1),
    cacheTemporaryObjectsSet_(false)
{}


Foam::objectRegistry::objectRegistry
(
    const IOobject& io,
    const label nIoObjects
)
:
    regIOobject(io),
    HashTable<regIOobject*>(nIoObjects),
    time_(io.time()),
    parent_(io.db()),
    dbDir_(parent_.dbDir()/local()/name()),
    event_(1),
    cacheTemporaryObjectsSet_(false)
{
    writeOpt() = IOobject::AUTO_WRITE;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::objectRegistry::~objectRegistry()
{
    cacheTemporaryObjects_.clear();
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::objectRegistry::names() const
{
    return HashTable<regIOobject*>::toc();
}


Foam::wordList Foam::objectRegistry::sortedNames() const
{
    return HashTable<regIOobject*>::sortedToc();
}


Foam::wordList Foam::objectRegistry::names(const word& ClassName) const
{
    wordList objectNames(size());

    label count=0;
    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        if (iter()->type() == ClassName)
        {
            objectNames[count++] = iter.key();
        }
    }

    objectNames.setSize(count);

    return objectNames;
}


Foam::wordList Foam::objectRegistry::sortedNames(const word& ClassName) const
{
    wordList sortedLst = names(ClassName);
    sort(sortedLst);

    return sortedLst;
}


const Foam::objectRegistry& Foam::objectRegistry::subRegistry
(
    const word& name,
    const bool forceCreate
) const
{
    if (forceCreate && !foundObject<objectRegistry>(name))
    {
        objectRegistry* fieldsCachePtr = new objectRegistry
        (
            IOobject
            (
                name,
                time().constant(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );
        fieldsCachePtr->store();
    }
    return lookupObject<objectRegistry>(name);
}


Foam::label Foam::objectRegistry::getEvent() const
{
    label curEvent = event_++;

    if (event_ == labelMax)
    {
        if (objectRegistry::debug)
        {
            WarningInFunction
                << "Event counter has overflowed. "
                << "Resetting counter on all dependent objects." << nl
                << "This might cause extra evaluations." << endl;
        }

        // Reset event counter
        curEvent = 1;
        event_ = 2;

        for (const_iterator iter = begin(); iter != end(); ++iter)
        {
            const regIOobject& io = *iter();

            if (objectRegistry::debug)
            {
                Pout<< "objectRegistry::getEvent() : "
                    << "resetting count on " << iter.key() << endl;
            }

            if (io.eventNo() != 0)
            {
                const_cast<regIOobject&>(io).eventNo() = curEvent;
            }
        }
    }

    return curEvent;
}


bool Foam::objectRegistry::checkIn(regIOobject& io) const
{
    if (objectRegistry::debug)
    {
        Pout<< "objectRegistry::checkIn(regIOobject&) : "
            << name() << " : checking in " << io.name()
            << " of type " << io.type()
            << endl;
    }

    // Delete cached object with the same name as io and if it is in the
    // cacheTemporaryObjects list
    if (cacheTemporaryObjects_.size())
    {
        HashTable<Pair<bool>>::iterator cacheIter
        (
            cacheTemporaryObjects_.find(io.name())
        );

        if (cacheIter != cacheTemporaryObjects_.end())
        {
            iterator iter = const_cast<objectRegistry&>(*this).find(io.name());

            if (iter != end() && iter() != &io && iter()->ownedByRegistry())
            {
                if (objectRegistry::debug)
                {
                    Pout<< "objectRegistry::checkIn(regIOobject&) : "
                        << name() << " : deleting cached object " << iter.key()
                        << endl;
                }

                cacheIter().first() = false;
                deleteCachedObject(*iter());
            }
        }
    }

    return const_cast<objectRegistry&>(*this).insert(io.name(), &io);
}


bool Foam::objectRegistry::checkOut(regIOobject& io) const
{
    iterator iter = const_cast<objectRegistry&>(*this).find(io.name());

    if (iter != end())
    {
        if (objectRegistry::debug)
        {
            Pout<< "objectRegistry::checkOut(regIOobject&) : "
                << name() << " : checking out " << iter.key()
                << endl;
        }

        if (iter() != &io)
        {
            if (objectRegistry::debug)
            {
                WarningInFunction
                    << name() << " : attempt to checkOut copy of "
                    << iter.key()
                    << endl;
            }

            return false;
        }
        else
        {
            regIOobject* object = iter();

            bool hasErased = const_cast<objectRegistry&>(*this).erase(iter);

            if (io.ownedByRegistry())
            {
                delete object;
            }

            return hasErased;
        }
    }
    else
    {
        if (objectRegistry::debug)
        {
            Pout<< "objectRegistry::checkOut(regIOobject&) : "
                << name() << " : could not find " << io.name()
                << " in registry " << name()
                << endl;
        }
    }

    return false;
}


void Foam::objectRegistry::clear()
{
    List<regIOobject*> myObjects(size());
    label nMyObjects = 0;

    for (iterator iter = begin(); iter != end(); ++iter)
    {
        if (iter()->ownedByRegistry())
        {
            myObjects[nMyObjects++] = iter();
        }
    }

    for (label i=0; i < nMyObjects; i++)
    {
        checkOut(*myObjects[i]);
    }
}


void Foam::objectRegistry::addTemporaryObject
(
    const word& name
) const
{
    if (!cacheTemporaryObjects_.found(name))
    {
        cacheTemporaryObjects_.insert(name, {false, false});
    }
}


bool Foam::objectRegistry::cacheTemporaryObject
(
    const word& name
) const
{
    return cacheTemporaryObjects_.found(name);
}


void Foam::objectRegistry::resetCacheTemporaryObject
(
    const regIOobject& ob
) const
{
    if (cacheTemporaryObjects_.size())
    {
        HashTable<Pair<bool>>::iterator iter
        (
            cacheTemporaryObjects_.find(ob.name())
        );

        // If object ob if is in the cacheTemporaryObjects list
        // and has been cached reset the cached flag
        if (iter != cacheTemporaryObjects_.end())
        {
            iter().first() = false;
        }
    }
}


bool Foam::objectRegistry::checkCacheTemporaryObjects() const
{
    bool enabled = cacheTemporaryObjects_.size();

    forAllConstIter(HashTable<regIOobject*>, *this, iter)
    {
        const objectRegistry* orPtr_ =
            dynamic_cast<const objectRegistry*>(iter());

        // Protect against re-searching the top-level registry
        if (orPtr_ && orPtr_ != this)
        {
            enabled = orPtr_->checkCacheTemporaryObjects() || enabled;
        }
    }

    if (enabled)
    {
        forAllIter
        (
            typename HashTable<Pair<bool>>,
            cacheTemporaryObjects_,
            iter
        )
        {
            if (!iter().second())
            {
                Warning
                    << "Could not find temporary object " << iter.key()
                    << " in registry " << name() << nl
                    << "Available temporary objects "
                    << temporaryObjects_
                    << endl;
            }
            else
            {
                iter().second() = false;
            }
        }

        temporaryObjects_.clear();
    }

    return enabled;
}


void Foam::objectRegistry::rename(const word& newName)
{
    regIOobject::rename(newName);

    // adjust dbDir_ as well
    string::size_type i = dbDir_.rfind('/');

    if (i == string::npos)
    {
        dbDir_ = newName;
    }
    else
    {
        dbDir_.replace(i+1, string::npos, newName);
    }
}


bool Foam::objectRegistry::modified() const
{
    forAllConstIter(HashTable<regIOobject*>, *this, iter)
    {
        if (iter()->modified())
        {
            return true;
        }
    }

    return false;
}


void Foam::objectRegistry::readModifiedObjects()
{
    for (iterator iter = begin(); iter != end(); ++iter)
    {
        if (objectRegistry::debug)
        {
            Pout<< "objectRegistry::readModifiedObjects() : "
                << name() << " : Considering reading object "
                << iter.key() << endl;
        }

        iter()->readIfModified();
    }
}


bool Foam::objectRegistry::readIfModified()
{
    readModifiedObjects();
    return true;
}


bool Foam::objectRegistry::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    bool ok = true;

    forAllConstIter(HashTable<regIOobject*>, *this, iter)
    {
        if (objectRegistry::debug)
        {
            Pout<< "objectRegistry::write() : "
                << name() << " : Considering writing object "
                << iter.key()
                << " of type " << iter()->type()
                << " with writeOpt " << iter()->writeOpt()
                << " to file " << iter()->objectPath()
                << endl;
        }

        if (iter()->writeOpt() != NO_WRITE)
        {
            ok = iter()->writeObject(fmt, ver, cmp, write) && ok;
        }
    }

    return ok;
}


// ************************************************************************* //
