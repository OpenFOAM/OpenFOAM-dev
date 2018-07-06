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

\*---------------------------------------------------------------------------*/

#include "IOobjectList.H"
#include "Time.H"
#include "OSspecific.H"
#include "IOList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IOobjectList::IOobjectList(const label nIoObjects)
:
    HashPtrTable<IOobject>(nIoObjects)
{}


Foam::IOobjectList::IOobjectList
(
    const objectRegistry& db,
    const fileName& instance,
    const fileName& local,
    IOobject::readOption r,
    IOobject::writeOption w,
    bool registerObject
)
:
    HashPtrTable<IOobject>()
{
    word newInstance;
    fileNameList ObjectNames = fileHandler().readObjects
    (
        db,
        instance,
        local,
        newInstance
    );


    forAll(ObjectNames, i)
    {
        IOobject* objectPtr = new IOobject
        (
            ObjectNames[i],
            newInstance,
            local,
            db,
            r,
            w,
            registerObject
        );

        // Use object with local scope
        if (objectPtr->typeHeaderOk<IOList<label>>(false))
        {
            insert(ObjectNames[i], objectPtr);
        }
        else
        {
            delete objectPtr;
        }
    }
}


Foam::IOobjectList::IOobjectList(const IOobjectList& ioOL)
:
    HashPtrTable<IOobject>(ioOL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::IOobjectList::~IOobjectList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::IOobjectList::add(IOobject& io)
{
    return insert(io.name(), &io);
}


bool Foam::IOobjectList::remove(IOobject& io)
{
    HashPtrTable<IOobject>::iterator iter =
        HashPtrTable<IOobject>::find(io.name());

    if (iter != end())
    {
        return erase(iter);
    }
    else
    {
        return false;
    }
}


Foam::IOobject* Foam::IOobjectList::lookup(const word& name) const
{
    HashPtrTable<IOobject>::const_iterator iter = find(name);

    if (iter != end())
    {
        if (IOobject::debug)
        {
            InfoInFunction << "Found " << name << endl;
        }

        return const_cast<IOobject*>(*iter);
    }
    else
    {
        if (IOobject::debug)
        {
            InfoInFunction << "Could not find " << name << endl;
        }

        return nullptr;
    }
}


Foam::IOobjectList Foam::IOobjectList::lookup(const wordRe& name) const
{
    IOobjectList objectsOfName(size());

    forAllConstIter(HashPtrTable<IOobject>, *this, iter)
    {
        if (name.match(iter()->name()))
        {
            if (IOobject::debug)
            {
                InfoInFunction << "Found " << iter.key() << endl;
            }

            objectsOfName.insert(iter.key(), new IOobject(*iter()));
        }
    }

    return objectsOfName;
}


Foam::IOobjectList Foam::IOobjectList::lookup(const wordReList& patterns) const
{
    wordReListMatcher names(patterns);

    IOobjectList objectsOfName(size());

    forAllConstIter(HashPtrTable<IOobject>, *this, iter)
    {
        if (names.match(iter()->name()))
        {
            if (IOobject::debug)
            {
                InfoInFunction << "Found " << iter.key() << endl;
            }

            objectsOfName.insert(iter.key(), new IOobject(*iter()));
        }
    }

    return objectsOfName;
}


Foam::IOobjectList Foam::IOobjectList::lookupClass(const word& ClassName) const
{
    IOobjectList objectsOfClass(size());

    forAllConstIter(HashPtrTable<IOobject>, *this, iter)
    {
        if (iter()->headerClassName() == ClassName)
        {
            if (IOobject::debug)
            {
                InfoInFunction << "Found " << iter.key() << endl;
            }

            objectsOfClass.insert(iter.key(), new IOobject(*iter()));
        }
    }

    return objectsOfClass;
}


Foam::wordList Foam::IOobjectList::names() const
{
    return HashPtrTable<IOobject>::toc();
}


Foam::wordList Foam::IOobjectList::sortedNames() const
{
    return HashPtrTable<IOobject>::sortedToc();
}


Foam::wordList Foam::IOobjectList::names(const word& ClassName) const
{
    wordList objectNames(size());

    label count = 0;
    forAllConstIter(HashPtrTable<IOobject>, *this, iter)
    {
        if (iter()->headerClassName() == ClassName)
        {
            objectNames[count++] = iter.key();
        }
    }

    objectNames.setSize(count);

    return objectNames;
}


Foam::wordList Foam::IOobjectList::sortedNames(const word& ClassName) const
{
    wordList sortedLst = names(ClassName);
    sort(sortedLst);

    return sortedLst;
}


// ************************************************************************* //
