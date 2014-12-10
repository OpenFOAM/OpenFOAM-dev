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

#include "functionObjectList.H"
#include "Time.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::functionObject* Foam::functionObjectList::remove
(
    const word& key,
    label& oldIndex
)
{
    functionObject* ptr = 0;

    // Find index of existing functionObject
    HashTable<label>::iterator fnd = indices_.find(key);

    if (fnd != indices_.end())
    {
        oldIndex = fnd();

        // retrieve the pointer and remove it from the old list
        ptr = this->set(oldIndex, 0).ptr();
        indices_.erase(fnd);
    }
    else
    {
        oldIndex = -1;
    }

    return ptr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjectList::functionObjectList
(
    const Time& t,
    const bool execution
)
:
    PtrList<functionObject>(),
    digests_(),
    indices_(),
    time_(t),
    parentDict_(t.controlDict()),
    execution_(execution),
    updated_(false)
{}


Foam::functionObjectList::functionObjectList
(
    const Time& t,
    const dictionary& parentDict,
    const bool execution
)
:
    PtrList<functionObject>(),
    digests_(),
    indices_(),
    time_(t),
    parentDict_(parentDict),
    execution_(execution),
    updated_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjectList::~functionObjectList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjectList::clear()
{
    PtrList<functionObject>::clear();
    digests_.clear();
    indices_.clear();
    updated_ = false;
}


Foam::label Foam::functionObjectList::findObjectID(const word& name) const
{
    forAll(*this, objectI)
    {
        if (operator[](objectI).name() == name)
        {
            return objectI;
        }
    }

    return -1;
}


void Foam::functionObjectList::on()
{
    execution_ = true;
}


void Foam::functionObjectList::off()
{
    // for safety, also force a read() when execution is turned back on
    updated_ = execution_ = false;
}


bool Foam::functionObjectList::status() const
{
    return execution_;
}


bool Foam::functionObjectList::start()
{
    return read();
}


bool Foam::functionObjectList::execute(const bool forceWrite)
{
    bool ok = true;

    if (execution_)
    {
        if (!updated_)
        {
            read();
        }

        forAll(*this, objectI)
        {
            ok = operator[](objectI).execute(forceWrite) && ok;
        }
    }

    return ok;
}


bool Foam::functionObjectList::end()
{
    bool ok = true;

    if (execution_)
    {
        if (!updated_)
        {
            read();
        }

        forAll(*this, objectI)
        {
            ok = operator[](objectI).end() && ok;
        }
    }

    return ok;
}


bool Foam::functionObjectList::timeSet()
{
    bool ok = true;

    if (execution_)
    {
        if (!updated_)
        {
            read();
        }

        forAll(*this, objectI)
        {
            ok = operator[](objectI).timeSet() && ok;
        }
    }

    return ok;
}


bool Foam::functionObjectList::adjustTimeStep()
{
    bool ok = true;

    if (execution_)
    {
        if (!updated_)
        {
            read();
        }

        forAll(*this, objectI)
        {
            ok = operator[](objectI).adjustTimeStep() && ok;
        }
    }

    return ok;
}


bool Foam::functionObjectList::read()
{
    bool ok = true;
    updated_ = execution_;

    // avoid reading/initializing if execution is off
    if (!execution_)
    {
        return ok;
    }

    // Update existing and add new functionObjects
    const entry* entryPtr = parentDict_.lookupEntryPtr
    (
        "functions",
        false,
        false
    );

    if (entryPtr)
    {
        PtrList<functionObject> newPtrs;
        List<SHA1Digest> newDigs;
        HashTable<label> newIndices;

        label nFunc = 0;

        if (entryPtr->isDict())
        {
            // a dictionary of functionObjects
            const dictionary& functionDicts = entryPtr->dict();

            newPtrs.setSize(functionDicts.size());
            newDigs.setSize(functionDicts.size());

            forAllConstIter(dictionary, functionDicts, iter)
            {
                // safety:
                if (!iter().isDict())
                {
                    continue;
                }
                const word& key = iter().keyword();
                const dictionary& dict = iter().dict();

                newDigs[nFunc] = dict.digest();

                label oldIndex;
                functionObject* objPtr = remove(key, oldIndex);
                if (objPtr)
                {
                    // an existing functionObject, and dictionary changed
                    if (newDigs[nFunc] != digests_[oldIndex])
                    {
                        ok = objPtr->read(dict) && ok;
                    }
                }
                else
                {
                    // new functionObject
                    objPtr = functionObject::New(key, time_, dict).ptr();
                    ok = objPtr->start() && ok;
                }

                newPtrs.set(nFunc, objPtr);
                newIndices.insert(key, nFunc);
                nFunc++;
            }
        }
        else
        {
            // a list of functionObjects
            PtrList<entry> functionDicts(entryPtr->stream());

            newPtrs.setSize(functionDicts.size());
            newDigs.setSize(functionDicts.size());

            forAllIter(PtrList<entry>, functionDicts, iter)
            {
                // safety:
                if (!iter().isDict())
                {
                    continue;
                }
                const word& key = iter().keyword();
                const dictionary& dict = iter().dict();

                newDigs[nFunc] = dict.digest();

                label oldIndex;
                functionObject* objPtr = remove(key, oldIndex);
                if (objPtr)
                {
                    // an existing functionObject, and dictionary changed
                    if (newDigs[nFunc] != digests_[oldIndex])
                    {
                        ok = objPtr->read(dict) && ok;
                    }
                }
                else
                {
                    // new functionObject
                    objPtr = functionObject::New(key, time_, dict).ptr();
                    ok = objPtr->start() && ok;
                }

                newPtrs.set(nFunc, objPtr);
                newIndices.insert(key, nFunc);
                nFunc++;
            }
        }

        // safety:
        newPtrs.setSize(nFunc);
        newDigs.setSize(nFunc);

        // updating the PtrList of functionObjects also deletes any existing,
        // but unused functionObjects
        PtrList<functionObject>::transfer(newPtrs);
        digests_.transfer(newDigs);
        indices_.transfer(newIndices);
    }
    else
    {
        PtrList<functionObject>::clear();
        digests_.clear();
        indices_.clear();
    }

    return ok;
}


void Foam::functionObjectList::updateMesh(const mapPolyMesh& mpm)
{
    if (execution_)
    {
        forAll(*this, objectI)
        {
            operator[](objectI).updateMesh(mpm);
        }
    }
}


void Foam::functionObjectList::movePoints(const polyMesh& mesh)
{
    if (execution_)
    {
        forAll(*this, objectI)
        {
            operator[](objectI).movePoints(mesh);
        }
    }
}


// ************************************************************************* //
