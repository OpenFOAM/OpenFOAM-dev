/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
#include "argList.H"
#include "timeControlFunctionObject.H"

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

        // Retrieve the pointer and remove it from the old list
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


Foam::autoPtr<Foam::functionObjectList> Foam::functionObjectList::New
(
    const argList& args,
    const Time& runTime,
    dictionary& functionObjectsDict
)
{
    autoPtr<functionObjectList> functionObjectsPtr;

    if (args.optionFound("dict"))
    {
        functionObjectsDict = IOdictionary
        (
            IOobject
            (
                args["dict"],
                runTime,
                IOobject::MUST_READ_IF_MODIFIED
            )
        );

        functionObjectsPtr.reset
        (
            new functionObjectList(runTime, functionObjectsDict)
        );
    }
    else
    {
        functionObjectsPtr.reset(new functionObjectList(runTime));
    }
    functionObjectsPtr->start();

    return functionObjectsPtr;
}


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
    // For safety, also force a read() when execution is turned back on
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


bool Foam::functionObjectList::execute(const bool postProcess)
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
            ok = operator[](objectI).execute(postProcess) && ok;
            ok = operator[](objectI).write(postProcess) && ok;
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

    // Avoid reading/initializing if execution is off
    if (!execution_)
    {
        return true;
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

        if (!entryPtr->isDict())
        {
            FatalIOErrorInFunction(parentDict_)
                << "'functions' entry is not a dictionary"
                << exit(FatalIOError);
        }

        const dictionary& functionsDict = entryPtr->dict();

        const_cast<Time&>(time_).libs().open
        (
            functionsDict,
            "libs",
            functionObject::dictionaryConstructorTablePtr_
        );

        newPtrs.setSize(functionsDict.size());
        newDigs.setSize(functionsDict.size());

        forAllConstIter(dictionary, functionsDict, iter)
        {
            const word& key = iter().keyword();

            if (!iter().isDict())
            {
                if (key != "libs")
                {
                    IOWarningInFunction(parentDict_)
                        << "Entry " << key << " is not a dictionary" << endl;
                }

                continue;
            }

            const dictionary& dict = iter().dict();
            bool enabled = dict.lookupOrDefault("enabled", true);

            newDigs[nFunc] = dict.digest();

            label oldIndex;
            functionObject* objPtr = remove(key, oldIndex);

            if (objPtr)
            {
                if (enabled)
                {
                    // Dictionary changed for an existing functionObject
                    if (newDigs[nFunc] != digests_[oldIndex])
                    {
                        ok = objPtr->read(dict) && ok;
                    }
                }
                else
                {
                    // Delete the disabled functionObject
                    delete objPtr;
                    objPtr = NULL;
                    continue;
                }
            }
            else if (enabled)
            {
                autoPtr<functionObject> foPtr;

                FatalError.throwExceptions();
                FatalIOError.throwExceptions();
                try
                {
                    if
                    (
                        dict.found("writeControl")
                     || dict.found("outputControl")
                    )
                    {
                        foPtr.set
                        (
                            new functionObjects::timeControl(key, time_, dict)
                        );
                    }
                    else
                    {
                        foPtr = functionObject::New(key, time_, dict);
                    }
                }
                catch (Foam::IOerror& ioErr)
                {
                    Info<< ioErr << nl << endl;
                    ::exit(1);
                }
                catch (Foam::error& err)
                {
                    WarningInFunction
                        << "Caught FatalError " << err << nl << endl;
                }
                FatalError.dontThrowExceptions();
                FatalIOError.dontThrowExceptions();

                if (foPtr.valid())
                {
                    objPtr = foPtr.ptr();
                }
                else
                {
                    ok = false;
                }
            }

            // Insert active functionObjects into the list
            if (objPtr)
            {
                newPtrs.set(nFunc, objPtr);
                newIndices.insert(key, nFunc);
                nFunc++;
            }
        }

        newPtrs.setSize(nFunc);
        newDigs.setSize(nFunc);

        // Updating the PtrList of functionObjects deletes any
        // existing unused functionObjects
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
