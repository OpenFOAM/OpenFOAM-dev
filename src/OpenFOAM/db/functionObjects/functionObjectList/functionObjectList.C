/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
#include "argList.H"
#include "timeControlFunctionObject.H"
#include "dictionaryEntry.H"

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


bool Foam::functionObjectList::writeData(Ostream&) const
{
    NotImplemented;
    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjectList::functionObjectList
(
    const Time& t,
    const bool execution
)
:
    regIOobject
    (
        IOobject
        (
            "functionObjectList",
            t.system(),
            t
        )
    ),
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
    regIOobject
    (
        IOobject
        (
            "functionObjectList",
            t.system(),
            t
        )
    ),
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
    dictionary& controlDict
)
{
    autoPtr<functionObjectList> functionsPtr;

    controlDict.add
    (
        dictionaryEntry("functions", controlDict, dictionary::null)
    );

    dictionary& functionsDict = controlDict.subDict("functions");

    word region = word::null;

    // Set the region name if specified
    if (args.optionFound("region"))
    {
        region = args["region"];
    }

    if
    (
        args.optionFound("dict")
     || args.optionFound("func")
     || args.optionFound("funcs")
    )
    {
        if (args.optionFound("dict"))
        {
            controlDict.merge
            (
                IOdictionary
                (
                    IOobject
                    (
                        args["dict"],
                        runTime,
                        IOobject::MUST_READ_IF_MODIFIED
                    )
                )
            );
        }

        if (args.optionFound("func"))
        {
            readConfigFile
            (
                "function",
                args["func"],
                functionsDict,
                functionEntries::includeFuncEntry::functionObjectDictPath,
                "system",
                {"command", args.commandLine()},
                region
            );
        }

        if (args.optionFound("funcs"))
        {
            wordList funcs(args.optionLookup("funcs")());

            forAll(funcs, i)
            {
                readConfigFile
                (
                    "function",
                    funcs[i],
                    functionsDict,
                    functionEntries::includeFuncEntry::functionObjectDictPath,
                    "system",
                    {"command", args.commandLine()},
                    region
                );
            }
        }

        functionsPtr.reset(new functionObjectList(runTime, controlDict));
    }
    else
    {
        functionsPtr.reset(new functionObjectList(runTime));
    }

    functionsPtr->read();

    return functionsPtr;
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
    forAll(*this, oi)
    {
        if (operator[](oi).name() == name)
        {
            return oi;
        }
    }

    return -1;
}


bool Foam::functionObjectList::start()
{
    bool ok = read();

    if (execution_)
    {
        forAll(*this, oi)
        {
            if (operator[](oi).executeAtStart())
            {
                ok = operator[](oi).execute() && ok;
                ok = operator[](oi).write() && ok;
            }
        }
    }

    return ok;
}


bool Foam::functionObjectList::execute()
{
    bool ok = true;

    if (execution_)
    {
        if (!updated_)
        {
            read();
        }

        forAll(*this, oi)
        {
            ok = operator[](oi).execute() && ok;
            ok = operator[](oi).write() && ok;
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

        forAll(*this, oi)
        {
            ok = operator[](oi).end() && ok;
        }
    }

    return ok;
}


Foam::scalar Foam::functionObjectList::timeToNextAction()
{
    scalar result = vGreat;

    if (execution_)
    {
        if (!updated_)
        {
            read();
        }

        forAll(*this, oi)
        {
            result = min(result, operator[](oi).timeToNextAction());
        }
    }

    return result;
}


Foam::scalar Foam::functionObjectList::maxDeltaT() const
{
    scalar result = vGreat;

    forAll(*this, oi)
    {
        result = min(result, operator[](oi).maxDeltaT());
    }

    return result;
}


bool Foam::functionObjectList::dependenciesModified() const
{
    if (&parentDict_ == &time_.controlDict())
    {
        return time_.controlDict().modified();
    }
    else
    {
        return false;
    }
}


bool Foam::functionObjectList::read()
{
    bool ok = true;
    updated_ = execution_;

    // Avoid reading/initialising if execution is off
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

        libs.open
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
            const bool enabled = dict.lookupOrDefault("enabled", true);

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
                    objPtr = nullptr;
                    continue;
                }
            }
            else if (enabled)
            {
                autoPtr<functionObject> foPtr;

                if
                (
                    dict.found("executeControl")
                 || dict.found("writeControl")
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


void Foam::functionObjectList::movePoints(const polyMesh& mesh)
{
    if (execution_)
    {
        forAll(*this, oi)
        {
            operator[](oi).movePoints(mesh);
        }
    }
}


void Foam::functionObjectList::topoChange(const polyTopoChangeMap& map)
{
    if (execution_)
    {
        forAll(*this, oi)
        {
            operator[](oi).topoChange(map);
        }
    }
}


void Foam::functionObjectList::mapMesh(const polyMeshMap& map)
{
    if (execution_)
    {
        forAll(*this, oi)
        {
            operator[](oi).mapMesh(map);
        }
    }
}


void Foam::functionObjectList::distribute(const polyDistributionMap& map)
{
    if (execution_)
    {
        forAll(*this, oi)
        {
            operator[](oi).distribute(map);
        }
    }
}


// ************************************************************************* //
