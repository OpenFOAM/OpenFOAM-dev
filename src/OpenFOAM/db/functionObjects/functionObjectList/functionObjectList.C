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

#include "functionObjectList.H"
#include "argList.H"
#include "timeControlFunctionObject.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::typeIOobject<Foam::IOdictionary>
Foam::functionObjectList::readFunctionsDict
(
    const Time& t,
    const bool execution
)
{
    if (execution)
    {
        typeIOobject<IOdictionary> functionsDict
        (
            "functions",
            t.system(),
            t,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            true
        );

        if (functionsDict.headerOk())
        {
            return functionsDict;
        }
    }

    return typeIOobject<IOdictionary>
    (
        "functions",
        t.system(),
        t,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    );
}


bool Foam::functionObjectList::readDict()
{
    bool ok = true;
    updated_ = execution_;

    // Avoid reading/initialising if execution is off
    if (!execution_)
    {
        return true;
    }

    if (IOdictionary::size())
    {
        PtrList<functionObject> newPtrs;
        List<SHA1Digest> newDigs;
        HashTable<label> newIndices;

        label nFunc = 0;

        const dictionary& functionsDict = *this;

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
                    IOWarningInFunction(*this)
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
        ptr = this->PtrList<functionObject>::set(oldIndex, 0).ptr();
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
    IOdictionary(readFunctionsDict(t, execution)),
    PtrList<functionObject>(),
    digests_(),
    indices_(),
    time_(t),
    execution_(execution),
    updated_(false)
{
    if (t.controlDict().found("functions"))
    {
        if (!headerOk())
        {
            merge(t.controlDict().subDict("functions"));
        }
        else
        {
            WarningInFunction
                << "Both " << relativeObjectPath()
                << " and " << t.controlDict().relativeObjectPath()/"functions"
                << " found, the latter will be ignored." << endl;
        }
    }
}


Foam::autoPtr<Foam::functionObjectList> Foam::functionObjectList::New
(
    const argList& args,
    const Time& runTime
)
{
    autoPtr<functionObjectList> functionsPtr(new functionObjectList(runTime));
    dictionary& functionsDict = *functionsPtr;

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
        // Remove functions specified in the functions dictionary
        functionsDict.clear();

        if (args.optionFound("dict"))
        {
            functionsDict.merge
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
    }

    functionsPtr->readDict();

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
    bool ok = readDict();

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
            readDict();
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
            readDict();
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
            readDict();
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


bool Foam::functionObjectList::read()
{
    if (regIOobject::read())
    {
        return readDict();
    }
    else
    {
        return false;
    }
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
