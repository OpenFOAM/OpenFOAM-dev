/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "functionObjectList.H"
#include "Time.H"
#include "mapPolyMesh.H"
#include "argList.H"
#include "timeControlFunctionObject.H"
#include "dictionaryEntry.H"
#include "stringOps.H"
#include "Tuple2.H"
#include "etcFiles.H"
#include "IOdictionary.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

Foam::fileName Foam::functionObjectList::functionObjectDictPath
(
    "caseDicts/postProcessing"
);


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


void Foam::functionObjectList::listDir
(
    const fileName& dir,
    HashSet<word>& foMap
)
{
    // Search specified directory for functionObject configuration files
    {
        fileNameList foFiles(fileHandler().readDir(dir));
        forAll(foFiles, f)
        {
            if (foFiles[f].ext().empty())
            {
                foMap.insert(foFiles[f]);
            }
        }
    }

    // Recurse into sub-directories
    {
        fileNameList foDirs(fileHandler().readDir(dir, fileName::DIRECTORY));
        forAll(foDirs, fd)
        {
            listDir(dir/foDirs[fd], foMap);
        }
    }
}


void Foam::functionObjectList::list()
{
    HashSet<word> foMap;

    fileNameList etcDirs(findEtcDirs(functionObjectDictPath));

    forAll(etcDirs, ed)
    {
        listDir(etcDirs[ed], foMap);
    }

    Info<< nl
        << "Available configured functionObjects:"
        << foMap.sortedToc()
        << nl;
}


Foam::fileName Foam::functionObjectList::findRegionDict
(
    const word& funcName,
    const word& region
)
{
    // First check if there is a functionObject dictionary file in the
    // case system directory
    fileName dictFile
    (
        stringOps::expand("$FOAM_CASE")/"system"/region/funcName
    );

    if (isFile(dictFile))
    {
        return dictFile;
    }
    else
    {
        fileNameList etcDirs(findEtcDirs(functionObjectDictPath));

        forAll(etcDirs, i)
        {
            dictFile = search(funcName, etcDirs[i]);
            if (!dictFile.empty())
            {
                return dictFile;
            }
        }
    }

    return fileName::null;
}


Foam::fileName Foam::functionObjectList::findDict
(
    const word& funcName,
    const word& region
)
{
    if (region == word::null)
    {
        return findRegionDict(funcName);
    }
    else
    {
        fileName dictFile(findRegionDict(funcName, region));

        if (dictFile != fileName::null)
        {
            return dictFile;
        }
        else
        {
            return findRegionDict(funcName);
        }
    }
}


bool Foam::functionObjectList::readFunctionObject
(
    const string& funcNameArgs,
    dictionary& functionsDict,
    HashSet<word>& requiredFields,
    const word& region
)
{
    // Parse the optional functionObject arguments:
    //     'Q(U)' -> funcName = Q; args = (U); field = U
    //
    // Supports named arguments:
    //     'patchAverage(patch=inlet, p)' -> funcName = patchAverage;
    //         args = (patch=inlet, p); field = p

    word funcName(funcNameArgs);

    int argLevel = 0;
    wordList args;

    List<Tuple2<word, string>> namedArgs;
    bool namedArg = false;
    word argName;

    word::size_type start = 0;
    word::size_type i = 0;

    for
    (
        word::const_iterator iter = funcNameArgs.begin();
        iter != funcNameArgs.end();
        ++iter
    )
    {
        char c = *iter;

        if (c == '(')
        {
            if (argLevel == 0)
            {
                funcName = funcNameArgs(start, i - start);
                start = i+1;
            }
            ++argLevel;
        }
        else if (c == ',' || c == ')')
        {
            if (argLevel == 1)
            {
                if (namedArg)
                {
                    namedArgs.append
                    (
                        Tuple2<word, string>
                        (
                            argName,
                            funcNameArgs(start, i - start)
                        )
                    );
                    namedArg = false;
                }
                else
                {
                    args.append
                    (
                        string::validate<word>(funcNameArgs(start, i - start))
                    );
                }
                start = i+1;
            }

            if (c == ')')
            {
                if (argLevel == 1)
                {
                    break;
                }
                --argLevel;
            }
        }
        else if (c == '=')
        {
            argName = string::validate<word>(funcNameArgs(start, i - start));
            start = i+1;
            namedArg = true;
        }

        ++i;
    }

    // Search for the functionObject dictionary
    fileName path = findDict(funcName, region);

    if (path == fileName::null)
    {
        WarningInFunction
            << "Cannot find functionObject file " << funcName << endl;
        return false;
    }

    // Read the functionObject dictionary
    //IFstream fileStream(path);
    autoPtr<ISstream> fileStreamPtr(fileHandler().NewIFstream(path));
    ISstream& fileStream = fileStreamPtr();

    dictionary funcsDict(fileStream);
    dictionary* funcDictPtr = &funcsDict;

    if (funcsDict.found(funcName) && funcsDict.isDict(funcName))
    {
        funcDictPtr = &funcsDict.subDict(funcName);
    }

    dictionary& funcDict = *funcDictPtr;

    // Insert the 'field' and/or 'fields' entry corresponding to the optional
    // arguments or read the 'field' or 'fields' entry and add the required
    // fields to requiredFields
    if (args.size() == 1)
    {
        funcDict.set("field", args[0]);
        funcDict.set("fields", args);
        requiredFields.insert(args[0]);
    }
    else if (args.size() > 1)
    {
        funcDict.set("fields", args);
        requiredFields.insert(args);
    }
    else if (funcDict.found("field"))
    {
        requiredFields.insert(word(funcDict.lookup("field")));
    }
    else if (funcDict.found("fields"))
    {
        requiredFields.insert(wordList(funcDict.lookup("fields")));
    }

    // Insert named arguments
    forAll(namedArgs, i)
    {
        IStringStream entryStream
        (
            namedArgs[i].first() + ' ' + namedArgs[i].second() + ';'
        );
        funcDict.set(entry::New(entryStream).ptr());
    }

    // Insert the region name if specified
    if (region != word::null)
    {
        funcDict.set("region", region);
    }

    // Merge this functionObject dictionary into functionsDict
    dictionary funcArgsDict;
    funcArgsDict.add(string::validate<word>(funcNameArgs), funcDict);
    functionsDict.merge(funcArgsDict);

    return true;
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
    dictionary& controlDict,
    HashSet<word>& requiredFields
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
            readFunctionObject
            (
                args["func"],
                functionsDict,
                requiredFields,
                region
            );
        }

        if (args.optionFound("funcs"))
        {
            wordList funcs(args.optionLookup("funcs")());

            forAll(funcs, i)
            {
                readFunctionObject
                (
                    funcs[i],
                    functionsDict,
                    requiredFields,
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

    functionsPtr->start();

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


bool Foam::functionObjectList::execute()
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
            ok = operator[](objectI).execute() && ok;
            ok = operator[](objectI).write() && ok;
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


bool Foam::functionObjectList::setTimeStep()
{
    bool set = true;

    if (execution_)
    {
        if (!updated_)
        {
            read();
        }

        wordList names;

        forAll(*this, objectI)
        {
            if (operator[](objectI).setTimeStep())
            {
                names.append(operator[](objectI).name());
                set = true;
            }
        }

        if (names.size() > 1)
        {
            WarningInFunction << "Multiple function objects (" << names[0];
            for (label i = 1; i < names.size(); ++ i)
            {
                WarningInFunction << ", " << names[i];
            }
            WarningInFunction << ") are setting the time step." << endl;
        }
    }

    return set;
}


Foam::scalar Foam::functionObjectList::timeToNextWrite()
{
    scalar result = vGreat;

    if (execution_)
    {
        if (!updated_)
        {
            read();
        }

        forAll(*this, objectI)
        {
            result = min(result, operator[](objectI).timeToNextWrite());
        }
    }

    return result;
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
                    objPtr = nullptr;
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
