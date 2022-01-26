/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
#include "stringOps.H"
#include "etcFiles.H"
#include "wordAndDictionary.H"


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
        fileNameList foDirs(fileHandler().readDir(dir, fileType::directory));
        forAll(foDirs, fd)
        {
            listDir(dir/foDirs[fd], foMap);
        }
    }
}


Foam::wordList Foam::functionObjectList::list()
{
    HashSet<word> foMap;

    fileNameList etcDirs(findEtcDirs(functionObjectDictPath));

    forAll(etcDirs, ed)
    {
        listDir(etcDirs[ed], foMap);
    }

    return foMap.sortedToc();
}


Foam::fileName Foam::functionObjectList::findDict
(
    const word& funcName,
    const word& region
)
{
    // First check if there is a functionObject dictionary file in the
    // region system directory
    {
        const fileName dictFile
        (
            stringOps::expand("$FOAM_CASE")/"system"/region/funcName
        );

        if (isFile(dictFile))
        {
            return dictFile;
        }
    }

    // Next, if the region is specified, check if there is a functionObject
    // dictionary file in the global system directory
    if (region != word::null)
    {
        const fileName dictFile
        (
            stringOps::expand("$FOAM_CASE")/"system"/funcName
        );

        if (isFile(dictFile))
        {
            return dictFile;
        }
    }

    // Finally, check etc directories
    {
        const fileNameList etcDirs(findEtcDirs(functionObjectDictPath));

        forAll(etcDirs, i)
        {
            const fileName dictFile(search(funcName, etcDirs[i]));

            if (!dictFile.empty())
            {
                return dictFile;
            }
        }
    }

    return fileName::null;
}


Foam::List<Foam::Tuple2<Foam::word, Foam::string>>
Foam::functionObjectList::unsetEntries(const dictionary& funcDict)
{
    const wordRe unsetPattern("<.*>");
    unsetPattern.compile();

    List<Tuple2<word, string>> unsetArgs;

    forAllConstIter(IDLList<entry>, funcDict, iter)
    {
        if (iter().isStream())
        {
            ITstream& its = iter().stream();
            OStringStream oss;
            bool isUnset = false;

            forAll(its, i)
            {
                oss << its[i];
                if (its[i].isWord() && unsetPattern.match(its[i].wordToken()))
                {
                    isUnset = true;
                }
            }

            if (isUnset)
            {
                unsetArgs.append
                (
                    Tuple2<word, string>
                    (
                        iter().keyword(),
                        oss.str()
                    )
                );
            }
        }
        else
        {
            List<Tuple2<word, string>> subUnsetArgs =
                unsetEntries(iter().dict());

            forAll(subUnsetArgs, i)
            {
                unsetArgs.append
                (
                    Tuple2<word, string>
                    (
                        iter().keyword() + '/' + subUnsetArgs[i].first(),
                        subUnsetArgs[i].second()
                    )
                );
            }
        }
    }

    return unsetArgs;
}


bool Foam::functionObjectList::readFunctionObject
(
    const string& funcArgs,
    dictionary& functionsDict,
    const Pair<string>& contextTypeAndValue,
    const word& region
)
{
    word funcType;
    wordReList args;
    List<Tuple2<word, string>> namedArgs;

    dictArgList(funcArgs, funcType, args, namedArgs);

    // Search for the functionObject dictionary
    fileName path = findDict(funcType, region);

    if (path == fileName::null)
    {
        FatalIOErrorInFunction(functionsDict)
            << "Cannot find functionObject configuration file "
            << funcType << nl << nl
            << "Available configured functionObjects:"
            << list()
            << exit(FatalIOError);
        return false;
    }

    // Read the functionObject dictionary
    // IFstream fileStream(path);
    autoPtr<ISstream> fileStreamPtr(fileHandler().NewIFstream(path));
    ISstream& fileStream = fileStreamPtr();

    // Delay processing the functionEntries
    // until after the function argument entries have been added
    entry::disableFunctionEntries = true;
    dictionary funcsDict(funcType, functionsDict, fileStream);
    entry::disableFunctionEntries = false;

    dictionary* funcDictPtr = &funcsDict;

    if (funcsDict.found(funcType) && funcsDict.isDict(funcType))
    {
        funcDictPtr = &funcsDict.subDict(funcType);
    }

    dictionary& funcDict = *funcDictPtr;

    // Store the funcDict as read for error reporting context
    const dictionary funcDict0(funcDict);

    // Insert the 'field' and/or 'fields' and 'objects' entries corresponding
    // to both the arguments and the named arguments
    DynamicList<wordAndDictionary> fieldArgs;
    forAll(args, i)
    {
        fieldArgs.append(wordAndDictionary(args[i], dictionary::null));
    }
    forAll(namedArgs, i)
    {
        if (namedArgs[i].first() == "field")
        {
            IStringStream iss(namedArgs[i].second());
            fieldArgs.append(wordAndDictionary(iss));
        }
        if
        (
            namedArgs[i].first() == "fields"
         || namedArgs[i].first() == "objects"
        )
        {
            IStringStream iss(namedArgs[i].second());
            fieldArgs.append(List<wordAndDictionary>(iss));
        }
    }
    if (fieldArgs.size() == 1)
    {
        funcDict.set("field", fieldArgs[0].first());
        funcDict.merge(fieldArgs[0].second());
    }
    if (fieldArgs.size() >= 1)
    {
        funcDict.set("fields", fieldArgs);
        funcDict.set("objects", fieldArgs);
    }

    // Insert non-field arguments
    forAll(namedArgs, i)
    {
        if
        (
            namedArgs[i].first() != "field"
         && namedArgs[i].first() != "fields"
         && namedArgs[i].first() != "objects"
        )
        {
            const Pair<word> dAk(dictAndKeyword(namedArgs[i].first()));
            dictionary& subDict(funcDict.scopedDict(dAk.first()));
            IStringStream entryStream
            (
                dAk.second() + ' ' + namedArgs[i].second() + ';'
            );
            subDict.set(entry::New(entryStream).ptr());
        }
    }

    // Insert the region name if specified
    if (region != word::null)
    {
        funcDict.set("region", region);
    }

    // Set the name of the function entry to that specified by the optional
    // funcName argument otherwise automatically generate a unique name
    // from the function type and arguments
    const word funcName
    (
        funcDict.lookupOrDefault("funcName", string::validate<word>(funcArgs))
    );

    // Check for anything in the configuration that has not been set
    List<Tuple2<word, string>> unsetArgs = unsetEntries(funcDict);
    bool hasUnsetError = false;
    forAll(unsetArgs, i)
    {
        if
        (
            unsetArgs[i].first() != "fields"
         && unsetArgs[i].first() != "objects"
        )
        {
            hasUnsetError = true;
        }
    }
    if (!hasUnsetError)
    {
        forAll(unsetArgs, i)
        {
            funcDict.set(unsetArgs[i].first(), wordList());
        }
    }
    else
    {
        FatalIOErrorInFunction(funcDict0)
            << nl;

        forAll(unsetArgs, i)
        {
            FatalIOErrorInFunction(funcDict0)
                << "Essential value for keyword '" << unsetArgs[i].first()
                << "' not set" << nl;
        }

        FatalIOErrorInFunction(funcDict0)
            << nl << "In function entry:" << nl
            << "    " << funcArgs.c_str() << nl
            << nl << "In " << contextTypeAndValue.first().c_str() << ":" << nl
            << "    " << contextTypeAndValue.second().c_str() << nl;

        word funcType;
        wordReList args;
        List<Tuple2<word, string>> namedArgs;
        dictArgList(funcArgs, funcType, args, namedArgs);

        string argList;
        forAll(args, i)
        {
            args[i].strip(" \n");
            argList += (argList.size() ? ", " : "") + args[i];
        }
        forAll(namedArgs, i)
        {
            namedArgs[i].second().strip(" \n");
            argList +=
                (argList.size() ? ", " : "")
              + namedArgs[i].first() + " = " + namedArgs[i].second();
        }
        forAll(unsetArgs, i)
        {
            unsetArgs[i].second().strip(" \n");
            argList +=
                (argList.size() ? ", " : "")
              + unsetArgs[i].first() + " = " + unsetArgs[i].second();
        }

        FatalIOErrorInFunction(funcDict0)
            << nl << "The function entry should be:" << nl
            << "    " << funcType << '(' << argList.c_str() << ')'
            << exit(FatalIOError);
    }

    // Re-parse the funcDict to execute the functionEntries
    // now that the function argument entries have been added
    dictionary funcArgsDict;
    funcArgsDict.add(funcName, funcDict);
    {
        OStringStream os;
        funcArgsDict.write(os);
        funcArgsDict = dictionary
        (
            funcType,
            functionsDict,
            IStringStream(os.str())()
        );
    }

    // Merge this functionObject dictionary into functionsDict
    functionsDict.merge(funcArgsDict);
    functionsDict.subDict(funcName).name() = funcDict.name();

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
            readFunctionObject
            (
                args["func"],
                functionsDict,
                {"command", args.commandLine()},
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


Foam::scalar Foam::functionObjectList::timeToNextWrite()
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
            result = min(result, operator[](oi).timeToNextWrite());
        }
    }

    return result;
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
        forAll(*this, oi)
        {
            operator[](oi).updateMesh(mpm);
        }
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


// ************************************************************************* //
