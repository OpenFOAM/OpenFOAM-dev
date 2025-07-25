/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "dictionary.H"
#include "IOobject.H"
#include "HashSet.H"
#include "inputModeEntry.H"
#include "calcIncludeEntry.H"
#include "stringOps.H"
#include "etcFiles.H"
#include "wordAndDictionary.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dictionary::dictionary
(
    const fileName& name,
    const dictionary& parentDict,
    Istream& is
)
:
    dictionaryName
    (
        parentDict.name().size()
      ? parentDict.name()/name
      : name
    ),
    parent_(parentDict),
    filePtr_(nullptr)
{
    read(is);
}


Foam::dictionary::dictionary(Istream& is, const bool keepHeader)
:
    dictionaryName(is.name()),
    parent_(dictionary::null),
    filePtr_(nullptr)
{
    // Reset input mode as this is a "top-level" dictionary
    functionEntries::inputModeEntry::clear();

    // Clear the cache of #calc include files
    functionEntries::calcIncludeEntry::clear();

    read(is, keepHeader);
}


Foam::dictionary::includedDictionary::includedDictionary
(
    const fileName& fName,
    const dictionary& parentDict
)
:
    dictionary(fName),
    global_(parentDict.topDict().global())
{
    autoPtr<ISstream> ifsPtr
    (
        fileHandler().NewIFstream(fName)
    );
    ISstream& ifs = ifsPtr();

    if (!ifs || !ifs.good())
    {
        FatalIOErrorInFunction(parentDict)
            << "Included dictionary file " << fName
            << " cannot be found for dictionary " << parentDict.name()
            << exit(FatalIOError);
    }

    read(ifs);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::dictionary> Foam::dictionary::New(Istream& is)
{
    return autoPtr<dictionary>(new dictionary(is));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::dictionary::read(Istream& is, const bool keepHeader)
{
    // Check for empty dictionary
    if (is.eof())
    {
        return true;
    }

    if (!is.good())
    {
        FatalIOErrorInFunction(is)
            << "Istream not OK for reading dictionary "
            << exit(FatalIOError);

        return false;
    }

    // Cache the current name and file/stream pointer
    const fileName name0(name());
    const Istream* filePtr0 = filePtr_;

    // Set the name and file/stream pointer to the given stream
    name() = is.name();
    filePtr_ = &is;

    token currToken(is);
    if (currToken != token::BEGIN_BLOCK)
    {
        is.putBack(currToken);
    }

    while (!is.eof() && entry::New(*this, is))
    {}

    // normally remove the FoamFile header entry if it exists
    if (!keepHeader)
    {
        remove(IOobject::foamFile);
    }

    if (is.bad())
    {
        InfoInFunction
            << "Istream not OK after reading dictionary " << name()
            << endl;

        return false;
    }

    // Reset the name and file/stream pointer to the original
    name() = name0;
    filePtr_ = filePtr0;

    return true;
}


bool Foam::dictionary::global() const
{
    if (&parent_ != &dictionary::null)
    {
        return parent_.global();
    }
    else
    {
        return false;
    }
}


bool Foam::dictionary::substituteKeyword(const word& keyword)
{
    word varName = keyword(1, keyword.size()-1);

    // lookup the variable name in the given dictionary
    const entry* ePtr = lookupEntryPtr(varName, true, true);

    // if defined insert its entries into this dictionary
    if (ePtr != nullptr)
    {
        const dictionary& addDict = ePtr->dict();

        forAllConstIter(IDLList<entry>, addDict, iter)
        {
            add(iter());
        }

        return true;
    }

    return false;
}



// * * * * * * * * * * * * * * Istream Operator  * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, dictionary& dict)
{
    // Reset input mode assuming this is a "top-level" dictionary
    functionEntries::inputModeEntry::clear();

    // Clear the cache of #calc include files
    functionEntries::calcIncludeEntry::clear();

    dict.clear();
    dict.name() = is.name();
    dict.read(is);

    return is;
}


// * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * * //

void Foam::dictionary::write(Ostream& os, bool subDict) const
{
    if (subDict)
    {
        os  << nl << indent << token::BEGIN_BLOCK << incrIndent << nl;
    }

    forAllConstIter(IDLList<entry>, *this, iter)
    {
        const entry& e = *iter;

        // Write entry
        os  << e;

        // Add extra new line between entries for "top-level" dictionaries
        if (!subDict && parent() == dictionary::null && (&e != last()))
        {
            os  << nl;
        }

        // Check stream before going to next entry.
        if (!os.good())
        {
            WarningInFunction
                << "Can't write entry " << iter().keyword()
                << " for dictionary " << name()
                << endl;
        }
    }

    if (subDict)
    {
        os  << decrIndent << indent << token::END_BLOCK << endl;
    }
}


Foam::Ostream& Foam::operator<<(Ostream& os, const dictionary& dict)
{
    dict.write(os, true);
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::List<Foam::Tuple2<Foam::word, Foam::string>>
unsetConfigEntries(const dictionary& configDict)
{
    const wordRe unsetPattern("<.*>");
    unsetPattern.compile();

    List<Tuple2<word, string>> unsetArgs;

    forAllConstIter(IDLList<entry>, configDict, iter)
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
                unsetConfigEntries(iter().dict());

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


void listConfigFiles
(
    const fileName& dir,
    HashSet<word>& foMap
)
{
    // Search specified directory for configuration files
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
            listConfigFiles(dir/foDirs[fd], foMap);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::fileName Foam::findConfigFile
(
    const word& configName,
    const fileName& configFilesPath,
    const word& configFilesDir,
    const word& region
)
{
    // First check if there is a configuration file in the
    // region configFilesDir directory
    {
        const fileName dictFile
        (
            stringOps::expandEnvVar("$FOAM_CASE")
           /configFilesDir/region/configName
        );

        if (isFile(dictFile))
        {
            return dictFile;
        }
    }

    // Next, if the region is specified, check if there is a configuration file
    // in the global configFilesDir directory
    if (region != word::null)
    {
        const fileName dictFile
        (
            stringOps::expandEnvVar("$FOAM_CASE")/configFilesDir/configName
        );

        if (isFile(dictFile))
        {
            return dictFile;
        }
    }

    // Finally, check etc directories
    {
        const fileNameList etcDirs(findEtcDirs(configFilesPath));

        forAll(etcDirs, i)
        {
            const fileName dictFile(search(configName, etcDirs[i]));

            if (!dictFile.empty())
            {
                return dictFile;
            }
        }
    }

    return fileName::null;
}


Foam::wordList Foam::listAllConfigFiles
(
    const fileName& configFilesPath
)
{
    HashSet<word> foMap;

    fileNameList etcDirs(findEtcDirs(configFilesPath));

    forAll(etcDirs, ed)
    {
        listConfigFiles(etcDirs[ed], foMap);
    }

    return foMap.sortedToc();
}


Foam::string Foam::expandArg
(
    const string& arg,
    dictionary& dict,
    const label lineNumber
)
{
    // Add a temporary dummy_ entry to set the arg lineNumber in dict
    dict.set(primitiveEntry("dummy_", token(word("<dummy>"), lineNumber)));

    string expandedArg(arg);
    stringOps::inplaceExpandEntry(expandedArg, dict, true, false);

    // Remove temporary dummy_ entry
    dict.remove("dummy_");

    return expandedArg;
}


void Foam::addArgEntry
(
    dictionary& dict,
    const word& keyword,
    const string& value,
    const label lineNumber
)
{
    IStringStream entryStream(dict.name(), keyword + ' ' + value + ';');
    entryStream.lineNumber() = lineNumber;
    autoPtr<entry> argEntry(entry::New(entryStream));
    if (argEntry.valid())
    {
        dict.set(argEntry.ptr());
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Cannot construct argument entry from string "
            << entryStream.str() << nl
            << "    on line " << lineNumber << " of dictionary " << dict.name()
            << exit(FatalIOError);
    }
}


bool Foam::readConfigFile
(
    const word& configType,
    const Tuple2<string, label>& argStringLine,
    dictionary& parentDict,
    const fileName& configFilesPath,
    const word& configFilesDir,
    const word& region
)
{
    word funcType;
    List<Tuple2<wordRe, label>> args;
    List<Tuple3<word, string, label>> namedArgs;

    dictArgList(argStringLine, funcType, args, namedArgs);

    // Search for the configuration file
    fileName path = findConfigFile
    (
        funcType,
        configFilesPath,
        configFilesDir,
        region
    );

    if (path == fileName::null)
    {
        if (funcType == word::null)
        {
            FatalIOErrorInFunction(parentDict)
                << "configuration file name not specified"
                << nl << nl
                << "Available configured objects:"
                << listAllConfigFiles(configFilesPath)
                << exit(FatalIOError);
        }
        else
        {
            FatalIOErrorInFunction(parentDict)
                << "Cannot find configuration file "
                << funcType << nl << nl
                << "Available configured objects:"
                << listAllConfigFiles(configFilesPath)
                << exit(FatalIOError);
        }

        return false;
    }

    // Read the configuration file
    autoPtr<ISstream> fileStreamPtr(fileHandler().NewIFstream(path));
    ISstream& fileStream = fileStreamPtr();

    // Delay processing the functionEntries
    // until after the argument entries have been added
    entry::disableFunctionEntries = true;
    dictionary funcDict(fileName(funcType), parentDict, fileStream);
    entry::disableFunctionEntries = false;

    // Store the funcDict as read for error reporting context
    const dictionary funcDict0(funcDict);

    // Insert the 'field' and/or 'fields' and 'objects' entries corresponding
    // to both the arguments and the named arguments
    DynamicList<wordAndDictionary> fieldArgs;
    forAll(args, i)
    {
        fieldArgs.append
        (
            wordAndDictionary
            (
                expandArg(args[i].first(), funcDict, args[i].second()),
                dictionary::null
            )
        );
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
         && namedArgs[i].first() != "funcName"
         && namedArgs[i].first() != "name"
        )
        {
            const Pair<word> dAk(dictAndKeyword(namedArgs[i].first()));
            dictionary& subDict(funcDict.scopedDict(dAk.first()));
            addArgEntry
            (
                subDict,
                dAk.second(),
                expandArg
                (
                    namedArgs[i].second(),
                    funcDict,
                    namedArgs[i].third()
                ),
                namedArgs[i].third()
            );
        }
    }

    // Insert the region name if specified
    if (region != word::null)
    {
        funcDict.set("region", region);
    }

    // Set the name of the entry to that specified by the optional
    // name argument otherwise automatically generate a unique name
    // from the type and arguments
    word entryName(funcType);
    if (args.size() || namedArgs.size())
    {
        bool named = false;
        forAll(namedArgs, i)
        {
            if
            (
                namedArgs[i].first() == "funcName"
             || namedArgs[i].first() == "name"
            )
            {
                entryName = expandArg
                (
                    namedArgs[i].second(),
                    funcDict,
                    namedArgs[i].third()
                );
                entryName.strip(" \n");
                named = true;
            }
        }

        if (!named)
        {
            entryName += '(';
            forAll(args, i)
            {
                if (i > 0)
                {
                    entryName += ',';
                }
                entryName += args[i].first();
            }
            forAll(namedArgs, i)
            {
                if (args.size() || i > 0)
                {
                    entryName += ',';
                }
                entryName += namedArgs[i].first();
                entryName += '=';
                entryName += expandArg
                (
                    namedArgs[i].second(),
                    funcDict,
                    namedArgs[i].third()
                );
            }
            entryName += ')';
            string::stripInvalid<word>(entryName);
        }
    }

    // Check for anything in the configuration that has not been set
    List<Tuple2<word, string>> unsetArgs = unsetConfigEntries(funcDict);
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
            << nl << "in " << configType << " entry:" << nl
            << argStringLine.first().c_str() << nl
            << nl << "in dictionary " << parentDict.name().c_str()
            << " starting at line " << argStringLine.second() << nl;

        word funcType;
        List<Tuple2<wordRe, label>> args;
        List<Tuple3<word, string, label>> namedArgs;
        dictArgList(argStringLine, funcType, args, namedArgs);

        string argList;
        forAll(args, i)
        {
            args[i].first().strip(" \n");
            argList += (argList.size() ? ", " : "") + args[i].first();
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
            << nl << "The " << configType << " entry should be:" << nl
            << "    " << funcType << '(' << argList.c_str() << ')'
            << exit(FatalIOError);
    }

    // Re-parse the funcDict to execute the functionEntries
    // now that the argument entries have been added
    dictionary funcArgsDict;
    funcArgsDict.add(entryName, funcDict);

    {
        OStringStream os;
        funcArgsDict.write(os);
        funcArgsDict = dictionary
        (
            funcType,
            funcDict,
            IStringStream(os.str())()
        );
    }

    // Merge this configuration dictionary into parentDict
    parentDict.merge(funcArgsDict);
    parentDict.subDict(entryName).name() = funcDict.name();

    return true;
}


// * * * * * * * * * * * * * * * IOstream Functions  * * * * * * * * * * * * //

void Foam::writeEntry(Ostream& os, const dictionary& value)
{
    os << value;
}


// ************************************************************************* //
