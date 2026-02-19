/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "dynamicCode.H"
#include "stringOps.H"
#include "IFstream.H"
#include "OFstream.H"
#include "OSHA1stream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::dynamicCode::allowSystemOperations
(
    Foam::debug::infoSwitch("allowSystemOperations", 0)
);

const Foam::fileName Foam::dynamicCode::codeTemplateDirName
(
    "codeTemplates/dynamicCode"
);

const Foam::word Foam::dynamicCode::topDirName
(
    "dynamicCode"
);

const char* const Foam::dynamicCode::libTargetRoot
(
    "LIB = $(PWD)/../platforms/$(WM_OPTIONS)/lib/lib"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::dynamicCode::addLineDirective
(
    string& code,
    const label lineNum,
    const fileName& name
)
{
    code = "#line " + Foam::name(lineNum) + " \"" + name + "\"\n" + code;
}


void Foam::dynamicCode::copyAndFilter
(
    ISstream& is,
    OSstream& os,
    const HashTable<string>& mapping
)
{
    if (!is.good())
    {
        FatalErrorInFunction
            << "Failed opening for reading " << is.name()
            << exit(FatalError);
    }

    if (!os.good())
    {
        FatalErrorInFunction
            << "Failed writing " << os.name()
            << exit(FatalError);
    }

    // Copy file while rewriting $VARS and ${VARS}
    string line;
    do
    {
        // Read the next line without continuation
        is.getLine(line, false);

        // Expand according to mapping.
        // Expanding according to env variables might cause too many
        // surprises
        stringOps::inplaceExpandCodeTemplate(line, mapping);
        os.writeQuoted(line, false) << nl;
    }
    while (is.good());
}


bool Foam::dynamicCode::resolveTemplates
(
    const wordList& templateNames,
    DynamicList<fileName>& resolvedFiles,
    DynamicList<fileName>& badFiles
)
{
    bool allOkay = true;
    forAll(templateNames, fileI)
    {
        const fileName& templateName = templateNames[fileI];

        const fileName file
        (
            findConfigFile
            (
                templateName,
                dynamicCode::codeTemplateDirName,
                "system"
            )
        );

        if (file.empty())
        {
            badFiles.append(templateName);
            allOkay = false;
        }
        else
        {
            resolvedFiles.append(file);
        }
    }

    return allOkay;
}


bool Foam::dynamicCode::createMakeFiles() const
{
    // Create Make/files
    if (compileFiles_.empty())
    {
        return false;
    }

    const fileName dstFile(codePath()/"Make/files");

    // Create dir
    mkDir(dstFile.path());

    OFstream os(dstFile);

    if (!os.good())
    {
        FatalErrorInFunction
                << "Failed writing " << dstFile
                << exit(FatalError);
    }

    // Write compile files
    forAll(compileFiles_, fileI)
    {
        os.writeQuoted(compileFiles_[fileI], false) << nl;
    }

    os  << nl << dynamicCode::libTargetRoot << codeSha1Name_ << nl;

    return true;
}


bool Foam::dynamicCode::createMakeOptions() const
{
    if (compileFiles_.empty())
    {
        return false;
    }

    // Read the options template file
    const fileName optionsFile
    (
        dynamicCode::resolveTemplate(optionsFileName_)
    );

    verbatimString options;
    verbatimString libs;

    if (!optionsFileName_.empty())
    {
        const fileName optionsFile
        (
            dynamicCode::resolveTemplate(optionsFileName_)
        );

        if (!optionsFile.empty())
        {
            IFstream is(optionsFile);
            if (!is.good())
            {
                FatalErrorInFunction
                    << "Failed opening " << optionsFile
                    << exit(FatalError);
            }

            dictionary optionsDict(is);

            options = optionsDict.lookupOrDefault<verbatimString>
            (
                "codeOptions",
                verbatimString::null
            );

            libs = optionsDict.lookupOrDefault<verbatimString>
            (
                "codeLibs",
                verbatimString::null
            );
        }
        else
        {
            FatalErrorInFunction
                << "Cannot find options file " << optionsFileName_
                << exit(FatalError);
        }
    }

    // Add the code specific options and libs
    options += options_;
    libs += libs_;

    if (options.empty() && libs.empty())
    {
        return false;
    }

    const fileName dstFile(codePath()/"Make/options");
    mkDir(dstFile.path());

    OFstream os(dstFile);

    if (!os.good())
    {
        FatalErrorInFunction
                << "Failed writing " << dstFile
                << exit(FatalError);
    }

    os.writeQuoted(options + "\n\n" + libs, false) << nl;

    return true;
}


bool Foam::dynamicCode::writeDigest() const
{
    const fileName file(digestFile());
    mkDir(file.path());

    OFstream os(file);
    os  << '_';
    os.writeQuoted(sha1_.str(), false) << nl;

    return os.good();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicCode::dynamicCode
(
    const dictionary& contextDict,
    const dictionary& codeDict,
    const word& codeName,
    const word& codeDirName,
    const wordList& codeKeys,
    const wordList& codeDictVars,
    const word& optionsFileName,
    const wordList& compileFiles,
    const wordList& copyFiles
)
:
    codeRoot_
    (
        stringOps::expandEnvVar("$FOAM_CASE")/topDirName
    ),
    libSubDir_(stringOps::expandEnvVar("platforms/$WM_OPTIONS/lib")),
    codeName_(codeName),
    codeKeys_(codeKeys),
    codeDictVars_(codeDictVars),
    optionsFileName_(optionsFileName),
    compileFiles_(compileFiles),
    copyFiles_(copyFiles),
    codeStrings_(codeKeys.size())
{
    if (isAdministrator())
    {
        FatalIOErrorInFunction(contextDict)
            << "This code should not be executed by someone with administrator"
            << " rights due to security reasons." << nl
            << "(it writes a shared library which then gets loaded "
            << "using dlopen)"
            << exit(FatalIOError);
    }

    if (!allowSystemOperations)
    {
        FatalIOErrorInFunction(contextDict)
            << "Loading a shared library using case-supplied code is not"
            << " enabled by default" << nl
            << "because of security issues. If you trust the code you can"
            << " enable this" << nl
            << "facility be adding to the InfoSwitches setting in the system"
            << " controlDict:" << nl << nl
            << "    allowSystemOperations 1" << nl << nl
            << "The system controlDict is either" << nl << nl
            << "    ~/.OpenFOAM/$WM_PROJECT_VERSION/controlDict" << nl << nl
            << "or" << nl << nl
            << "    $WM_PROJECT_DIR/etc/controlDict" << nl
            << endl
            << exit(FatalIOError);
    }

    read(contextDict, codeDict);

    const word sha1Str(sha1_.str());

    codeSha1Name_ = codeName_ + '_' + sha1Str;

    codeDirName_ =
    (
        codeDirName.empty()
      ? word('_' + sha1Str)
      : codeDirName
    );

    varSubstitutions_.set("typeName", codeName_);
    varSubstitutions_.set("uniqueFunctionName", codeSha1Name_);
    varSubstitutions_.set("SHA1sum", sha1Str);
}


Foam::dynamicCode::dynamicCode
(
    const dictionary& contextDict,
    const word& codeName,
    const word& codeDirName,
    const wordList& codeKeys,
    const wordList& codeDictVars,
    const word& codeOptionsFileName,
    const wordList& compileFiles,
    const wordList& copyFiles
)
:
    dynamicCode
    (
        contextDict,
        contextDict,
        codeName,
        codeDirName,
        codeKeys,
        codeDictVars,
        codeOptionsFileName,
        compileFiles,
        copyFiles
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dynamicCode::read
(
    const dictionary& contextDict,
    const dictionary& codeDict
)
{
    // Expand all dictionary entries. Note that this removes any leading or
    // trailing whitespace, which is necessary for compilation options, and
    // doesn't hurt for everything else
    List<const entry*> codePtrs(codeKeys_.size(), nullptr);
    codeKeySubstitutions_.clear();
    forAll(codeKeys_, i)
    {
        const word& key = codeKeys_[i];
        codePtrs[i] = codeDict.lookupEntryPtr(key, false, false);
        if (codePtrs[i])
        {
            codeStrings_[i] = verbatimString(codePtrs[i]->stream());
            stringOps::inplaceExpandCodeString
            (
                codeStrings_[i],
                contextDict, // Lookup variables from the context dictionary
                codeDictVars_[i]
            );
            codeKeySubstitutions_.insert(key, codeStrings_[i]);
        }
        else
        {
            codeKeySubstitutions_.insert(key, "");
        }
    }

    // Code options
    const entry* optionsPtr =
        codeDict.lookupEntryPtr("codeOptions", false, false);
    if (optionsPtr)
    {
        optionsString_ = verbatimString(optionsPtr->stream());
        stringOps::inplaceExpandCodeString
        (
            optionsString_,
            contextDict,
            word::null
        );
        options_ = stringOps::trim(optionsString_);
    }

    // Code libs
    const entry* libsPtr = codeDict.lookupEntryPtr("codeLibs", false, false);
    if (libsPtr)
    {
        libsString_ = verbatimString(libsPtr->stream());
        stringOps::inplaceExpandCodeString
        (
            libsString_,
            contextDict,
            word::null
        );
        libs_ = stringOps::trim(libsString_);
    }

    // Calculate SHA1 digest from all entries
    OSHA1stream os;
    forAllConstIter(HashTable<string>, codeKeySubstitutions_, iter)
    {
        os << iter();
    }
    os << options_ << libs_;
    sha1_ = os.digest();

    // Add line directives after calculating SHA1
    forAll(codeKeys_, i)
    {
        if (codePtrs[i])
        {
            const word& key = codeKeys_[i];
            addLineDirective
            (
                codeKeySubstitutions_[key],
                codePtrs[i]->startLineNumber(),
                codeDict.name()
            );
        }
    }
}


Foam::word Foam::dynamicCode::libraryBaseName(const fileName& libPath)
{
    word libName(libPath.name(true));
    libName.erase(0, 3);    // Remove leading 'lib' from name
    return libName;
}


Foam::fileName Foam::dynamicCode::resolveTemplate
(
    const fileName& templateName
)
{
    return findConfigFile
    (
        templateName,
        codeTemplateDirName,
        "system"
    );
}


bool Foam::dynamicCode::copyOrCreateFiles(const bool verbose) const
{
    if (verbose)
    {
        Info<< "Creating new library in " << libRelPath() << endl;
    }

    HashTable<string> filterVars(varSubstitutions_);

    // Collect all the filter variables
    forAllConstIter(HashTable<string>, codeKeySubstitutions_, iter)
    {
        filterVars.set(iter.key(), iter());
    }

    const label nFiles =
        compileFiles_.size() + copyFiles_.size();

    DynamicList<fileName> resolvedFiles(nFiles);
    DynamicList<fileName> badFiles(nFiles);

    // Resolve template, or add to bad-files
    dynamicCode::resolveTemplates
    (
        compileFiles_,
        resolvedFiles,
        badFiles
    );
    dynamicCode::resolveTemplates
    (
        copyFiles_,
        resolvedFiles,
        badFiles
    );

    if (!badFiles.empty())
    {
        FatalErrorInFunction
            << "Could not find the code template(s): "
            << badFiles << nl
            << exit(FatalError);
    }

    // Create dir
    const fileName outputDir(codePath());

    // Create dir
    mkDir(outputDir);

    // Copy/filter files
    forAll(resolvedFiles, fileI)
    {
        const fileName& srcFile = resolvedFiles[fileI];
        const fileName dstFile(outputDir/srcFile.name());

        if (verbose)
        {
            Info << "    Copying " << srcFile << " to " << dstFile << endl;
        }

        IFstream is(srcFile);
        if (!is.good())
        {
            FatalErrorInFunction
                << "Failed opening " << srcFile
                << exit(FatalError);
        }

        OFstream os(dstFile);
        if (!os.good())
        {
            FatalErrorInFunction
                << "Failed writing " << dstFile
                << exit(FatalError);
        }

        // Copy lines while expanding variables
        dynamicCode::copyAndFilter(is, os, filterVars);
    }


    // Create Make/files + Make/options
    createMakeFiles();
    createMakeOptions();

    writeDigest();

    return true;
}


bool Foam::dynamicCode::wmakeLibso() const
{
    const string wmakeCmd("wmake -s libso " + codePath());

    if (system(wmakeCmd))
    {
        return false;
    }
    else
    {
        return true;
    }
}


bool Foam::dynamicCode::upToDate() const
{
    const fileName file(digestFile());

    if (!exists(file, false, true) || SHA1Digest(IFstream(file)()) != sha1_)
    {
        return false;
    }

    return true;
}


void Foam::dynamicCode::read(const dictionary& contextDict)
{
    read(contextDict, contextDict);
}


void Foam::dynamicCode::write(Ostream& os) const
{
    forAll(codeStrings_, i)
    {
        if (codeStrings_[i] != verbatimString::null)
        {
            writeEntry(os, codeKeys_[i], codeStrings_[i]);
        }
    }

    if (optionsString_ != verbatimString::null)
    {
        writeEntry(os, "codeOptions", optionsString_);
    }

    if (libsString_ != verbatimString::null)
    {
        writeEntry(os, "codeLibs", libsString_);
    }
}


// ************************************************************************* //
