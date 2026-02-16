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
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::fileName Foam::dynamicCode::codeTemplateDirName
(
    "codeTemplates/dynamicCode"
);

const char* const Foam::dynamicCode::libTargetRoot
(
    "LIB = $(PWD)/../platforms/$(WM_OPTIONS)/lib/lib"
);

const Foam::word Foam::dynamicCode::topDirName("dynamicCode");


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::word Foam::dynamicCode::libraryBaseName(const fileName& libPath)
{
    word libName(libPath.name(true));
    libName.erase(0, 3);    // Remove leading 'lib' from name
    return libName;
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

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


bool Foam::dynamicCode::resolveTemplates
(
    const UList<fileName>& templateNames,
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
                codeTemplateDirName,
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


bool Foam::dynamicCode::writeCommentSHA1(Ostream& os) const
{
    const bool hasSHA1 = filterVars_.found("SHA1sum");

    if (hasSHA1)
    {
        os  << "# dynamicCode:\n# SHA1 = ";
        os.writeQuoted(filterVars_["SHA1sum"], false) << "\n\n";
    }

    return hasSHA1;
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

    writeCommentSHA1(os);

    // Write compile files
    forAll(compileFiles_, fileI)
    {
        os.writeQuoted(compileFiles_[fileI].name(), false) << nl;
    }

    os  << nl
        << libTargetRoot << codeName_.c_str() << nl;

    return true;
}


bool Foam::dynamicCode::createMakeOptions() const
{
    if (compileFiles_.empty())
    {
        return false;
    }

    const fileName optionsFile(resolveTemplate(codeOptionsFileName_));

    verbatimString codeOptions;
    verbatimString codeLibs;

    if (!codeOptionsFileName_.empty())
    {
        const fileName optionsFile(resolveTemplate(codeOptionsFileName_));

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

            codeOptions = optionsDict.lookupOrDefault<verbatimString>
            (
                "codeOptions",
                verbatimString::null
            );

            codeLibs = optionsDict.lookupOrDefault<verbatimString>
            (
                "codeLibs",
                verbatimString::null
            );
        }
        else
        {
            FatalErrorInFunction
                << "Cannot find options file " << codeOptionsFileName_
                << exit(FatalError);
        }
    }

    if
    (
        makeOptions_.empty() && codeOptions.empty()
     && makeLibs_.empty() && codeLibs.empty()
    )
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

    writeCommentSHA1(os);

    os.writeQuoted
    (
        codeOptions + makeOptions_ + "\n\n" + codeLibs + makeLibs_,
        false
    ) << nl;

    return true;
}


bool Foam::dynamicCode::writeDigest(const std::string& sha1) const
{
    const fileName file(digestFile());
    mkDir(file.path());

    OFstream os(file);
    os  << '_';
    os.writeQuoted(sha1, false) << nl;

    return os.good();
}


bool Foam::dynamicCode::upToDate(const SHA1Digest& sha1) const
{
    const fileName file(digestFile());

    if (!exists(file, false, true) || SHA1Digest(IFstream(file)()) != sha1)
    {
        return false;
    }

    return true;
}


Foam::fileName Foam::dynamicCode::codeRelPath() const
{
    return topDirName/codeDirName_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicCode::dynamicCode
(
    const dynamicCodeContext& context,
    const word& codeName,
    const word& codeDirName
)
:
    context_(context),
    codeRoot_(stringOps::expandEnvVar("$FOAM_CASE")/topDirName),
    libSubDir_(stringOps::expandEnvVar("platforms/$WM_OPTIONS/lib")),
    codeName_(codeName + context.sha1().str(true)),
    codeDirName_
    (
        codeDirName.empty()
      ? word(context.sha1().str(true))
      : codeDirName
    ),
    filterVars_
    {
        {"typeName", codeName_},
        {"SHA1sum", SHA1Digest().str()}
    }
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::dynamicCode::libRelPath() const
{
    return codeRelPath()/libSubDir_/"lib" + codeName_ + ".so";
}


void Foam::dynamicCode::filter()
{
    compileFiles_.clear();
    copyFiles_.clear();

    forAllConstIter(HashTable<string>, context_.code(), iter)
    {
        setFilterVariable(iter.key(), iter());
    }

    codeOptionsFileName_ = context_.codeOptionsFileName_;
    makeOptions_ = context_.options_;
    makeLibs_ = context_.libs_;

    setFilterVariable("SHA1sum", context_.sha1().str());
}


void Foam::dynamicCode::addCompileFile(const fileName& name)
{
    compileFiles_.append(name);
}


void Foam::dynamicCode::addCopyFile(const fileName& name)
{
    copyFiles_.append(name);
}


void Foam::dynamicCode::setFilterVariable
(
    const word& key,
    const std::string& value
)
{
    filterVars_.set(key, value);
}


bool Foam::dynamicCode::copyOrCreateFiles(const bool verbose) const
{
    if (verbose)
    {
        Info<< "Creating new library in " << libRelPath() << endl;
    }

    const label nFiles = compileFiles_.size() + copyFiles_.size();

    DynamicList<fileName> resolvedFiles(nFiles);
    DynamicList<fileName> badFiles(nFiles);

    // Resolve template, or add to bad-files
    resolveTemplates(compileFiles_, resolvedFiles, badFiles);
    resolveTemplates(copyFiles_, resolvedFiles, badFiles);

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
        copyAndFilter(is, os, filterVars_);
    }


    // Create Make/files + Make/options
    createMakeFiles();
    createMakeOptions();

    writeDigest(filterVars_["SHA1sum"]);

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
    return upToDate(context_.sha1());
}


// ************************************************************************* //
