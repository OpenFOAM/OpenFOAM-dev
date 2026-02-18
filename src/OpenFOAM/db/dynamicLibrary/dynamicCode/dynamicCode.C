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

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

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
    if (context_.compileFiles_.empty())
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
    forAll(context_.compileFiles_, fileI)
    {
        os.writeQuoted(context_.compileFiles_[fileI], false) << nl;
    }

    os  << nl << dynamicCodeContext::libTargetRoot << codeName_ << nl;

    return true;
}


bool Foam::dynamicCode::createMakeOptions() const
{
    if (context_.compileFiles_.empty())
    {
        return false;
    }

    // Read the options template file
    const fileName optionsFile
    (
        dynamicCodeContext::resolveTemplate(context_.optionsFileName_)
    );

    verbatimString options;
    verbatimString libs;

    if (!context_.optionsFileName_.empty())
    {
        const fileName optionsFile
        (
            dynamicCodeContext::resolveTemplate(context_.optionsFileName_)
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
                << "Cannot find options file " << context_.optionsFileName_
                << exit(FatalError);
        }
    }

    // Add the code specific options and libs
    options += context_.options_;
    libs += context_.libs_;

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

    writeCommentSHA1(os);

    os.writeQuoted(options + "\n\n" + libs, false) << nl;

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
    return dynamicCodeContext::topDirName/codeDirName_;
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
    codeRoot_
    (
        stringOps::expandEnvVar("$FOAM_CASE")/dynamicCodeContext::topDirName
    ),
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

    HashTable<string> filterVars(filterVars_);

    // Collect all the filter variables
    forAllConstIter(HashTable<string>, context_.filterVars(), iter)
    {
        filterVars.set(iter.key(), iter());
    }

    filterVars.set("SHA1sum", context_.sha1().str());


    const label nFiles =
        context_.compileFiles_.size() + context_.copyFiles_.size();

    DynamicList<fileName> resolvedFiles(nFiles);
    DynamicList<fileName> badFiles(nFiles);

    // Resolve template, or add to bad-files
    dynamicCodeContext::resolveTemplates
    (
        context_.compileFiles_,
        resolvedFiles,
        badFiles
    );
    dynamicCodeContext::resolveTemplates
    (
        context_.copyFiles_,
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
        dynamicCodeContext::copyAndFilter(is, os, filterVars);
    }


    // Create Make/files + Make/options
    createMakeFiles();
    createMakeOptions();

    writeDigest(filterVars["SHA1sum"]);

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
