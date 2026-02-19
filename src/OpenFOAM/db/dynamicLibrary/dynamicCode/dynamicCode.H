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

Class
    Foam::dynamicCode

Description
    Encapsulation of dynamic code dictionaries and functionality

SourceFiles
    dynamicCode.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicCode_H
#define dynamicCode_H

#include "dictionary.H"
#include "HashTable.H"
#include "SHA1Digest.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class dynamicCode Declaration
\*---------------------------------------------------------------------------*/

class dynamicCode
{
    // Private Data

        //- Top-level directory name for copy/compiling
        static const word topDirName;

        //- Name of the code template sub-directory
        //  Used when locating the codeTemplateName via Foam::findEtcFile
        static const fileName codeTemplateDirName;

        //- Root of the LIB target for Make/files
        static const char* const libTargetRoot;

        //- Root for dynamic code compilation
        fileName codeRoot_;

        //- Subdirectory name for loading libraries
        const fileName libSubDir_;

        //- Name for code
        word codeName_;

        //- Unique name for code with the SHA1 appendedm
        word codeSha1Name_;

        //- Name for code subdirectory
        word codeDirName_;

        //- Keywords associated with source code
        const wordList& codeKeys_;

        //- Name of the dictionary variables in the source code
        const wordList& codeDictVars_;

        //- Options file name
        const word& optionsFileName_;

        //- Files to copy and filter and compile
        const wordList& compileFiles_;

        //- Files to copy and filter
        const wordList& copyFiles_;

        //- Code strings
        //  Stored to write for restart
        List<verbatimString> codeStrings_;

        //- Variable substitutions
        HashTable<string> varSubstitutions_;

        //- Code key substitutions
        HashTable<string> codeKeySubstitutions_;

        //- Optional "codeOptions" string
        //  Stored to write for restart
        verbatimString optionsString_;

        //- Optional "codeOptions" entry
        string options_;

        //- Optional "codeLibs" string
        //  Stored to write for restart
        verbatimString libsString_;

        //- Optional "codeLibs" entry
        string libs_;

        //- Calculated SHA1Digest
        SHA1Digest sha1_;


    // Private Member Functions

        //- Path for SHA1Digest
        //  Corresponds to codePath()/Make/SHA1Digest
        fileName digestFile() const
        {
            return codeRoot_/codeDirName_/"Make/SHA1Digest";
        }

        //- Path for specified code name
        //  Corresponds to codeRoot()/codeDirName()
        fileName codePath() const
        {
            return codeRoot_/codeDirName_;
        }

        //- Path for specified code name relative to \$FOAM_CASE
        //  Corresponds to topDirName/codeDirName()
        fileName codeRelPath() const
        {
            return topDirName/codeDirName_;
        }

        //- Add a \#line directive to the start of the given source string that
        //  compilation messages are meaningful
        static void addLineDirective
        (
            string&,
            const label lineNum,
            const fileName& name
        );

        //- Copy lines while expanding variables
        static void copyAndFilter
        (
            ISstream&,
            OSstream&,
            const HashTable<string>& mapping
        );

        //- Resolve code-templates via Foam::findConfigFile
        static bool resolveTemplates
        (
            const wordList& templateNames,
            DynamicList<fileName>& resolvedFiles,
            DynamicList<fileName>& badFiles
        );

        //- Copy/create Make/files prior to compilation
        bool createMakeFiles() const;

        //- Copy/create Make/options prior to compilation
        bool createMakeOptions() const;

        //- Write digest to Make/SHA1Digest
        bool writeDigest() const;


public:

    // Static Data Members

        //- Flag if system operations are allowed
        static int allowSystemOperations;


    // Constructors

        //- Construct from the context and code dictionaries
        //  and lists of which entries correspond to code
        dynamicCode
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
        );

        //- Construct from the context dictionary also containing the code
        //  and lists of which entries correspond
        dynamicCode
        (
            const dictionary& contextDict,
            const word& codeName,
            const word& codeDirName,
            const wordList& codeKeys,
            const wordList& codeDictVars,
            const word& optionsFileName,
            const wordList& compileFiles,
            const wordList& copyFiles
        );


    // Member Functions

        //- Return the code-name
        const word& codeName() const
        {
            return codeName_;
        }

        //- Return the unique code name
        const word& codeSha1Name() const
        {
            return codeSha1Name_;
        }

        //- Library path for specified code name
        //  Corresponds to codeRoot()/libSubDir()/lib\<codeSha1Name\>.so
        fileName libPath() const
        {
            return codeRoot_/libSubDir_/"lib" + codeSha1Name_ + ".so";
        }

        //- Library path for specified code name relative to \$FOAM_CASE
        //  Corresponds to
        //  dynamicCode/codeDirName()/libSubDir()/lib\<codeSha1Name\>.so
        fileName libRelPath() const
        {
            return codeRelPath()/libSubDir_/"lib" + codeSha1Name_ + ".so";
        }

        HashTable<string>& varSubstitutions()
        {
            return varSubstitutions_;
        }

        //- Return the library basename without leading 'lib' or trailing '.so'
        static word libraryBaseName(const fileName& libPath);

        //- Resolve code-template via Foam::findConfigFile
        static fileName resolveTemplate(const fileName& templateName);

        //- Verify if the copied code is up-to-date, based on Make/SHA1Digest
        bool upToDate() const;

        //- Copy/create files prior to compilation
        bool copyOrCreateFiles(const bool verbose = false) const;

        //- Compile a libso
        bool wmakeLibso() const;

        void read(const dictionary& contextDict, const dictionary& codeDict);

        void read(const dictionary& contextDict);

        //- Write the code for restart
        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
