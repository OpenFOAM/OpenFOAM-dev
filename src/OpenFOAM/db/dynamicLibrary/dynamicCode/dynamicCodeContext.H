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
    Foam::dynamicCodeContext

Description
    Encapsulation of dynamic code dictionaries

SourceFiles
    dynamicCodeContext.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicCodeContext_H
#define dynamicCodeContext_H

#include "dictionary.H"
#include "HashTable.H"
#include "SHA1Digest.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dynamicCode;

/*---------------------------------------------------------------------------*\
                     Class dynamicCodeContext Declaration
\*---------------------------------------------------------------------------*/

class dynamicCodeContext
{
    // Private Data

        //- Top-level directory name for copy/compiling
        static const word topDirName;

        //- Name of the code template sub-directory
        //  Used when locating the codeTemplateName via Foam::findEtcFile
        static const fileName codeTemplateDirName;

        //- Root of the LIB target for Make/files
        static const char* const libTargetRoot;

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

        //- Variables to use during filtering
        HashTable<string> filterVars_;

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


public:

    friend class dynamicCode;


    // Static Data Members

        //- Flag if system operations are allowed
        static int allowSystemOperations;


    // Constructors

        //- Construct from the context and code dictionaries
        //  and lists of which entries correspond to code
        dynamicCodeContext
        (
            const dictionary& contextDict,
            const dictionary& codeDict,
            const wordList& codeKeys,
            const wordList& codeDictVars,
            const word& optionsFileName,
            const wordList& compileFiles,
            const wordList& copyFiles
        );

        //- Construct from the context dictionary also containing the code
        //  and lists of which entries correspond
        dynamicCodeContext
        (
            const dictionary& contextDict,
            const wordList& codeKeys,
            const wordList& codeDictVars,
            const word& optionsFileName,
            const wordList& compileFiles,
            const wordList& copyFiles
        );


    // Static Member functions

        //- Return the library basename without leading 'lib' or trailing '.so'
        static word libraryBaseName(const fileName& libPath);


    // Member Functions

        //- Variables to use during filtering
        const HashTable<string>& filterVars() const
        {
            return filterVars_;
        }

        //- Return SHA1 digest calculated from include, options, code
        const SHA1Digest& sha1() const
        {
            return sha1_;
        }

        //- Resolve code-template via Foam::findConfigFile
        static fileName resolveTemplate(const fileName& templateName);

        void read
        (
            const dictionary& contextDict,
            const dictionary& codeDict
        );

        void read(const dictionary& contextDict);

        //- Write the code for restart
        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
