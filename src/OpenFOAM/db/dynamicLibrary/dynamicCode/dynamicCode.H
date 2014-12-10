/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Tools for handling dynamic code compilation

SourceFiles
    dynamicCode.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicCode_H
#define dynamicCode_H

#include "Tuple2.H"
#include "HashTable.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class dynamicCodeContext;
class ISstream;
class OSstream;
class SHA1Digest;

/*---------------------------------------------------------------------------*\
                        Class dynamicCode Declaration
\*---------------------------------------------------------------------------*/

class dynamicCode
{
public:
    typedef Tuple2<fileName, string> fileAndContent;

private:
    // Private data

        //- Root for dynamic code compilation
        fileName codeRoot_;

        //- Subdirectory name for loading libraries
        const fileName libSubDir_;

        //- Name for code
        word codeName_;

        //- Name for code subdirectory
        word codeDirName_;

        //- Files to copy and filter
        DynamicList<fileName> compileFiles_;

        //- Files to copy and filter
        DynamicList<fileName> copyFiles_;

        //- Direct contents for files
        DynamicList<fileAndContent> createFiles_;

        //- Variables to use during filtering
        HashTable<string> filterVars_;

        //- Contents for Make/options
        std::string makeOptions_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        dynamicCode(const dynamicCode&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicCode&);


protected:

    // Static data members

        //- Root of the LIB target for Make/files
        static const char* const libTargetRoot;

        //- Top-level directory name for copy/compiling
        static const char* const topDirName;


    // Protected Member Functions

        //- Copy lines while expanding variables
        static void copyAndFilter
        (
            ISstream&,
            OSstream&,
            const HashTable<string>& mapping
        );

        //- Resolve code-templates via the codeTemplateEnvName
        //  alternatively in the codeTemplateDirName via Foam::findEtcFile
        static bool resolveTemplates
        (
            const UList<fileName>& templateNames,
            DynamicList<fileName>& resolvedFiles,
            DynamicList<fileName>& badFiles
        );

        //- Write SHA1 value as C-comment
        bool writeCommentSHA1(Ostream&) const;

        //- Copy/create Make/files prior to compilation
        bool createMakeFiles() const;

        //- Copy/create Make/options prior to compilation
        bool createMakeOptions() const;


        //- Write digest to Make/SHA1Digest
        bool writeDigest(const SHA1Digest&) const;

        //- Write digest to Make/SHA1Digest
        bool writeDigest(const std::string&) const;


public:

    // Static data members

        //- Name of the code template environment variable
        //  Used to located the codeTemplateName
        static const word codeTemplateEnvName;

        //- Name of the code template sub-directory
        //  Used when locating the codeTemplateName via Foam::findEtcFile
        static const fileName codeTemplateDirName;

        //- Flag if system operations are allowed
        static int allowSystemOperations;


    // Static Member functions

        //- Check security for creating dynamic code
        static void checkSecurity(const char* title, const dictionary&);

        //- Return the library basename without leading 'lib' or trailing '.so'
        static word libraryBaseName(const fileName& libPath);


    // Constructors

        //- Construct for a specified code name and code directory name
        //  Defaults to using the code name for the code directory name
        dynamicCode
        (
            const word& codeName,
            const word& codeDirName = ""
        );


    // Member functions

        //- Return the code-name
        const word& codeName() const
        {
            return codeName_;
        }

        //- Return the code-dirname
        const word& codeDirName() const
        {
            return codeDirName_;
        }

        //- Root for dynamic code compilation
        //  Expanded from \$FOAM_CASE/dynamicCode
        const fileName& codeRoot() const
        {
            return codeRoot_;
        }

        //- Subdirectory name for loading libraries
        //  Expanded from platforms/\$WM_OPTIONS/lib
        fileName libSubDir() const
        {
            return libSubDir_;
        }

        //- Path for specified code name
        //  Corresponds to codeRoot()/codeDirName()
        fileName codePath() const
        {
            return codeRoot_/codeDirName_;
        }

        //- Library path for specified code name
        //  Corresponds to codeRoot()/libSubDir()/lib\<codeName\>.so
        fileName libPath() const
        {
            return codeRoot_/libSubDir_/"lib" + codeName_ + ".so";
        }

        //- Path for specified code name relative to \$FOAM_CASE
        //  Corresponds to topDirName/codeDirName()
        fileName codeRelPath() const;


        //- Library path for specified code name relative to \$FOAM_CASE
        //  Corresponds to
        //  dynamicCode/codeDirName()/libSubDir()/lib\<codeName\>.so
        fileName libRelPath() const;


        //- Path for SHA1Digest
        //  Corresponds to codePath()/Make/SHA1Digest
        fileName digestFile() const
        {
            return codeRoot_/codeDirName_/"Make/SHA1Digest";
        }


        //- Clear files and variables
        void clear();

        //- Clear files and reset variables to specified context
        void reset(const dynamicCodeContext&);


        //- Add a file template name, which will be found and filtered
        void addCompileFile(const fileName& name);

        //- Add a file template name, which will be found and filtered
        void addCopyFile(const fileName& name);

        //- Add a file to create with its contents. Will not be filtered
        void addCreateFile(const fileName& name, const string& contents);

        //- Define filter variables for code, codeInclude, SHA1sum
        void setFilterContext(const dynamicCodeContext&);

        //- Define a filter variable
        void setFilterVariable(const word& key, const std::string& value);

        //- Define contents for Make/options
        void setMakeOptions(const std::string& content);


        //- Verify if the copied code is up-to-date, based on Make/SHA1Digest
        bool upToDate(const dynamicCodeContext& context) const;

        //- Verify if the copied code is up-to-date, based on Make/SHA1Digest
        bool upToDate(const SHA1Digest& sha1) const;

        //- Copy/create files prior to compilation
        bool copyOrCreateFiles(const bool verbose = false) const;

        //- Compile a libso
        bool wmakeLibso() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
