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
    Extension of dynamicCodeContext providing code compilation

SourceFiles
    dynamicCode.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicCode_H
#define dynamicCode_H

#include "HashTable.H"
#include "DynamicList.H"
#include "dynamicCodeContext.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class ISstream;
class OSstream;

/*---------------------------------------------------------------------------*\
                         Class dynamicCode Declaration
\*---------------------------------------------------------------------------*/

class dynamicCode
{
    // Private Data

        const dynamicCodeContext& context_;

        //- Root for dynamic code compilation
        fileName codeRoot_;

        //- Subdirectory name for loading libraries
        const fileName libSubDir_;

        //- Name for code
        word codeName_;

        //- Name for code subdirectory
        word codeDirName_;

        //- Variables to use during filtering
        HashTable<string> filterVars_;


    // Private Member Functions

        //- Write SHA1 value as C-comment
        bool writeCommentSHA1(Ostream&) const;

        //- Copy/create Make/files prior to compilation
        bool createMakeFiles() const;

        //- Copy/create Make/options prior to compilation
        bool createMakeOptions() const;

        //- Write digest to Make/SHA1Digest
        bool writeDigest(const std::string&) const;

        //- Verify if the copied code is up-to-date, based on Make/SHA1Digest
        bool upToDate(const SHA1Digest& sha1) const;

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
        fileName codeRelPath() const;


public:

    // Constructors

        //- Construct for a specified code name and code directory name
        //  Defaults to using the sha1 for the code directory name
        dynamicCode
        (
            const dynamicCodeContext& context,
            const word& codeName,
            const word& codeDirName = word::null
        );

        //- Disallow default bitwise copy construction
        dynamicCode(const dynamicCode&) = delete;


    // Member Functions

        //- Return the code-name
        const word& codeName() const
        {
            return codeName_;
        }

        //- Library path for specified code name
        //  Corresponds to codeRoot()/libSubDir()/lib\<codeName\>.so
        fileName libPath() const
        {
            return codeRoot_/libSubDir_/"lib" + codeName_ + ".so";
        }

        //- Library path for specified code name relative to \$FOAM_CASE
        //  Corresponds to
        //  dynamicCode/codeDirName()/libSubDir()/lib\<codeName\>.so
        fileName libRelPath() const;

        //- Define a filter variable
        void setFilterVariable(const word& key, const std::string& value);

        //- Return SHA1 digest calculated from include, options, code
        const SHA1Digest& sha1() const
        {
            return context_.sha1();
        }

        //- Verify if the copied code is up-to-date, based on Make/SHA1Digest
        bool upToDate() const;

        //- Copy/create files prior to compilation
        bool copyOrCreateFiles(const bool verbose = false) const;

        //- Compile a libso
        bool wmakeLibso() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const dynamicCode&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
