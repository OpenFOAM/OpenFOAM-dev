/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::functionEntries::codeDict

Description
    Compiles and executes C++ OpenFOAM code string expressions

    and provides non-const access to the dictionary into which generated
    entries can be inserted or manipulated in any other way.

    \c \#codeDict reads three entries: \c code, \c codeInclude (optional),
    \c codeOptions (optional) to generate the library source code stored in the
    local \c dynamicCode directory with a subdirectory name corresponding to the
    SHA1 of the code.  The code is then compiled into a dynamically loaded
    library libcodeDict_<SHA1>.so stored in the
    \c dynamicCode/platforms/$WM_OPTIONS/lib directory using 'wmake libso'.
    The resulting library is loaded in executed with arguments
    \code
        (dictionary& dict, Istream& is)
    \endcode
    where the dict is the current dictionary and is the stream from which the
    dictionary is being read.

    The verbatim string format \c \#{ ... \c \#} is used to allow multi-line
    input without the need to escape the newlines.

    Dictionary entries constructed with \c \#codeDict can conveniently access
    and use typed variables.  This means calculations involving vectors and
    tensors and list etc. are possible.  To access a variable and construct it
    as a given type within a \c \#codeDict entry, put the type immediately
    after the $ symbol inside angled brackets <>. So, $<vector>var or
    $<vector>{var} substitutes a variable named var as a vector.

    The variable values are looked-up at run-time making the code corresponding
    to each \c \#codeDict independent of the values in the dictionary and of
    each other.  Hence the \c \#codeDict code does not need to be recompiled
    when the values in the dictionary are changed, only if the code is changed.

Usage
    Example to set a zoneGenerator in createZonesDict:
    \verbatim
        inner
        {
            type        union;

            // Create inner baffles along the length of the pipe
            #codeDict
            {
                code
                #{
                    const label nBaffles = $<label>blockMeshDict!nBaffles;
                    for (label i = 1; i < nBaffles; i += 2)
                    {
                        const scalar y =
                            scalar(2*i + 1)/(2*nBaffles)
                           *$<scalar>blockMeshDict!height;

                        dict.set
                        (
                            word("baffle" + name(i)),
                            dictionary::entries
                            (
                                "type", "plane",
                                "point", vector(0, y, 0),
                                "normal", vector(0, 1, 0)
                            )
                        );
                    }
                #};
            }
        }
    \endverbatim

See also
    Foam::functionEntries::codeStream

SourceFiles
    codeDict.C

\*---------------------------------------------------------------------------*/

#ifndef codeDict_H
#define codeDict_H

#include "codeStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dlLibraryTable;

namespace functionEntries
{

// Forward declaration of friend classes
class streamEntry;
class calcEntry;
class codeBlockEntry;
class codeBlockDictEntry;

/*---------------------------------------------------------------------------*\
                         Class codeDict Declaration
\*---------------------------------------------------------------------------*/

class codeDict
:
    public codeStream
{
    //- Interpreter function type
    typedef void (*codeDictFunctionType)(dictionary&, Istream& is);


    // Private Member Functions

        //- Return the code string
        static string codeString
        (
            const label index,
            const dictionary& codeDict,
            Istream&
        );

        //- Construct, compile, load and return streaming function
        codeDictFunctionType getFunction
        (
            const dictionary& contextDict,
            const dictionary& codeDict
        );

        //- Compile and execute the code to manipulate the contextDict
        bool resultStream(dictionary& contextDict, Istream& is);


public:

    // Static Data Members

        //- Name of the C code template to be used
        static const word codeTemplateC;


    // Related types

        //- Declare friendship with the codeBlockEntry class
        friend class codeBlockEntry;
        friend class codeBlockDictEntry;


    //- Runtime type information
    FunctionTypeName("#codeDict");


    // Constructors

        //- Construct from line number, dictionary and Istream
        codeDict
        (
            const label lineNumber,
            const dictionary& parentDict,
            Istream& is
        );

        //- Copy construct
        codeDict(const codeDict&) = default;

        //- Clone
        virtual autoPtr<entry> clone(const dictionary&) const
        {
            return autoPtr<entry>(new codeDict(*this));
        }


    // Member Functions

        //- Expand the functionEntry into the contextDict
        virtual bool execute(dictionary& contextDict, Istream&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const codeDict&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
