/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2025 OpenFOAM Foundation
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
    Foam::functionEntries::ifeqEntry

Description
    Conditional parsing of and within dictionary entries.

    E.g.
    \verbatim
        a #calc "0.123";
        b 1.23e-1;

        #ifeq ($a $b)
            ..
        #else
            ..
        #endif

        c #ifeq ($a $b) 5.53 #else 7.89 #endif;
    \endverbatim

    \verbatim
    ddtSchemes
    {
    #ifeq (${FOAM_APPLICATION} foamRun)
        default         steadyState;
    #else
        default         Euler;
    #endif
    }
    \endverbatim

    Note:
    - supports dictionary variables and environment variables
    - the two arguments should be tokens
    - the equality test supports any word/string/variable, integers
      and floating point numbers
    - parsing of (non)matching \c \#else, \c \#endif is not very sophisticated

See also
    Foam::functionEntries::ifEntry

SourceFiles
    ifeqEntry.C

\*---------------------------------------------------------------------------*/

#ifndef ifeqEntry_H
#define ifeqEntry_H

#include "functionEntry.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                          Class ifeqEntry Declaration
\*---------------------------------------------------------------------------*/

class ifeqEntry
:
    public functionEntry
{
protected:

    typedef Tuple2<fileName, label> filePos;


private:

    // Private Member Functions

        tokenList readArgList(Istream& is);

        //- Read tokens. Skip dummy tokens
        static void readToken(token& t, Istream& is);

        //- Expand token if it is a variable
        static token expand(const dictionary& dict, const token& t);

        static bool equalToken(const token& t1, const token& t2);

        //- Consume tokens until reached a specific word
        static void skipUntil
        (
            DynamicList<filePos>& stack,
            const dictionary& parentDict,
            const functionName& endWord,
            Istream& is
        );

        template<class Context>
        static bool evaluate
        (
            const bool doIf,
            DynamicList<filePos>& stack,
            const dictionary& contextDict,
            Context& context,
            Istream& is
        );

        //- Main driver: depending on 'equal' starts evaluating or
        //  skips forward to #else
        template<class Context>
        static bool execute
        (
            DynamicList<filePos>& stack,
            const dictionary& contextDict,
            Context& context,
            Istream& is
        );


protected:

    // Protected Member Functions

        //- Main driver: depending on 'equal' starts evaluating or
        //  skips forward to #else
        template<class Context>
        static bool execute
        (
            const bool equal,
            DynamicList<filePos>& stack,
            const dictionary& contextDict,
            Context& context,
            Istream& is
        );


public:

    //- Runtime type information
    FunctionTypeName("#ifeq");


    // Constructors

        //- Construct from function type, dictionary and tokenList
        ifeqEntry
        (
            const functionName& functionType,
            const dictionary& parentDict,
            const tokenList& tokens
        );

        //- Construct from dictionary and Istream
        ifeqEntry(const dictionary& parentDict, Istream& is);

        //- Disallow default bitwise copy construction
        ifeqEntry(const ifeqEntry&) = delete;


    // Member Functions

        //- Expand the functionEntry into the contextDict
        static bool execute(dictionary& contextDict, Istream& is);

        //- Expand the functionEntry adding tokens to the contextEntry
        static bool execute
        (
            const dictionary& contextDict,
            primitiveEntry& contextEntry,
            Istream&
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ifeqEntry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
