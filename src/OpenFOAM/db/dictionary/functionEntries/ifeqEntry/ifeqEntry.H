/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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
    Conditional parsing of dictionary entries.

    E.g.
    \verbatim
        a #calc "0.123";
        b 1.23e-1;

        #ifeq $a $b
            ..
        #else
            ..
        #endif
    \endverbatim

    \verbatim
    ddtSchemes
    {
    #ifeq ${FOAM_APPLICATION} simpleFoam
        default         steadyState;
    #else
        default         Euler;
    #endif
    }
    \endverbatim

    Note:
    - supports dictionary variables and environment variables
    - the two arguments should be two tokens
    - the comparison is a string comparison for any word/string/variable,
      integer comparison for two integers and floating point comparison for
      any floating point number.
    - parsing of (non)matching \c #else, \c #endif is not very sophisticated

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

        static bool evaluate
        (
            const bool doIf,
            DynamicList<filePos>& stack,
            dictionary& parentDict,
            Istream& is
        );

        //- Main driver: depending on 'equal' starts evaluating or
        //  skips forward to #else
        static bool execute
        (
            DynamicList<filePos>& stack,
            dictionary& parentDict,
            Istream& is
        );


protected:

    // Protected Member Functions

        //- Main driver: depending on 'equal' starts evaluating or
        //  skips forward to #else
        static bool execute
        (
            const bool equal,
            DynamicList<filePos>& stack,
            dictionary& parentDict,
            Istream& is
        );


public:

    //- Runtime type information
    ClassName("ifeq");

    // Function names for the individual elements of the conditional
    static const functionName ifName;
    static const functionName ifeqName;
    static const functionName elifName;
    static const functionName elseName;
    static const functionName endifName;


    // Constructors

        //- Disallow default bitwise copy construction
        ifeqEntry(const ifeqEntry&) = delete;


    // Member Functions

        //- Execute the functionEntry in a sub-dict context
        static bool execute(dictionary& parentDict, Istream& is);


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
