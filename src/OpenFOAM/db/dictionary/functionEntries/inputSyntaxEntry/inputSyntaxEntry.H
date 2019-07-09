/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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
    Foam::functionEntries::inputSyntaxEntry

Description
    Specify the dictionary scoping syntax, expects a single word to follow.

    An example of \c \#inputSyntax directive:
    \verbatim
        #inputSyntax dot
    \endverbatim

    The possible input syntax:
      - \par slash      use '/' as the scope operator
                        and '../' to go to the parent dictionary
      - \par dot        use '.' as the scope operator
                        and '..' to go to the parent dictionary
      - \par default    currently identical to dot

    The default dictionary syntax is \c dot but changed to \c slash in
    etc/controlDict
    \verbatim
        OptimisationSwitches
        {
            .
            .
            .
            inputSyntax slash;
        }
    \endverbatim

SourceFiles
    inputSyntaxEntry.C

\*---------------------------------------------------------------------------*/

#ifndef inputSyntaxEntry_H
#define inputSyntaxEntry_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                        Class inputSyntaxEntry Declaration
\*---------------------------------------------------------------------------*/

class inputSyntaxEntry
:
    public functionEntry
{
        //- The input syntax options
        enum inputSyntax
        {
            SLASH,
            DOT
        };

        //- The current input syntax
        static inputSyntax syntax_;


    // Private Member Functions

        //- Read the syntax as a word and return the corresponding enum
        static inputSyntax syntax(Istream&);


public:

    //- Runtime type information
    ClassName("inputSyntax");


    // Constructors

        //- Disallow default bitwise copy construction
        inputSyntaxEntry(const inputSyntaxEntry&) = delete;


    // Member Functions

        //- Execute the functionEntry in a sub-dict context
        static bool execute(dictionary& parentDict, Istream&);

        //- Reset the inputSyntax to %default (ie, %merge)
        static void clear();

        //- Return true if the inputSyntax is %slash
        static bool slash();

        //- Return true if the inputSyntax is %dot
        static bool dot();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const inputSyntaxEntry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
