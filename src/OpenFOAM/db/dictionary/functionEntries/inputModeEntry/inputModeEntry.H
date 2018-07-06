/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::functionEntries::inputModeEntry

Description
    Specify the input mode when reading dictionaries, expects
    a single word to follow.

    An example of \c \#inputMode directive:
    \verbatim
        #inputMode merge
    \endverbatim

    The possible input modes:
      - \par merge      merge sub-dictionaries when possible
      - \par overwrite  keep last entry and silently remove previous ones
      - \par protect    keep initial entry and silently ignore subsequent ones
      - \par warn       keep initial entry and warn about subsequent ones
      - \par error      issue a FatalError for duplicate entries
      - \par default    currently identical to merge

SourceFiles
    inputModeEntry.C

\*---------------------------------------------------------------------------*/

#ifndef inputModeEntry_H
#define inputModeEntry_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                        Class inputModeEntry Declaration
\*---------------------------------------------------------------------------*/

class inputModeEntry
:
    public functionEntry
{
        //- The input mode options
        enum inputMode
        {
            MERGE,
            OVERWRITE,
            PROTECT,
            WARN,
            ERROR
        };

        //- The current input mode
        static inputMode mode_;


    // Private Member Functions

        //- Read the mode as a word and set enum appropriately
        static void setMode(Istream&);

        //- Disallow default bitwise copy construct
        inputModeEntry(const inputModeEntry&);

        //- Disallow default bitwise assignment
        void operator=(const inputModeEntry&);


public:

    //- Runtime type information
    ClassName("inputMode");


    // Member Functions

        //- Execute the functionEntry in a sub-dict context
        static bool execute(dictionary& parentDict, Istream&);

        //- Reset the inputMode to %default (ie, %merge)
        static void clear();

        //- Return true if the inputMode is %merge
        static bool merge();

        //- Return true if the inputMode is %overwrite
        static bool overwrite();

        //- Return true if the inputMode is %protect
        static bool protect();

        //- Return true if the inputMode is %error
        static bool error();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
