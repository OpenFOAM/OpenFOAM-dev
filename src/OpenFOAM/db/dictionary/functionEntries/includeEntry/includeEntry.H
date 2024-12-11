/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::functionEntries::includeEntry

Description
    Specify an include file when reading dictionaries.

    The file name is provided as a string in order to support the '/' directory
    separator and unusual characters as well as the usual expansion of
    environment variables and other constructs
    (eg, the \c ~OpenFOAM/ expansion).

    e.g. to include the file \c includeFile from the \c include directory:
    \verbatim
        #include "include/includeFile"
    \endverbatim

    Additionally optional arguments may be specified which are available for
    substitution in the included file e.g.
    \verbatim
        #include "regionFunctions"(region = heatedWall)
    \endverbatim

    will substitute \c heatedWall for \c $region in the \c regionFunctions file

    \verbatim
        #includeFunc sample(region = $region)
    \endverbatim

See also
    fileName

SourceFiles
    includeEntry.C

\*---------------------------------------------------------------------------*/

#ifndef includeEntry_H
#define includeEntry_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                        Class includeEntry Declaration
\*---------------------------------------------------------------------------*/

class includeEntry
:
    public functionEntry
{
    // Private Member Functions

        //- Return the optional named arguments
        //  temporarily inserted into parentDict
        static List<Tuple3<word, string, label>> insertNamedArgs
        (
            dictionary& parentDict,
            Istream& is
        );

        //- Remove named argument entries from parentDict
        //  renaming any existing entries which had the same name
        static void removeInsertNamedArgs
        (
            dictionary& parentDict,
            const List<Tuple3<word, string, label>>& namedArgs
        );


protected:

    // Protected Member Functions

        //- Read the include fileName from Istream, expand and return
        static fileName includeFileName(Istream&, const dictionary&);

        //- Expand include fileName and return
        static fileName includeFileName
        (
            const fileName& dir,
            const fileName&,
            const dictionary&
        );


public:

    // Static Data Members

        //- Report which file is included to stdout
        static bool log;


    //- Runtime type information
    ClassName("include");


    // Constructors

        //- Disallow default bitwise copy construction
        includeEntry(const includeEntry&) = delete;


    // Member Functions

        //- Execute the functionEntry in a sub-dict context
        static bool execute(dictionary& parentDict, Istream&);

        //- Execute the functionEntry in a primitiveEntry context
        static bool execute
        (
            const dictionary& parentDict,
            primitiveEntry&,
            Istream&
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const includeEntry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
