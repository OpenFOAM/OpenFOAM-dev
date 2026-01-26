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
        //  temporarily inserted into contextDict
        List<Tuple3<word, string, label>> insertNamedArgs
        (
            dictionary& contextDict,
            const Tuple2<string, label>& args
        ) const;

        //- Remove named argument entries from contextDict
        //  renaming any existing entries which had the same name
        void removeInsertNamedArgs
        (
            dictionary& contextDict,
            const List<Tuple3<word, string, label>>& namedArgs
        ) const;



protected:

    // Protected Member Functions

        //- Return the file name
        const string& fName() const
        {
            return operator[](0).stringToken();
        }

        //- Return the arguments string for processing by insertNamedArgs
        Tuple2<string, label> args() const
        {
            if (size() > 1)
            {
                return Tuple2<string, label>
                (
                    char(operator[](1).pToken())
                  + operator[](2).stringToken()
                  + char(operator[](3).pToken()),
                    operator[](1).lineNumber()
                );
            }
            else
            {
                return Tuple2<string, label>(string::null, 0);
            }
        }

        //- Expand fileName and return
        static fileName includeFileName
        (
            const Istream& is,
            const fileName&,
            const dictionary&
        );

        //- Expand include fileName and return
        static fileName includeFileName
        (
            const fileName& dir,
            const fileName&,
            const dictionary&
        );

        //- Expand the functionEntry into the contextDict
        //  with support for includeIfPresent
        bool execute
        (
            dictionary& contextDict,
            Istream&,
            const bool ifPresent
        );

        //- Expand the functionEntry into the contextEntry
        //  with support for includeIfPresent
        static bool execute
        (
            const dictionary& contextDict,
            primitiveEntry& contextEntry,
            Istream&,
            const bool ifPresent
        );


public:

    // Static Data Members

        //- Report which file is included to stdout
        static bool log;

        //- Runtime type information
        FunctionTypeName("#include");


    // Constructors

        //- Construct from functionType, dictionary and Istream
        includeEntry
        (
            const functionName& functionType,
            const label lineNumber,
            const dictionary& parentDict,
            Istream& is
        );

        //- Construct from line number, dictionary and Istream
        includeEntry
        (
            const label lineNumber,
            const dictionary& parentDict,
            Istream& is
        );

        //- Copy construct
        includeEntry(const includeEntry&) = default;

        //- Clone
        virtual autoPtr<entry> clone(const dictionary&) const
        {
            return autoPtr<entry>(new includeEntry(*this));
        }


    // Member Functions

        //- Expand the functionEntry into the contextDict
        virtual bool execute(dictionary& contextDict, Istream&);

        //- Expand the functionEntry into the contextEntry
        static bool execute
        (
            const dictionary& contextDict,
            primitiveEntry& contextEntry,
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
