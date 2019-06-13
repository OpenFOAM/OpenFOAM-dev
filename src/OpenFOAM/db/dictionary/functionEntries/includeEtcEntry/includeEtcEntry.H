/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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
    Foam::functionEntries::includeEtcEntry

Description
    Specify an etc file to include when reading dictionaries, expects a
    single string to follow.

    Searches for files from user/group/shipped directories.
    The search scheme allows for version-specific and
    version-independent files using the following hierarchy:
    - \b user settings:
      - ~/.OpenFOAM/\<VERSION\>
      - ~/.OpenFOAM/
    - \b group (site) settings (when $WM_PROJECT_SITE is set):
      - $WM_PROJECT_SITE/\<VERSION\>/etc/
      - $WM_PROJECT_SITE/etc/
    - \b group (site) settings (when $WM_PROJECT_SITE is not set):
      - $WM_PROJECT_INST_DIR/site/\<VERSION\>/etc/
      - $WM_PROJECT_INST_DIR/site/etc/
    - \b other (shipped) settings:
      - $WM_PROJECT_DIR/etc/

    An example of the \c \#includeEtc directive:
    \verbatim
        #includeEtc "etcFile"
    \endverbatim

    The usual expansion of environment variables and other constructs is
    retained.

See also
    findEtcFile, fileName, string::expand()

SourceFiles
    includeEtcEntry.C

\*---------------------------------------------------------------------------*/

#ifndef includeEtcEntry_H
#define includeEtcEntry_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                        Class includeEtcEntry Declaration
\*---------------------------------------------------------------------------*/

class includeEtcEntry
:
    public functionEntry
{
    // Private Member Functions

        //- Expand include fileName and return
        static fileName includeEtcFileName
        (
            const fileName&,
            const dictionary&
        );


public:

    // Static Data Members

        //- Report which file is included to stdout
        static bool log;


    //- Runtime type information
    ClassName("includeEtc");


    // Constructors

        //- Disallow default bitwise copy construction
        includeEtcEntry(const includeEtcEntry&) = delete;


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
        void operator=(const includeEtcEntry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
