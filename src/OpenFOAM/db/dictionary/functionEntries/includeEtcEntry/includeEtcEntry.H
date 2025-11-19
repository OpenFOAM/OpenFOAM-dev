/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2025 OpenFOAM Foundation
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
        fileName includeEtcFileName
        (
            const fileName&,
            const dictionary&
        ) const;

        //- Return the file name
        const string& fName() const
        {
            return operator[](0).stringToken();
        }


public:

    // Static Data Members

        //- Report which file is included to stdout
        static bool log;


    //- Runtime type information
    FunctionTypeName("#includeEtc");


    // Constructors

        //- Construct from line number, dictionary and Istream
        includeEtcEntry
        (
            const label lineNumber,
            const dictionary& parentDict,
            Istream& is
        );

        //- Copy construct
        includeEtcEntry(const includeEtcEntry&) = default;

        //- Clone
        virtual autoPtr<entry> clone(const dictionary&) const
        {
            return autoPtr<entry>(new includeEtcEntry(*this));
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
        void operator=(const includeEtcEntry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
