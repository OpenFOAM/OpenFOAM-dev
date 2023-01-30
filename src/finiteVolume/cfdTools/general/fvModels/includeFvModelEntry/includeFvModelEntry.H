/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::functionEntries::includeFvModelEntry

Description
    Specify a fvModel dictionary file to include, expects the
    fvModel name to follow with option arguments (without quotes).

    Searches for fvModel dictionary file in user/group/shipped
    directories allowing for version-specific and version-independent files
    using the following hierarchy:
    - \b user settings:
      - ~/.OpenFOAM/\<VERSION\>/caseDicts/fvModels
      - ~/.OpenFOAM/caseDicts/fvModels
    - \b group (site) settings (when $WM_PROJECT_SITE is set):
      - $WM_PROJECT_SITE/\<VERSION\>/etc/caseDicts/fvModels
      - $WM_PROJECT_SITE/etc/caseDicts/fvModels
    - \b group (site) settings (when $WM_PROJECT_SITE is not set):
      - $WM_PROJECT_INST_DIR/site/\<VERSION\>/etc/caseDicts/fvModels
      - $WM_PROJECT_INST_DIR/site/etc/caseDicts/fvModels
    - \b other (shipped) settings:
      - $WM_PROJECT_DIR/etc/caseDicts/fvModels

    The optional field arguments included in the name are inserted in 'field' or
    'fields' entries in the fvModel dictionary and included in the name
    of the fvModel entry to avoid conflict.

    Examples:
    \verbatim
        #includeModel clouds
        #includeModel surfaceFilms
    \endverbatim

    Other dictionary entries may also be specified using named arguments.

See also
    Foam::includeFvConstraintEntry
    Foam::includeFuncEntry

SourceFiles
    includeFvModelEntry.C

\*---------------------------------------------------------------------------*/

#ifndef includeFvModelEntry_H
#define includeFvModelEntry_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                      Class includeFvModelEntry Declaration
\*---------------------------------------------------------------------------*/

class includeFvModelEntry
:
    public functionEntry
{

public:

    //- Runtime type information
    ClassName("includeModel");


    // Static Data Members

        //- Default relative path to the directory structure
        //  containing the fvModel dictionary files
        static fileName fvModelDictPath;


    // Member Functions

        //- Execute the functionEntry in a sub-dict context
        static bool execute(dictionary& parentDict, Istream&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
