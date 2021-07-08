/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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
    Foam::functionEntries::includeFuncEntry

Description
    Specify a functionObject dictionary file to include, expects the
    functionObject name to follow with option arguments (without quotes).

    Searches for functionObject dictionary file in user/group/shipped
    directories allowing for version-specific and version-independent files
    using the following hierarchy:
    - \b user settings:
      - ~/.OpenFOAM/\<VERSION\>/caseDicts/postProcessing
      - ~/.OpenFOAM/caseDicts/postProcessing
    - \b group (site) settings (when $WM_PROJECT_SITE is set):
      - $WM_PROJECT_SITE/\<VERSION\>/etc/caseDicts/postProcessing
      - $WM_PROJECT_SITE/etc/caseDicts/postProcessing
    - \b group (site) settings (when $WM_PROJECT_SITE is not set):
      - $WM_PROJECT_INST_DIR/site/\<VERSION\>/etc/caseDicts/postProcessing
      - $WM_PROJECT_INST_DIR/site/etc/caseDicts/postProcessing
    - \b other (shipped) settings:
      - $WM_PROJECT_DIR/etc/caseDicts/postProcessing

    The optional field arguments included in the name are inserted in 'field' or
    'fields' entries in the functionObject dictionary and included in the name
    of the functionObject entry to avoid conflict.

    Examples:
    \verbatim
        #includeFunc Q
        #includeFunc components(U)
        #includeFunc mag(Ux)
        #includeFunc mag(p)
    \endverbatim

    Other dictionary entries may also be specified using named arguments.

See also
    Foam::functionObjectList

SourceFiles
    includeFuncEntry.C

\*---------------------------------------------------------------------------*/

#ifndef includeFuncEntry_H
#define includeFuncEntry_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                      Class includeFuncEntry Declaration
\*---------------------------------------------------------------------------*/

class includeFuncEntry
:
    public functionEntry
{

public:

    //- Runtime type information
    ClassName("includeFunc");


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
