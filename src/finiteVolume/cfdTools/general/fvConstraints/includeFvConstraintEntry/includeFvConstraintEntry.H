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
    Foam::functionEntries::includeFvConstraintEntry

Description
    Specify a fvConstraint dictionary file to include, expects the
    fvConstraint name to follow with option arguments (without quotes).

    Searches for fvConstraint dictionary file in user/group/shipped
    directories allowing for version-specific and version-independent files
    using the following hierarchy:
    - \b user settings:
      - ~/.OpenFOAM/\<VERSION\>/caseDicts/fvConstraints
      - ~/.OpenFOAM/caseDicts/fvConstraints
    - \b group (site) settings (when $WM_PROJECT_SITE is set):
      - $WM_PROJECT_SITE/\<VERSION\>/etc/caseDicts/fvConstraints
      - $WM_PROJECT_SITE/etc/caseDicts/fvConstraints
    - \b group (site) settings (when $WM_PROJECT_SITE is not set):
      - $WM_PROJECT_INST_DIR/site/\<VERSION\>/etc/caseDicts/fvConstraints
      - $WM_PROJECT_INST_DIR/site/etc/caseDicts/fvConstraints
    - \b other (shipped) settings:
      - $WM_PROJECT_DIR/etc/caseDicts/fvConstraints

    The optional field arguments included in the name are inserted in 'field' or
    'fields' entries in the fvConstraint dictionary and included in the name
    of the fvConstraint entry to avoid conflict.

    Examples:
    \verbatim
        #includeConstraint limitPressure(minFactor=0.1, maxFactor=2)
        #includeConstraint limitTemperature(min=101, max=1000)
    \endverbatim
    or for a multiphase case:
    \verbatim
        #includeConstraint limitLowPressure(min=1e4)
        #includeConstraint limitTemperature(phase=steam, min=270, max=2000)
        #includeConstraint limitTemperature(phase=water, min=270, max=2000)
    \endverbatim

    Other dictionary entries may also be specified using named arguments.

See also
    Foam::includeFvModelEntry
    Foam::includeFuncEntry

SourceFiles
    includeFvConstraintEntry.C

\*---------------------------------------------------------------------------*/

#ifndef includeFvConstraintEntry_H
#define includeFvConstraintEntry_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                      Class includeFvConstraintEntry Declaration
\*---------------------------------------------------------------------------*/

class includeFvConstraintEntry
:
    public functionEntry
{

public:

    //- Runtime type information
    ClassName("includeConstraint");


    // Static Data Members

        //- Default relative path to the directory structure
        //  containing the fvConstraint dictionary files
        static fileName fvConstraintDictPath;


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
