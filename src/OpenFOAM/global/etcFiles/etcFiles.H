/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

InNamespace
    Foam

Description
    Functions to search 'etc' directories for configuration files etc.

SourceFiles
    etcFiles.C

\*---------------------------------------------------------------------------*/

#ifndef etcFiles_H
#define etcFiles_H

#include "fileNameList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Search for directories from user/group/shipped directories.
//  The search scheme allows for version-specific and
//  version-independent files using the following hierarchy:
//  - \b user settings:
//    - ~/.OpenFOAM/\<VERSION\>
//    - ~/.OpenFOAM/
//  - \b group (site) settings (when $WM_PROJECT_SITE is set):
//    - $WM_PROJECT_SITE/\<VERSION\>
//    - $WM_PROJECT_SITE
//  - \b group (site) settings (when $WM_PROJECT_SITE is not set):
//    - $WM_PROJECT_INST_DIR/site/\<VERSION\>
//    - $WM_PROJECT_INST_DIR/site/
//  - \b other (shipped) settings:
//    - $WM_PROJECT_DIR/etc/
//
//  \return The list of full paths of all the matching directories or
//  an empty list if the name cannot be found.
fileNameList findEtcDirs(const fileName& local = fileName::null);

//- Search for files from user/group/shipped directories.
//  The search scheme allows for version-specific and
//  version-independent files using the following hierarchy:
//  - \b user settings:
//    - ~/.OpenFOAM/\<VERSION\>
//    - ~/.OpenFOAM/
//  - \b group (site) settings (when $WM_PROJECT_SITE is set):
//    - $WM_PROJECT_SITE/\<VERSION\>
//    - $WM_PROJECT_SITE
//  - \b group (site) settings (when $WM_PROJECT_SITE is not set):
//    - $WM_PROJECT_INST_DIR/site/\<VERSION\>
//    - $WM_PROJECT_INST_DIR/site/
//  - \b other (shipped) settings:
//    - $WM_PROJECT_DIR/etc/
//
//  \return The list of full paths of all the matching files or
//  an empty list if the name cannot be found.
//  Optionally abort if the file cannot be found.
//  Optionally stop search after the first file has been found.
fileNameList findEtcFiles
(
    const fileName&,
    bool mandatory=false,
    bool findFirst=false
);

//- Search for a file using findEtcFiles.
//  \return The full path name of the first file found in the
//  search hierarchy or an empty fileName if the name cannot be found.
//  Optionally abort if the file cannot be found.
fileName findEtcFile(const fileName&, bool mandatory=false);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
