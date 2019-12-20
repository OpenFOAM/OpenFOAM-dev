/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "etcFiles.H"
#include "OSspecific.H"
#include "foamVersion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::fileNameList Foam::findEtcDirs(const fileName& local)
{
    fileNameList dirs;

    // Search for user directories in
    // * ~/.OpenFOAM/VERSION
    // * ~/.OpenFOAM
    //
    fileName searchDir = home()/".OpenFOAM";
    if (isDir(searchDir))
    {
        fileName fullName = searchDir/FOAMversion/local;
        if (isDir(fullName))
        {
            dirs.append(fullName);
        }

        fullName = searchDir/local;
        if (isDir(fullName))
        {
            dirs.append(fullName);
        }
    }

    // Search for group (site) directories in
    // * $WM_PROJECT_SITE/VERSION/etc/
    // * $WM_PROJECT_SITE/etc/
    //
    searchDir = getEnv("WM_PROJECT_SITE");
    if (searchDir.size())
    {
        if (isDir(searchDir))
        {
            fileName fullName = searchDir/FOAMversion/"etc"/local;
            if (isDir(fullName))
            {
                dirs.append(fullName);
            }

            fullName = searchDir/"etc"/local;
            if (isDir(fullName))
            {
                dirs.append(fullName);
            }
        }
    }
    else
    {
        // Or search for group (site) files in
        // * $WM_PROJECT_INST_DIR/site/VERSION/etc/
        // * $WM_PROJECT_INST_DIR/site/etc/
        //
        searchDir = getEnv("WM_PROJECT_INST_DIR");
        if (isDir(searchDir))
        {
            fileName fullName = searchDir/"site"/FOAMversion/"etc"/local;
            if (isDir(fullName))
            {
                dirs.append(fullName);
            }

            fullName = searchDir/"site/etc"/local;
            if (isDir(fullName))
            {
                dirs.append(fullName);
            }
        }
    }

    // Search for other (shipped) files in
    // * $WM_PROJECT_DIR/etc
    //
    searchDir = getEnv("WM_PROJECT_DIR");
    if (isDir(searchDir))
    {
        fileName fullName = searchDir/"etc"/local;
        if (isDir(fullName))
        {
            dirs.append(fullName);
        }
    }

    return dirs;
}


Foam::fileNameList Foam::findEtcFiles
(
    const fileName& name,
    bool mandatory,
    bool findFirst
)
{
    fileNameList results;

    // Search for user files in
    // * ~/.OpenFOAM/VERSION
    // * ~/.OpenFOAM
    //
    fileName searchDir = home()/".OpenFOAM";
    if (isDir(searchDir))
    {
        fileName fullName = searchDir/FOAMversion/name;
        if (isFile(fullName))
        {
            results.append(fullName);
            if (findFirst)
            {
                return results;
            }
        }

        fullName = searchDir/name;
        if (isFile(fullName))
        {
            results.append(fullName);
            if (findFirst)
            {
                return results;
            }
        }
    }

    // Search for group (site) files in
    // * $WM_PROJECT_SITE/VERSION/etc/
    // * $WM_PROJECT_SITE/etc/
    //
    searchDir = getEnv("WM_PROJECT_SITE");
    if (searchDir.size())
    {
        if (isDir(searchDir))
        {
            fileName fullName = searchDir/FOAMversion/"etc"/name;
            if (isFile(fullName))
            {
                results.append(fullName);
                if (findFirst)
                {
                    return results;
                }
            }

            fullName = searchDir/"etc"/name;
            if (isFile(fullName))
            {
                results.append(fullName);
                if (findFirst)
                {
                    return results;
                }
            }
        }
    }
    else
    {
        // Or search for group (site) files in
        // * $WM_PROJECT_INST_DIR/site/VERSION/etc/
        // * $WM_PROJECT_INST_DIR/site/etc/
        //
        searchDir = getEnv("WM_PROJECT_INST_DIR");
        if (isDir(searchDir))
        {
            fileName fullName = searchDir/"site"/FOAMversion/"etc"/name;
            if (isFile(fullName))
            {
                results.append(fullName);
                if (findFirst)
                {
                    return results;
                }
            }

            fullName = searchDir/"site/etc"/name;
            if (isFile(fullName))
            {
                results.append(fullName);
                if (findFirst)
                {
                    return results;
                }
            }
        }
    }

    // Search for other (shipped) files in
    // * $WM_PROJECT_DIR/etc
    //
    searchDir = getEnv("WM_PROJECT_DIR");
    if (isDir(searchDir))
    {
        fileName fullName = searchDir/"etc"/name;
        if (isFile(fullName))
        {
            results.append(fullName);
            if (findFirst)
            {
                return results;
            }
        }
    }

    // Not found
    if (results.empty())
    {
        // Abort if the file is mandatory, otherwise return null
        if (mandatory)
        {
            std::cerr
                << "--> FOAM FATAL ERROR in Foam::findEtcFiles() :"
                   " could not find mandatory file\n    '"
                << name.c_str() << "'\n\n" << std::endl;
            ::exit(1);
        }
    }

    // Return list of matching paths or empty list if none found
    return results;
}


Foam::fileName Foam::findEtcFile(const fileName& name, bool mandatory)
{
    fileNameList results(findEtcFiles(name, mandatory, true));

    if (results.size())
    {
        return results[0];
    }
    else
    {
        return fileName();
    }
}


// ************************************************************************* //
