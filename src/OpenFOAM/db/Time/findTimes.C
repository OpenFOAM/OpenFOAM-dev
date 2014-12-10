/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Description
    Searches the current case directory for valid times
    and sets the time list to these.
    This is done if a times File does not exist.

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "OSspecific.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::instantList Foam::Time::findTimes
(
    const fileName& directory,
    const word& constantName
)
{
    if (debug)
    {
        Info<< "Time::findTimes(const fileName&): finding times in directory "
            << directory << endl;
    }

    // Read directory entries into a list
    fileNameList dirEntries(readDir(directory, fileName::DIRECTORY));

    // Initialise instant list
    instantList Times(dirEntries.size() + 1);
    label nTimes = 0;

    // Check for "constant"
    bool haveConstant = false;
    forAll(dirEntries, i)
    {
        if (dirEntries[i] == constantName)
        {
            Times[nTimes].value() = 0;
            Times[nTimes].name() = dirEntries[i];
            nTimes++;
            haveConstant = true;
            break;
        }
    }

    // Read and parse all the entries in the directory
    forAll(dirEntries, i)
    {
        IStringStream timeStream(dirEntries[i]);
        token timeToken(timeStream);

        if (timeToken.isNumber() && timeStream.eof())
        {
            Times[nTimes].value() = timeToken.number();
            Times[nTimes].name() = dirEntries[i];
            nTimes++;
        }
    }

    // Reset the length of the times list
    Times.setSize(nTimes);

    if (haveConstant)
    {
        if (nTimes > 2)
        {
            std::sort(&Times[1], Times.end(), instant::less());
        }
    }
    else if (nTimes > 1)
    {
        std::sort(&Times[0], Times.end(), instant::less());
    }

    return Times;
}

// ************************************************************************* //
