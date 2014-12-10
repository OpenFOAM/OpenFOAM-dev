/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "TimePaths.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TimePaths::TimePaths
(
    const fileName& rootPath,
    const fileName& caseName,
    const word& systemName,
    const word& constantName
)
:
    processorCase_(false),
    rootPath_(rootPath),
    case_(caseName),
    system_(systemName),
    constant_(constantName)
{
    // Find out from case name whether a processor directory
    std::string::size_type pos = caseName.find("processor");
    if (pos != string::npos)
    {
        processorCase_ = true;

        if (pos == 0)
        {
            globalCaseName_ = ".";
        }
        else
        {
            globalCaseName_ = caseName(pos-1);
        }
    }
    else
    {
        globalCaseName_ = caseName;
    }
}


Foam::TimePaths::TimePaths
(
    const bool processorCase,
    const fileName& rootPath,
    const fileName& globalCaseName,
    const fileName& caseName,
    const word& systemName,
    const word& constantName
)
:
    processorCase_(processorCase),
    rootPath_(rootPath),
    globalCaseName_(globalCaseName),
    case_(caseName),
    system_(systemName),
    constant_(constantName)
{
    if (!processorCase)
    {
        // For convenience: find out from case name whether it is a
        // processor directory and set processorCase flag so file searching
        // goes up one level.
        std::string::size_type pos = caseName.find("processor");

        if (pos != string::npos)
        {
            processorCase_ = true;

            if (pos == 0)
            {
                globalCaseName_ = ".";
            }
            else
            {
                globalCaseName_ = caseName(pos-1);
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::TimePaths::caseSystem() const
{
    if (processorCase_)
    {
        return ".."/system();
    }
    else
    {
        return system();
    }
}


Foam::fileName Foam::TimePaths::caseConstant() const
{
    if (processorCase_)
    {
        return ".."/constant();
    }
    else
    {
        return constant();
    }
}



// ************************************************************************* //
