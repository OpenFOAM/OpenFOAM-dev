/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "templateName.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::wordList Foam::templateName::split(const word& templateName)
{
    DynamicList<word> dynCmpts;

    string::size_type beg=0, end=0, endb=0, endc=0;

    while
    (
        (endb = templateName.find('<', beg)) != string::npos
     || (endc = templateName.find(',', beg)) != string::npos
    )
    {
        if (endb == string::npos)
        {
            end = endc;
        }
        else if ((endc = templateName.find(',', beg)) != string::npos)
        {
            end = min(endb, endc);
        }
        else
        {
            end = endb;
        }

        if (beg < end)
        {
            dynCmpts.append(templateName.substr(beg, end - beg));
            dynCmpts.last().replaceAll(">","");
        }

        beg = end + 1;
    }

    if (beg < templateName.size())
    {
        dynCmpts.append(templateName.substr(beg, string::npos));
        dynCmpts.last().replaceAll(">","");
    }

    wordList cmpts;
    cmpts.transfer(dynCmpts);
    return cmpts;
}


// ************************************************************************* //
