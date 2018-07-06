/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "foamVersion.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Name>
inline Foam::word Foam::IOobject::groupName(Name name, const word& group)
{
    if (group != word::null)
    {
        return name + ('.' + group);
    }
    else
    {
        return name;
    }
}


template<class Stream>
inline Stream& Foam::IOobject::writeBanner(Stream& os, bool noHint)
{
    static bool spacesSet(false);
    static char spaces[40];

    if (!spacesSet)
    {
        memset(spaces, ' ', 40);

        size_t len = strlen(Foam::FOAMversion);
        if (len < 38)
        {
            spaces[38 - len] = '\0';
        }
        else
        {
            spaces[0] = '\0';
        }
        spacesSet = true;
    }

    if (noHint)
    {
        os  <<
            "/*--------------------------------------"
            "-------------------------------------*\\\n";
    }
    else
    {
        os  <<
            "/*--------------------------------*- C++ "
            "-*----------------------------------*\\\n";
    }

    os  <<
        "| =========                 |"
        "                                                 |\n"
        "| \\\\      /  F ield         |"
        " OpenFOAM: The Open Source CFD Toolbox           |\n"
        "|  \\\\    /   O peration     |"
        " Version:  " << FOAMversion << spaces << "|\n"
        "|   \\\\  /    A nd           |"
        " Website:  https://openfoam.org                  |\n"
        "|    \\\\/     M anipulation  |"
        "                                                 |\n"
        "\\*-----------------------------------------"
        "----------------------------------*/\n";

    return os;
}


template<class Stream>
inline Stream& Foam::IOobject::writeDivider(Stream& os)
{
    os  <<
        "// * * * * * * * * * * * * * * * * * "
        "* * * * * * * * * * * * * * * * * * * * //\n";

    return os;
}

template<class Stream>
inline Stream& Foam::IOobject::writeEndDivider(Stream& os)
{
    os  << "\n\n"
        "// *****************************************"
        "******************************** //\n";

    return os;
}


// ************************************************************************* //
