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

#include "STARCDsurfaceFormatCore.H"
#include "clock.H"
#include "regExp.H"
#include "IStringStream.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// parse things like this:
//     CTNAME  1  someName
// don't bother with the older comma-delimited format

Foam::Map<Foam::word>
Foam::fileFormats::STARCDsurfaceFormatCore::readInpCellTable
(
    IFstream& is
)
{
    Map<word> lookup;

    regExp ctnameRE
    (
        " *CTNA[^ ]*"        // keyword - min 4 chars
        "[[:space:]]+"       // space delimited
        "([0-9]+)"           // 1: <digits>
        "[[:space:]]+"       // space delimited
        "([^,[:space:]].*)", // 2: <name>
        true                 // ignore case
    );

    string line;
    List<string> groups;
    while (is.good() && is.getLine(line).good())
    {
        if (ctnameRE.match(line, groups))
        {
            const label tableId = atoi(groups[0].c_str());

            // strip bad chars
            string::stripInvalid<word>(groups[1]);

            if (!groups[1].empty())
            {
                lookup.insert(tableId, groups[1]);
            }
        }
    }

    return lookup;
}


void Foam::fileFormats::STARCDsurfaceFormatCore::writeCase
(
    Ostream& os,
    const pointField& pointLst,
    const label nFaces,
    const UList<surfZone>& zoneLst
)
{
    word caseName = os.name().lessExt().name();

    os  << "! STAR-CD file written " << clock::dateTime().c_str() << nl
        << "! " << pointLst.size() << " points, " << nFaces << " faces" << nl
        << "! case " << caseName << nl
        << "! ------------------------------" << nl;

    forAll(zoneLst, zoneI)
    {
        os  << "ctable " << zoneI + 1 << " shell" << " ,,,,,," << nl
            << "ctname " << zoneI + 1 << " "
            << zoneLst[zoneI].name() << nl;
    }

    os  << "! ------------------------------" << nl
        << "*set icvo mxv - 1" << nl
        << "vread " << caseName << ".vrt icvo,,,coded" << nl
        << "cread " << caseName << ".cel icvo,,,add,coded" << nl
        << "*set icvo" << nl
        << "! end" << nl;

    os.flush();
}


// ************************************************************************* //
