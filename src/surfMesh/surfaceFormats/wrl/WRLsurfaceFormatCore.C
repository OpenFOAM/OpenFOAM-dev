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

#include "WRLsurfaceFormatCore.H"
#include "clock.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fileFormats::WRLsurfaceFormatCore::writeHeader
(
    Ostream& os,
    const pointField& pointLst,
    const label nFaces,
    const UList<surfZone>& zoneLst
)
{
    os  << "#VRML V2.0 utf8" << nl
        << nl
        << "# written " << clock::dateTime().c_str() << nl
        << "# points : " << pointLst.size() << nl
        << "# faces  : " << nFaces << nl
        << "# zones  : " << zoneLst.size() << nl;

    // Print zone names as comment
    forAll(zoneLst, zoneI)
    {
        os  << "#   " << zoneI << "  " << zoneLst[zoneI].name()
            << "  (nFaces: " << zoneLst[zoneI].size() << ")" << nl;
    }
}


void Foam::fileFormats::WRLsurfaceFormatCore::writeAppearance
(
    Ostream& os
)
{
    os  <<
        "   appearance Appearance {\n"
        "    material Material {\n"
        "     diffuseColor   0.8 0.8 0.8\n"
        "     specularColor  1.0 1.0 1.0\n"
        "     shininess      0.5\n"
        "     transparency   0.0\n"
        "    }\n"                      // end material
        "   }\n";                      // end appearance

}


// ************************************************************************* //
