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

#include "VTKsurfaceFormatCore.H"
#include "clock.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fileFormats::VTKsurfaceFormatCore::writeHeader
(
    Ostream& os,
    const pointField& pointLst
)
{
    // Write header
    os  << "# vtk DataFile Version 2.0" << nl
        << "surface written " << clock::dateTime().c_str() << nl
        << "ASCII" << nl
        << nl
        << "DATASET POLYDATA" << nl;

    // Write vertex coords
    os  << "POINTS " << pointLst.size() << " float" << nl;
    forAll(pointLst, ptI)
    {
        const point& pt = pointLst[ptI];

        os  << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
    }
}


void Foam::fileFormats::VTKsurfaceFormatCore::writeTail
(
    Ostream& os,
    const UList<surfZone>& zoneLst
)
{
    label nFaces = 0;
    forAll(zoneLst, zoneI)
    {
        nFaces += zoneLst[zoneI].size();
    }

    // Print zone numbers
    os  << nl
        << "CELL_DATA " << nFaces << nl
        << "FIELD attributes 1" << nl
        << "region 1 " << nFaces << " float" << nl;


    forAll(zoneLst, zoneI)
    {
        forAll(zoneLst[zoneI], localFacei)
        {
            if (localFacei)
            {
                if (localFacei % 20)
                {
                    os << ' ';
                }
                else
                {
                    os << nl;
                }
            }
            os  << zoneI + 1;
        }
        os  << nl;
    }
}


void Foam::fileFormats::VTKsurfaceFormatCore::writeTail
(
    Ostream& os,
    const labelUList& zoneIds
)
{
    // Print zone numbers
    os  << nl
        << "CELL_DATA " << zoneIds.size() << nl
        << "FIELD attributes 1" << nl
        << "region 1 " << zoneIds.size() << " float" << nl;

    forAll(zoneIds, facei)
    {
        if (facei)
        {
            if (facei % 20)
            {
                os << ' ';
            }
            else
            {
                os << nl;
            }
        }
        os  << zoneIds[facei] + 1;
    }
    os  << nl;
}


// ************************************************************************* //
