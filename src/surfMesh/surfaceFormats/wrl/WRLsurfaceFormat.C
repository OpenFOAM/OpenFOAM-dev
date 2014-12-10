/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "WRLsurfaceFormat.H"

#include "Ostream.H"
#include "OFstream.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::WRLsurfaceFormat<Face>::WRLsurfaceFormat()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::fileFormats::WRLsurfaceFormat<Face>::write
(
    const fileName& filename,
    const MeshedSurfaceProxy<Face>& surf
)
{
    const pointField& pointLst = surf.points();
    const List<Face>&  faceLst = surf.faces();
    const List<label>& faceMap = surf.faceMap();

    // for no zones, suppress the group name
    const List<surfZone>& zones =
    (
        surf.surfZones().empty()
      ? surfaceFormatsCore::oneZone(faceLst, "")
      : surf.surfZones()
    );

    const bool useFaceMap = (surf.useFaceMap() && zones.size() > 1);

    OFstream os(filename);
    if (!os.good())
    {
        FatalErrorIn
        (
            "fileFormats::WRLsurfaceFormat::write"
            "(const fileName&, const MeshedSurfaceProxy<Face>&)"
        )
            << "Cannot open file for writing " << filename
            << exit(FatalError);
    }

    writeHeader(os, pointLst, faceLst.size(), zones);

    os  << "\n"
        "Group {\n"
        " children [\n"
        "  Shape {\n";

   writeAppearance(os);

   os  <<
        "   geometry IndexedFaceSet {\n"
        "    coord Coordinate {\n"
        "     point [\n";

    // Write vertex coords
    forAll(pointLst, ptI)
    {
        const point& pt = pointLst[ptI];

        os  << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
    }

    os  <<
        "     ]\n"                     // end point
        "    }\n"                      // end coord Coordinate
        "    coordIndex [\n";

    label faceIndex = 0;
    forAll(zones, zoneI)
    {
        const surfZone& zone = zones[zoneI];

        if (useFaceMap)
        {
            forAll(zone, localFaceI)
            {
                const Face& f = faceLst[faceMap[faceIndex++]];

                forAll(f, fp)
                {
                    os << f[fp] << ' ';
                }
                os << "-1,\n";
            }
        }
        else
        {
            forAll(zone, localFaceI)
            {
                const Face& f = faceLst[faceIndex++];

                forAll(f, fp)
                {
                    os << ' ' << f[fp];
                }
                os << " -1,\n";
            }
        }
    }

    os  <<
        "    ]\n"                      // end coordIndex
        "   }\n"                       // end geometry IndexedFaceSet
        "  }\n"                        // end Shape
        " ]\n"                         // end children
        "}\n";                         // end Group
}


// ************************************************************************* //
