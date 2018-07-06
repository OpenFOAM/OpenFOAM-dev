/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "triSurface.H"
#include "VTKsurfaceFormat.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::triSurface::readVTK(const fileName& fName)
{
    // Read (and triangulate) point, faces, zone info
    fileFormats::VTKsurfaceFormat<triFace> surf(fName);

    List<labelledTri> tris(surf.faces().size());
    forAll(tris, i)
    {
        const triFace& f = surf[i];
        tris[i] = labelledTri(f[0], f[1], f[2], 0);
    }

    // Add regions from zone
    const List<surfZone>& surfZones = surf.surfZones();

    geometricSurfacePatchList patches;

    if (surfZones.size())
    {
        patches.setSize(surfZones.size());
        forAll(surfZones, zoneI)
        {
            const surfZone& zone = surfZones[zoneI];

            // Add patch. Convert synthetic 'zone' name into 'patch' for now.
            // (vtk format does not contain region names)
            word regionName = zone.name();
            if (regionName != (string("zone") + name(zoneI)))
            {
                regionName = string("patch") + name(zoneI);
            }

            patches[zoneI] = geometricSurfacePatch
            (
                (
                    zone.geometricType() != word::null
                  ? zone.geometricType()
                  : "empty"
                ),
                regionName,
                zoneI
            );

            // Set triangle regions
            for (label i = zone.start(); i < zone.start()+zone.size(); i++)
            {
                tris[i].region() = zoneI;
            }
        }
    }
    else
    {
        // Add single patch
        patches[0] = geometricSurfacePatch("empty", "patch0", 0);

        // Triangle regions already set to 0
    }


    // Create triSurface
    *this = triSurface
    (
        tris.xfer(),
        patches,
        xferCopy<List<point>>(surf.points())
    );

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
