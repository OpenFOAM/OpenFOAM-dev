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

#include "triSurface.H"
#include "IOmanip.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::triSurface::writeAC(Ostream& os) const
{
    // Write with patches as separate objects under "world" object.
    // Header is taken over from sample file.
    // Defines separate materials for all patches. Recycle colours.

    // Define 8 standard colours as r,g,b components
    static scalar colourMap[] =
    {
        1, 1, 1,
        1, 0, 0,
        0, 1, 0,
        0, 0, 1,
        1, 1, 0,
        0, 1, 1,
        1, 0, 1,
        0.5, 0.5, 1
    };

    // Calculate patch face indexing

    labelList faceMap;

    surfacePatchList patches(calcPatches(faceMap));


    // Write header. Define materials.

    os  << "AC3Db" << endl;

    forAll(patches, patchi)
    {
        const word& pName = patches[patchi].name();

        label colourI = patchi % 8;
        label colourCompI = 3 * colourI;

        os  << "MATERIAL \"" << pName << "Mat\" rgb "
            << colourMap[colourCompI] << ' ' << colourMap[colourCompI+1]
            << ' ' << colourMap[colourCompI+2]
            << "  amb 0.2 0.2 0.2  emis 0 0 0  spec 0.5 0.5 0.5  shi 10"
            << "  trans 0"
            << endl;
    }

    os  << "OBJECT world" << endl
        << "kids " << patches.size() << endl;


    // Write patch points & faces.

    label faceIndex = 0;

    forAll(patches, patchi)
    {
        const surfacePatch& sp = patches[patchi];

        os  << "OBJECT poly" << endl
            << "name \"" << sp.name() << '"' << endl;

        // Create patch with only patch faces included for ease of addressing

        boolList include(size(), false);

        forAll(sp, patchFacei)
        {
            const label facei = faceMap[faceIndex++];

            include[facei] = true;
        }

        labelList pointMap;
        labelList faceMap;

        triSurface patch = subsetMesh(include, pointMap, faceMap);

        // Now we have triSurface for this patch alone. Write it.

        os << "numvert " << patch.nPoints() << endl;

        forAll(patch.localPoints(), ptI)
        {
            const point& pt = patch.localPoints()[ptI];

            os << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
        }

        os << "numsurf " << patch.localFaces().size() << endl;

        forAll(patch.localFaces(), facei)
        {
            const labelledTri& f = patch.localFaces()[facei];

            os  << "SURF 0x20" << endl          // polygon
                << "mat " << patchi << endl
                << "refs " << f.size() << endl;

            os << f[0] << " 0 0" << endl;
            os << f[1] << " 0 0" << endl;
            os << f[2] << " 0 0" << endl;
        }

        os << "kids 0" << endl;
    }
}


// ************************************************************************* //
