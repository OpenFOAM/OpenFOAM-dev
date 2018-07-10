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

Description
    Lightwave OBJ format.

    Note: Java obj loader does not support '#' on line

\*---------------------------------------------------------------------------*/

#include "triSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void triSurface::writeOBJ(const bool writeSorted, Ostream& os) const
{
    // Write header
    os  << "# Wavefront OBJ file" << nl
        << "# Regions:" << nl;

    labelList faceMap;

    surfacePatchList patches(calcPatches(faceMap));

    const pointField& ps = points();

    // Print patch names as comment
    forAll(patches, patchi)
    {
        os  << "#     " << patchi << "    "
            << patches[patchi].name() << nl;
    }
    os  << "#" << nl;

    os  << "# points    : " << ps.size() << nl
        << "# triangles : " << size() << nl
        << "#" << nl;


    // Write vertex coords
    forAll(ps, pointi)
    {
        os  << "v "
            << ps[pointi].x() << ' '
            << ps[pointi].y() << ' '
            << ps[pointi].z() << nl;
    }

    if (writeSorted)
    {
        label faceIndex = 0;

        forAll(patches, patchi)
        {
            // Print all faces belonging to this patch

            os  << "g " << patches[patchi].name() << nl;

            for
            (
                label patchFacei = 0;
                patchFacei < patches[patchi].size();
                patchFacei++
            )
            {
                const label facei = faceMap[faceIndex++];

                os  << "f "
                    << operator[](facei)[0] + 1 << ' '
                    << operator[](facei)[1] + 1 << ' '
                    << operator[](facei)[2] + 1
                    //<< "  # " << operator[](facei).region()
                    << nl;
            }
        }
    }
    else
    {
        // Get patch (=compact region) per face
        labelList patchIDs(size());
        forAll(patches, patchi)
        {
            label facei = patches[patchi].start();

            forAll(patches[patchi], i)
            {
                patchIDs[faceMap[facei++]] = patchi;
            }
        }


        label prevPatchi = -1;

        forAll(*this, facei)
        {
            if (prevPatchi != patchIDs[facei])
            {
                prevPatchi = patchIDs[facei];
                os  << "g " << patches[patchIDs[facei]].name() << nl;
            }
            os  << "f "
                << operator[](facei)[0] + 1 << ' '
                << operator[](facei)[1] + 1 << ' '
                << operator[](facei)[2] + 1
                << nl;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
