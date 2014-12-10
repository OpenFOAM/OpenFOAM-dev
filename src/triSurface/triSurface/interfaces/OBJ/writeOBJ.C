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

    surfacePatchList myPatches(calcPatches(faceMap));

    const pointField& ps = points();

    // Print patch names as comment
    forAll(myPatches, patchI)
    {
        os  << "#     " << patchI << "    "
            << myPatches[patchI].name() << nl;
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

        forAll(myPatches, patchI)
        {
            // Print all faces belonging to this patch

            os  << "g " << myPatches[patchI].name() << nl;

            for
            (
                label patchFaceI = 0;
                patchFaceI < myPatches[patchI].size();
                patchFaceI++
            )
            {
                const label faceI = faceMap[faceIndex++];

                os  << "f "
                    << operator[](faceI)[0] + 1 << ' '
                    << operator[](faceI)[1] + 1 << ' '
                    << operator[](faceI)[2] + 1
                    //<< "  # " << operator[](faceI).region()
                    << nl;
            }
        }
    }
    else
    {
        // Get patch (=compact region) per face
        labelList patchIDs(size());
        forAll(myPatches, patchI)
        {
            label faceI = myPatches[patchI].start();

            forAll(myPatches[patchI], i)
            {
                patchIDs[faceMap[faceI++]] = patchI;
            }
        }


        label prevPatchI = -1;

        forAll(*this, faceI)
        {
            if (prevPatchI != patchIDs[faceI])
            {
                prevPatchI = patchIDs[faceI];
                os  << "g " << myPatches[patchIDs[faceI]].name() << nl;
            }
            os  << "f "
                << operator[](faceI)[0] + 1 << ' '
                << operator[](faceI)[1] + 1 << ' '
                << operator[](faceI)[2] + 1
                << nl;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
