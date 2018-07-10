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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void triSurface::writeSMESH(const bool writeSorted, Ostream& os) const
{
    const pointField& ps = points();

    // Write header
    os  << "# tetgen .smesh file" << endl
        << ps.size() << " 3" << endl;   // 3 dimensions

    // Write vertex coords
    forAll(ps, pointi)
    {
        os  << pointi << ' '
            << ps[pointi].x() << ' '
            << ps[pointi].y() << ' '
            << ps[pointi].z() << endl;
    }

    if (writeSorted)
    {
        labelList faceMap;

        surfacePatchList patches(calcPatches(faceMap));

        os  << size() << " 1" << endl;   // 1 attribute: region number

        label faceIndex = 0;

        forAll(patches, patchi)
        {
            // Print all faces belonging to this patch

            for
            (
                label patchFacei = 0;
                patchFacei < patches[patchi].size();
                patchFacei++
            )
            {
                const label facei = faceMap[faceIndex++];

                os  << "3 " // triangles
                    << operator[](facei)[0] << ' '
                    << operator[](facei)[1] << ' '
                    << operator[](facei)[2] << ' '
                    << operator[](facei).region()   // region number
                    << endl;
            }
        }

        os  << '0' << endl      // holes
            << '0' << endl;     // regions
    }
    else
    {
        os  << size() << " 1" << endl;   // 1 attribute: region number

        forAll(*this, facei)
        {
            os  << "3 "
                << operator[](facei)[0] << ' '
                << operator[](facei)[1] << ' '
                << operator[](facei)[2] << ' '
                << operator[](facei).region()       // region number
                << endl;
        }

        os  << '0' << endl      // holes
            << '0' << endl;     // regions
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
