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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void triSurface::writeTRI(const bool writeSorted, Ostream& os) const
{
    const pointField& ps = points();

    // Write as cloud of triangles

    labelList faceMap;

    surfacePatchList patches(calcPatches(faceMap));

    if (writeSorted)
    {
        label faceIndex = 0;

        forAll(patches, patchi)
        {
            for
            (
                label patchFacei = 0;
                patchFacei < patches[patchi].size();
                patchFacei++
            )
            {
                const label facei = faceMap[faceIndex++];

                const point& p1 = ps[operator[](facei)[0]];
                const point& p2 = ps[operator[](facei)[1]];
                const point& p3 = ps[operator[](facei)[2]];

                os  << p1.x() << token::SPACE
                    << p1.y() << token::SPACE
                    << p1.z() << token::SPACE

                    << p2.x() << token::SPACE
                    << p2.y() << token::SPACE
                    << p2.z() << token::SPACE

                    << p3.x() << token::SPACE
                    << p3.y() << token::SPACE
                    << p3.z() << token::SPACE

                    << "0x" << hex << operator[](facei).region() << dec
                    << endl;
            }
        }
    }
    else
    {
        forAll(*this, facei)
        {
            const point& p1 = ps[operator[](facei)[0]];
            const point& p2 = ps[operator[](facei)[1]];
            const point& p3 = ps[operator[](facei)[2]];

            os  << p1.x() << token::SPACE
                << p1.y() << token::SPACE
                << p1.z() << token::SPACE

                << p2.x() << token::SPACE
                << p2.y() << token::SPACE
                << p2.z() << token::SPACE

                << p3.x() << token::SPACE
                << p3.y() << token::SPACE
                << p3.z() << token::SPACE

                << "0x" << hex << operator[](facei).region() << dec
                << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
