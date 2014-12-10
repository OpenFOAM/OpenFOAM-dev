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

Description
    OpenDX format. Both data-only and scalar/vector data.

\*---------------------------------------------------------------------------*/

#include "triSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Geometry (positions + connections)
// writeSorted: sort acc. to patch
void triSurface::writeDXGeometry
(
    const bool writeSorted,
    Ostream& os
) const
{
    labelList faceMap;
    surfacePatchList myPatches(calcPatches(faceMap));

    // Print patch names as comment
    os  << "# Patches:" << endl;
    forAll(myPatches, patchI)
    {
        os  << "#     " << patchI << "    "
            << myPatches[patchI].name() << endl;
    }
    os  << nl << endl;

    // Write vertex coordinates

    os  << "# The irregular positions" << endl
        << "object 1 class array type float rank 1 shape 3 items "
        << nPoints() << " data follows" << endl;
    forAll(localPoints(), pointI)
    {
        const point& pt = localPoints()[pointI];
        os  << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
    }
    os  << endl;

    os  << "# The irregular connections (triangles)" << endl
        << "object 2 class array type int rank 1 shape 3 items "
        << size() << " data follows" << endl;

    if (writeSorted)
    {
        label faceIndex = 0;

        forAll(myPatches, patchI)
        {
            // Print all faces belonging to this patch

            for
            (
                label patchFaceI = 0;
                patchFaceI < myPatches[patchI].size();
                patchFaceI++
            )
            {
                const label faceI = faceMap[faceIndex++];
                const labelledTri& f = localFaces()[faceI];

                os  << f[0] << ' ' << f[1] << ' ' << f[2] << endl;
            }
        }
    }
    else
    {
        forAll(*this, faceI)
        {
            const labelledTri& f = localFaces()[faceI];

            os  << f[0] << ' ' << f[1] << ' ' << f[2] << endl;
        }
    }
    os  << "attribute \"element type\" string \"triangles\"" << endl
        << "attribute \"ref\" string \"positions\"" << endl << endl;
}


// Standard trailer
void triSurface::writeDXTrailer(Ostream& os) const
{
    os  << "# the field, with three components: \"positions\", \"connections\""
        << ", and \"data\"" << endl
        << "object \"irregular positions irregular connections\" class field"
        << endl
        << "component \"positions\" value 1" << endl
        << "component \"connections\" value 2" << endl
        << "component \"data\" value 3" << endl;
}


// Geometry only (data field is either faceIndex or patchIndex)
void triSurface::writeDX(const bool writeSorted, Ostream& os) const
{
    writeDXGeometry(writeSorted, os);

    os  << "object 3 class array type float rank 0 items " << size()
        << " data follows" << endl;
    if (writeSorted)
    {
        // Write patch number as data

        labelList faceMap;
        surfacePatchList myPatches(calcPatches(faceMap));

        forAll(myPatches, patchI)
        {
            forAll(myPatches[patchI], patchFaceI)
            {
                os  << patchI << endl;
            }
        }
    }
    else
    {
        // Write face number as data

        forAll(*this, faceI)
        {
            os  << faceI << endl;
        }
    }

    os  << endl << "attribute \"dep\" string \"connections\"" << endl << endl;

    writeDXTrailer(os);

    os  << "end" << endl;
}


// Geometry + scalar data
void triSurface::writeDX(const scalarField& field, Ostream& os) const
{
    writeDXGeometry(false, os);

    if (field.size() == size())
    {
        // Connections dependent data
        os  << "object 3 class array type float rank 0 items " << field.size()
            << " data follows" << endl;
        forAll(field, faceI)
        {
            os  << field[faceI] << endl;
        }
        os  << endl
            << "attribute \"dep\" string \"connections\"" << endl << endl;
    }
    else if (field.size() == nPoints())
    {
        // Positions dependent data
        os  << "object 3 class array type float rank 0 items " << field.size()
            << " data follows" << endl;
        forAll(field, pointI)
        {
            os  << field[pointI] << endl;
        }
        os  << endl
            << "attribute \"dep\" string \"positions\"" << endl << endl;
    }
    else
    {
        FatalErrorIn
        (
            "writeDX(const scalarField&, Ostream&)"
        )   << "Illegal field size " << field.size() << " is not equal "
            << " to number of faces " << size() << " or to number "
            << " of points " << nPoints() << exit(FatalError);
    }

    writeDXTrailer(os);

    os  << "end" << endl;
}


// Geometry + vector data
void triSurface::writeDX(const vectorField& field, Ostream& os) const
{
    writeDXGeometry(false, os);

    if (field.size() == size())
    {
        // Connections dependent data
        os  << "object 3 class array type float rank 1 shape 3 items "
            << field.size() << " data follows" << endl;
        forAll(field, faceI)
        {
            os  << field[faceI].x() << ' '
                << field[faceI].y() << ' '
                << field[faceI].z() << endl;
        }
        os  << endl
            << "attribute \"dep\" string \"connections\"" << endl << endl;
    }
    else if (field.size() == nPoints())
    {
        // Positions dependent data
        os  << "object 3 class array type float rank 1 shape 3 items "
            << field.size() << " data follows" << endl;
        forAll(field, pointI)
        {
            os  << field[pointI].x() << ' '
                << field[pointI].y() << ' '
                << field[pointI].z() << endl;
        }
        os  << endl
            << "attribute \"dep\" string \"positions\"" << endl << endl;
    }
    else
    {
        FatalErrorIn
        (
            "writeDX(const vectorField&, Ostream&)"
        )   << "Illegal field size " << field.size() << " is not equal "
            << " to number of faces " << size() << " or to number "
            << " of points " << nPoints() << exit(FatalError);
    }

    writeDXTrailer(os);

    os  << "end" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
