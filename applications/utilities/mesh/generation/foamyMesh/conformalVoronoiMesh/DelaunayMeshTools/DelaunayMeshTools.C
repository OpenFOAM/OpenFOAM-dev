/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "DelaunayMeshTools.H"
#include "meshTools.H"
#include "OFstream.H"
#include "pointConversion.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::DelaunayMeshTools::writeOBJ
(
    const fileName& fName,
    const List<Foam::point>& points
)
{
    if (points.size())
    {
        OFstream str(fName);

        Pout<< nl
            << "Writing " << points.size() << " points from pointList to "
            << str.name() << endl;

        forAll(points, p)
        {
            meshTools::writeOBJ(str, points[p]);
        }
    }
}


void Foam::DelaunayMeshTools::writeOBJ
(
    const fileName& fName,
    const List<Vb>& points
)
{
    if (points.size())
    {
        OFstream str(fName);

        Pout<< nl
            << "Writing " << points.size() << " points from pointList to "
            << str.name() << endl;

        forAll(points, p)
        {
            meshTools::writeOBJ(str, topoint(points[p].point()));
        }
    }
}


void Foam::DelaunayMeshTools::writeObjMesh
(
    const fileName& fName,
    const pointField& points,
    const faceList& faces
)
{
    OFstream str(fName);

    Pout<< nl
        << "Writing points and faces to " << str.name() << endl;

    forAll(points, p)
    {
        meshTools::writeOBJ(str, points[p]);
    }

    forAll(faces, f)
    {
        str<< 'f';

        const face& fP = faces[f];

        forAll(fP, p)
        {
            str<< ' ' << fP[p] + 1;
        }

        str<< nl;
    }
}


// ************************************************************************* //
