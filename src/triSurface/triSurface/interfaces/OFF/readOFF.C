/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Geomview OFF polyList format. Does triangulation.

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "polygonTriangulate.H"
#include "IFstream.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::triSurface::readOFF(const fileName& OFFfileName)
{
    IFstream OFFfile(OFFfileName);

    if (!OFFfile.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << OFFfileName
            << exit(FatalError);
    }

    // Read header
    string hdr = getLineNoComment(OFFfile);
    if (hdr != "OFF")
    {
        FatalErrorInFunction
            << "OFF file " << OFFfileName
            << " does not start with 'OFF'"
            << exit(FatalError);
    }


    label nPoints, nEdges, nElems;

    string line = getLineNoComment(OFFfile);
    {
        IStringStream lineStream(line);
        lineStream >> nPoints >> nElems >> nEdges;
    }

    // Read points
    pointField points(nPoints);

    forAll(points, pointi)
    {
        scalar x, y, z;
        line = getLineNoComment(OFFfile);
        {
            IStringStream lineStream(line);
            lineStream >> x >> y >> z;
        }
        points[pointi] = point(x, y, z);
    }

    // Create a triangulation engine
    polygonTriangulate triEngine;

    // Read faces & triangulate them,
    DynamicList<labelledTri> tris(nElems);

    for (label facei = 0; facei < nElems; facei++)
    {
        line = getLineNoComment(OFFfile);
        {
            IStringStream lineStream(line);

            label nVerts;
            lineStream >> nVerts;

            face f(nVerts);

            forAll(f, fp)
            {
                lineStream >> f[fp];
            }

            // Triangulate.
            if (nVerts == 3)
            {
                tris.append(labelledTri(f[0], f[1], f[2], 0));
            }
            else if (nVerts == 4)
            {
                tris.append(labelledTri(f[0], f[1], f[2], 0));
                tris.append(labelledTri(f[2], f[3], f[0], 0));
            }
            else
            {
                triEngine.triangulate(UIndirectList<point>(points, f));

                forAll(triEngine.triPoints(), trii)
                {
                    tris.append(labelledTri(triEngine.triPoints(trii, f), 0));
                }
            }
        }
    }

    tris.shrink();

    *this = triSurface(tris, points);

    return true;
}


// ************************************************************************* //
