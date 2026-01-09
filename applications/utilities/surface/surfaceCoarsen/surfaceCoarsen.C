/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

Application
    surfaceCoarsen

Description
    Surface coarsening using 'bunnylod':

        Polygon Reduction Demo
        By Stan Melax (c) 1998
        mailto:melax@cs.ualberta.ca
        http://www.cs.ualberta.ca/~melax

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fileName.H"
#include "triSurface.H"
#include "OFstream.H"
#include "triFace.H"
#include "triFaceList.H"

// From bunnylod
#include "progmesh.h"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int mapVertex(::List<int>& collapse_map, int a, int mx)
{
    if (mx <= 0)
    {
        return 0;
    }
    while (a >= mx)
    {
        a = collapse_map[a];
    }
    return a;
}



int main(int argc, char *argv[])
{
    #include "removeCaseOptions.H"

    argList::validArgs.append("surface file");
    argList::validArgs.append("output surface file");
    argList::validArgs.append("reductionFactor");
    argList args(argc, argv);

    const fileName inFileName = args[1];
    const fileName outFileName = args[2];
    const scalar reduction = args.argRead<scalar>(3);

    if (reduction <= 0 || reduction > 1)
    {
        FatalErrorInFunction
            << "Reduction factor " << reduction
            << " should be within 0..1" << endl
            << "(it is the reduction in number of vertices)"
            << exit(FatalError);
    }

    Info<< nl << "Input surface    : " << inFileName
        << nl << "Reduction factor : " << reduction
        << nl << "Output surface   : " << outFileName << nl << endl;

    const triSurface surf(inFileName);

    Info<< "Surface:" << incrIndent << endl;
    surf.writeStats(Info);
    Info<< decrIndent << endl;

    ::List<::Vector> vertices;     // global list of vertices
    ::List<::tridata> triangles;     // global list of triangles

    // Convert triSurface to progmesh format. Note: can use global point
    // numbering since surface read in from file.
    const pointField& pts = surf.points();
    forAll(pts, ptI)
    {
        const point& pt = pts[ptI];
        vertices.Add(::Vector(pt.x(), pt.y(), pt.z()));
    }
    forAll(surf, facei)
    {
        const labelledTri& f = surf[facei];
        tridata td;
        td.v[0] = f[0];
        td.v[1] = f[1];
        td.v[2] = f[2];
        triangles.Add(td);
    }

    ::List<int> collapse_map;   // to which neighbor each vertex collapses
    ::List<int> permutation;
    ::ProgressiveMesh(vertices, triangles, collapse_map, permutation);

    // rearrange the vertex list
    ::List<::Vector> temp_list;
    for (int i = 0; i < vertices.num; i ++)
    {
        temp_list.Add(vertices[i]);
    }
    for (int i = 0; i < vertices.num; i ++)
    {
        vertices[permutation[i]] = temp_list[i];
    }

    // update the changes in the entries in the triangle list
    for (int i=0; i < triangles.num; i++)
    {
        for (int j=0;j<3;j++)
        {
            triangles[i].v[j] = permutation[triangles[i].v[j]];
        }
    }

    // Only get triangles with non-collapsed edges.
    int render_num = int(reduction * surf.nPoints());

    Info<< "Reducing to " << render_num << " vertices" << endl << endl;

    // Convert triangles into labelledTris
    Foam::List<labelledTri> newTriangles(surf.size());
    label newTrianglei = 0;
    for (int i = 0; i < triangles.num; i ++)
    {
        int p0 = mapVertex(collapse_map, triangles[i].v[0], render_num);
        int p1 = mapVertex(collapse_map, triangles[i].v[1], render_num);
        int p2 = mapVertex(collapse_map, triangles[i].v[2], render_num);

        if (p0 == p1 || p1 == p2 || p2 == p0) continue;

        newTriangles[newTrianglei++] = labelledTri(p0, p1, p2, 0);
    }
    newTriangles.setSize(newTrianglei);

    // Convert vertices into points
    pointField newPoints(vertices.num);
    for (int i = 0; i < vertices.num; i ++)
    {
        newPoints[i] = point(vertices[i].x, vertices[i].y, vertices[i].z);
    }

    // Construct the new surface
    triSurface newSurf(newTriangles, newPoints);

    // Remove any unused points
    newSurf = triSurface(newSurf.localFaces(), newSurf.localPoints());

    Info<< "Coarsened surface:" << incrIndent << endl;
    newSurf.writeStats(Info);
    Info<< decrIndent << endl;

    Info<< "Writing to file " << outFileName << endl << endl;

    newSurf.write(outFileName);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
