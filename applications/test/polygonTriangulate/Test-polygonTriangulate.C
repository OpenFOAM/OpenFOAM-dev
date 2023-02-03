/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

#include "argList.H"
#include "clock.H"
#include "OBJstream.H"
#include "polygonTriangulate.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("number of edges");
    argList::addOption
    (
        "error",
        "value",
        "polygon error - default is 0"
    );
    argList::addOption
    (
        "seed",
        "value",
        "random number generator seed - default is clock time"
    );
    argList::addBoolOption
    (
        "simple",
        "assume polygon as no self intersections"
    );
    argList::addBoolOption
    (
        "nonOptimal",
        "do not optimise the triangulation quality"
    );

    Foam::argList args(argc, argv);
    Info<< nl;

    // Get the size of the polygon
    const label n = args.argRead<label>(1);

    // Initialise the random number generator
    label seed;
    if (args.optionFound("seed"))
    {
        seed = args.optionRead<label>("seed");
    }
    else
    {
        seed = clock::getTime();
        Info<< "Seeding random number generator with value " << seed
            << nl << endl;
    }
    Random rndGen(seed);

    // Get controls
    const bool simple = args.optionFound("simple");
    const bool nonOptimal = args.optionFound("nonOptimal");

    // Generate a random polygon
    const List<point> polygon =
        polygonTriangulate::randomPolygon
        (
            rndGen,
            n,
            args.optionLookupOrDefault<scalar>("error", 0)
        );

    // Write the polygon
    {
        OBJstream os(args[0] + "_polygon.obj");
        os.write(face(identityMap(polygon.size())), polygon, false);
    }

    // Triangulate the polygon
    const List<triFace> triPoints =
        polygonTriangulate().triangulate(polygon, simple, !nonOptimal);

    // Write the triangulation
    faceList triFacePoints(triPoints.size());
    forAll(triPoints, trii)
    {
        triFacePoints[trii] = triPoints[trii];
    }
    {
        OBJstream os(args[0] + "_triangulation.obj");
        os.write(triFacePoints, pointField(polygon), false);
    }

    Info<< "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
