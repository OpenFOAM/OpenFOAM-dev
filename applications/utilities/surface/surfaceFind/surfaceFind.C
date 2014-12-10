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

Application
    surfaceFind

Description
    Finds nearest face and vertex.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "OFstream.H"

#include "MeshedSurfaces.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("surfaceFile");
    argList::addOption("x", "X", "The point x-coordinate (if non-zero)");
    argList::addOption("y", "Y", "The point y-coordinate (if non-zero)");
    argList::addOption("z", "Z", "The point y-coordinate (if non-zero)");

    argList args(argc, argv);

    const point samplePt
    (
        args.optionLookupOrDefault<scalar>("x", 0),
        args.optionLookupOrDefault<scalar>("y", 0),
        args.optionLookupOrDefault<scalar>("z", 0)
    );
    Info<< "Looking for nearest face/vertex to " << samplePt << endl;


    Info<< "Reading surf ..." << endl;
    meshedSurface surf1(args[1]);

    //
    // Nearest vertex
    //

    const pointField& localPoints = surf1.localPoints();

    label minIndex = -1;
    scalar minDist = GREAT;

    forAll(localPoints, pointI)
    {
        const scalar dist = mag(localPoints[pointI] - samplePt);
        if (dist < minDist)
        {
            minDist = dist;
            minIndex = pointI;
        }
    }

    Info<< "Nearest vertex:" << nl
        << "    index      :" << minIndex << " (in localPoints)" << nl
        << "    index      :" << surf1.meshPoints()[minIndex]
        << " (in points)" << nl
        << "    coordinates:" << localPoints[minIndex] << nl
        << endl;

    //
    // Nearest face
    //

    const pointField& points = surf1.points();

    minIndex = -1;
    minDist = GREAT;

    forAll(surf1, faceI)
    {
        const point centre = surf1[faceI].centre(points);

        const scalar dist = mag(centre - samplePt);
        if (dist < minDist)
        {
            minDist = dist;
            minIndex = faceI;
        }
    }

    const face& f = surf1[minIndex];

    Info<< "Face with nearest centre:" << nl
        << "    index        :" << minIndex << nl
        << "    centre       :" << f.centre(points) << nl
        << "    face         :" << f << nl
        << "    vertex coords:\n";
    forAll(f, fp)
    {
        Info<< "        " << points[f[fp]] << "\n";
    }

    Info<< endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
