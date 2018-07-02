/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

Application
    surfaceOrient

Description
    Set normal consistent with respect to a user provided 'outside' point.
    If the -inside option is used the point is considered inside.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "triSurfaceSearch.H"
#include "orientedSurface.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    #include "removeCaseOptions.H"

    argList::addNote
    (
        "set face normals consistent with a user-provided 'outside' point"
    );

    argList::validArgs.append("surface file");
    argList::validArgs.append("output surface file");
    argList::validArgs.append("visiblePoint");
    argList::addBoolOption
    (
        "inside",
        "treat provided point as being inside"
    );
    argList::addBoolOption
    (
        "usePierceTest",
        "determine orientation by counting number of intersections"
    );

    argList args(argc, argv);

    const fileName surfFileName = args[1];
    const fileName outFileName  = args[2];
    const point visiblePoint    = args.argRead<point>(3);

    const bool orientInside = args.optionFound("inside");
    const bool usePierceTest = args.optionFound("usePierceTest");

    Info<< "Reading surface from " << surfFileName << nl
        << "Orienting surface such that visiblePoint " << visiblePoint
        << " is ";

    if (orientInside)
    {
        Info<< "inside" << endl;
    }
    else
    {
        Info<< "outside" << endl;
    }



    // Load surface
    triSurface surf(surfFileName);


    bool anyFlipped = false;

    if (usePierceTest)
    {
        triSurfaceSearch surfSearches(surf);

        anyFlipped = orientedSurface::orient
        (
            surf,
            surfSearches,
            visiblePoint,
           !orientInside
        );
    }
    else
    {
        anyFlipped = orientedSurface::orient
        (
            surf,
            visiblePoint,
           !orientInside
        );
    }

    if (anyFlipped)
    {
        Info<< "Flipped orientation of (part of) surface." << endl;
    }
    else
    {
        Info<< "Did not flip orientation of any triangle of surface." << endl;
    }

    Info<< "Writing new surface to " << outFileName << endl;

    surf.write(outFileName);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
