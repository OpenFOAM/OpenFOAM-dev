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
    surfaceClean

Description
    - removes baffles
    - collapses small edges, removing triangles.
    - converts sliver triangles into split edges by projecting point onto
      base of triangle.

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "argList.H"
#include "OFstream.H"

#include "collapseBase.H"
#include "collapseEdge.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    #include "removeCaseOptions.H"

    argList::validArgs.append("surface file");
    argList::validArgs.append("output surface file");
    argList::validArgs.append("min length");
    argList::validArgs.append("min quality");
    argList::addBoolOption
    (
        "noClean",
        "perform some surface checking/cleanup on the input surface"
    );
    argList args(argc, argv);

    const fileName inFileName = args[1];
    const fileName outFileName = args[2];
    const scalar minLen = args.argRead<scalar>(3);
    const scalar minQuality = args.argRead<scalar>(4);

    Info<< "Reading surface " << inFileName << nl
        << "Collapsing all triangles with" << nl
        << "    edges or heights < " << minLen << nl
        << "    quality          < " << minQuality << nl
        << "Writing result to " << outFileName << nl << endl;


    Info<< "Reading surface from " << inFileName << " ..." << nl << endl;
    triSurface surf(inFileName);
    surf.writeStats(Info);

    if (!args.optionFound("noClean"))
    {
        Info<< "Removing duplicate and illegal triangles ..." << nl << endl;
        surf.cleanup(true);
    }

    Info<< "Collapsing triangles to edges ..." << nl << endl;

    while (true)
    {
        label nEdgeCollapse = collapseEdge(surf, minLen);

        if (nEdgeCollapse == 0)
        {
            break;
        }
    }
    while (true)
    {
        label nSplitEdge = collapseBase(surf, minLen, minQuality);

        if (nSplitEdge == 0)
        {
            break;
        }
    }

    Info<< nl
        << "Resulting surface:" << endl;
    surf.writeStats(Info);

    Info<< nl
        << "Writing refined surface to " << outFileName << " ..." << endl;
    surf.write(outFileName);

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
