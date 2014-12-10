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
    surfacePointMerge

Description
    Merges points on surface if they are within absolute distance.
    Since absolute distance use with care!

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "triSurfaceTools.H"
#include "argList.H"
#include "OFstream.H"
#include "boundBox.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("surfaceFile");
    argList::validArgs.append("merge distance");
    argList::validArgs.append("output surfaceFile");
    argList args(argc, argv);

    const fileName surfFileName = args[1];
    const scalar   mergeTol = args.argRead<scalar>(2);
    const fileName outFileName = args[3];

    Info<< "Reading surface from " << surfFileName << " ..." << endl;
    Info<< "Merging points within " << mergeTol << " metre." << endl;

    triSurface surf1(surfFileName);

    Info<< "Original surface:" << endl;

    surf1.writeStats(Info);


    triSurface cleanSurf(surf1);

    while (true)
    {
        label nOldVert = cleanSurf.nPoints();

        cleanSurf = triSurfaceTools::mergePoints(cleanSurf, mergeTol);

        Info<< "After merging points:" << endl;

        cleanSurf.writeStats(Info);

        if (nOldVert == cleanSurf.nPoints())
        {
            break;
        }
    }

    cleanSurf.write(outFileName);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
