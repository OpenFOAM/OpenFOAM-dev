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
    surfaceRefineRedGreen

Description
    Refine by splitting all three edges of triangle ('red' refinement).
    Neighbouring triangles (which are not marked for refinement get split
    in half ('green' refinement). (R. Verfuerth, "A review of a posteriori
    error estimation and adaptive mesh refinement techniques",
    Wiley-Teubner, 1996)

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "triSurfaceTools.H"
#include "argList.H"
#include "OFstream.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("surfaceFile");
    argList::validArgs.append("output surfaceFile");
    argList args(argc, argv);

    const fileName surfFileName = args[1];
    const fileName outFileName  = args[2];

    Info<< "Reading surface from " << surfFileName << " ..." << endl;

    triSurface surf1(surfFileName);

    // Refine
    triSurface surf2 = triSurfaceTools::redGreenRefine
    (
        surf1,
        identity(surf1.size())  //Hack: refine all
    );

    Info<< "Original surface:" << endl
        << "    triangles     :" << surf1.size() << endl
        << "    vertices(used):" << surf1.nPoints() << endl
        << "Refined surface:" << endl
        << "    triangles     :" << surf2.size() << endl
        << "    vertices(used):" << surf2.nPoints() << endl << endl;


    Info<< "Writing refined surface to " << outFileName << " ..." << endl;

    surf2.write(outFileName);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
