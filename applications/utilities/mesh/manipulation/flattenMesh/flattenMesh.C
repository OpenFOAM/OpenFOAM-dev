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
    flattenMesh

Description
    Flattens the front and back planes of a 2D cartesian mesh.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMesh.H"
#include "twoDPointCorrector.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "addMeshOption.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createSpecifiedPolyMesh.H"

    pointIOField& points = mesh.lookupObjectRef<pointIOField>("points");

    const point midPoint = gAverage(points);

    const twoDPointCorrector& twoDCorr = twoDPointCorrector::New(mesh);

    const direction planeNormalCmpt = twoDCorr.normalDir();

    const scalar midCmptVal = midPoint[planeNormalCmpt];
    const scalar minCmptVal = mesh.bounds().min()[planeNormalCmpt];
    const scalar maxCmptVal = mesh.bounds().max()[planeNormalCmpt];

    forAll(points, pointi)
    {
        if (points[pointi][planeNormalCmpt] < midCmptVal)
        {
            points[pointi][planeNormalCmpt] = minCmptVal;
        }
        else
        {
            points[pointi][planeNormalCmpt] = maxCmptVal;
        }
    }

    twoDCorr.correctPoints(points);

    // Ensure the points are written to a sufficient precision
    IOstream::defaultPrecision(IOstream::highPrecision());

    Info<< "Writing points into directory " << points.path() << nl << endl;

    points.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
