/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

Description
    Test pointEdgeWave.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "pointMesh.H"
#include "OSspecific.H"
#include "IFstream.H"
#include "pointEdgePoint.H"
#include "PointEdgeWave.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("(patches)");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    labelList patchIDs
    (
        pbm.patchSet(wordReList(IStringStream(args[1])())).sortedToc()
    );

    Info<< "Starting walk from patches "
        << UIndirectList<word>(pbm.names(), patchIDs)
        << nl
        << endl;

    label nPoints = 0;
    forAll(patchIDs, i)
    {
        nPoints += pbm[patchIDs[i]].nPoints();
    }

    Info<< "Seeding " << returnReduce(nPoints, sumOp<label>())
        << " patch points" << nl << endl;


    // Set initial changed points to all the patch points(if patch present)
    List<pointEdgePoint> wallInfo(nPoints);
    labelList wallPoints(nPoints);
    nPoints = 0;

    forAll(patchIDs, i)
    {
        // Retrieve the patch now we have its index in patches.
        const polyPatch& pp = pbm[patchIDs[i]];

        forAll(pp.meshPoints(), ppI)
        {
            label meshPointi = pp.meshPoints()[ppI];
            wallPoints[nPoints] = meshPointi;
            wallInfo[nPoints] = pointEdgePoint(mesh.points()[meshPointi], 0.0);
            nPoints++;
        }
    }

    // Current info on points
    List<pointEdgePoint> allPointInfo(mesh.nPoints());

    // Current info on edges
    List<pointEdgePoint> allEdgeInfo(mesh.nEdges());

    PointEdgeWave<pointEdgePoint> wallCalc
    (
        mesh,
        wallPoints,
        wallInfo,

        allPointInfo,
        allEdgeInfo,
        mesh.nPoints()  // max iterations
    );


    pointScalarField psf
    (
        IOobject
        (
            "wallDist",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(mesh),
        dimensionedScalar(dimLength, 0)
    );

    forAll(allPointInfo, pointi)
    {
        psf[pointi] = Foam::sqrt(allPointInfo[pointi].distSqr());
    }

    Info<< "Writing wallDist pointScalarField to " << runTime.value()
        << endl;

    psf.write();

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
