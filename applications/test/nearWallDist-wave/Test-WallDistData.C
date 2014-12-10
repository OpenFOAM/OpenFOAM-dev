/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Calculate distance to wall.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wallDistData.H"
#include "wallPointData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    Info<< "Mesh read in = "
        << runTime.cpuTimeIncrement()
        << " s\n" << endl << endl;


    Info<< "Time now = " << runTime.timeName() << endl;

    // wall distance and reflection vectors

    volVectorField n
    (
        IOobject
        (
            "n",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedVector("n", dimLength, point::max)
    );

    // Fill wall patches with unit normal
    forAll(mesh.boundary(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];

        if (isA<wallFvPatch>(patch))
        {
            fvPatchVectorField& wallData = n.boundaryField()[patchI];

            forAll(patch.Cf(), patchFaceI)
            {
                wallData[patchFaceI] = patch.Cf()[patchFaceI];
                wallData[patchFaceI] /= Foam::mag(wallData[patchFaceI]);
            }
        }
    }


    // Do distance calculation, transporting values of n.
    wallDistData<wallPointData<vector> > y(mesh, n, true);

    if (y.nUnset() != 0)
    {
        WarningIn(args.executable())
            << "There are " << y.nUnset()
            << " remaining unset cells and/or boundary values" << endl;
    }


    y.write();
    y.data().write();

    runTime++;

    Info<< "Time now = " << runTime.timeName() << endl;


    // Move points

    boundBox meshBb(mesh.points());

    pointField newPoints(mesh.points());

    const point half(0.5*(meshBb.min() + meshBb.max()));

    forAll(newPoints, pointI)
    {
        point& pt = newPoints[pointI];

        // expand around half
        pt.y() +=  pt.y() - half.y();
    }

    mesh.movePoints(newPoints);

    mesh.write();

    y.correct();

    y.write();
    y.data().write();


    Info<< "End\n" << endl;

    return 0;



}


// ************************************************************************* //
