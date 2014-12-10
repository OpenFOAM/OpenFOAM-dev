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
    Wave propagation of nearwall distance through grid. Every iteration
    information goes through one layer of cells.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wallFvPatch.H"
#include "FaceCellWave.H"
#include "wallPoint.H"


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

    Info<< "Creating field wDistNC\n" << endl;
    volScalarField wallDistUncorrected
    (
        IOobject
        (
            "wDistNC",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "wallDist",
            dimensionSet(0, 1, 0, 0, 0),
            0.0
        )
    );

    //
    // Set initial changed faces: set wallPoint for wall faces to wall centre
    //

    // Count walls
    label nWalls = 0;
    forAll(mesh.boundary(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];

        if (isA<wallFvPatch>(patch))
        {
            nWalls += patch.size();
        }
    }

    List<wallPoint> faceDist(nWalls);
    labelList changedFaces(nWalls);

    label nChangedFaces = 0;
    forAll(mesh.boundary(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];

        if (isA<wallFvPatch>(patch))
        {
            forAll(patch.Cf(), patchFaceI)
            {
                const polyPatch& polyPatch = mesh.boundaryMesh()[patchI];

                label meshFaceI = polyPatch.start() + patchFaceI;

                changedFaces[nChangedFaces] = meshFaceI;

                faceDist[nChangedFaces] =
                    wallPoint(patch.Cf()[patchFaceI], 0.0);

                nChangedFaces++;
            }
        }
    }

    List<wallPoint> allFaceInfo(mesh.nFaces());
    List<wallPoint> allCellInfo(mesh.nCells());

    FaceCellWave<wallPoint> wallDistCalc
    (
        mesh,
        changedFaces,
        faceDist,
        allFaceInfo,
        allCellInfo,
        0             // max iterations
    );

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << endl;


        label nCells = wallDistCalc.faceToCell();

        Info<< "    Total changed cells   : " << nCells << endl;

        if (nCells == 0)
        {
            break;
        }


        label nFaces = wallDistCalc.cellToFace();

        Info<< "    Total changed faces   : " << nFaces << endl;

        if (nFaces == 0)
        {
            break;
        }


        //
        // Copy face and cell values into field
        //

        label nIllegal = 0;

        // Copy cell values
        forAll(allCellInfo, cellI)
        {
            scalar dist = allCellInfo[cellI].distSqr();
            if (allCellInfo[cellI].valid(wallDistCalc.data()))
            {
                wallDistUncorrected[cellI] = Foam::sqrt(dist);
            }
            else
            {
                wallDistUncorrected[cellI] = -1;
                nIllegal++;
            }
        }

        // Copy boundary values
        forAll(wallDistUncorrected.boundaryField(), patchI)
        {
            fvPatchScalarField& patchField =
                wallDistUncorrected.boundaryField()[patchI];

            forAll(patchField, patchFaceI)
            {
                const label meshFaceI = patchField.patch().start() + patchFaceI;

                scalar dist = allFaceInfo[meshFaceI].distSqr();
                if (allFaceInfo[meshFaceI].valid(wallDistCalc.data()))
                {
                    patchField[patchFaceI] = Foam::sqrt(dist);
                }
                else
                {
                    patchField[patchFaceI] = dist;
                    nIllegal++;
                }
            }
        }

        Info<< "nIllegal:" << nIllegal << endl;


        //
        // Write it
        //

        wallDistUncorrected.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n" << endl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
