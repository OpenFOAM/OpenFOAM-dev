/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Test-hexRef8

Description
    Test app for refinement and unrefinement. Runs a few iterations refining
    and unrefining.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "hexRef8.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "Random.H"
#include "zeroGradientFvPatchFields.H"
#include "calculatedPointPatchFields.H"
#include "pointConstraints.H"
#include "fvcDiv.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool notEqual(const scalar s1, const scalar s2, const scalar tol)
{
    return mag(s1-s2) > tol;
}


// Main program:
int main(int argc, char *argv[])
{
    #include "addTimeOptions.H"
    argList::validArgs.append("inflate (true|false)");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"


    const pointConstraints& pc = pointConstraints::New(pointMesh::New(mesh));

    const Switch inflate(args.args()[1]);

    if (inflate)
    {
        Info<< "Splitting/deleting cells using inflation/deflation" << nl
            << endl;
    }
    else
    {
        Info<< "Splitting/deleting cells, introducing points at new position"
            << nl << endl;
    }


    Random rndGen(0);


    // Force generation of V()
    (void)mesh.V();


    // Test mapping
    // ------------

    // 1. uniform field stays uniform
    volScalarField one
    (
        IOobject
        (
            "one",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1.0),
        zeroGradientFvPatchScalarField::typeName
    );
    Info<< "Writing one field "
        << one.name() << " in " << runTime.timeName() << endl;
    one.write();


    // 2. linear profile gets preserved
    volScalarField ccX
    (
        IOobject
        (
            "ccX",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh.C().component(0)
    );
    Info<< "Writing x component of cell centres to "
        << ccX.name()
        << " in " << runTime.timeName() << endl;
    ccX.write();


    // Uniform surface field
    surfaceScalarField surfaceOne
    (
        IOobject
        (
            "surfaceOne",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1.0),
        calculatedFvsPatchScalarField::typeName
    );
    Info<< "Writing surface one field "
        << surfaceOne.name() << " in " << runTime.timeName() << endl;
    surfaceOne.write();


    // Uniform point field
    pointScalarField pointX
    (
        IOobject
        (
            "pointX",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(mesh),
        dimensionedScalar("one", dimless, 1.0),
        calculatedPointPatchScalarField::typeName
    );
    pointX.primitiveFieldRef() = mesh.points().component(0);
    pointX.correctBoundaryConditions();
    Info<< "Writing x-component field "
        << pointX.name() << " in " << runTime.timeName() << endl;
    pointX.write();



    // Force allocation of V. Important for any mesh changes since otherwise
    // old time volumes are not stored
    const scalar totalVol = gSum(mesh.V());


    // Construct refiner. Read initial cell and point levels.
    hexRef8 meshCutter(mesh);


    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (mesh.globalData().nTotalCells() == 0)
        {
            break;
        }


        mesh.moving(false);
        mesh.topoChanging(false);


        label action = rndGen.sampleAB<label>(0, 6);


        if (action == 0)
        {
            Info<< nl << "-- moving only" << endl;
            mesh.movePoints(pointField(mesh.points()));
        }
        else if (action == 1 || action == 2)
        {
            // Mesh changing engine.
            polyTopoChange meshMod(mesh);

            if (action == 1)
            {
                // Refine
                label nRefine = mesh.nCells()/20;
                DynamicList<label> refineCandidates(nRefine);

                for (label i=0; i<nRefine; i++)
                {
                    refineCandidates.append
                    (
                        rndGen.sampleAB<label>(0, mesh.nCells())
                    );
                }

                labelList cellsToRefine
                (
                    meshCutter.consistentRefinement
                    (
                        refineCandidates,
                        true                  // buffer layer
                    )
                );
                Info<< nl << "-- selected "
                    << returnReduce(cellsToRefine.size(), sumOp<label>())
                    << " cells out of " << mesh.globalData().nTotalCells()
                    << " for refinement" << endl;

                // Play refinement commands into mesh changer.
                meshCutter.setRefinement(cellsToRefine, meshMod);
            }
            else
            {
                // Unrefine
                labelList allSplitPoints(meshCutter.getSplitPoints());

                label nUnrefine = allSplitPoints.size()/20;
                labelHashSet candidates(2*nUnrefine);

                for (label i=0; i<nUnrefine; i++)
                {
                    label index =
                        rndGen.sampleAB<label>(0, allSplitPoints.size());
                    candidates.insert(allSplitPoints[index]);
                }

                labelList splitPoints = meshCutter.consistentUnrefinement
                (
                    candidates.toc(),
                    false
                );
                Info<< nl << "-- selected "
                    << returnReduce(splitPoints.size(), sumOp<label>())
                    << " points out of "
                    << returnReduce(allSplitPoints.size(), sumOp<label>())
                    << " for unrefinement" << endl;

                // Play refinement commands into mesh changer.
                meshCutter.setUnrefinement(splitPoints, meshMod);
            }




            // Create mesh, return map from old to new mesh.
            Info<< nl << "-- actually changing mesh" << endl;
            autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, inflate);

            // Update fields
            Info<< nl << "-- mapping mesh data" << endl;
            mesh.updateMesh(map);

            // Inflate mesh
            if (map().hasMotionPoints())
            {
                Info<< nl << "-- moving mesh" << endl;
                mesh.movePoints(map().preMotionPoints());
            }

            // Update numbering of cells/vertices.
            Info<< nl << "-- mapping hexRef8 data" << endl;
            meshCutter.updateMesh(map);
        }


        Info<< nl<< "-- Mesh : moving:" << mesh.moving()
            << " topoChanging:" << mesh.topoChanging()
            << " changing:" << mesh.changing()
            << endl;



        Info<< "Writing fields" << nl << endl;
        runTime.write();



        // Check mesh volume conservation
        if (mesh.moving())
        {
            #include "volContinuity.H"
        }
        else
        {
            if (mesh.V().size() != mesh.nCells())
            {
                FatalErrorInFunction
                    << "Volume not mapped. V:" << mesh.V().size()
                    << " nCells:" << mesh.nCells()
                    << exit(FatalError);
            }

            const scalar newVol = gSum(mesh.V());
            Info<< "Initial volume = " << totalVol
                << "  New volume = " << newVol
                << endl;

            if (mag(newVol-totalVol)/totalVol > 1e-10)
            {
                FatalErrorInFunction
                    << "Volume loss: old volume:" << totalVol
                    << "  new volume:" << newVol
                    << exit(FatalError);
            }
            else
            {
                Info<< "Volume check OK" << nl << endl;
            }
        }


        // Check constant profile
        {
            const scalar max = gMax(one);
            const scalar min = gMin(one);

            Info<< "Uniform one field min = " << min
                << "  max = " << max << endl;

            if (notEqual(max, 1.0, 1e-10) || notEqual(min, 1.0, 1e-10))
            {
                FatalErrorInFunction
                    << "Uniform volVectorField not preserved."
                    << " Min and max should both be 1.0. min:" << min
                    << " max:" << max
                    << exit(FatalError);
            }
            else
            {
                Info<< "Uniform field mapping check OK" << nl << endl;
            }
        }

        // Check linear profile
        {
            const scalarField diff = ccX-mesh.C().component(0);

            const scalar max = gMax(diff);
            const scalar min = gMin(diff);

            Info<< "Linear profile field min = " << min
                << "  max = " << max << endl;

            if (notEqual(max, 0.0, 1e-10) || notEqual(min, 0.0, 1e-10))
            {
                Info<< "Linear profile not preserved."
                    << " Min and max should both be 0.0. min:" << min
                    << " max:" << max << nl << endl;
            }
            else
            {
                Info<< "Linear profile mapping check OK" << nl << endl;
            }
        }

        // Check face field mapping
        if (surfaceOne.size())
        {
            const scalar max = gMax(surfaceOne.primitiveField());
            const scalar min = gMin(surfaceOne.primitiveField());

            Info<< "Uniform surface field min = " << min
                << "  max = " << max << endl;

            if (notEqual(max, 1.0, 1e-10) || notEqual(min, 1.0, 1e-10))
            {
                FatalErrorInFunction
                    << "Uniform surfaceScalarField not preserved."
                    << " Min and max should both be 1.0. min:" << min
                    << " max:" << max
                    << exit(FatalError);
            }
            else
            {
                Info<< "Uniform surfaceScalarField mapping check OK" << nl
                    << endl;
            }
        }


        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "pc:" << pc.patchPatchPointConstraintPoints().size() << endl;


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
