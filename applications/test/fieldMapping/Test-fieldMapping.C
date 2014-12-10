/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Test-fieldMapping

Description
    Test app for mapping of fields.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "Time.H"
#include "OFstream.H"
#include "meshTools.H"
#include "removeFaces.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "fvcDiv.H"
#include "Random.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool notEqual(const scalar s1, const scalar s2, const scalar tol)
{
    return mag(s1-s2) > tol;
}


// Main program:

int main(int argc, char *argv[])
{
#   include "addTimeOptions.H"
    argList::validArgs.append("inflate (true|false)");
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    const Switch inflate(args.args()[1]);

    if (inflate)
    {
        Info<< "Deleting cells using inflation/deflation" << nl << endl;
    }
    else
    {
        Info<< "Deleting cells, introducing points at new position" << nl
            << endl;
    }


    Random rndGen(0);



    // Test mapping
    // ------------
    // Mapping is volume averaged


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




    // Force allocation of V. Important for any mesh changes since otherwise
    // old time volumes are not stored
    const scalar totalVol = gSum(mesh.V());

    // Face removal engine. No checking for not merging boundary faces.
    removeFaces faceRemover(mesh, GREAT);


    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (!mesh.nInternalFaces())
        {
            break;
        }

        // Remove face
        label candidateFaceI = rndGen.integer(0, mesh.nInternalFaces()-1);
        Info<< "Wanting to delete face " << mesh.faceCentres()[candidateFaceI]
            << nl << endl;

        labelList candidates(1, candidateFaceI);


        // Get compatible set of faces and connected sets of cells.
        labelList cellRegion;
        labelList cellRegionMaster;
        labelList facesToRemove;

        faceRemover.compatibleRemoves
        (
            candidates,
            cellRegion,
            cellRegionMaster,
            facesToRemove
        );

        // Topo changes container
        polyTopoChange meshMod(mesh);

        // Insert mesh refinement into polyTopoChange.
        faceRemover.setRefinement
        (
            facesToRemove,
            cellRegion,
            cellRegionMaster,
            meshMod
        );

        // Change mesh and inflate
        Info<< "Actually changing mesh" << nl << endl;
        autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, inflate);

        Info<< "Mapping fields" << nl << endl;
        mesh.updateMesh(morphMap);

        // Move mesh (since morphing does not do this)
        if (morphMap().hasMotionPoints())
        {
            Info<< "Moving mesh" << nl << endl;
            mesh.movePoints(morphMap().preMotionPoints());
        }

        // Update numbering of cells/vertices.
        faceRemover.updateMesh(morphMap);


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
                FatalErrorIn(args.executable())
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
                FatalErrorIn(args.executable())
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
                FatalErrorIn(args.executable())
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
                FatalErrorIn(args.executable())
                    << "Linear profile not preserved."
                    << " Min and max should both be 0.0. min:" << min
                    << " max:" << max
                    << exit(FatalError);
            }
            else
            {
                Info<< "Linear profile mapping check OK" << nl << endl;
            }
        }

        // Check face field mapping
        if (surfaceOne.size())
        {
            const scalar max = gMax(surfaceOne.internalField());
            const scalar min = gMin(surfaceOne.internalField());

            Info<< "Uniform surface field min = " << min
                << "  max = " << max << endl;

            if (notEqual(max, 1.0, 1e-10) || notEqual(min, 1.0, 1e-10))
            {
                FatalErrorIn(args.executable())
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

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
