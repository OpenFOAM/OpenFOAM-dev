/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    moveMesh

Description
    Mesh motion and topological mesh change utility.

    Executes the mover, topoChanger and distributor specified in the
    dynamicMeshDict in a time-loop.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "pimpleControl.H"
#include "vtkSurfaceWriter.H"
#include "cyclicAMIPolyPatch.H"
#include "PatchTools.H"
#include "checkGeometry.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    argList::addBoolOption
    (
        "checkAMI",
        "check AMI weights"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    Foam::word regionName;

    if (args.optionReadIfPresent("region", regionName))
    {
        Foam::Info
            << "Create mesh " << regionName << " for time = "
            << runTime.timeName() << Foam::nl << Foam::endl;
    }
    else
    {
        regionName = Foam::fvMesh::defaultRegion;
        Foam::Info
            << "Create mesh for time = "
            << runTime.timeName() << Foam::nl << Foam::endl;
    }

    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            regionName,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

    const bool checkAMI  = args.optionFound("checkAMI");

    if (checkAMI)
    {
        Info<< "Writing VTK files with weights of AMI patches." << nl << endl;
    }

    pimpleControl pimple(mesh);

    bool moveMeshOuterCorrectors
    (
        pimple.dict().lookupOrDefault<Switch>("moveMeshOuterCorrectors", false)
    );

    while (runTime.run())
    {
        // Update the mesh for topology change, mesh to mesh mapping
        mesh.update();

        runTime++;

        Info<< "Time = " << runTime.userTimeName() << endl;

        while (pimple.loop())
        {
            if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
            {
                // Move the mesh
                mesh.move();
            }
        }

        mesh.checkMesh(true);

        if (checkAMI)
        {
            writeAMIWeightsSums(mesh);
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
