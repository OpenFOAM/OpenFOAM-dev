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
    foamyHexMesh

Description
    Conformal Voronoi automatic mesh generator

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IOdictionary.H"
#include "searchableSurfaces.H"
#include "conformalVoronoiMesh.H"
#include "vtkSetWriter.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::argList::addBoolOption
    (
        "checkGeometry",
        "check all surface geometry for quality"
    );

    Foam::argList::addBoolOption
    (
        "conformationOnly",
        "conform to the initial points without any point motion"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    runTime.functionObjects().off();

    const bool checkGeometry = args.optionFound("checkGeometry");
    const bool conformationOnly = args.optionFound("conformationOnly");

    IOdictionary foamyHexMeshDict
    (
        IOobject
        (
            args.executable() + "Dict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );


    if (checkGeometry)
    {
        const searchableSurfaces allGeometry
        (
            IOobject
            (
                "cvSearchableSurfaces",
                runTime.constant(),
                "triSurface",
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            foamyHexMeshDict.subDict("geometry"),
            foamyHexMeshDict.lookupOrDefault("singleRegionName", true)
        );

        // Write some stats
        allGeometry.writeStats(List<wordList>(0), Info);
        // Check topology
        allGeometry.checkTopology(true);
        // Check geometry
        allGeometry.checkGeometry
        (
            100.0,      // max size ratio
            1e-9,       // intersection tolerance
            autoPtr<writer<scalar>>(new vtkSetWriter<scalar>()),
            0.01,       // min triangle quality
            true
        );

        return 0;
    }


    conformalVoronoiMesh::debug = true;

    Info<< "Create mesh for time = " << runTime.timeName() << nl << endl;

    conformalVoronoiMesh mesh(runTime, foamyHexMeshDict);


    if (conformationOnly)
    {
        mesh.initialiseForConformation();

        runTime++;

        mesh.writeMesh(runTime.timeName());
    }
    else
    {
        mesh.initialiseForMotion();

        while (runTime.run())
        {
            runTime++;

            Info<< nl << "Time = " << runTime.timeName() << endl;

            mesh.move();

            Info<< nl
                << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;
        }
    }


    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
