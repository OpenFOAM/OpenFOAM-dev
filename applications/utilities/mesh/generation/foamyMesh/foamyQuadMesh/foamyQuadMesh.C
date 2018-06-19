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
    foamyQuadMesh

Description
    Conformal-Voronoi 2D extruding automatic mesher with grid or read
    initial points and point position relaxation with optional
    "squarification".

\*---------------------------------------------------------------------------*/

#include "CV2D.H"
#include "argList.H"

#include "MeshedSurfaces.H"
#include "shortEdgeFilter2D.H"
#include "extrude2DMesh.H"
#include "polyMesh.H"
#include "patchToPoly2DMesh.H"
#include "extrudeModel.H"
#include "polyTopoChange.H"
#include "edgeCollapser.H"
#include "globalIndex.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validOptions.insert
    (
        "pointsFile",
        "filename"
    );

    #include "addOverwriteOption.H"

    #include "setRootCase.H"
    #include "createTime.H"

    // Read control dictionary
    // ~~~~~~~~~~~~~~~~~~~~~~~
    IOdictionary controlDict
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

    const dictionary& shortEdgeFilterDict
    (
        controlDict.subDict("shortEdgeFilter")
    );
    const dictionary& extrusionDict(controlDict.subDict("extrusion"));

    Switch extrude(extrusionDict.lookup("extrude"));
    const bool overwrite = args.optionFound("overwrite");

    // Read and triangulation
    // ~~~~~~~~~~~~~~~~~~~~~~
    CV2D mesh(runTime, controlDict);

    if (args.options().found("pointsFile"))
    {
        fileName pointFileName(IStringStream(args.options()["pointsFile"])());
        mesh.insertPoints(pointFileName);
    }
    else
    {
        mesh.insertGrid();
    }

    mesh.insertSurfacePointPairs();
    mesh.boundaryConform();

    while (runTime.loop())
    {
        Info<< nl << "Time = " << runTime.timeName() << endl;

        mesh.newPoints();
    }

    mesh.write();

    Info<< "Finished Delaunay in = " << runTime.cpuTimeIncrement() << " s."
        << endl;

    Info<< "Begin filtering short edges:" << endl;
    shortEdgeFilter2D sef(mesh, shortEdgeFilterDict);

    sef.filter();

    Info<< "Meshed surface after edge filtering :" << endl;
    sef.fMesh().writeStats(Info);

    if (mesh.meshControls().meshedSurfaceOutput())
    {
        Info<< "Write .obj file of the 2D mesh: MeshedSurface.obj" << endl;
        sef.fMesh().write("MeshedSurface.obj");
    }

    Info<< "Finished filtering in = " << runTime.cpuTimeIncrement() << " s."
        << endl;

    Info<< "Begin constructing a polyMesh:" << endl;

    patchToPoly2DMesh poly2DMesh
    (
        sef.fMesh(),
        sef.patchNames(),
        sef.patchSizes(),
        sef.mapEdgesRegion()
    );

    poly2DMesh.createMesh();

    polyMesh pMesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        xferMove(poly2DMesh.points()),
        xferMove(poly2DMesh.faces()),
        xferMove(poly2DMesh.owner()),
        xferMove(poly2DMesh.neighbour())
    );

    Info<< "Constructing patches." << endl;
    List<polyPatch*> patches(poly2DMesh.patchNames().size());
    label countPatches = 0;

    forAll(patches, patchi)
    {
        if (poly2DMesh.patchSizes()[patchi] != 0)
        {
            patches[countPatches] = new polyPatch
            (
                poly2DMesh.patchNames()[patchi],
                poly2DMesh.patchSizes()[patchi],
                poly2DMesh.patchStarts()[patchi],
                countPatches,
                pMesh.boundaryMesh(),
                word::null
            );

            countPatches++;
        }
    }
    patches.setSize(countPatches);
    pMesh.addPatches(patches);

    if (extrude)
    {
        Info<< "Begin extruding the polyMesh:" << endl;

        {
            // Point generator
            autoPtr<extrudeModel> model(extrudeModel::New(extrusionDict));

            extrude2DMesh extruder(pMesh, extrusionDict, model());

            extruder.addFrontBackPatches();

            polyTopoChange meshMod(pMesh.boundaryMesh().size());

            extruder.setRefinement(meshMod);

            autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(pMesh, false);

            pMesh.updateMesh(morphMap);
        }
    }

    if (!overwrite)
    {
        runTime++;
    }
    else
    {
        pMesh.setInstance("constant");
    }

    pMesh.write();

    Info<< "Finished extruding in = "
        << runTime.cpuTimeIncrement() << " s." << endl;

    Info<< nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
