/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    splitBaffles

Description
    Detects faces that share points (baffles) and duplicates the points to
    separate them

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyTopoChange.H"
#include "localPointRegion.H"
#include "duplicatePoints.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Detect faces that share points (baffles)\n"
        "and duplicate the points to separate them."
    );

    #include "addOverwriteOption.H"
    #include "addRegionOption.H"
    argList::addBoolOption
    (
        "fields",
        "update fields"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();
    #include "createNamedMesh.H"

    const bool overwrite  = args.optionFound("overwrite");
    const bool fields     = args.optionFound("fields");

    const word oldInstance = mesh.pointsInstance();

    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    if (fields) Info<< "Reading geometric fields" << nl << endl;

    #include "readVolFields.H"
    #include "readSurfaceFields.H"
    #include "readPointFields.H"

    // Mesh change engine
    polyTopoChange meshMod(mesh);

    // Analyse which points need to be duplicated
    localPointRegion regionSide(mesh);

    // Point duplication engine
    duplicatePoints pointDuplicator(mesh);

    // Insert topo changes
    pointDuplicator.setRefinement(regionSide, meshMod);

    if (!overwrite)
    {
        runTime++;
    }

    // Change the mesh. No inflation.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

    // Update fields
    mesh.updateMesh(map);

    // Move mesh (since morphing does not do this)
    if (map().hasMotionPoints())
    {
        mesh.movePoints(map().preMotionPoints());
    }

    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }

    Info<< "Writing mesh to time " << runTime.timeName() << endl;
    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
