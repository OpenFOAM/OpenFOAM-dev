/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    singleCellMesh

Description
    Reads all fields and maps them to a mesh with all internal faces removed
    (singleCellFvMesh) which gets written to region "singleCell".

    Used to generate mesh and fields that can be used for boundary-only data.
    Might easily result in illegal mesh though so only look at boundaries
    in paraview.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "Time.H"
#include "ReadFields.H"
#include "singleCellFvMesh.H"
#include "timeSelector.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Name of region to create
const string singleCellName = "singleCell";


template<class GeoField>
void interpolateFields
(
    const singleCellFvMesh& scMesh,
    const PtrList<GeoField>& flds
)
{
    forAll(flds, i)
    {
        tmp<GeoField> scFld = scMesh.interpolate(flds[i]);
        GeoField* scFldPtr = scFld.ptr();
        scFldPtr->writeOpt() = IOobject::AUTO_WRITE;
        scFldPtr->store();
    }
}



int main(int argc, char *argv[])
{
    // constant, not false
    timeSelector::addOptions(true, false);
    argList::addBoolOption
    (
        "noFields",
        "do not update fields"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    const bool fields = !args.optionFound("noFields");

    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createNamedMesh.H"

    if (regionName == singleCellName)
    {
        FatalErrorInFunction
            << "Cannot convert region " << singleCellName
            << " since result would overwrite it. Please rename your region."
            << exit(FatalError);
    }

    // Create the mesh
    Info<< "Creating singleCell mesh" << nl << endl;
    autoPtr<singleCellFvMesh> scMesh
    (
        new singleCellFvMesh
        (
            IOobject
            (
                singleCellName,
                mesh.polyMesh::instance(),
                runTime,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
    // For convenience create any fv* files
    if (!exists(scMesh().fvSolution::objectPath()))
    {
        mkDir(scMesh().fvSolution::path());
        ln("../fvSolution", scMesh().fvSolution::objectPath());
    }
    if (!exists(scMesh().fvSchemes::objectPath()))
    {
        mkDir(scMesh().fvSolution::path());
        ln("../fvSchemes", scMesh().fvSchemes::objectPath());
    }


    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< nl << "Time = " << runTime.timeName() << endl;


        // Check for new mesh
        if (mesh.readUpdate() != polyMesh::UNCHANGED)
        {
            Info<< "Detected changed mesh. Recreating singleCell mesh." << endl;
            scMesh.clear(); // remove any registered objects
            scMesh.reset
            (
                new singleCellFvMesh
                (
                    IOobject
                    (
                        singleCellName,
                        mesh.polyMesh::instance(),
                        runTime,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        }


        // Read objects in time directory
        IOobjectList objects(mesh, runTime.timeName());

        if (fields) Info<< "Reading geometric fields" << nl << endl;

        #include "readVolFields.H"

        // Map and store the fields on the scMesh.
        if (fields) interpolateFields(scMesh(), vsFlds);
        if (fields) interpolateFields(scMesh(), vvFlds);
        if (fields) interpolateFields(scMesh(), vstFlds);
        if (fields) interpolateFields(scMesh(), vsymtFlds);
        if (fields) interpolateFields(scMesh(), vtFlds);


        // Write
        Info<< "Writing mesh to time " << runTime.timeName() << endl;
        scMesh().write();
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
