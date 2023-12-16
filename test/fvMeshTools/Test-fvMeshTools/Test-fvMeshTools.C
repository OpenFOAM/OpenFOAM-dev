/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2023 OpenFOAM Foundation
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
    Test-fvMeshTools

Description
    Testing of adding and removing patches.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "fvMeshTools.H"
#include "wallPolyPatch.H"
#include "processorFvPatchField.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote("Test patch manipulation");

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTimeNoFunctionObjects.H"
    #include "createNamedMesh.H"

    // Read objects in time directory
    IOobjectList objects(mesh, runTime.name());

    Info<< "Reading geometric fields" << nl << endl;

    const bool fields = true;
    #include "readVolFields.H"
    #include "readSurfaceFields.H"
    #include "readPointFields.H"

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // Add/insert a (global) wall patch
    {
        wallPolyPatch pp
        (
            "myWall",
            0,          // dummy
            0,          // dummy
            0,          // dummy
            pbm,
            wallPolyPatch::typeName
        );

        label newPatchi = fvMeshTools::addPatch
        (
            mesh,
            pp,
            dictionary(),   // no specialised patch fields
            calculatedFvPatchField<scalar>::typeName,
            true            // parallel sync'ed addition
        );

        Info<< "Inserted patch " << mesh.boundaryMesh()[newPatchi].name()
            << " type " << mesh.boundaryMesh()[newPatchi].type()
            << " at index " << newPatchi << endl;

        runTime++;
        mesh.setInstance(runTime.name());
        Info<< "Writing mesh with added patch to " << runTime.name()
            << endl;
        mesh.write();
    }

    // Remove a (zero-sized!) patch everywhere
    const label removei = 0;
    if (!isA<processorPolyPatch>(pbm[removei]) && pbm[removei].size() == 0)
    {
        Info<< "Removing patch " << pbm[removei].name() << endl;

        labelList oldToNew(pbm.size());
        for (label i = 0; i < removei; i++)
        {
            oldToNew[i] = i;
        }
        oldToNew[removei] = pbm.size()-1;
        for (label i = removei+1; i < oldToNew.size(); i++)
        {
            oldToNew[i] = i-1;
        }
        fvMeshTools::reorderPatches(mesh, oldToNew, pbm.size()-1, true);

        runTime++;
        mesh.setInstance(runTime.name());
        Info<< "Writing mesh with removed patch to " << runTime.name()
            << endl;
        mesh.write();
    }

    // Add a pair of processor patches
    if (Pstream::parRun())
    {
        word newPatchName;

        if (Pstream::myProcNo() == 0 || Pstream::myProcNo() == 1)
        {
            const label nbrProcNo = (1-Pstream::myProcNo());
            newPatchName =
                processorPolyPatch::newName(Pstream::myProcNo(), nbrProcNo)
              + "_extra";

            dictionary dict;
            dict.add("myProcNo", Pstream::myProcNo());
            dict.add("neighbProcNo", nbrProcNo);
            dict.add("startFace", 0);
            dict.add("nFaces", 0);

            processorPolyPatch pp
            (
                newPatchName,
                dict,
                0,          // dummy index
                pbm,
                processorPolyPatch::typeName
            );

            label newPatchi = fvMeshTools::addPatch
            (
                mesh,
                pp,
                dictionary(),   // no specialised patch fields
                processorFvPatchField<scalar>::typeName,
                false            // parallel sync'ed addition
            );

            Pout<< "Inserted patch " << mesh.boundaryMesh()[newPatchi].name()
                << " type " << mesh.boundaryMesh()[newPatchi].type()
                << " at index " << newPatchi << endl;
        }

        runTime++;
        mesh.setInstance(runTime.name());
        Info<< "Writing mesh with added (local) patch to "
            << runTime.name() << endl;
        mesh.write();

        // Remove the added patch
        if (newPatchName.size())
        {
            label removei = pbm.findIndex(newPatchName);
            if (removei == -1)
            {
                FatalErrorInFunction << "Problem" << exit(FatalError);
            }
            Pout<< "Removing patch " << pbm[removei].name() << endl;

            labelList oldToNew(pbm.size());
            for (label i = 0; i < removei; i++)
            {
                oldToNew[i] = i;
            }
            oldToNew[removei] = pbm.size()-1;
            for (label i = removei+1; i < oldToNew.size(); i++)
            {
                oldToNew[i] = i-1;
            }
            fvMeshTools::reorderPatches(mesh, oldToNew, pbm.size()-1, false);
        }

        runTime++;
        mesh.setInstance(runTime.name());
        Info<< "Writing mesh with removed (local) patch to "
            << runTime.name() << endl;
        mesh.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
