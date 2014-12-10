/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    removeFaces

Description
    Utility to remove faces (combines cells on both sides).

    Takes faceSet of candidates for removal and writes faceSet with faces that
    will actually be removed. (because e.g. would cause two faces between the
    same cells). See removeFaces in dynamicMesh library for constraints.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyTopoChange.H"
#include "faceSet.H"
#include "removeFaces.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
#   include "addOverwriteOption.H"
    argList::validArgs.append("faceSet");

#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
#   include "createMesh.H"
    const word oldInstance = mesh.pointsInstance();

    const word setName = args[1];
    const bool overwrite = args.optionFound("overwrite");

    // Read faces
    faceSet candidateSet(mesh, setName);

    Pout<< "Read " << candidateSet.size() << " faces to remove" << nl
        << endl;


    labelList candidates(candidateSet.toc());

    // Face removal engine. No checking for not merging boundary faces.
    removeFaces faceRemover(mesh, 2);

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

    {
        faceSet compatibleRemoves(mesh, "compatibleRemoves", facesToRemove);

        Pout<< "Original faces to be removed:" << candidateSet.size() << nl
            << "New faces to be removed:" << compatibleRemoves.size() << nl
            << endl;

        Pout<< "Writing new faces to be removed to faceSet "
            << compatibleRemoves.instance()
              /compatibleRemoves.local()
              /compatibleRemoves.name()
            << endl;

        compatibleRemoves.write();
    }


    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    // Read vol fields.
    PtrList<volScalarField> vsFlds;
    ReadFields(mesh, objects, vsFlds);

    PtrList<volVectorField> vvFlds;
    ReadFields(mesh, objects, vvFlds);

    PtrList<volSphericalTensorField> vstFlds;
    ReadFields(mesh, objects, vstFlds);

    PtrList<volSymmTensorField> vsymtFlds;
    ReadFields(mesh, objects, vsymtFlds);

    PtrList<volTensorField> vtFlds;
    ReadFields(mesh, objects, vtFlds);

    // Read surface fields.
    PtrList<surfaceScalarField> ssFlds;
    ReadFields(mesh, objects, ssFlds);

    PtrList<surfaceVectorField> svFlds;
    ReadFields(mesh, objects, svFlds);

    PtrList<surfaceSphericalTensorField> sstFlds;
    ReadFields(mesh, objects, sstFlds);

    PtrList<surfaceSymmTensorField> ssymtFlds;
    ReadFields(mesh, objects, ssymtFlds);

    PtrList<surfaceTensorField> stFlds;
    ReadFields(mesh, objects, stFlds);


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

    autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

    mesh.updateMesh(morphMap);

    // Move mesh (since morphing does not do this)
    if (morphMap().hasMotionPoints())
    {
        mesh.movePoints(morphMap().preMotionPoints());
    }

    // Update numbering of cells/vertices.
    faceRemover.updateMesh(morphMap);

    if (!overwrite)
    {
        runTime++;
    }
    else
    {
        mesh.setInstance(oldInstance);
    }

    // Take over refinement levels and write to new time directory.
    Pout<< "Writing mesh to time " << runTime.timeName() << endl;
    mesh.write();

    Pout<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
