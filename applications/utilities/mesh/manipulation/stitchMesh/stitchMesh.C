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
    stitchMesh

Description
    'Stitches' a mesh.

    Takes a mesh and two patches and merges the faces on the two patches
    (if geometrically possible) so the faces become internal.

    Can do
    - 'perfect' match: faces and points on patches align exactly. Order might
    be different though.
    - 'integral' match: where the surfaces on both patches exactly
    match but the individual faces not
    - 'partial' match: where the non-overlapping part of the surface remains
    in the respective patch.

    Note : Is just a front-end to perfectInterface/slidingInterface.

    Comparable to running a meshModifier of the form
    (if masterPatch is called "M" and slavePatch "S"):

    \verbatim
    couple
    {
        type                    slidingInterface;
        masterFaceZoneName      MSMasterZone
        slaveFaceZoneName       MSSlaveZone
        cutPointZoneName        MSCutPointZone
        cutFaceZoneName         MSCutFaceZone
        masterPatchName         M;
        slavePatchName          S;
        typeOfMatch             partial or integral
    }
    \endverbatim


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "polyTopoChanger.H"
#include "mapPolyMesh.H"
#include "ListOps.H"
#include "slidingInterface.H"
#include "perfectInterface.H"
#include "IOobjectList.H"
#include "ReadFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label addPointZone(const polyMesh& mesh, const word& name)
{
    label zoneID = mesh.pointZones().findZoneID(name);

    if (zoneID != -1)
    {
        Info<< "Reusing existing pointZone "
            << mesh.pointZones()[zoneID].name()
            << " at index " << zoneID << endl;
    }
    else
    {
        pointZoneMesh& pointZones = const_cast<polyMesh&>(mesh).pointZones();
        zoneID = pointZones.size();
        Info<< "Adding pointZone " << name << " at index " << zoneID << endl;

        pointZones.setSize(zoneID+1);
        pointZones.set
        (
            zoneID,
            new pointZone
            (
                name,
                labelList(0),
                zoneID,
                pointZones
            )
        );
    }
    return zoneID;
}


label addFaceZone(const polyMesh& mesh, const word& name)
{
    label zoneID = mesh.faceZones().findZoneID(name);

    if (zoneID != -1)
    {
        Info<< "Reusing existing faceZone " << mesh.faceZones()[zoneID].name()
            << " at index " << zoneID << endl;
    }
    else
    {
        faceZoneMesh& faceZones = const_cast<polyMesh&>(mesh).faceZones();
        zoneID = faceZones.size();
        Info<< "Adding faceZone " << name << " at index " << zoneID << endl;

        faceZones.setSize(zoneID+1);
        faceZones.set
        (
            zoneID,
            new faceZone
            (
                name,
                labelList(0),
                boolList(),
                zoneID,
                faceZones
            )
        );
    }
    return zoneID;
}


label addCellZone(const polyMesh& mesh, const word& name)
{
    label zoneID = mesh.cellZones().findZoneID(name);

    if (zoneID != -1)
    {
        Info<< "Reusing existing cellZone " << mesh.cellZones()[zoneID].name()
            << " at index " << zoneID << endl;
    }
    else
    {
        cellZoneMesh& cellZones = const_cast<polyMesh&>(mesh).cellZones();
        zoneID = cellZones.size();
        Info<< "Adding cellZone " << name << " at index " << zoneID << endl;

        cellZones.setSize(zoneID+1);
        cellZones.set
        (
            zoneID,
            new cellZone
            (
                name,
                labelList(0),
                zoneID,
                cellZones
            )
        );
    }
    return zoneID;
}


// Checks whether patch present
void checkPatch(const polyBoundaryMesh& bMesh, const word& name)
{
    const label patchI = bMesh.findPatchID(name);

    if (patchI == -1)
    {
        FatalErrorIn("checkPatch(const polyBoundaryMesh&, const word&)")
            << "Cannot find patch " << name << endl
            << "It should be present and of non-zero size" << endl
            << "Valid patches are " << bMesh.names()
            << exit(FatalError);
    }

    if (bMesh[patchI].empty())
    {
        FatalErrorIn("checkPatch(const polyBoundaryMesh&, const word&)")
            << "Patch " << name << " is present but zero size"
            << exit(FatalError);
    }
}



int main(int argc, char *argv[])
{
    argList::addNote
    (
        "merge the faces on the specified patches (if geometrically possible)\n"
        "so the faces become internal"
    );

    argList::noParallel();
    #include "addOverwriteOption.H"
    #include "addRegionOption.H"

    argList::validArgs.append("masterPatch");
    argList::validArgs.append("slavePatch");

    argList::addBoolOption
    (
        "partial",
        "couple partially overlapping patches"
    );
    argList::addBoolOption
    (
        "perfect",
        "couple perfectly aligned patches"
    );
    argList::addOption
    (
        "toleranceDict",
        "file",
        "dictionary file with tolerances"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();
    #include "createNamedMesh.H"

    const word oldInstance = mesh.pointsInstance();

    const word masterPatchName = args[1];
    const word slavePatchName  = args[2];

    const bool partialCover = args.optionFound("partial");
    const bool perfectCover = args.optionFound("perfect");
    const bool overwrite    = args.optionFound("overwrite");

    if (partialCover && perfectCover)
    {
        FatalErrorIn(args.executable())
            << "Cannot supply both partial and perfect." << endl
            << "Use perfect match option if the patches perfectly align"
            << " (both vertex positions and face centres)" << endl
            << exit(FatalError);
    }


    const word mergePatchName(masterPatchName + slavePatchName);
    const word cutZoneName(mergePatchName + "CutFaceZone");

    slidingInterface::typeOfMatch tom = slidingInterface::INTEGRAL;

    if (partialCover)
    {
        Info<< "Coupling partially overlapping patches "
            << masterPatchName << " and " << slavePatchName << nl
            << "Resulting internal faces will be in faceZone " << cutZoneName
            << nl
            << "Any uncovered faces will remain in their patch"
            << endl;

        tom = slidingInterface::PARTIAL;
    }
    else if (perfectCover)
    {
        Info<< "Coupling perfectly aligned patches "
            << masterPatchName << " and " << slavePatchName << nl
            << "Resulting (internal) faces will be in faceZone " << cutZoneName
            << nl << nl
            << "Note: both patches need to align perfectly." << nl
            << "Both the vertex"
            << " positions and the face centres need to align to within" << nl
            << "a tolerance given by the minimum edge length on the patch"
            << endl;
    }
    else
    {
        Info<< "Coupling patches " << masterPatchName << " and "
            << slavePatchName << nl
            << "Resulting (internal) faces will be in faceZone " << cutZoneName
            << nl << nl
            << "Note: the overall area covered by both patches should be"
            << " identical (\"integral\" interface)." << endl
            << "If this is not the case use the -partial option" << nl << endl;
    }

    // set up the tolerances for the sliding mesh
    dictionary slidingTolerances;
    if (args.options().found("toleranceDict"))
    {
        IOdictionary toleranceFile
        (
            IOobject
            (
                args.options()["toleranceDict"],
                runTime.constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );
        slidingTolerances += toleranceFile;
    }

    // Check for non-empty master and slave patches
    checkPatch(mesh.boundaryMesh(), masterPatchName);
    checkPatch(mesh.boundaryMesh(), slavePatchName);

    // Create and add face zones and mesh modifiers

    // Master patch
    const polyPatch& masterPatch = mesh.boundaryMesh()[masterPatchName];

    // Make list of masterPatch faces
    labelList isf(masterPatch.size());

    forAll(isf, i)
    {
        isf[i] = masterPatch.start() + i;
    }

    polyTopoChanger stitcher(mesh);
    stitcher.setSize(1);

    mesh.pointZones().clearAddressing();
    mesh.faceZones().clearAddressing();
    mesh.cellZones().clearAddressing();

    if (perfectCover)
    {
        // Add empty zone for resulting internal faces
        label cutZoneID = addFaceZone(mesh, cutZoneName);

        mesh.faceZones()[cutZoneID].resetAddressing
        (
            isf,
            boolList(masterPatch.size(), false)
        );

        // Add the perfect interface mesh modifier
        stitcher.set
        (
            0,
            new perfectInterface
            (
                "couple",
                0,
                stitcher,
                cutZoneName,
                masterPatchName,
                slavePatchName
            )
        );
    }
    else
    {
        label pointZoneID = addPointZone(mesh, mergePatchName + "CutPointZone");
        mesh.pointZones()[pointZoneID] = labelList(0);

        label masterZoneID = addFaceZone(mesh, mergePatchName + "MasterZone");

        mesh.faceZones()[masterZoneID].resetAddressing
        (
            isf,
            boolList(masterPatch.size(), false)
        );

        // Slave patch
        const polyPatch& slavePatch = mesh.boundaryMesh()[slavePatchName];

        labelList osf(slavePatch.size());

        forAll(osf, i)
        {
            osf[i] = slavePatch.start() + i;
        }

        label slaveZoneID = addFaceZone(mesh, mergePatchName + "SlaveZone");
        mesh.faceZones()[slaveZoneID].resetAddressing
        (
            osf,
            boolList(slavePatch.size(), false)
        );

        // Add empty zone for cut faces
        label cutZoneID = addFaceZone(mesh, cutZoneName);
        mesh.faceZones()[cutZoneID].resetAddressing
        (
            labelList(0),
            boolList(0, false)
        );


        // Add the sliding interface mesh modifier
        stitcher.set
        (
            0,
            new slidingInterface
            (
                "couple",
                0,
                stitcher,
                mergePatchName + "MasterZone",
                mergePatchName + "SlaveZone",
                mergePatchName + "CutPointZone",
                cutZoneName,
                masterPatchName,
                slavePatchName,
                tom,                    // integral or partial
                true                    // couple/decouple mode
            )
        );
        static_cast<slidingInterface&>(stitcher[0]).setTolerances
        (
            slidingTolerances,
            true
        );
    }


    // Search for list of objects for this time
    IOobjectList objects(mesh, runTime.timeName());

    // Read all current fvFields so they will get mapped
    Info<< "Reading all current volfields" << endl;
    PtrList<volScalarField> volScalarFields;
    ReadFields(mesh, objects, volScalarFields);

    PtrList<volVectorField> volVectorFields;
    ReadFields(mesh, objects, volVectorFields);

    PtrList<volSphericalTensorField> volSphericalTensorFields;
    ReadFields(mesh, objects, volSphericalTensorFields);

    PtrList<volSymmTensorField> volSymmTensorFields;
    ReadFields(mesh, objects, volSymmTensorFields);

    PtrList<volTensorField> volTensorFields;
    ReadFields(mesh, objects, volTensorFields);

    //- uncomment if you want to interpolate surface fields (usually bad idea)
    //Info<< "Reading all current surfaceFields" << endl;
    //PtrList<surfaceScalarField> surfaceScalarFields;
    //ReadFields(mesh, objects, surfaceScalarFields);
    //
    //PtrList<surfaceVectorField> surfaceVectorFields;
    //ReadFields(mesh, objects, surfaceVectorFields);
    //
    //PtrList<surfaceTensorField> surfaceTensorFields;
    //ReadFields(mesh, objects, surfaceTensorFields);

    if (!overwrite)
    {
        runTime++;
    }

    // Execute all polyMeshModifiers
    autoPtr<mapPolyMesh> morphMap = stitcher.changeMesh(true);

    mesh.movePoints(morphMap->preMotionPoints());

    // Write mesh
    if (overwrite)
    {
        mesh.setInstance(oldInstance);
        stitcher.instance() = oldInstance;
    }
    Info<< nl << "Writing polyMesh to time " << runTime.timeName() << endl;

    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    // Bypass runTime write (since only writes at outputTime)
    if
    (
       !runTime.objectRegistry::writeObject
        (
            runTime.writeFormat(),
            IOstream::currentVersion,
            runTime.writeCompression()
        )
    )
    {
        FatalErrorIn(args.executable())
            << "Failed writing polyMesh."
            << exit(FatalError);
    }

    mesh.faceZones().write();
    mesh.pointZones().write();
    mesh.cellZones().write();

    // Write fields
    runTime.write();

    Info<< nl << "end" << endl;

    return 0;
}


// ************************************************************************* //
