/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    reconstructPar

Description
    Reconstructs fields of a case that is decomposed for parallel
    execution of OpenFOAM.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "IOobjectList.H"
#include "processorRunTimes.H"
#include "domainDecomposition.H"
#include "fvFieldReconstructor.H"
#include "pointFieldReconstructor.H"
#include "reconstructLagrangian.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

bool haveAllTimes
(
    const HashSet<word>& masterTimeDirSet,
    const instantList& timeDirs
)
{
    // Loop over all times
    forAll(timeDirs, timei)
    {
        if (!masterTimeDirSet.found(timeDirs[timei].name()))
        {
            return false;
        }
    }
    return true;
}


void writeDecomposition(const domainDecomposition& meshes)
{
    // Write as volScalarField::Internal for postprocessing.
    volScalarField::Internal cellProc
    (
        IOobject
        (
            "cellProc",
            meshes.completeMesh().time().name(),
            meshes.completeMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        meshes.completeMesh(),
        dimless,
        scalarField(scalarList(meshes.cellProc()))
    );

    cellProc.write();

    Info<< "Wrote decomposition as volScalarField::Internal to "
        << cellProc.name() << " for use in postprocessing."
        << nl << endl;
}


}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Reconstruct fields of a parallel case"
    );

    // Enable -constant ... if someone really wants it
    // Enable -withZero to prevent accidentally trashing the initial fields
    timeSelector::addOptions(true, true);
    argList::noParallel();
    #include "addRegionOption.H"
    #include "addAllRegionsOption.H"
    argList::addBoolOption
    (
        "cellProc",
        "write cell processor indices as a volScalarField::Internal for "
        "post-processing."
    );
    argList::addOption
    (
        "fields",
        "list",
        "specify a list of fields to be reconstructed. Eg, '(U T p)' - "
        "regular expressions not currently supported"
    );
    argList::addBoolOption
    (
        "noFields",
        "skip reconstructing fields"
    );
    argList::addOption
    (
        "lagrangianFields",
        "list",
        "specify a list of lagrangian fields to be reconstructed. Eg, '(U d)' -"
        "regular expressions not currently supported, "
        "positions always included."
    );
    argList::addBoolOption
    (
        "noLagrangian",
        "skip reconstructing lagrangian positions and fields"
    );
    argList::addBoolOption
    (
        "noSets",
        "skip reconstructing cellSets, faceSets, pointSets"
    );
    argList::addBoolOption
    (
        "newTimes",
        "only reconstruct new times (i.e. that do not exist already)"
    );

    #include "setRootCase.H"

    const bool writeCellProc = args.optionFound("cellProc");

    HashSet<word> selectedFields;
    if (args.optionFound("fields"))
    {
        args.optionLookup("fields")() >> selectedFields;
    }

    const bool noFields = args.optionFound("noFields");

    if (noFields)
    {
        Info<< "Skipping reconstructing fields"
            << nl << endl;
    }

    const bool noLagrangian = args.optionFound("noLagrangian");

    if (noLagrangian)
    {
        Info<< "Skipping reconstructing lagrangian positions and fields"
            << nl << endl;
    }

    const bool noReconstructSets = args.optionFound("noSets");

    if (noReconstructSets)
    {
        Info<< "Skipping reconstructing cellSets, faceSets and pointSets"
            << nl << endl;
    }

    HashSet<word> selectedLagrangianFields;
    if (args.optionFound("lagrangianFields"))
    {
        if (noLagrangian)
        {
            FatalErrorInFunction
                << "Cannot specify noLagrangian and lagrangianFields "
                << "options together."
                << exit(FatalError);
        }

        args.optionLookup("lagrangianFields")() >> selectedLagrangianFields;
    }

    // Set time from database
    Info<< "Create time\n" << endl;
    processorRunTimes runTimes(Foam::Time::controlDictName, args);

    // Allow override of time
    const instantList times = runTimes.selectProc(args);

    const Time& runTime = runTimes.procTimes()[0];

    #include "setRegionNames.H"

    // Determine the processor count
    const label nProcs = fileHandler().nProcs
    (
        args.path(),
        regionNames[0] == polyMesh::defaultRegion
      ? word::null
      : regionNames[0]
    );

    if (!nProcs)
    {
        FatalErrorInFunction
            << "No processor* directories found"
            << exit(FatalError);
    }

    // Warn fileHandler of number of processors
    const_cast<fileOperation&>(fileHandler()).setNProcs(nProcs);

    // Note that we do not set the runTime time so it is still the
    // one set through the controlDict. The -time option
    // only affects the selected set of times from processor0.
    // - can be illogical
    // + any point motion handled through mesh.readUpdate
    if (times.empty())
    {
        WarningInFunction << "No times selected" << endl;
        exit(1);
    }

    // Get current times if -newTimes
    const bool newTimes   = args.optionFound("newTimes");
    instantList masterTimeDirs;
    if (newTimes)
    {
        masterTimeDirs = runTimes.completeTime().times();
    }
    HashSet<word> masterTimeDirSet(2*masterTimeDirs.size());
    forAll(masterTimeDirs, i)
    {
        masterTimeDirSet.insert(masterTimeDirs[i].name());
    }
    if
    (
        newTimes
     && regionNames.size() == 1
     && regionNames[0] == fvMesh::defaultRegion
     && haveAllTimes(masterTimeDirSet, times)
    )
    {
        Info<< "All times already reconstructed.\n\nEnd\n" << endl;
        return 0;
    }

    // Reconstruct all regions
    forAll(regionNames, regioni)
    {
        const word& regionName = regionNames[regioni];

        const word& regionDir =
            regionName == polyMesh::defaultRegion
          ? word::null
          : regionName;

        // Create meshes
        Info<< "\n\nReconstructing mesh " << regionName << nl << endl;
        domainDecomposition meshes(runTimes, regionName);
        if (meshes.readReconstruct(!noReconstructSets) && writeCellProc)
        {
            writeDecomposition(meshes);
            fileHandler().flush();
        }

        // Loop over all times
        forAll(times, timei)
        {
            if (newTimes && masterTimeDirSet.found(times[timei].name()))
            {
                Info<< "Skipping time " << times[timei].name()
                    << endl << endl;
                continue;
            }

            // Set the time
            runTimes.setTime(times[timei], timei);

            Info<< "Time = " << runTimes.completeTime().userTimeName()
                << nl << endl;

            // Update the meshes
            const fvMesh::readUpdateState state =
                meshes.readUpdateReconstruct();

            // Write the mesh out, if necessary
            if (state != fvMesh::UNCHANGED)
            {
                meshes.writeComplete(!noReconstructSets);
            }

            // Write the decomposition, if necessary
            if
            (
                writeCellProc
             && meshes.completeMesh().facesInstance()
             == runTimes.completeTime().name()
            )
            {
                writeDecomposition(meshes);
                fileHandler().flush();
            }

            // Get list of objects from processor0 database
            IOobjectList objects
            (
                meshes.procMeshes()[0],
                runTimes.procTimes()[0].name()
            );

            if (!noFields)
            {
                // If there are any FV fields, reconstruct them
                Info<< "Reconstructing FV fields" << nl << endl;

                fvFieldReconstructor fvReconstructor
                (
                    meshes.completeMesh(),
                    meshes.procMeshes(),
                    meshes.procFaceAddressing(),
                    meshes.procCellAddressing(),
                    meshes.procFaceAddressingBf()
                );

                fvReconstructor.reconstructFvVolumeInternalFields<scalar>
                (
                    objects,
                    selectedFields
                );
                fvReconstructor.reconstructFvVolumeInternalFields<vector>
                (
                    objects,
                    selectedFields
                );
                fvReconstructor.reconstructFvVolumeInternalFields
                <sphericalTensor>
                (
                    objects,
                    selectedFields
                );
                fvReconstructor.reconstructFvVolumeInternalFields<symmTensor>
                (
                    objects,
                    selectedFields
                );
                fvReconstructor.reconstructFvVolumeInternalFields<tensor>
                (
                    objects,
                    selectedFields
                );

                fvReconstructor.reconstructFvVolumeFields<scalar>
                (
                    objects,
                    selectedFields
                );
                fvReconstructor.reconstructFvVolumeFields<vector>
                (
                    objects,
                    selectedFields
                );
                fvReconstructor.reconstructFvVolumeFields<sphericalTensor>
                (
                    objects,
                    selectedFields
                );
                fvReconstructor.reconstructFvVolumeFields<symmTensor>
                (
                    objects,
                    selectedFields
                );
                fvReconstructor.reconstructFvVolumeFields<tensor>
                (
                    objects,
                    selectedFields
                );

                fvReconstructor.reconstructFvSurfaceFields<scalar>
                (
                    objects,
                    selectedFields
                );
                fvReconstructor.reconstructFvSurfaceFields<vector>
                (
                    objects,
                    selectedFields
                );
                fvReconstructor.reconstructFvSurfaceFields<sphericalTensor>
                (
                    objects,
                    selectedFields
                );
                fvReconstructor.reconstructFvSurfaceFields<symmTensor>
                (
                    objects,
                    selectedFields
                );
                fvReconstructor.reconstructFvSurfaceFields<tensor>
                (
                    objects,
                    selectedFields
                );

                if (fvReconstructor.nReconstructed() == 0)
                {
                    Info<< "No FV fields" << nl << endl;
                }
            }

            if (!noFields)
            {
                Info<< "Reconstructing point fields" << nl << endl;

                const pointMesh& completePMesh =
                    pointMesh::New(meshes.completeMesh());

                pointFieldReconstructor pointReconstructor
                (
                    completePMesh,
                    meshes.procMeshes(),
                    meshes.procPointAddressing()
                );

                pointReconstructor.reconstructFields<scalar>
                (
                    objects,
                    selectedFields
                );
                pointReconstructor.reconstructFields<vector>
                (
                    objects,
                    selectedFields
                );
                pointReconstructor.reconstructFields<sphericalTensor>
                (
                    objects,
                    selectedFields
                );
                pointReconstructor.reconstructFields<symmTensor>
                (
                    objects,
                    selectedFields
                );
                pointReconstructor.reconstructFields<tensor>
                (
                    objects,
                    selectedFields
                );

                if (pointReconstructor.nReconstructed() == 0)
                {
                    Info<< "No point fields" << nl << endl;
                }
            }


            // If there are any clouds, reconstruct them.
            // The problem is that a cloud of size zero will not get written so
            // in pass 1 we determine the cloud names and per cloud name the
            // fields. Note that the fields are stored as IOobjectList from
            // the first processor that has them. They are in pass2 only used
            // for name and type (scalar, vector etc).

            if (!noLagrangian)
            {
                HashTable<IOobjectList> cloudObjects;

                forAll(runTimes.procTimes(), proci)
                {
                    fileName lagrangianDir
                    (
                        fileHandler().filePath
                        (
                            runTimes.procTimes()[proci].timePath()
                           /regionDir
                           /cloud::prefix
                        )
                    );

                    fileNameList cloudDirs;
                    if (!lagrangianDir.empty())
                    {
                        cloudDirs = fileHandler().readDir
                        (
                            lagrangianDir,
                            fileType::directory
                        );
                    }

                    forAll(cloudDirs, i)
                    {
                        // Check if we already have cloud objects for this
                        // cloudname
                        HashTable<IOobjectList>::const_iterator iter =
                            cloudObjects.find(cloudDirs[i]);

                        if (iter == cloudObjects.end())
                        {
                            // Do local scan for valid cloud objects
                            IOobjectList sprayObjs
                            (
                                meshes.procMeshes()[proci],
                                runTimes.procTimes()[proci].name(),
                                cloud::prefix/cloudDirs[i]
                            );

                            IOobject* positionsPtr =
                                sprayObjs.lookup(word("positions"));

                            if (positionsPtr)
                            {
                                cloudObjects.insert(cloudDirs[i], sprayObjs);
                            }
                        }
                    }
                }

                if (cloudObjects.size())
                {
                    // Pass2: reconstruct the cloud
                    forAllConstIter(HashTable<IOobjectList>, cloudObjects, iter)
                    {
                        const word cloudName =
                            string::validate<word>(iter.key());

                        // Objects (on arbitrary processor)
                        const IOobjectList& sprayObjs = iter();

                        Info<< "Reconstructing lagrangian fields for cloud "
                            << cloudName << nl << endl;

                        reconstructLagrangianPositions
                        (
                            meshes.completeMesh(),
                            cloudName,
                            meshes.procMeshes(),
                            meshes.procFaceAddressing(),
                            meshes.procCellAddressing()
                        );
                        reconstructLagrangianFields<label>
                        (
                            cloudName,
                            meshes.completeMesh(),
                            meshes.procMeshes(),
                            sprayObjs,
                            selectedLagrangianFields
                        );
                        reconstructLagrangianFieldFields<label>
                        (
                            cloudName,
                            meshes.completeMesh(),
                            meshes.procMeshes(),
                            sprayObjs,
                            selectedLagrangianFields
                        );
                        reconstructLagrangianFields<scalar>
                        (
                            cloudName,
                            meshes.completeMesh(),
                            meshes.procMeshes(),
                            sprayObjs,
                            selectedLagrangianFields
                        );
                        reconstructLagrangianFieldFields<scalar>
                        (
                            cloudName,
                            meshes.completeMesh(),
                            meshes.procMeshes(),
                            sprayObjs,
                            selectedLagrangianFields
                        );
                        reconstructLagrangianFields<vector>
                        (
                            cloudName,
                            meshes.completeMesh(),
                            meshes.procMeshes(),
                            sprayObjs,
                            selectedLagrangianFields
                        );
                        reconstructLagrangianFieldFields<vector>
                        (
                            cloudName,
                            meshes.completeMesh(),
                            meshes.procMeshes(),
                            sprayObjs,
                            selectedLagrangianFields
                        );
                        reconstructLagrangianFields<sphericalTensor>
                        (
                            cloudName,
                            meshes.completeMesh(),
                            meshes.procMeshes(),
                            sprayObjs,
                            selectedLagrangianFields
                        );
                        reconstructLagrangianFieldFields<sphericalTensor>
                        (
                            cloudName,
                            meshes.completeMesh(),
                            meshes.procMeshes(),
                            sprayObjs,
                            selectedLagrangianFields
                        );
                        reconstructLagrangianFields<symmTensor>
                        (
                            cloudName,
                            meshes.completeMesh(),
                            meshes.procMeshes(),
                            sprayObjs,
                            selectedLagrangianFields
                        );
                        reconstructLagrangianFieldFields<symmTensor>
                        (
                            cloudName,
                            meshes.completeMesh(),
                            meshes.procMeshes(),
                            sprayObjs,
                            selectedLagrangianFields
                        );
                        reconstructLagrangianFields<tensor>
                        (
                            cloudName,
                            meshes.completeMesh(),
                            meshes.procMeshes(),
                            sprayObjs,
                            selectedLagrangianFields
                        );
                        reconstructLagrangianFieldFields<tensor>
                        (
                            cloudName,
                            meshes.completeMesh(),
                            meshes.procMeshes(),
                            sprayObjs,
                            selectedLagrangianFields
                        );
                    }
                }
                else
                {
                    Info<< "No lagrangian fields" << nl << endl;
                }
            }

            // If there is a "uniform" directory in the time region
            // directory copy from the master processor
            {
                fileName uniformDir0
                (
                    fileHandler().filePath
                    (
                        runTimes.procTimes()[0].timePath()/regionDir/"uniform"
                    )
                );

                if (!uniformDir0.empty() && fileHandler().isDir(uniformDir0))
                {
                    fileHandler().cp
                    (
                        uniformDir0,
                        runTimes.completeTime().timePath()/regionDir
                    );
                }
            }

            // For the first region of a multi-region case additionally
            // copy the "uniform" directory in the time directory
            if (regioni == 0 && regionDir != word::null)
            {
                fileName uniformDir0
                (
                    fileHandler().filePath
                    (
                        runTimes.procTimes()[0].timePath()/"uniform"
                    )
                );

                if (!uniformDir0.empty() && fileHandler().isDir(uniformDir0))
                {
                    fileHandler().cp
                    (
                        uniformDir0,
                        runTimes.completeTime().timePath()
                    );
                }
            }
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
