/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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
    reconstructPar

Description
    Reconstructs fields of a case that is decomposed for parallel
    execution of OpenFOAM.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"

#include "fvCFD.H"
#include "IOobjectList.H"
#include "processorMeshes.H"
#include "regionProperties.H"
#include "fvFieldReconstructor.H"
#include "pointFieldReconstructor.H"
#include "reconstructLagrangian.H"

#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"

#include "hexRef8Data.H"

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
}


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
    #include "createTime.H"

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

    const wordList regionNames(selectRegionNames(args, runTime));

    // Determine the processor count
    const label nProcs =
        fileHandler().nProcs(args.path(), regionDir(regionNames[0]));

    if (!nProcs)
    {
        FatalErrorInFunction
            << "No processor* directories found"
            << exit(FatalError);
    }

    // Warn fileHandler of number of processors
    const_cast<fileOperation&>(fileHandler()).setNProcs(nProcs);

    // Create the processor databases
    PtrList<Time> databases(nProcs);

    forAll(databases, proci)
    {
        databases.set
        (
            proci,
            new Time
            (
                Time::controlDictName,
                args.rootPath(),
                args.caseName()/fileName(word("processor") + name(proci))
            )
        );
    }

    // Use the times list from the master processor
    // and select a subset based on the command-line options
    instantList timeDirs = timeSelector::select
    (
        databases[0].times(),
        args
    );

    // Note that we do not set the runTime time so it is still the
    // one set through the controlDict. The -time option
    // only affects the selected set of times from processor0.
    // - can be illogical
    // + any point motion handled through mesh.readUpdate

    if (timeDirs.empty())
    {
        WarningInFunction << "No times selected" << endl;
        exit(1);
    }

    // Get current times if -newTimes
    const bool newTimes   = args.optionFound("newTimes");
    instantList masterTimeDirs;
    if (newTimes)
    {
        masterTimeDirs = runTime.times();
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
     && haveAllTimes(masterTimeDirSet, timeDirs)
    )
    {
        Info<< "All times already reconstructed.\n\nEnd\n" << endl;
        return 0;
    }

    // Set all times on processor meshes equal to reconstructed mesh
    forAll(databases, proci)
    {
        databases[proci].setTime(runTime);
    }

    forAll(regionNames, regioni)
    {
        const word& regionName = regionNames[regioni];
        const word& regionDir = Foam::regionDir(regionName);

        Info<< "\n\nReconstructing fields for mesh " << regionName << nl
            << endl;

        fvMesh mesh
        (
            IOobject
            (
                regionName,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        );


        // Read all meshes and addressing to reconstructed mesh
        processorMeshes procMeshes(databases, regionName);


        // Check face addressing for meshes that have been decomposed
        // with a very old foam version
        #include "checkFaceAddressingComp.H"

        // Loop over all times
        forAll(timeDirs, timei)
        {
            if (newTimes && masterTimeDirSet.found(timeDirs[timei].name()))
            {
                Info<< "Skipping time " << timeDirs[timei].name()
                    << endl << endl;
                continue;
            }


            // Set time for global database
            runTime.setTime(timeDirs[timei], timei);

            Info<< "Time = " << runTime.timeName() << endl << endl;

            // Set time for all databases
            forAll(databases, proci)
            {
                databases[proci].setTime(timeDirs[timei], timei);
            }

            // Check if any new meshes need to be read.
            fvMesh::readUpdateState meshStat = mesh.readUpdate();

            fvMesh::readUpdateState procStat = procMeshes.readUpdate();

            if (procStat == fvMesh::POINTS_MOVED)
            {
                // Reconstruct the points for moving mesh cases and write
                // them out
                procMeshes.reconstructPoints(mesh);
            }
            else if (meshStat != procStat)
            {
                WarningInFunction
                    << "readUpdate for the reconstructed mesh:"
                    << meshStat << nl
                    << "readUpdate for the processor meshes  :"
                    << procStat << nl
                    << "These should be equal or your addressing"
                    << " might be incorrect."
                    << " Please check your time directories for any "
                    << "mesh directories." << endl;
            }


            // Get list of objects from processor0 database
            IOobjectList objects
            (
                procMeshes.meshes()[0],
                databases[0].timeName()
            );

            if (!noFields)
            {
                // If there are any FV fields, reconstruct them
                Info<< "Reconstructing FV fields" << nl << endl;

                fvFieldReconstructor fvReconstructor
                (
                    mesh,
                    procMeshes.meshes(),
                    procMeshes.faceProcAddressing(),
                    procMeshes.cellProcAddressing(),
                    procMeshes.boundaryProcAddressing()
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

                const pointMesh& pMesh = pointMesh::New(mesh);
                PtrList<pointMesh> pMeshes(procMeshes.meshes().size());

                forAll(pMeshes, proci)
                {
                    pMeshes.set
                    (
                        proci,
                        new pointMesh(procMeshes.meshes()[proci])
                    );
                }

                pointFieldReconstructor pointReconstructor
                (
                    pMesh,
                    pMeshes,
                    procMeshes.pointProcAddressing(),
                    procMeshes.boundaryProcAddressing()
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

                forAll(databases, proci)
                {
                    fileName lagrangianDir
                    (
                        fileHandler().filePath
                        (
                            databases[proci].timePath()
                          / regionDir
                          / cloud::prefix
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
                                procMeshes.meshes()[proci],
                                databases[proci].timeName(),
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
                        const word cloudName = string::validate<word>
                        (
                            iter.key()
                        );

                        // Objects (on arbitrary processor)
                        const IOobjectList& sprayObjs = iter();

                        Info<< "Reconstructing lagrangian fields for cloud "
                            << cloudName << nl << endl;

                        reconstructLagrangianPositions
                        (
                            mesh,
                            cloudName,
                            procMeshes.meshes(),
                            procMeshes.faceProcAddressing(),
                            procMeshes.cellProcAddressing()
                        );
                        reconstructLagrangianFields<label>
                        (
                            cloudName,
                            mesh,
                            procMeshes.meshes(),
                            sprayObjs,
                            selectedLagrangianFields
                        );
                        reconstructLagrangianFieldFields<label>
                        (
                            cloudName,
                            mesh,
                            procMeshes.meshes(),
                            sprayObjs,
                            selectedLagrangianFields
                        );
                        reconstructLagrangianFields<scalar>
                        (
                            cloudName,
                            mesh,
                            procMeshes.meshes(),
                            sprayObjs,
                            selectedLagrangianFields
                        );
                        reconstructLagrangianFieldFields<scalar>
                        (
                            cloudName,
                            mesh,
                            procMeshes.meshes(),
                            sprayObjs,
                            selectedLagrangianFields
                        );
                        reconstructLagrangianFields<vector>
                        (
                            cloudName,
                            mesh,
                            procMeshes.meshes(),
                            sprayObjs,
                            selectedLagrangianFields
                        );
                        reconstructLagrangianFieldFields<vector>
                        (
                            cloudName,
                            mesh,
                            procMeshes.meshes(),
                            sprayObjs,
                            selectedLagrangianFields
                        );
                        reconstructLagrangianFields<sphericalTensor>
                        (
                            cloudName,
                            mesh,
                            procMeshes.meshes(),
                            sprayObjs,
                            selectedLagrangianFields
                        );
                        reconstructLagrangianFieldFields<sphericalTensor>
                        (
                            cloudName,
                            mesh,
                            procMeshes.meshes(),
                            sprayObjs,
                            selectedLagrangianFields
                        );
                        reconstructLagrangianFields<symmTensor>
                        (
                            cloudName,
                            mesh,
                            procMeshes.meshes(),
                            sprayObjs,
                            selectedLagrangianFields
                        );
                        reconstructLagrangianFieldFields<symmTensor>
                        (
                            cloudName,
                            mesh,
                            procMeshes.meshes(),
                            sprayObjs,
                            selectedLagrangianFields
                        );
                        reconstructLagrangianFields<tensor>
                        (
                            cloudName,
                            mesh,
                            procMeshes.meshes(),
                            sprayObjs,
                            selectedLagrangianFields
                        );
                        reconstructLagrangianFieldFields<tensor>
                        (
                            cloudName,
                            mesh,
                            procMeshes.meshes(),
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


            if (!noReconstructSets)
            {
                // Scan to find all sets
                HashTable<label> cSetNames;
                HashTable<label> fSetNames;
                HashTable<label> pSetNames;

                forAll(procMeshes.meshes(), proci)
                {
                    const fvMesh& procMesh = procMeshes.meshes()[proci];

                    // Note: look at sets in current time only or between
                    // mesh and current time?. For now current time. This will
                    // miss out on sets in intermediate times that have not
                    // been reconstructed.
                    IOobjectList objects
                    (
                        procMesh,
                        databases[0].timeName(),    // procMesh.facesInstance()
                        polyMesh::meshSubDir/"sets"
                    );

                    IOobjectList cSets(objects.lookupClass(cellSet::typeName));
                    forAllConstIter(IOobjectList, cSets, iter)
                    {
                        cSetNames.insert(iter.key(), cSetNames.size());
                    }

                    IOobjectList fSets(objects.lookupClass(faceSet::typeName));
                    forAllConstIter(IOobjectList, fSets, iter)
                    {
                        fSetNames.insert(iter.key(), fSetNames.size());
                    }
                    IOobjectList pSets(objects.lookupClass(pointSet::typeName));
                    forAllConstIter(IOobjectList, pSets, iter)
                    {
                        pSetNames.insert(iter.key(), pSetNames.size());
                    }
                }

                if (cSetNames.size() || fSetNames.size() || pSetNames.size())
                {
                    // Construct all sets
                    PtrList<cellSet> cellSets(cSetNames.size());
                    PtrList<faceSet> faceSets(fSetNames.size());
                    PtrList<pointSet> pointSets(pSetNames.size());

                    Info<< "Reconstructing sets:" << endl;
                    if (cSetNames.size())
                    {
                        Info<< "    cellSets "
                            << cSetNames.sortedToc() << endl;
                    }
                    if (fSetNames.size())
                    {
                        Info<< "    faceSets "
                            << fSetNames.sortedToc() << endl;
                    }
                    if (pSetNames.size())
                    {
                        Info<< "    pointSets "
                            << pSetNames.sortedToc() << endl;
                    }

                    // Load sets
                    forAll(procMeshes.meshes(), proci)
                    {
                        const fvMesh& procMesh = procMeshes.meshes()[proci];

                        IOobjectList objects
                        (
                            procMesh,
                            databases[0].timeName(),
                            polyMesh::meshSubDir/"sets"
                        );

                        // cellSets
                        const labelList& cellMap =
                            procMeshes.cellProcAddressing()[proci];

                        IOobjectList cSets
                        (
                            objects.lookupClass(cellSet::typeName)
                        );

                        forAllConstIter(IOobjectList, cSets, iter)
                        {
                            // Load cellSet
                            const cellSet procSet(*iter());
                            label setI = cSetNames[iter.key()];
                            if (!cellSets.set(setI))
                            {
                                cellSets.set
                                (
                                    setI,
                                    new cellSet
                                    (
                                        mesh,
                                        iter.key(),
                                        procSet.size()
                                    )
                                );
                            }
                            cellSet& cSet = cellSets[setI];
                            cSet.instance() = runTime.timeName();

                            forAllConstIter(cellSet, procSet, iter)
                            {
                                cSet.insert(cellMap[iter.key()]);
                            }
                        }

                        // faceSets
                        const labelList& faceMap =
                        procMeshes.faceProcAddressing()[proci];

                        IOobjectList fSets
                        (
                            objects.lookupClass(faceSet::typeName)
                        );

                        forAllConstIter(IOobjectList, fSets, iter)
                        {
                            // Load faceSet
                            const faceSet procSet(*iter());
                            label setI = fSetNames[iter.key()];
                            if (!faceSets.set(setI))
                            {
                                faceSets.set
                                (
                                    setI,
                                    new faceSet
                                    (
                                        mesh,
                                        iter.key(),
                                        procSet.size()
                                    )
                                );
                            }
                            faceSet& fSet = faceSets[setI];
                            fSet.instance() = runTime.timeName();

                            forAllConstIter(faceSet, procSet, iter)
                            {
                                fSet.insert(mag(faceMap[iter.key()])-1);
                            }
                        }
                        // pointSets
                        const labelList& pointMap =
                            procMeshes.pointProcAddressing()[proci];

                        IOobjectList pSets
                        (
                            objects.lookupClass(pointSet::typeName)
                        );
                        forAllConstIter(IOobjectList, pSets, iter)
                        {
                            // Load pointSet
                            const pointSet propSet(*iter());
                            label setI = pSetNames[iter.key()];
                            if (!pointSets.set(setI))
                            {
                                pointSets.set
                                (
                                    setI,
                                    new pointSet
                                    (
                                        mesh,
                                        iter.key(),
                                        propSet.size()
                                    )
                                );
                            }
                            pointSet& pSet = pointSets[setI];
                            pSet.instance() = runTime.timeName();

                            forAllConstIter(pointSet, propSet, iter)
                            {
                                pSet.insert(pointMap[iter.key()]);
                            }
                        }
                    }

                    // Write sets
                    forAll(cellSets, i)
                    {
                        cellSets[i].write();
                    }
                    forAll(faceSets, i)
                    {
                        faceSets[i].write();
                    }
                    forAll(pointSets, i)
                    {
                        pointSets[i].write();
                    }
                }
            }


            // Reconstruct refinement data
            {
                PtrList<hexRef8Data> procData(procMeshes.meshes().size());

                forAll(procMeshes.meshes(), procI)
                {
                    const fvMesh& procMesh = procMeshes.meshes()[procI];

                    procData.set
                    (
                        procI,
                        new hexRef8Data
                        (
                            IOobject
                            (
                                "dummy",
                                procMesh.time().timeName(),
                                polyMesh::meshSubDir,
                                procMesh,
                                IOobject::READ_IF_PRESENT,
                                IOobject::NO_WRITE,
                                false
                            )
                        )
                    );
                }

                // Combine individual parts

                const PtrList<labelIOList>& cellAddr =
                    procMeshes.cellProcAddressing();

                UPtrList<const labelList> cellMaps(cellAddr.size());
                forAll(cellAddr, i)
                {
                    cellMaps.set(i, &cellAddr[i]);
                }

                const PtrList<labelIOList>& pointAddr =
                    procMeshes.pointProcAddressing();

                UPtrList<const labelList> pointMaps(pointAddr.size());
                forAll(pointAddr, i)
                {
                    pointMaps.set(i, &pointAddr[i]);
                }

                UPtrList<const hexRef8Data> procRefs(procData.size());
                forAll(procData, i)
                {
                    procRefs.set(i, &procData[i]);
                }

                hexRef8Data
                (
                    IOobject
                    (
                        "dummy",
                        mesh.time().timeName(),
                        polyMesh::meshSubDir,
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    cellMaps,
                    pointMaps,
                    procRefs
                ).write();
            }

            // If there is a "uniform" directory in the time region
            // directory copy from the master processor
            {
                fileName uniformDir0
                (
                    fileHandler().filePath
                    (
                        databases[0].timePath()/regionDir/"uniform"
                    )
                );

                if (!uniformDir0.empty() && fileHandler().isDir(uniformDir0))
                {
                    fileHandler().cp(uniformDir0, runTime.timePath()/regionDir);
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
                        databases[0].timePath()/"uniform"
                    )
                );

                if (!uniformDir0.empty() && fileHandler().isDir(uniformDir0))
                {
                    fileHandler().cp(uniformDir0, runTime.timePath());
                }
            }
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
