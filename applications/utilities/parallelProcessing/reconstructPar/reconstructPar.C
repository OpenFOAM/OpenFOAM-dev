/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
#include "multiDomainDecomposition.H"
#include "fvFieldReconstructor.H"
#include "pointFieldReconstructor.H"
#include "lagrangianFieldReconstructor.H"
#include "LagrangianFieldReconstructor.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

bool haveUniform
(
    const processorRunTimes& runTimes,
    const word& regionDir = word::null
)
{
    return
        fileHandler().isDir
        (
            fileHandler().filePath
            (
                runTimes.procTimes()[0].timePath()/regionDir/"uniform"
            )
        );
}


void reconstructUniform
(
    const processorRunTimes& runTimes,
    const word& regionDir = word::null
)
{
    fileHandler().cp
    (
        fileHandler().filePath
        (
            runTimes.procTimes()[0].timePath()/regionDir/"uniform"
        ),
        runTimes.completeTime().timePath()/regionDir
    );
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
        << cellProc.name() << " for use in postprocessing"
        << endl;
}

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class delayedNewLine
{
    mutable bool first_;

public:

    delayedNewLine()
    :
        first_(true)
    {}

    friend Ostream& operator<<(Ostream& os, const delayedNewLine& dnl)
    {
        if (!dnl.first_) os << nl;
        dnl.first_ = false;
        return os;
    }
};

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Reconstruct fields of a parallel case"
    );

    argList::noParallel();
    #include "addMeshOption.H"
    #include "addRegionOption.H"
    #include "addAllRegionsOption.H"
    argList::addBoolOption
    (
        "cellProc",
        "write cell processor indices as a volScalarField::Internal for "
        "post-processing"
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
        "positions always included"
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
    argList::addBoolOption
    (
        "rm",
        "remove processor time directories after reconstruction"
    );

    // Include explicit constant options, and explicit zero option (to prevent
    // the user accidentally trashing the initial fields)
    timeSelector::addOptions(true, true);

    #include "setRootCase.H"
    #include "setMeshPath.H"

    const bool writeCellProc = args.optionFound("cellProc");

    HashSet<word> selectedFields;
    if (args.optionFound("fields"))
    {
        args.optionLookup("fields")() >> selectedFields;
    }

    const bool noFields = args.optionFound("noFields");

    if (noFields)
    {
        Info<< "Skipping reconstructing fields" << nl << endl;
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
                << "options together"
                << exit(FatalError);
        }

        args.optionLookup("lagrangianFields")() >> selectedLagrangianFields;
    }

    // Set time from database
    Info<< "Create time" << nl << endl;
    processorRunTimes runTimes(Foam::Time::controlDictName, args);

    // Get the times to reconstruct
    instantList times = runTimes.selectProc(args);

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

    // Quit if no times
    if (times.empty())
    {
        WarningInFunction << "No times selected" << nl << endl;
        exit(1);
    }

    // If only reconstructing new times then filter out existing times
    if (args.optionFound("newTimes"))
    {
        // Get all existing times
        const instantList existingTimes = runTimes.completeTime().times();

        // Put into a set
        HashSet<word> existingTimesSet;
        existingTimesSet.resize(2*existingTimes.size());
        forAll(existingTimes, i)
        {
            existingTimesSet.insert(existingTimes[i].name());
        }

        // Remove times from the existing time set by shuffling up
        label timei = 0;
        forAll(times, timej)
        {
            if (!existingTimesSet.found(times[timej].name()))
            {
                times[timei ++] = times[timej];
            }
        }
        times.resize(timei);
    }

    // Quit if no times
    if (times.empty())
    {
        Info<< "All times already reconstructed" << nl << nl
            << "End" << nl << endl;
        return 0;
    }

    // Create meshes
    multiDomainDecomposition regionMeshes(runTimes, meshPath, regionNames);
    if (regionMeshes.readReconstruct(!noReconstructSets))
    {
        Info<< endl;

        if (writeCellProc)
        {
            forAll(regionNames, regioni)
            {
                writeDecomposition(regionMeshes[regioni]());
                Info<< endl;
                fileHandler().flush();
            }
        }
    }

    // Loop over all times
    forAll(times, timei)
    {
        // Set the time
        runTimes.setTime(times[timei], timei);

        Info<< "Time = " << runTimes.completeTime().userTimeName()
            << nl << endl;

        // Update the meshes
        const fvMesh::readUpdateState stat =
            regionMeshes.readUpdateReconstruct();
        if (stat >= fvMesh::TOPO_CHANGE) Info<< endl;

        // Write the mesh out (if anything has changed)
        regionMeshes.writeComplete(!noReconstructSets);

        // Write the decomposition, if necessary
        forAll(regionNames, regioni)
        {
            if (writeCellProc && stat >= fvMesh::TOPO_CHANGE)
            {
                writeDecomposition(regionMeshes[regioni]());
                Info<< endl;
                fileHandler().flush();
            }
        }

        // Do a region-by-region reconstruction of all the available fields
        forAll(regionNames, regioni)
        {
            const word& regionName = regionNames[regioni];
            const word regionDir =
                regionName == polyMesh::defaultRegion ? word::null : regionName;

            const delayedNewLine dnl;

            // Prefixed scope
            {
                const RegionRef<domainDecomposition> meshes =
                    regionMeshes[regioni];

                // Search for objects at this time
                IOobjectList objects
                (
                    meshes().procMeshes()[0],
                    runTimes.procTimes()[0].name()
                );

                if (!noFields)
                {
                    Info<< dnl << "Reconstructing FV fields" << endl;

                    if
                    (
                        fvFieldReconstructor::reconstructs
                        (
                            objects,
                            selectedFields
                        )
                    )
                    {
                        fvFieldReconstructor fvReconstructor
                        (
                            meshes().completeMesh(),
                            meshes().procMeshes(),
                            meshes().procFaceAddressing(),
                            meshes().procCellAddressing(),
                            meshes().procFaceAddressingBf()
                        );

                        #define DO_FV_VOL_INTERNAL_FIELDS_TYPE(Type, nullArg)  \
                            fvReconstructor.reconstructVolInternalFields<Type> \
                            (objects, selectedFields);
                        FOR_ALL_FIELD_TYPES(DO_FV_VOL_INTERNAL_FIELDS_TYPE)
                        #undef DO_FV_VOL_INTERNAL_FIELDS_TYPE

                        #define DO_FV_VOL_FIELDS_TYPE(Type, nullArg)           \
                            fvReconstructor.reconstructVolFields<Type>         \
                            (objects, selectedFields);
                        FOR_ALL_FIELD_TYPES(DO_FV_VOL_FIELDS_TYPE)
                        #undef DO_FV_VOL_FIELDS_TYPE

                        #define DO_FV_SURFACE_FIELDS_TYPE(Type, nullArg)       \
                            fvReconstructor.reconstructFvSurfaceFields<Type>   \
                            (objects, selectedFields);
                        FOR_ALL_FIELD_TYPES(DO_FV_SURFACE_FIELDS_TYPE)
                        #undef DO_FV_SURFACE_FIELDS_TYPE
                    }
                    else
                    {
                        Info<< dnl << "    (no FV fields)" << endl;
                    }
                }

                if (!noFields)
                {
                    Info<< dnl << "Reconstructing point fields" << endl;

                    if
                    (
                        pointFieldReconstructor::reconstructs
                        (
                            objects,
                            selectedFields
                        )
                    )
                    {
                        pointFieldReconstructor pointReconstructor
                        (
                            pointMesh::New(meshes().completeMesh()),
                            meshes().procMeshes(),
                            meshes().procPointAddressing()
                        );

                        #define DO_POINT_FIELDS_TYPE(Type, nullArg)            \
                            pointReconstructor.reconstructFields<Type>         \
                            (objects, selectedFields);
                        FOR_ALL_FIELD_TYPES(DO_POINT_FIELDS_TYPE)
                        #undef DO_POINT_FIELDS_TYPE
                    }
                    else
                    {
                        Info<< dnl << "    (no point fields)" << endl;
                    }
                }

                if (!noLagrangian)
                {
                    // Search for clouds that exist on any processor and add
                    // them into this table of cloud objects
                    HashTable<IOobjectList> cloudsObjects;
                    forAll(runTimes.procTimes(), proci)
                    {
                        // Find cloud directories
                        fileNameList cloudDirs
                        (
                            fileHandler().readDir
                            (
                                fileHandler().filePath
                                (
                                    runTimes.procTimes()[proci].timePath()
                                   /regionDir
                                   /lagrangian::cloud::prefix
                                ),
                                fileType::directory
                            )
                        );

                        // Add objects in any found cloud directories
                        forAll(cloudDirs, i)
                        {
                            // Pass if we already have an objects for this name
                            HashTable<IOobjectList>::const_iterator iter =
                                cloudsObjects.find(cloudDirs[i]);
                            if (iter != cloudsObjects.end()) continue;

                            // Do local scan for valid cloud objects
                            IOobjectList cloudObjs
                            (
                                meshes().procMeshes()[proci],
                                runTimes.procTimes()[proci].name(),
                                lagrangian::cloud::prefix/cloudDirs[i],
                                IOobject::MUST_READ,
                                IOobject::NO_WRITE,
                                false
                            );

                            // If "positions" is present, then add to the table
                            if (cloudObjs.lookup(word("positions")))
                            {
                                cloudsObjects.insert(cloudDirs[i], cloudObjs);
                            }
                        }
                    }

                    // Reconstruct the objects found above
                    if (cloudsObjects.size())
                    {
                        forAllConstIter
                        (
                            HashTable<IOobjectList>,
                            cloudsObjects,
                            iter
                        )
                        {
                            const word cloudName =
                                string::validate<word>(iter.key());

                            const IOobjectList& cloudObjects = iter();

                            Info<< dnl << "Reconstructing lagrangian fields "
                                << "for cloud " << cloudName << endl;

                            if
                            (
                                lagrangianFieldReconstructor::reconstructs
                                (
                                    cloudObjects,
                                    selectedLagrangianFields
                                )
                            )
                            {
                                lagrangianFieldReconstructor
                                    lagrangianReconstructor
                                    (
                                        meshes().completeMesh(),
                                        meshes().procMeshes(),
                                        meshes().procFaceAddressing(),
                                        meshes().procCellAddressing(),
                                        cloudName
                                    );

                                #define DO_CLOUD_FIELDS_TYPE(Type, nullArg)    \
                                    lagrangianReconstructor                    \
                                   .reconstructFields<Type>                    \
                                    (cloudObjects, selectedLagrangianFields);
                                DO_CLOUD_FIELDS_TYPE(label, )
                                FOR_ALL_FIELD_TYPES(DO_CLOUD_FIELDS_TYPE)
                                #undef DO_CLOUD_FIELDS_TYPE
                            }
                            else
                            {
                                Info<< dnl << "    (no lagrangian fields)"
                                    << endl;
                            }
                        }
                    }
                }

                if (!noLagrangian)
                {
                    // Search for Lagrangian meshes that exist on any processor
                    // and add them into this table of objects
                    HashTable<IOobjectList> LagrangianObjects;
                    forAll(runTimes.procTimes(), proci)
                    {
                        // Find Lagrangian directories
                        fileNameList LagrangianDirs
                        (
                            fileHandler().readDir
                            (
                                fileHandler().filePath
                                (
                                    runTimes.procTimes()[proci].timePath()
                                   /regionDir
                                   /LagrangianMesh::prefix
                                ),
                                fileType::directory
                            )
                        );

                        // Add objects in any found Lagrangian directories
                        forAll(LagrangianDirs, i)
                        {
                            // Pass if we already have an objects for this name
                            if
                            (
                                LagrangianObjects.find(LagrangianDirs[i])
                             != LagrangianObjects.end()
                            ) continue;

                            // Do local scan for valid Lagrangian objects
                            IOobjectList objects
                            (
                                meshes().procMeshes()[proci],
                                runTimes.procTimes()[proci].name(),
                                LagrangianMesh::prefix/LagrangianDirs[i],
                                IOobject::MUST_READ,
                                IOobject::NO_WRITE,
                                false
                            );

                            // If coordinates or fields are present then add
                            // this set of objects to the table
                            if
                            (
                                objects.found(LagrangianMesh::coordinatesName)
                             || LagrangianFieldReconstructor::reconstructs
                                (
                                    objects,
                                    selectedLagrangianFields
                                )
                            )
                            {
                                LagrangianObjects.insert
                                (
                                    LagrangianDirs[i],
                                    objects
                                );
                            }
                        }
                    }

                    // Reconstruct the objects found above
                    if (LagrangianObjects.size())
                    {
                        forAllConstIter
                        (
                            HashTable<IOobjectList>,
                            LagrangianObjects,
                            iter
                        )
                        {
                            const word LagrangianName =
                                string::validate<word>(iter.key());

                            Info<< dnl << "Reconstructing Lagrangian fields "
                                << "for " << LagrangianName << endl;

                            const LagrangianFieldReconstructor
                                LagrangianReconstructor
                                (
                                    meshes().completeMesh(),
                                    meshes().procMeshes(),
                                    meshes().procFaceAddressing(),
                                    meshes().procCellAddressing(),
                                    LagrangianName
                                );

                            if
                            (
                                LagrangianFieldReconstructor::reconstructs
                                (
                                    iter(),
                                    selectedLagrangianFields
                                )
                            )
                            {
                                #define DO_LAGRANGIAN_FIELDS_TYPE(             \
                                    Type, GeoField)                            \
                                    LagrangianReconstructor                    \
                                   .reconstructFields<GeoField<Type>>          \
                                    (iter(), selectedLagrangianFields);
                                DO_LAGRANGIAN_FIELDS_TYPE
                                (
                                    label,
                                    LagrangianField
                                )
                                FOR_ALL_FIELD_TYPES
                                (
                                    DO_LAGRANGIAN_FIELDS_TYPE,
                                    LagrangianField
                                );
                                DO_LAGRANGIAN_FIELDS_TYPE
                                (
                                    label,
                                    LagrangianInternalField
                                )
                                FOR_ALL_FIELD_TYPES
                                (
                                    DO_LAGRANGIAN_FIELDS_TYPE,
                                    LagrangianInternalField
                                );
                                #undef DO_LAGRANGIAN_FIELDS_TYPE

                                // --> Note we don't have to explicitly
                                // reconstruct the dynamic variants of these
                                // fields as they are IO compatible with the
                                // non-dynamic fields
                            }
                            else
                            {
                                Info<< dnl << "    (no Lagrangian fields)"
                                    << endl;
                            }
                        }
                    }
                }
            }

            Info<< dnl;
        }

        // Collect the uniform directory
        if (haveUniform(runTimes))
        {
            Info<< "Collecting uniform files" << endl;

            reconstructUniform(runTimes);

            Info<< endl;
        }

        if (regionNames != wordList(1, polyMesh::defaultRegion))
        {
            // Collect the region uniform directories
            forAll(regionNames, regioni)
            {
                const word& regionName = regionNames[regioni];
                const word regionDir =
                    regionName == polyMesh::defaultRegion
                  ? word::null
                  : regionName;

                if (haveUniform(runTimes, regionDir))
                {
                    // Prefixed scope
                    {
                        const RegionRef<domainDecomposition> meshes =
                            regionMeshes[regioni];

                        Info<< "Collecting uniform files" << endl;

                        reconstructUniform(runTimes, regionDir);
                    }

                    Info<< endl;
                }
            }
        }

        if (args.optionFound("rm") && times[timei].name() != Time::constantName)
        {
            const bool allRegions = args.optionFound("allRegions");

            forAll(regionNames, regioni)
            {
                const word& regionName = regionNames[regioni];
                const word regionDir =
                    allRegions || regionName == polyMesh::defaultRegion
                  ? word::null
                  : regionName;

                const RegionRef<domainDecomposition> meshes =
                    regionMeshes[regioni];

                Info<< "Removing processors time directory" << endl;

                for (label proci=0; proci<nProcs; proci++)
                {
                    const fileName procTimePath = fileHandler().filePath
                    (
                        runTimes.procTimes()[proci].timePath()/regionDir
                    );

                    if (isDir(procTimePath))
                    {
                        rmDir(procTimePath);
                    }
                }
            }

            Info<< endl;
        }
    }

    Info<< "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
