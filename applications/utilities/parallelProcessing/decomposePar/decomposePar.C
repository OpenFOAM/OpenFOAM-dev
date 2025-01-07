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
    decomposePar

Description
    Automatically decomposes a mesh and fields of a case for parallel
    execution of OpenFOAM.

Usage
    \b decomposePar [OPTION]

    Options:
      - \par -cellProc
        Write cell processor indices as a volScalarField::Internal for
        post-processing.

      - \par -region \<regionName\> \n
        Decompose named region. Does not check for existence of processor*.

      - \par -allRegions \n
        Decompose all regions in regionSolvers. Does not check for
        existence of processor*.

      - \par -copyZero \n
        Copy \a 0 directory to processor* rather than decompose the fields.

      - \par -copyUniform \n
        Copy any \a uniform directories too.

      - \par -constant
        Decompose mesh and fields in the constant directory.

      - \par -time xxx:yyy \n
        Override controlDict settings and decompose selected times.

      - \par -fields \n
        Use existing geometry decomposition and convert fields only.

      - \par -noSets \n
        Skip decomposing cellSets, faceSets, pointSets.

      - \par -force \n
        Remove any existing \a processor subdirectories before decomposing the
        geometry.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "IOobjectList.H"
#include "processorRunTimes.H"
#include "multiDomainDecomposition.H"
#include "decompositionMethod.H"
#include "fvFieldDecomposer.H"
#include "pointFieldDecomposer.H"
#include "lagrangianFieldDecomposer.H"
#include "LagrangianFieldDecomposer.H"

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
            runTimes.completeTime().timePath()/regionDir/"uniform"
        );
}


void decomposeUniform
(
    const bool copyUniform,
    const processorRunTimes& runTimes,
    const word& regionDir = word::null
)
{
    const fileName uniformDir(regionDir/"uniform");

    forAll(runTimes.procTimes(), proci)
    {
        const fileName procTimePath =
            fileHandler().filePath(runTimes.procTimes()[proci].timePath());

        if (!fileHandler().isDir(procTimePath))
        {
            fileHandler().mkDir(procTimePath);
        }

        if (copyUniform)
        {
            if (!fileHandler().exists(procTimePath/uniformDir))
            {
                fileHandler().cp
                (
                    runTimes.completeTime().timePath()/uniformDir,
                    procTimePath/uniformDir
                );
            }
        }
        else
        {
            // Link with relative paths
            string parentPath = string("..")/"..";

            if (regionDir != word::null)
            {
                parentPath = parentPath/"..";
            }

            fileName currentDir(cwd());

            chDir(procTimePath);

            if (!fileHandler().exists(uniformDir))
            {
                fileHandler().ln
                (
                    parentPath/runTimes.completeTime().name()/uniformDir,
                    uniformDir
                );
            }

            chDir(currentDir);
        }
    }
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
        "decompose a mesh and fields of a case for parallel execution"
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
    argList::addBoolOption
    (
        "copyZero",
        "Copy \a 0 directory to processor* rather than decompose the fields"
    );
    argList::addBoolOption
    (
        "copyUniform",
        "copy any uniform/ directories too"
    );
    argList::addBoolOption
    (
        "fields",
        "use existing geometry decomposition and convert fields only"
    );
    argList::addBoolOption
    (
        "noFields",
        "opposite of -fields; only decompose geometry"
    );
    argList::addBoolOption
    (
        "noSets",
        "skip decomposing cellSets, faceSets, pointSets"
    );
    argList::addBoolOption
    (
        "force",
        "remove existing processor*/ subdirs before decomposing the geometry"
    );

    // Include explicit constant option, execute from zero by default
    timeSelector::addOptions(true, false);

    #include "setRootCase.H"
    #include "setMeshPath.H"

    const bool region              = args.optionFound("region");
    const bool writeCellProc       = args.optionFound("cellProc");
    const bool copyZero            = args.optionFound("copyZero");
    const bool copyUniform         = args.optionFound("copyUniform");
    const bool decomposeFieldsOnly = args.optionFound("fields");
    const bool decomposeGeomOnly   = args.optionFound("noFields");
    const bool decomposeSets       = !args.optionFound("noSets");
    const bool forceOverwrite      = args.optionFound("force");

    if (decomposeGeomOnly)
    {
        Info<< "Skipping decomposing fields" << nl << endl;

        if (decomposeFieldsOnly || copyZero)
        {
            FatalErrorInFunction
                << "Cannot combine geometry-only decomposition (-noFields)"
                << " with field decomposition (-fields or -copyZero)"
                << exit(FatalError);
        }
    }

    // Set time from database
    Info<< "Create time" << nl << endl;
    processorRunTimes runTimes(Foam::Time::controlDictName, args);
    const Time& runTime = runTimes.completeTime();

    // Allow override of time
    const instantList times = runTimes.selectComplete(args);

    #include "setRegionNames.H"

    // Remove existing processor directories if requested
    if (forceOverwrite)
    {
        if (region)
        {
            FatalErrorInFunction
                << "Cannot force the decomposition of a single region"
                << exit(FatalError);
        }

        const label nProcs0 =
            fileHandler().nProcs(runTimes.completeTime().path());

        Info<< "Removing " << nProcs0
            << " existing processor directories" << nl << endl;

        // Remove existing processor directories
        const fileNameList dirs
        (
            fileHandler().readDir
            (
                runTimes.completeTime().path(),
                fileType::directory
            )
        );
        forAllReverse(dirs, diri)
        {
            const fileName& d = dirs[diri];

            // Starts with 'processors'
            if (d.find("processors") == 0)
            {
                if (fileHandler().exists(d))
                {
                    fileHandler().rmDir(d);
                }
            }

            // Starts with 'processor'
            if (d.find("processor") == 0)
            {
                // Check that integer after processor
                fileName num(d.substr(9));
                label proci = -1;
                if (Foam::read(num.c_str(), proci))
                {
                    if (fileHandler().exists(d))
                    {
                        fileHandler().rmDir(d);
                    }
                }
            }
        }

        // Flush file handler to clear any detected processor directories
        fileHandler().flush();
    }

    // Check the specified number of processes is consistent with any existing
    // processor directories
    {
        const label nProcs0 =
            fileHandler().nProcs(runTimes.completeTime().path());

        if (nProcs0 && nProcs0 != runTimes.nProcs())
        {
            FatalErrorInFunction
                << "Case is already decomposed with " << nProcs0
                << " domains, use the -force option or manually" << nl
                << "remove processor directories before decomposing. e.g.,"
                << nl
                << "    rm -rf " << runTimes.completeTime().path().c_str()
                << "/processor*"
                << nl
                << exit(FatalError);
        }
    }

    // Get the decomposition dictionary
    const dictionary decomposeParDict =
        decompositionMethod::decomposeParDict(runTimes.completeTime());

    // Check existing decomposition
    forAll(regionNames, regioni)
    {
        const word& regionName = regionNames[regioni];
        const word regionDir =
            regionName == polyMesh::defaultRegion ? word::null : regionName;

        // Determine the existing processor count directly
        const label nProcs =
            fileHandler().nProcs(runTimes.completeTime().path(), regionDir);

        // Get requested numberOfSubdomains
        const label nDomains =
            decomposeParDict.lookup<label>("numberOfSubdomains");

        // Give file handler a chance to determine the output directory
        const_cast<fileOperation&>(fileHandler()).setNProcs(nDomains);

        // Sanity check number of processors in a previously decomposed case
        if (decomposeFieldsOnly && nProcs != nDomains)
        {
            FatalErrorInFunction
                << "Specified -fields, but the case was decomposed with "
                << nProcs << " domains" << nl << "instead of " << nDomains
                << " domains as specified in decomposeParDict" << nl
                << exit(FatalError);
        }
    }

    // Create meshes
    multiDomainDecomposition regionMeshes(runTimes, meshPath, regionNames);
    if
    (
        !(decomposeFieldsOnly && copyZero)
     && regionMeshes.readDecompose(decomposeSets)
    )
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

    // Get flag to determine whether or not to distribute uniform data
    const bool distributed =
        decomposeParDict.lookupOrDefault<bool>("distributed", false);

    // Loop over all times
    forAll(times, timei)
    {
        // Set the time
        runTimes.setTime(times[timei], timei);

        Info<< "Time = " << runTimes.completeTime().userTimeName()
            << nl << endl;

        // Update the meshes, if necessary
        const fvMesh::readUpdateState stat =
            !(decomposeFieldsOnly && copyZero)
          ? regionMeshes.readUpdateDecompose()
          : fvMesh::UNCHANGED;
        if (stat >= fvMesh::TOPO_CHANGE) Info<< endl;

        // Write the mesh out (if anything has changed), if necessary
        if (!decomposeFieldsOnly)
        {
            regionMeshes.writeProcs(decomposeSets);
        }

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

        // If only decomposing geometry then there is no more to do
        if (decomposeGeomOnly)
        {
            continue;
        }

        // If copying from zero then just copy everything from the <time>
        // directory to the processor*/<time> directories without altering them
        if (copyZero)
        {
            const fileName completeTimePath =
                runTimes.completeTime().timePath();

            fileName prevProcTimePath;
            for (label proci = 0; proci < runTimes.nProcs(); proci++)
            {
                const Time& procRunTime = runTimes.procTimes()[proci];

                if (fileHandler().isDir(completeTimePath))
                {
                    const fileName procTimePath
                    (
                        fileHandler().objectPath
                        (
                            IOobject(word::null, fileName::null, procRunTime)
                        )
                    );

                    if (procTimePath != prevProcTimePath)
                    {
                        Info<< "Processor " << proci
                            << ": copying " << completeTimePath << nl
                            << " to " << procTimePath << endl;
                        fileHandler().cp(completeTimePath, procTimePath);
                        prevProcTimePath = procTimePath;
                    }
                }
            }

            Info<< endl;

            continue;
        }

        // Otherwise we are doing full region-by-region decomposition of all
        // the available fields
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
                    meshes().completeMesh(),
                    runTimes.completeTime().name()
                );

                {
                    Info<< dnl << "Decomposing FV fields" << endl;

                    if (fvFieldDecomposer::decomposes(objects))
                    {
                        fvFieldDecomposer fvDecomposer
                        (
                            meshes().completeMesh(),
                            meshes().procMeshes(),
                            meshes().procFaceAddressing(),
                            meshes().procCellAddressing(),
                            meshes().procFaceAddressingBf()
                        );

                        #define DO_FV_VOL_INTERNAL_FIELDS_TYPE(Type, nullArg)  \
                            fvDecomposer.decomposeVolInternalFields<Type>      \
                            (objects);
                        FOR_ALL_FIELD_TYPES(DO_FV_VOL_INTERNAL_FIELDS_TYPE)
                        #undef DO_FV_VOL_INTERNAL_FIELDS_TYPE

                        #define DO_FV_VOL_FIELDS_TYPE(Type, nullArg)           \
                            fvDecomposer.decomposeVolFields<Type>              \
                            (objects);
                        FOR_ALL_FIELD_TYPES(DO_FV_VOL_FIELDS_TYPE)
                        #undef DO_FV_VOL_FIELDS_TYPE

                        #define DO_FV_SURFACE_FIELDS_TYPE(Type, nullArg)       \
                            fvDecomposer.decomposeFvSurfaceFields<Type>        \
                            (objects);
                        FOR_ALL_FIELD_TYPES(DO_FV_SURFACE_FIELDS_TYPE)
                        #undef DO_FV_SURFACE_FIELDS_TYPE
                    }
                    else
                    {
                        Info<< dnl << "    (no FV fields)" << endl;
                    }
                }

                {
                    Info<< dnl << "Decomposing point fields" << endl;

                    if (pointFieldDecomposer::decomposes(objects))
                    {
                        pointFieldDecomposer pointDecomposer
                        (
                            pointMesh::New(meshes().completeMesh()),
                            meshes().procMeshes(),
                            meshes().procPointAddressing()
                        );

                        #define DO_POINT_FIELDS_TYPE(Type, nullArg)            \
                            pointDecomposer.decomposeFields<Type>              \
                            (objects);
                        FOR_ALL_FIELD_TYPES(DO_POINT_FIELDS_TYPE)
                        #undef DO_POINT_FIELDS_TYPE
                    }
                    else
                    {
                        Info<< dnl << "    (no point fields)" << endl;
                    }
                }

                {
                    // Find cloud directories
                    fileNameList cloudDirs
                    (
                        fileHandler().readDir
                        (
                            runTimes.completeTime().timePath()
                           /regionDir
                           /lagrangian::cloud::prefix,
                            fileType::directory
                        )
                    );

                    // Add objects in any found cloud directories
                    HashTable<IOobjectList> cloudsObjects;
                    forAll(cloudDirs, i)
                    {
                        // Do local scan for valid cloud objects
                        IOobjectList cloudObjs
                        (
                            meshes().completeMesh(),
                            runTimes.completeTime().name(),
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

                    // Decompose the objects found above
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

                            Info<< dnl << "Decomposing lagrangian fields for "
                                << "cloud " << cloudName << endl;

                            if
                            (
                                lagrangianFieldDecomposer::decomposes
                                (
                                    cloudObjects
                                )
                            )
                            {
                                lagrangianFieldDecomposer
                                    lagrangianDecomposer
                                    (
                                        meshes().completeMesh(),
                                        meshes().procMeshes(),
                                        meshes().procFaceAddressing(),
                                        meshes().procCellAddressing(),
                                        cloudName
                                    );

                                #define DO_CLOUD_FIELDS_TYPE(Type, nullArg)    \
                                    lagrangianDecomposer.decomposeFields<Type> \
                                    (cloudObjects);
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

                {
                    // Find Lagrangian directories
                    fileNameList LagrangianDirs
                    (
                        fileHandler().readDir
                        (
                            runTimes.completeTime().timePath()
                           /regionDir
                           /LagrangianMesh::prefix,
                            fileType::directory
                        )
                    );

                    // Add objects in any found Lagrangian directories
                    HashTable<IOobjectList> LagrangianObjects;
                    forAll(LagrangianDirs, i)
                    {
                        // Do local scan for valid Lagrangian objects
                        IOobjectList objects
                        (
                            meshes().completeMesh(),
                            runTimes.completeTime().name(),
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
                         || LagrangianFieldDecomposer::decomposes
                            (
                                objects
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

                    // Decompose the objects found above
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

                            Info<< dnl << "Decomposing Lagrangian fields "
                                << "for " << LagrangianName << endl;

                            const LagrangianFieldDecomposer
                                LagrangianDecomposer
                                (
                                    meshes().completeMesh(),
                                    meshes().procMeshes(),
                                    meshes().procFaceAddressing(),
                                    meshes().procCellAddressing(),
                                    LagrangianName
                                );

                            if
                            (
                                LagrangianFieldDecomposer::decomposes
                                (
                                    iter()
                                )
                            )
                            {
                                #define DO_LAGRANGIAN_FIELDS_TYPE(             \
                                    Type, GeoField)                            \
                                    LagrangianDecomposer                       \
                                   .decomposeFields<GeoField<Type>>            \
                                    (iter());
                                DO_LAGRANGIAN_FIELDS_TYPE
                                (
                                    label,
                                    LagrangianField
                                )
                                FOR_ALL_FIELD_TYPES
                                (
                                    DO_LAGRANGIAN_FIELDS_TYPE,
                                    LagrangianField
                                )
                                DO_LAGRANGIAN_FIELDS_TYPE
                                (
                                    label,
                                    LagrangianInternalField
                                )
                                FOR_ALL_FIELD_TYPES
                                (
                                    DO_LAGRANGIAN_FIELDS_TYPE,
                                    LagrangianInternalField
                                )
                                #undef DO_LAGRANGIAN_FIELDS_TYPE

                                // --> Note we don't have to explicitly
                                // decompose the dynamic variants of these
                                // fields as they are IO compatible with the
                                // non-dynamic fields
                            }
                            else
                            {
                                Info<< dnl << "    (no lagrangian fields)"
                                    << endl;
                            }
                        }
                    }
                }
            }

            Info<< dnl;
        }

        // Distribute the uniform directory
        if (haveUniform(runTimes))
        {
            Info<< "Distributing uniform files" << endl;

            decomposeUniform
            (
                copyUniform || distributed,
                runTimes
            );

            Info<< endl;
        }

        if (regionNames == wordList(1, polyMesh::defaultRegion)) continue;

        // Distribute the region uniform directories
        forAll(regionNames, regioni)
        {
            const word& regionName = regionNames[regioni];
            const word regionDir =
                regionName == polyMesh::defaultRegion ? word::null : regionName;

            if (haveUniform(runTimes, regionDir))
            {
                // Prefixed scope
                {
                    const RegionRef<domainDecomposition> meshes =
                        regionMeshes[regioni];

                    Info<< "Distributing uniform files" << endl;

                    decomposeUniform
                    (
                        copyUniform || distributed,
                        runTimes,
                        regionDir
                    );
                }

                Info<< endl;
            }
        }
    }

    Info<< "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
