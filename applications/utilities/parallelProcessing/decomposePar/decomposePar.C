/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
        Decompose all regions in regionProperties. Does not check for
        existence of processor*.

      - \par -copyZero \n
        Copy \a 0 directory to processor* rather than decompose the fields.

      - \par -copyUniform \n
        Copy any \a uniform directories too.

      - \par -constant

      - \par -time xxx:yyy \n
        Override controlDict settings and decompose selected times. Does not
        re-decompose the mesh i.e. does not handle moving mesh or changing
        mesh cases.

      - \par -fields \n
        Use existing geometry decomposition and convert fields only.

      - \par -noSets \n
        Skip decomposing cellSets, faceSets, pointSets.

      - \par -force \n
        Remove any existing \a processor subdirectories before decomposing the
        geometry.

      - \par -dict \<filename\>
        Specify alternative dictionary for the decomposition.

\*---------------------------------------------------------------------------*/

#include "processorRunTimes.H"
#include "domainDecomposition.H"
#include "decompositionMethod.H"
#include "argList.H"
#include "timeSelector.H"
#include "regionProperties.H"

#include "labelIOField.H"
#include "labelFieldIOField.H"
#include "scalarIOField.H"
#include "scalarFieldIOField.H"
#include "vectorIOField.H"
#include "vectorFieldIOField.H"
#include "sphericalTensorIOField.H"
#include "sphericalTensorFieldIOField.H"
#include "symmTensorIOField.H"
#include "symmTensorFieldIOField.H"
#include "tensorIOField.H"
#include "tensorFieldIOField.H"

#include "readFields.H"
#include "dimFieldDecomposer.H"
#include "fvFieldDecomposer.H"
#include "pointFieldDecomposer.H"
#include "lagrangianFieldDecomposer.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void decomposeUniform
(
    const bool copyUniform,
    const bool distributeUniform,
    const Time& runTime,
    const Time& procRunTime,
    const word& regionDir = word::null
)
{
    const fileName uniformDir(regionDir/"uniform");

    if (fileHandler().isDir(runTime.timePath()/uniformDir))
    {
        Info<< "Detected additional non-decomposed files in "
            << runTime.timePath()/uniformDir
            << endl;

        const fileName timePath =
            fileHandler().filePath(procRunTime.timePath());

        if (copyUniform || distributeUniform)
        {
            if (!fileHandler().exists(timePath/uniformDir))
            {
                fileHandler().cp
                (
                    runTime.timePath()/uniformDir,
                    timePath/uniformDir
                );
            }
        }
        else
        {
            // link with relative paths
            string parentPath = string("..")/"..";

            if (regionDir != word::null)
            {
                parentPath = parentPath/"..";
            }

            fileName currentDir(cwd());
            chDir(timePath);

            if (!fileHandler().exists(uniformDir))
            {
                fileHandler().ln
                (
                    parentPath/runTime.timeName()/uniformDir,
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
            meshes.completeMesh().time().timeName(),
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
        "decompose a mesh and fields of a case for parallel execution"
    );

    argList::noParallel();
    #include "addDictOption.H"
    #include "addRegionOption.H"
    #include "addAllRegionsOption.H"
    argList::addBoolOption
    (
        "cellProc",
        "write cell processor indices as a volScalarField::Internal for "
        "post-processing."
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

    // Include explicit constant options, have zero from time range
    timeSelector::addOptions(true, false);

    #include "setRootCase.H"

    bool region                  = args.optionFound("region");
    bool writeCellProc           = args.optionFound("cellProc");
    bool copyZero                = args.optionFound("copyZero");
    bool copyUniform             = args.optionFound("copyUniform");
    bool decomposeFieldsOnly     = args.optionFound("fields");
    bool decomposeGeomOnly       = args.optionFound("noFields");
    bool decomposeSets           = !args.optionFound("noSets");
    bool forceOverwrite          = args.optionFound("force");

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
    Info<< "Create time\n" << endl;
    processorRunTimes runTimes(Foam::Time::controlDictName, args);

    // Allow override of time
    const instantList times = runTimes.selectComplete(args);

    // Get region names
    const wordList regionNames =
        selectRegionNames(args, runTimes.completeTime());

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
            << " existing processor directories" << endl;

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

    // Decompose all regions
    forAll(regionNames, regioni)
    {
        const word& regionName = regionNames[regioni];
        const word& regionDir = Foam::regionDir(regionName);

        Info<< "\n\nDecomposing mesh " << regionName << nl << endl;

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

        // Get flag to determine whether or not to distribute uniform data
        const label distributeUniform =
            decomposeParDict.lookupOrDefault<bool>("distributed", false);

        // Create meshes
        Info<< "Create mesh" << endl;
        domainDecomposition meshes(runTimes, regionName);
        if (!decomposeFieldsOnly || !copyZero)
        {
            if (meshes.readDecompose(decomposeSets) && writeCellProc)
            {
                writeDecomposition(meshes);
                fileHandler().flush();
            }
        }

        // Field maps. These are preserved if decomposing multiple times.
        PtrList<fvFieldDecomposer> fieldDecomposerList
        (
            meshes.nProcs()
        );
        PtrList<dimFieldDecomposer> dimFieldDecomposerList
        (
            meshes.nProcs()
        );
        PtrList<pointFieldDecomposer> pointFieldDecomposerList
        (
            meshes.nProcs()
        );

        // Loop over all times
        forAll(times, timei)
        {
            // Set the time
            runTimes.setTime(times[timei], timei);

            Info<< "Time = " << runTimes.completeTime().userTimeName() << endl;

            // Update the meshes, if necessary
            fvMesh::readUpdateState state = fvMesh::UNCHANGED;
            if (!decomposeFieldsOnly || !copyZero)
            {
                state = meshes.readUpdateDecompose();
            }

            // Write the mesh out, if necessary
            if (decomposeFieldsOnly)
            {
                // Nothing to do
            }
            else if (state != fvMesh::UNCHANGED)
            {
                meshes.writeProcs(decomposeSets);
            }

            // Write the decomposition, if necessary
            if
            (
                writeCellProc
             && meshes.completeMesh().facesInstance()
             == runTimes.completeTime().timeName()
            )
            {
                writeDecomposition(meshes);
                fileHandler().flush();
            }

            // Clear the field maps if there has been topology change
            if
            (
                state == fvMesh::TOPO_CHANGE
             || state == fvMesh::TOPO_PATCH_CHANGE
            )
            {
                for (label proci = 0; proci < meshes.nProcs(); proci++)
                {
                    fieldDecomposerList.set(proci, nullptr);
                    dimFieldDecomposerList.set(proci, nullptr);
                    pointFieldDecomposerList.set(proci, nullptr);
                }
            }

            // Decompose the fields at this time
            if (decomposeGeomOnly)
            {
                // Do nothing
            }
            else if (copyZero)
            {
                // Copy the field files from the <time> directory to the
                // processor*/<time> directories without altering them
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
                                IOobject
                                (
                                    "",
                                    procRunTime.timeName(),
                                    procRunTime
                                ),
                                word::null
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
            }
            else
            {
                // Decompose the fields

                // Search for list of objects for this time
                IOobjectList objects
                (
                    meshes.completeMesh(),
                    runTimes.completeTime().timeName()
                );

                // Construct the vol fields
                PtrList<volScalarField> volScalarFields;
                readFields(meshes.completeMesh(), objects, volScalarFields);
                PtrList<volVectorField> volVectorFields;
                readFields(meshes.completeMesh(), objects, volVectorFields);
                PtrList<volSphericalTensorField> volSphericalTensorFields;
                readFields
                (
                    meshes.completeMesh(),
                    objects,
                    volSphericalTensorFields
                );
                PtrList<volSymmTensorField> volSymmTensorFields;
                readFields(meshes.completeMesh(), objects, volSymmTensorFields);
                PtrList<volTensorField> volTensorFields;
                readFields(meshes.completeMesh(), objects, volTensorFields);

                // Construct the dimensioned fields
                PtrList<DimensionedField<scalar, volMesh>> dimScalarFields;
                readFields(meshes.completeMesh(), objects, dimScalarFields);
                PtrList<DimensionedField<vector, volMesh>> dimVectorFields;
                readFields(meshes.completeMesh(), objects, dimVectorFields);
                PtrList<DimensionedField<sphericalTensor, volMesh>>
                    dimSphericalTensorFields;
                readFields
                (
                    meshes.completeMesh(),
                    objects,
                    dimSphericalTensorFields
                );
                PtrList<DimensionedField<symmTensor, volMesh>>
                    dimSymmTensorFields;
                readFields(meshes.completeMesh(), objects, dimSymmTensorFields);
                PtrList<DimensionedField<tensor, volMesh>> dimTensorFields;
                readFields(meshes.completeMesh(), objects, dimTensorFields);

                // Construct the surface fields
                PtrList<surfaceScalarField> surfaceScalarFields;
                readFields(meshes.completeMesh(), objects, surfaceScalarFields);
                PtrList<surfaceVectorField> surfaceVectorFields;
                readFields(meshes.completeMesh(), objects, surfaceVectorFields);
                PtrList<surfaceSphericalTensorField>
                    surfaceSphericalTensorFields;
                readFields
                (
                    meshes.completeMesh(),
                    objects,
                    surfaceSphericalTensorFields
                );
                PtrList<surfaceSymmTensorField> surfaceSymmTensorFields;
                readFields
                (
                    meshes.completeMesh(),
                    objects,
                    surfaceSymmTensorFields
                );
                PtrList<surfaceTensorField> surfaceTensorFields;
                readFields(meshes.completeMesh(), objects, surfaceTensorFields);

                // Construct the point fields
                const pointMesh& pMesh = pointMesh::New(meshes.completeMesh());
                PtrList<pointScalarField> pointScalarFields;
                readFields(pMesh, objects, pointScalarFields);
                PtrList<pointVectorField> pointVectorFields;
                readFields(pMesh, objects, pointVectorFields);
                PtrList<pointSphericalTensorField> pointSphericalTensorFields;
                readFields(pMesh, objects, pointSphericalTensorFields);
                PtrList<pointSymmTensorField> pointSymmTensorFields;
                readFields(pMesh, objects, pointSymmTensorFields);
                PtrList<pointTensorField> pointTensorFields;
                readFields(pMesh, objects, pointTensorFields);

                // Construct the Lagrangian fields
                fileNameList cloudDirs
                (
                    fileHandler().readDir
                    (
                        runTimes.completeTime().timePath()/cloud::prefix,
                        fileType::directory
                    )
                );
                PtrList<Cloud<indexedParticle>>
                    lagrangianPositions(cloudDirs.size());
                PtrList<List<SLList<indexedParticle*>*>>
                    cellParticles(cloudDirs.size());
                PtrList<PtrList<labelIOField>>
                    lagrangianLabelFields(cloudDirs.size());
                PtrList<PtrList<labelFieldCompactIOField>>
                    lagrangianLabelFieldFields(cloudDirs.size());
                PtrList<PtrList<scalarIOField>>
                    lagrangianScalarFields(cloudDirs.size());
                PtrList<PtrList<scalarFieldCompactIOField>>
                    lagrangianScalarFieldFields(cloudDirs.size());
                PtrList<PtrList<vectorIOField>>
                    lagrangianVectorFields(cloudDirs.size());
                PtrList<PtrList<vectorFieldCompactIOField>>
                    lagrangianVectorFieldFields(cloudDirs.size());
                PtrList<PtrList<sphericalTensorIOField>>
                    lagrangianSphericalTensorFields(cloudDirs.size());
                PtrList<PtrList<sphericalTensorFieldCompactIOField>>
                    lagrangianSphericalTensorFieldFields(cloudDirs.size());
                PtrList<PtrList<symmTensorIOField>>
                    lagrangianSymmTensorFields(cloudDirs.size());
                PtrList<PtrList<symmTensorFieldCompactIOField>>
                    lagrangianSymmTensorFieldFields(cloudDirs.size());
                PtrList<PtrList<tensorIOField>>
                    lagrangianTensorFields(cloudDirs.size());
                PtrList<PtrList<tensorFieldCompactIOField>>
                    lagrangianTensorFieldFields(cloudDirs.size());

                label cloudI = 0;

                forAll(cloudDirs, i)
                {
                    IOobjectList sprayObjs
                    (
                        meshes.completeMesh(),
                        runTimes.completeTime().timeName(),
                        cloud::prefix/cloudDirs[i],
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    );

                    IOobject* positionsPtr = sprayObjs.lookup
                    (
                        word("positions")
                    );

                    if (positionsPtr)
                    {
                        // Read lagrangian particles
                        Info<< "Identified lagrangian data set: "
                            << cloudDirs[i] << endl;
                        lagrangianPositions.set
                        (
                            cloudI,
                            new Cloud<indexedParticle>
                            (
                                meshes.completeMesh(),
                                cloudDirs[i],
                                false
                            )
                        );

                        // Sort particles per cell
                        cellParticles.set
                        (
                            cloudI,
                            new List<SLList<indexedParticle*>*>
                            (
                                meshes.completeMesh().nCells(),
                                static_cast<SLList<indexedParticle*>*>(nullptr)
                            )
                        );

                        // Populate the cloud
                        label index = 0;
                        forAllIter
                        (
                            Cloud<indexedParticle>,
                            lagrangianPositions[cloudI],
                            iter
                        )
                        {
                            iter().index() = index ++;

                            label celli = iter().cell();

                            // Check
                            if
                            (
                                celli < 0
                             || celli >= meshes.completeMesh().nCells()
                            )
                            {
                                FatalErrorInFunction
                                    << "Illegal cell number " << celli
                                    << " for particle with index "
                                    << iter().index()
                                    << " at position "
                                    << iter().position() << nl
                                    << "Cell number should be between 0 and "
                                    << meshes.completeMesh().nCells()-1 << nl
                                    << "On this mesh the particle should"
                                    << " be in cell "
                                    << meshes.completeMesh().findCell
                                       (iter().position())
                                    << exit(FatalError);
                            }

                            if (!cellParticles[cloudI][celli])
                            {
                                cellParticles[cloudI][celli] =
                                    new SLList<indexedParticle*>();
                            }

                            cellParticles[cloudI][celli]->append(&iter());
                        }

                        // Read fields
                        IOobjectList lagrangianObjects
                        (
                            meshes.completeMesh(),
                            runTimes.completeTime().timeName(),
                            cloud::prefix/cloudDirs[cloudI],
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE,
                            false
                        );
                        lagrangianFieldDecomposer::readFields
                        (
                            cloudI,
                            lagrangianObjects,
                            lagrangianLabelFields
                        );
                        lagrangianFieldDecomposer::readFieldFields
                        (
                            cloudI,
                            lagrangianObjects,
                            lagrangianLabelFieldFields
                        );
                        lagrangianFieldDecomposer::readFields
                        (
                            cloudI,
                            lagrangianObjects,
                            lagrangianScalarFields
                        );
                        lagrangianFieldDecomposer::readFieldFields
                        (
                            cloudI,
                            lagrangianObjects,
                            lagrangianScalarFieldFields
                        );
                        lagrangianFieldDecomposer::readFields
                        (
                            cloudI,
                            lagrangianObjects,
                            lagrangianVectorFields
                        );
                        lagrangianFieldDecomposer::readFieldFields
                        (
                            cloudI,
                            lagrangianObjects,
                            lagrangianVectorFieldFields
                        );
                        lagrangianFieldDecomposer::readFields
                        (
                            cloudI,
                            lagrangianObjects,
                            lagrangianSphericalTensorFields
                        );
                        lagrangianFieldDecomposer::readFieldFields
                        (
                            cloudI,
                            lagrangianObjects,
                            lagrangianSphericalTensorFieldFields
                        );
                        lagrangianFieldDecomposer::readFields
                        (
                            cloudI,
                            lagrangianObjects,
                            lagrangianSymmTensorFields
                        );
                        lagrangianFieldDecomposer::readFieldFields
                        (
                            cloudI,
                            lagrangianObjects,
                            lagrangianSymmTensorFieldFields
                        );
                        lagrangianFieldDecomposer::readFields
                        (
                            cloudI,
                            lagrangianObjects,
                            lagrangianTensorFields
                        );
                        lagrangianFieldDecomposer::readFieldFields
                        (
                            cloudI,
                            lagrangianObjects,
                            lagrangianTensorFieldFields
                        );

                        cloudI++;
                    }
                }

                lagrangianPositions.setSize(cloudI);
                cellParticles.setSize(cloudI);
                lagrangianLabelFields.setSize(cloudI);
                lagrangianLabelFieldFields.setSize(cloudI);
                lagrangianScalarFields.setSize(cloudI);
                lagrangianScalarFieldFields.setSize(cloudI);
                lagrangianVectorFields.setSize(cloudI);
                lagrangianVectorFieldFields.setSize(cloudI);
                lagrangianSphericalTensorFields.setSize(cloudI);
                lagrangianSphericalTensorFieldFields.setSize(cloudI);
                lagrangianSymmTensorFields.setSize(cloudI);
                lagrangianSymmTensorFieldFields.setSize(cloudI);
                lagrangianTensorFields.setSize(cloudI);
                lagrangianTensorFieldFields.setSize(cloudI);

                Info<< endl;

                // split the fields over processors
                for (label proci = 0; proci < meshes.nProcs(); proci++)
                {
                    Info<< "Processor " << proci << ": field transfer" << endl;

                    // FV fields
                    {
                        if (!fieldDecomposerList.set(proci))
                        {
                            fieldDecomposerList.set
                            (
                                proci,
                                new fvFieldDecomposer
                                (
                                    meshes.completeMesh(),
                                    meshes.procMeshes()[proci],
                                    meshes.procFaceAddressing()[proci],
                                    meshes.procCellAddressing()[proci],
                                    meshes.procFaceAddressingBf()[proci]
                                )
                            );
                        }
                        const fvFieldDecomposer& fieldDecomposer =
                            fieldDecomposerList[proci];

                        fieldDecomposer.decomposeFields(volScalarFields);
                        fieldDecomposer.decomposeFields(volVectorFields);
                        fieldDecomposer.decomposeFields
                        (
                            volSphericalTensorFields
                        );
                        fieldDecomposer.decomposeFields(volSymmTensorFields);
                        fieldDecomposer.decomposeFields(volTensorFields);

                        fieldDecomposer.decomposeFields(surfaceScalarFields);
                        fieldDecomposer.decomposeFields(surfaceVectorFields);
                        fieldDecomposer.decomposeFields
                        (
                            surfaceSphericalTensorFields
                        );
                        fieldDecomposer.decomposeFields
                        (
                            surfaceSymmTensorFields
                        );
                        fieldDecomposer.decomposeFields(surfaceTensorFields);

                        if (times.size() == 1)
                        {
                            // Clear cached decomposer
                            fieldDecomposerList.set(proci, nullptr);
                        }
                    }

                    // Dimensioned fields
                    {
                        if (!dimFieldDecomposerList.set(proci))
                        {
                            dimFieldDecomposerList.set
                            (
                                proci,
                                new dimFieldDecomposer
                                (
                                    meshes.completeMesh(),
                                    meshes.procMeshes()[proci],
                                    meshes.procFaceAddressing()[proci],
                                    meshes.procCellAddressing()[proci]
                                )
                            );
                        }
                        const dimFieldDecomposer& dimDecomposer =
                            dimFieldDecomposerList[proci];

                        dimDecomposer.decomposeFields(dimScalarFields);
                        dimDecomposer.decomposeFields(dimVectorFields);
                        dimDecomposer.decomposeFields(dimSphericalTensorFields);
                        dimDecomposer.decomposeFields(dimSymmTensorFields);
                        dimDecomposer.decomposeFields(dimTensorFields);

                        if (times.size() == 1)
                        {
                            dimFieldDecomposerList.set(proci, nullptr);
                        }
                    }

                    // Point fields
                    if
                    (
                        pointScalarFields.size()
                     || pointVectorFields.size()
                     || pointSphericalTensorFields.size()
                     || pointSymmTensorFields.size()
                     || pointTensorFields.size()
                    )
                    {
                        const pointMesh& procPMesh =
                            pointMesh::New(meshes.procMeshes()[proci]);

                        if (!pointFieldDecomposerList.set(proci))
                        {
                            pointFieldDecomposerList.set
                            (
                                proci,
                                new pointFieldDecomposer
                                (
                                    pMesh,
                                    procPMesh,
                                    meshes.procPointAddressing()[proci]
                                )
                            );
                        }
                        const pointFieldDecomposer& pointDecomposer =
                            pointFieldDecomposerList[proci];

                        pointDecomposer.decomposeFields(pointScalarFields);
                        pointDecomposer.decomposeFields(pointVectorFields);
                        pointDecomposer.decomposeFields
                        (
                            pointSphericalTensorFields
                        );
                        pointDecomposer.decomposeFields(pointSymmTensorFields);
                        pointDecomposer.decomposeFields(pointTensorFields);

                        if (times.size() == 1)
                        {
                            pointFieldDecomposerList.set(proci, nullptr);
                        }
                    }

                    // If there is lagrangian data write it out
                    forAll(lagrangianPositions, cloudI)
                    {
                        if (lagrangianPositions[cloudI].size())
                        {
                            lagrangianFieldDecomposer fieldDecomposer
                            (
                                meshes.completeMesh(),
                                meshes.procMeshes()[proci],
                                meshes.procFaceAddressing()[proci],
                                meshes.procCellAddressing()[proci],
                                cloudDirs[cloudI],
                                lagrangianPositions[cloudI],
                                cellParticles[cloudI]
                            );

                            // Lagrangian fields
                            {
                                fieldDecomposer.decomposeFields
                                (
                                    cloudDirs[cloudI],
                                    lagrangianLabelFields[cloudI]
                                );
                                fieldDecomposer.decomposeFieldFields
                                (
                                    cloudDirs[cloudI],
                                    lagrangianLabelFieldFields[cloudI]
                                );
                                fieldDecomposer.decomposeFields
                                (
                                    cloudDirs[cloudI],
                                    lagrangianScalarFields[cloudI]
                                );
                                fieldDecomposer.decomposeFieldFields
                                (
                                    cloudDirs[cloudI],
                                    lagrangianScalarFieldFields[cloudI]
                                );
                                fieldDecomposer.decomposeFields
                                (
                                    cloudDirs[cloudI],
                                    lagrangianVectorFields[cloudI]
                                );
                                fieldDecomposer.decomposeFieldFields
                                (
                                    cloudDirs[cloudI],
                                    lagrangianVectorFieldFields[cloudI]
                                );
                                fieldDecomposer.decomposeFields
                                (
                                    cloudDirs[cloudI],
                                    lagrangianSphericalTensorFields[cloudI]
                                );
                                fieldDecomposer.decomposeFieldFields
                                (
                                    cloudDirs[cloudI],
                                    lagrangianSphericalTensorFieldFields[cloudI]
                                );
                                fieldDecomposer.decomposeFields
                                (
                                    cloudDirs[cloudI],
                                    lagrangianSymmTensorFields[cloudI]
                                );
                                fieldDecomposer.decomposeFieldFields
                                (
                                    cloudDirs[cloudI],
                                    lagrangianSymmTensorFieldFields[cloudI]
                                );
                                fieldDecomposer.decomposeFields
                                (
                                    cloudDirs[cloudI],
                                    lagrangianTensorFields[cloudI]
                                );
                                fieldDecomposer.decomposeFieldFields
                                (
                                    cloudDirs[cloudI],
                                    lagrangianTensorFieldFields[cloudI]
                                );
                            }
                        }
                    }

                    // Decompose the "uniform" directory in the region time
                    // directory
                    decomposeUniform
                    (
                        copyUniform,
                        distributeUniform,
                        runTimes.completeTime(),
                        meshes.procMeshes()[proci].time(),
                        regionDir
                    );

                    // For the first region of a multi-region case additionally
                    // decompose the "uniform" directory in the no-region time
                    // directory
                    if (regionNames.size() > 1 && regioni == 0)
                    {
                        decomposeUniform
                        (
                            copyUniform,
                            distributeUniform,
                            runTimes.completeTime(),
                            meshes.procMeshes()[proci].time()
                        );
                    }
                }
            }
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
