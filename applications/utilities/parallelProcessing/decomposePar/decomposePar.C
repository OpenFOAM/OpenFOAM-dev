/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

    - decomposePar [OPTION]

    \param -cellDist \n
    Write the cell distribution as a labelList, for use with 'manual'
    decomposition method or as a volScalarField for post-processing.

    \param -region regionName \n
    Decompose named region. Does not check for existence of processor*.

    \param -allRegions \n
    Decompose all regions in regionProperties. Does not check for
    existence of processor*.

    \param -copyUniform \n
    Copy any \a uniform directories too.

    \param -constant \n
    \param -time xxx:yyy \n
    Override controlDict settings and decompose selected times. Does not
    re-decompose the mesh i.e. does not handle moving mesh or changing
    mesh cases.

    \param -fields \n
    Use existing geometry decomposition and convert fields only.

    \param -noSets \n
    Skip decomposing cellSets, faceSets, pointSets.

    \param -force \n
    Remove any existing \a processor subdirectories before decomposing the
    geometry.

    \param -ifRequired \n
    Only decompose the geometry if the number of domains has changed from a
    previous decomposition. No \a processor subdirectories will be removed
    unless the \a -force option is also specified. This option can be used
    to avoid redundant geometry decomposition (eg, in scripts), but should
    be used with caution when the underlying (serial) geometry or the
    decomposition method etc. have been changed between decompositions.

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "fvCFD.H"
#include "IOobjectList.H"
#include "domainDecomposition.H"
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
#include "pointFields.H"
#include "regionProperties.H"

#include "readFields.H"
#include "dimFieldDecomposer.H"
#include "fvFieldDecomposer.H"
#include "pointFieldDecomposer.H"
#include "lagrangianFieldDecomposer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const labelIOList& procAddressing
(
    const PtrList<fvMesh>& procMeshList,
    const label procI,
    const word& name,
    PtrList<labelIOList>& procAddressingList
)
{
    const fvMesh& procMesh = procMeshList[procI];

    if (!procAddressingList.set(procI))
    {
        procAddressingList.set
        (
            procI,
            new labelIOList
            (
                IOobject
                (
                    name,
                    procMesh.facesInstance(),
                    procMesh.meshSubDir,
                    procMesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );
    }
    return procAddressingList[procI];
}



int main(int argc, char *argv[])
{
    argList::addNote
    (
        "decompose a mesh and fields of a case for parallel execution"
    );

    argList::noParallel();
    #include "addRegionOption.H"
    argList::addBoolOption
    (
        "allRegions",
        "operate on all regions in regionProperties"
    );
    argList::addBoolOption
    (
        "cellDist",
        "write cell distribution as a labelList - for use with 'manual' "
        "decomposition method or as a volScalarField for post-processing."
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
        "noSets",
        "skip decomposing cellSets, faceSets, pointSets"
    );
    argList::addBoolOption
    (
        "force",
        "remove existing processor*/ subdirs before decomposing the geometry"
    );
    argList::addBoolOption
    (
        "ifRequired",
        "only decompose geometry if the number of domains has changed"
    );

    // Include explicit constant options, have zero from time range
    timeSelector::addOptions(true, false);

    #include "setRootCase.H"

    bool allRegions              = args.optionFound("allRegions");
    bool writeCellDist           = args.optionFound("cellDist");
    bool copyUniform             = args.optionFound("copyUniform");
    bool decomposeFieldsOnly     = args.optionFound("fields");
    bool decomposeSets           = !args.optionFound("noSets");
    bool forceOverwrite          = args.optionFound("force");
    bool ifRequiredDecomposition = args.optionFound("ifRequired");

    // Set time from database
    #include "createTime.H"
    // Allow override of time
    instantList times = timeSelector::selectIfPresent(runTime, args);


    wordList regionNames;
    wordList regionDirs;
    if (allRegions)
    {
        Info<< "Decomposing all regions in regionProperties" << nl << endl;
        regionProperties rp(runTime);
        forAllConstIter(HashTable<wordList>, rp, iter)
        {
            const wordList& regions = iter();
            forAll(regions, i)
            {
                if (findIndex(regionNames, regions[i]) == -1)
                {
                    regionNames.append(regions[i]);
                }
            }
        }
        regionDirs = regionNames;
    }
    else
    {
        word regionName;
        if (args.optionReadIfPresent("region", regionName))
        {
            regionNames = wordList(1, regionName);
            regionDirs = regionNames;
        }
        else
        {
            regionNames = wordList(1, fvMesh::defaultRegion);
            regionDirs = wordList(1, word::null);
        }
    }



    forAll(regionNames, regionI)
    {
        const word& regionName = regionNames[regionI];
        const word& regionDir = regionDirs[regionI];

        Info<< "\n\nDecomposing mesh " << regionName << nl << endl;


        // determine the existing processor count directly
        label nProcs = 0;
        while
        (
            isDir
            (
                runTime.path()
              / (word("processor") + name(nProcs))
              / runTime.constant()
              / regionDir
              / polyMesh::meshSubDir
            )
        )
        {
            ++nProcs;
        }

        // get requested numberOfSubdomains
        const label nDomains = readLabel
        (
            IOdictionary
            (
                IOobject
                (
                    "decomposeParDict",
                    runTime.time().system(),
                    regionDir,          // use region if non-standard
                    runTime,
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE,
                    false
                )
            ).lookup("numberOfSubdomains")
        );

        if (decomposeFieldsOnly)
        {
            // Sanity check on previously decomposed case
            if (nProcs != nDomains)
            {
                FatalErrorIn(args.executable())
                    << "Specified -fields, but the case was decomposed with "
                    << nProcs << " domains"
                    << nl
                    << "instead of " << nDomains
                    << " domains as specified in decomposeParDict"
                    << nl
                    << exit(FatalError);
            }
        }
        else if (nProcs)
        {
            bool procDirsProblem = true;

            if (ifRequiredDecomposition && nProcs == nDomains)
            {
                // we can reuse the decomposition
                decomposeFieldsOnly = true;
                procDirsProblem = false;
                forceOverwrite = false;

                Info<< "Using existing processor directories" << nl;
            }

            if (forceOverwrite)
            {
                Info<< "Removing " << nProcs
                    << " existing processor directories" << endl;

                // remove existing processor dirs
                // reverse order to avoid gaps if someone interrupts the process
                for (label procI = nProcs-1; procI >= 0; --procI)
                {
                    fileName procDir
                    (
                        runTime.path()/(word("processor") + name(procI))
                    );

                    rmDir(procDir);
                }

                procDirsProblem = false;
            }

            if (procDirsProblem)
            {
                FatalErrorIn(args.executable())
                    << "Case is already decomposed with " << nProcs
                    << " domains, use the -force option or manually" << nl
                    << "remove processor directories before decomposing. e.g.,"
                    << nl
                    << "    rm -rf " << runTime.path().c_str() << "/processor*"
                    << nl
                    << exit(FatalError);
            }
        }

        Info<< "Create mesh" << endl;
        domainDecomposition mesh
        (
            IOobject
            (
                regionName,
                runTime.timeName(),
                runTime
            )
        );

        // Decompose the mesh
        if (!decomposeFieldsOnly)
        {
            mesh.decomposeMesh();

            mesh.writeDecomposition(decomposeSets);

            if (writeCellDist)
            {
                const labelList& procIds = mesh.cellToProc();

                // Write the decomposition as labelList for use with 'manual'
                // decomposition method.
                labelIOList cellDecomposition
                (
                    IOobject
                    (
                        "cellDecomposition",
                        mesh.facesInstance(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    procIds
                );
                cellDecomposition.write();

                Info<< nl << "Wrote decomposition to "
                    << cellDecomposition.objectPath()
                    << " for use in manual decomposition." << endl;

                // Write as volScalarField for postprocessing.
                volScalarField cellDist
                (
                    IOobject
                    (
                        "cellDist",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimensionedScalar("cellDist", dimless, 0),
                    zeroGradientFvPatchScalarField::typeName
                );

                forAll(procIds, celli)
                {
                   cellDist[celli] = procIds[celli];
                }

                cellDist.write();

                Info<< nl << "Wrote decomposition as volScalarField to "
                    << cellDist.name() << " for use in postprocessing."
                    << endl;
            }
        }



        // Caches
        // ~~~~~~
        // Cached processor meshes and maps. These are only preserved if running
        // with multiple times.
        PtrList<Time> processorDbList(mesh.nProcs());
        PtrList<fvMesh> procMeshList(mesh.nProcs());
        PtrList<labelIOList> faceProcAddressingList(mesh.nProcs());
        PtrList<labelIOList> cellProcAddressingList(mesh.nProcs());
        PtrList<labelIOList> boundaryProcAddressingList(mesh.nProcs());
        PtrList<fvFieldDecomposer> fieldDecomposerList(mesh.nProcs());
        PtrList<dimFieldDecomposer> dimFieldDecomposerList(mesh.nProcs());
        PtrList<labelIOList> pointProcAddressingList(mesh.nProcs());
        PtrList<pointFieldDecomposer> pointFieldDecomposerList(mesh.nProcs());



        // Loop over all times
        forAll(times, timeI)
        {
            runTime.setTime(times[timeI], timeI);

            Info<< "Time = " << runTime.timeName() << endl;

            // Search for list of objects for this time
            IOobjectList objects(mesh, runTime.timeName());


            // Construct the vol fields
            // ~~~~~~~~~~~~~~~~~~~~~~~~
            PtrList<volScalarField> volScalarFields;
            readFields(mesh, objects, volScalarFields);
            PtrList<volVectorField> volVectorFields;
            readFields(mesh, objects, volVectorFields);
            PtrList<volSphericalTensorField> volSphericalTensorFields;
            readFields(mesh, objects, volSphericalTensorFields);
            PtrList<volSymmTensorField> volSymmTensorFields;
            readFields(mesh, objects, volSymmTensorFields);
            PtrList<volTensorField> volTensorFields;
            readFields(mesh, objects, volTensorFields);


            // Construct the dimensioned fields
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            PtrList<DimensionedField<scalar, volMesh> > dimScalarFields;
            readFields(mesh, objects, dimScalarFields);
            PtrList<DimensionedField<vector, volMesh> > dimVectorFields;
            readFields(mesh, objects, dimVectorFields);
            PtrList<DimensionedField<sphericalTensor, volMesh> >
                dimSphericalTensorFields;
            readFields(mesh, objects, dimSphericalTensorFields);
            PtrList<DimensionedField<symmTensor, volMesh> > dimSymmTensorFields;
            readFields(mesh, objects, dimSymmTensorFields);
            PtrList<DimensionedField<tensor, volMesh> > dimTensorFields;
            readFields(mesh, objects, dimTensorFields);


            // Construct the surface fields
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            PtrList<surfaceScalarField> surfaceScalarFields;
            readFields(mesh, objects, surfaceScalarFields);
            PtrList<surfaceVectorField> surfaceVectorFields;
            readFields(mesh, objects, surfaceVectorFields);
            PtrList<surfaceSphericalTensorField> surfaceSphericalTensorFields;
            readFields(mesh, objects, surfaceSphericalTensorFields);
            PtrList<surfaceSymmTensorField> surfaceSymmTensorFields;
            readFields(mesh, objects, surfaceSymmTensorFields);
            PtrList<surfaceTensorField> surfaceTensorFields;
            readFields(mesh, objects, surfaceTensorFields);


            // Construct the point fields
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~
            const pointMesh& pMesh = pointMesh::New(mesh);

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
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            fileNameList cloudDirs
            (
                readDir(runTime.timePath()/cloud::prefix, fileName::DIRECTORY)
            );

            // Particles
            PtrList<Cloud<indexedParticle> > lagrangianPositions
            (
                cloudDirs.size()
            );
            // Particles per cell
            PtrList< List<SLList<indexedParticle*>*> > cellParticles
            (
                cloudDirs.size()
            );

            PtrList<PtrList<labelIOField> > lagrangianLabelFields
            (
                cloudDirs.size()
            );
            PtrList<PtrList<labelFieldCompactIOField> >
            lagrangianLabelFieldFields
            (
                cloudDirs.size()
            );
            PtrList<PtrList<scalarIOField> > lagrangianScalarFields
            (
                cloudDirs.size()
            );
            PtrList<PtrList<scalarFieldCompactIOField> >
            lagrangianScalarFieldFields
            (
                cloudDirs.size()
            );
            PtrList<PtrList<vectorIOField> > lagrangianVectorFields
            (
                cloudDirs.size()
            );
            PtrList<PtrList<vectorFieldCompactIOField> >
            lagrangianVectorFieldFields
            (
                cloudDirs.size()
            );
            PtrList<PtrList<sphericalTensorIOField> >
            lagrangianSphericalTensorFields
            (
                cloudDirs.size()
            );
            PtrList<PtrList<sphericalTensorFieldCompactIOField> >
                lagrangianSphericalTensorFieldFields(cloudDirs.size());
            PtrList<PtrList<symmTensorIOField> > lagrangianSymmTensorFields
            (
                cloudDirs.size()
            );
            PtrList<PtrList<symmTensorFieldCompactIOField> >
            lagrangianSymmTensorFieldFields
            (
                cloudDirs.size()
            );
            PtrList<PtrList<tensorIOField> > lagrangianTensorFields
            (
                cloudDirs.size()
            );
            PtrList<PtrList<tensorFieldCompactIOField> >
            lagrangianTensorFieldFields
            (
                cloudDirs.size()
            );

            label cloudI = 0;

            forAll(cloudDirs, i)
            {
                IOobjectList sprayObjs
                (
                    mesh,
                    runTime.timeName(),
                    cloud::prefix/cloudDirs[i]
                );

                IOobject* positionsPtr = sprayObjs.lookup(word("positions"));

                if (positionsPtr)
                {
                    // Read lagrangian particles
                    // ~~~~~~~~~~~~~~~~~~~~~~~~~

                    Info<< "Identified lagrangian data set: " << cloudDirs[i]
                        << endl;

                    lagrangianPositions.set
                    (
                        cloudI,
                        new Cloud<indexedParticle>
                        (
                            mesh,
                            cloudDirs[i],
                            false
                        )
                    );


                    // Sort particles per cell
                    // ~~~~~~~~~~~~~~~~~~~~~~~

                    cellParticles.set
                    (
                        cloudI,
                        new List<SLList<indexedParticle*>*>
                        (
                            mesh.nCells(),
                            static_cast<SLList<indexedParticle*>*>(NULL)
                        )
                    );

                    label i = 0;

                    forAllIter
                    (
                        Cloud<indexedParticle>,
                        lagrangianPositions[cloudI],
                        iter
                    )
                    {
                        iter().index() = i++;

                        label celli = iter().cell();

                        // Check
                        if (celli < 0 || celli >= mesh.nCells())
                        {
                            FatalErrorIn(args.executable())
                                << "Illegal cell number " << celli
                                << " for particle with index " << iter().index()
                                << " at position " << iter().position() << nl
                                << "Cell number should be between 0 and "
                                << mesh.nCells()-1 << nl
                                << "On this mesh the particle should"
                                << " be in cell "
                                << mesh.findCell(iter().position())
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
                    // ~~~~~~~~~~~

                    IOobjectList lagrangianObjects
                    (
                        mesh,
                        runTime.timeName(),
                        cloud::prefix/cloudDirs[cloudI]
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


            // Any uniform data to copy/link?
            fileName uniformDir("uniform");

            if (isDir(runTime.timePath()/uniformDir))
            {
                Info<< "Detected additional non-decomposed files in "
                    << runTime.timePath()/uniformDir
                    << endl;
            }
            else
            {
                uniformDir.clear();
            }

            Info<< endl;

            // split the fields over processors
            for (label procI = 0; procI < mesh.nProcs(); procI++)
            {
                Info<< "Processor " << procI << ": field transfer" << endl;


                // open the database
                if (!processorDbList.set(procI))
                {
                    processorDbList.set
                    (
                        procI,
                        new Time
                        (
                            Time::controlDictName,
                            args.rootPath(),
                            args.caseName()
                           /fileName(word("processor") + name(procI))
                        )
                    );
                }
                Time& processorDb = processorDbList[procI];


                processorDb.setTime(runTime);

                // remove files remnants that can cause horrible problems
                // - mut and nut are used to mark the new turbulence models,
                //   their existence prevents old models from being upgraded
                {
                    fileName timeDir(processorDb.path()/processorDb.timeName());

                    rm(timeDir/"mut");
                    rm(timeDir/"nut");
                }

                // read the mesh
                if (!procMeshList.set(procI))
                {
                    procMeshList.set
                    (
                        procI,
                        new fvMesh
                        (
                            IOobject
                            (
                                regionName,
                                processorDb.timeName(),
                                processorDb
                            )
                        )
                    );
                }
                const fvMesh& procMesh = procMeshList[procI];

                const labelIOList& faceProcAddressing = procAddressing
                (
                    procMeshList,
                    procI,
                    "faceProcAddressing",
                    faceProcAddressingList
                );

                const labelIOList& cellProcAddressing = procAddressing
                (
                    procMeshList,
                    procI,
                    "cellProcAddressing",
                    cellProcAddressingList
                );

                const labelIOList& boundaryProcAddressing = procAddressing
                (
                    procMeshList,
                    procI,
                    "boundaryProcAddressing",
                    boundaryProcAddressingList
                );


                // FV fields
                {
                    if (!fieldDecomposerList.set(procI))
                    {
                        fieldDecomposerList.set
                        (
                            procI,
                            new fvFieldDecomposer
                            (
                                mesh,
                                procMesh,
                                faceProcAddressing,
                                cellProcAddressing,
                                boundaryProcAddressing
                            )
                        );
                    }
                    const fvFieldDecomposer& fieldDecomposer =
                        fieldDecomposerList[procI];

                    fieldDecomposer.decomposeFields(volScalarFields);
                    fieldDecomposer.decomposeFields(volVectorFields);
                    fieldDecomposer.decomposeFields(volSphericalTensorFields);
                    fieldDecomposer.decomposeFields(volSymmTensorFields);
                    fieldDecomposer.decomposeFields(volTensorFields);

                    fieldDecomposer.decomposeFields(surfaceScalarFields);
                    fieldDecomposer.decomposeFields(surfaceVectorFields);
                    fieldDecomposer.decomposeFields
                    (
                        surfaceSphericalTensorFields
                    );
                    fieldDecomposer.decomposeFields(surfaceSymmTensorFields);
                    fieldDecomposer.decomposeFields(surfaceTensorFields);

                    if (times.size() == 1)
                    {
                        // Clear cached decomposer
                        fieldDecomposerList.set(procI, NULL);
                    }
                }

                // Dimensioned fields
                {
                    if (!dimFieldDecomposerList.set(procI))
                    {
                        dimFieldDecomposerList.set
                        (
                            procI,
                            new dimFieldDecomposer
                            (
                                mesh,
                                procMesh,
                                faceProcAddressing,
                                cellProcAddressing
                            )
                        );
                    }
                    const dimFieldDecomposer& dimDecomposer =
                        dimFieldDecomposerList[procI];

                    dimDecomposer.decomposeFields(dimScalarFields);
                    dimDecomposer.decomposeFields(dimVectorFields);
                    dimDecomposer.decomposeFields(dimSphericalTensorFields);
                    dimDecomposer.decomposeFields(dimSymmTensorFields);
                    dimDecomposer.decomposeFields(dimTensorFields);

                    if (times.size() == 1)
                    {
                        dimFieldDecomposerList.set(procI, NULL);
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
                    const labelIOList& pointProcAddressing = procAddressing
                    (
                        procMeshList,
                        procI,
                        "pointProcAddressing",
                        pointProcAddressingList
                    );

                    const pointMesh& procPMesh = pointMesh::New(procMesh);

                    if (!pointFieldDecomposerList.set(procI))
                    {
                        pointFieldDecomposerList.set
                        (
                            procI,
                            new pointFieldDecomposer
                            (
                                pMesh,
                                procPMesh,
                                pointProcAddressing,
                                boundaryProcAddressing
                            )
                        );
                    }
                    const pointFieldDecomposer& pointDecomposer =
                        pointFieldDecomposerList[procI];

                    pointDecomposer.decomposeFields(pointScalarFields);
                    pointDecomposer.decomposeFields(pointVectorFields);
                    pointDecomposer.decomposeFields(pointSphericalTensorFields);
                    pointDecomposer.decomposeFields(pointSymmTensorFields);
                    pointDecomposer.decomposeFields(pointTensorFields);


                    if (times.size() == 1)
                    {
                        pointProcAddressingList.set(procI, NULL);
                        pointFieldDecomposerList.set(procI, NULL);
                    }
                }


                // If there is lagrangian data write it out
                forAll(lagrangianPositions, cloudI)
                {
                    if (lagrangianPositions[cloudI].size())
                    {
                        lagrangianFieldDecomposer fieldDecomposer
                        (
                            mesh,
                            procMesh,
                            faceProcAddressing,
                            cellProcAddressing,
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


                // Any non-decomposed data to copy?
                if (uniformDir.size())
                {
                    const fileName timePath = processorDb.timePath();

                    if (copyUniform || mesh.distributed())
                    {
                        cp
                        (
                            runTime.timePath()/uniformDir,
                            timePath/uniformDir
                        );
                    }
                    else
                    {
                        // link with relative paths
                        const string parentPath = string("..")/"..";

                        fileName currentDir(cwd());
                        chDir(timePath);
                        ln
                        (
                            parentPath/runTime.timeName()/uniformDir,
                            uniformDir
                        );
                        chDir(currentDir);
                    }
                }



                // We have cached all the constant mesh data for the current
                // processor. This is only important if running with multiple
                // times, otherwise it is just extra storage.
                if (times.size() == 1)
                {
                    boundaryProcAddressingList.set(procI, NULL);
                    cellProcAddressingList.set(procI, NULL);
                    faceProcAddressingList.set(procI, NULL);
                    procMeshList.set(procI, NULL);
                    processorDbList.set(procI, NULL);
                }
            }
        }
    }

    Info<< "\nEnd.\n" << endl;

    return 0;
}


// ************************************************************************* //
