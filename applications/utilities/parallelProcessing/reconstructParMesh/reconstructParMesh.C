/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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
    reconstructParMesh

Description
    Reconstructs a mesh using geometric information only.

    Writes point/face/cell procAddressing so afterwards reconstructPar can be
    used to reconstruct fields.

    Note:
    - uses geometric matching tolerance (set with -mergeTol (at your option)

    If the parallel case does not have correct procBoundaries use the
    -fullMatch option which will check all boundary faces (bit slower).

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"

#include "IOobjectList.H"
#include "labelIOList.H"
#include "processorPolyPatch.H"
#include "mapAddedPolyMesh.H"
#include "polyMeshAdder.H"
#include "faceCoupleInfo.H"
#include "fvMeshAdder.H"
#include "polyTopoChange.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "regionProperties.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Tolerance (as fraction of the bounding box). Needs to be fairly lax since
// usually meshes get written with limited precision (6 digits)
static const scalar defaultMergeTol = 1e-7;


static void renumber
(
    const labelList& map,
    labelList& elems
)
{
    forAll(elems, i)
    {
        if (elems[i] >= 0)
        {
            elems[i] = map[elems[i]];
        }
    }
}


// Determine which faces are coupled. Uses geometric merge distance.
// Looks either at all boundaryFaces (fullMatch) or only at the
// procBoundaries for proci. Assumes that masterMesh contains already merged
// all the processors < proci.
autoPtr<faceCoupleInfo> determineCoupledFaces
(
    const bool fullMatch,
    const label masterMeshProcStart,
    const label masterMeshProcEnd,
    const polyMesh& masterMesh,
    const label meshToAddProcStart,
    const label meshToAddProcEnd,
    const polyMesh& meshToAdd,
    const scalar mergeDist
)
{
    if (fullMatch || masterMesh.nCells() == 0)
    {
        return autoPtr<faceCoupleInfo>
        (
            new faceCoupleInfo
            (
                masterMesh,
                meshToAdd,
                mergeDist,      // Absolute merging distance
                true            // Matching faces identical
            )
        );
    }
    else
    {
        // Pick up all patches on masterMesh ending in "toDDD" where DDD is
        // the processor number proci.

        const polyBoundaryMesh& masterPatches = masterMesh.boundaryMesh();


        DynamicList<label> masterFaces
        (
            masterMesh.nFaces()
          - masterMesh.nInternalFaces()
        );


        forAll(masterPatches, patchi)
        {
            const polyPatch& pp = masterPatches[patchi];

            if (isA<processorPolyPatch>(pp))
            {
                for
                (
                    label proci=meshToAddProcStart;
                    proci<meshToAddProcEnd;
                    proci++
                )
                {
                    const string toProcString("to" + name(proci));
                    if (
                        pp.name().rfind(toProcString)
                     == (pp.name().size()-toProcString.size())
                    )
                    {
                        label meshFacei = pp.start();
                        forAll(pp, i)
                        {
                            masterFaces.append(meshFacei++);
                        }
                        break;
                    }
                }

            }
        }
        masterFaces.shrink();


        // Pick up all patches on meshToAdd ending in "procBoundaryDDDtoYYY"
        // where DDD is the processor number proci and YYY is < proci.

        const polyBoundaryMesh& addPatches = meshToAdd.boundaryMesh();

        DynamicList<label> addFaces
        (
            meshToAdd.nFaces()
          - meshToAdd.nInternalFaces()
        );

        forAll(addPatches, patchi)
        {
            const polyPatch& pp = addPatches[patchi];

            if (isA<processorPolyPatch>(pp))
            {
                bool isConnected = false;

                for
                (
                    label mergedProci=masterMeshProcStart;
                    !isConnected && (mergedProci < masterMeshProcEnd);
                    mergedProci++
                )
                {
                    for
                    (
                        label proci = meshToAddProcStart;
                        proci < meshToAddProcEnd;
                        proci++
                    )
                    {
                        const word fromProcString
                        (
                            processorPolyPatch::newName(proci, mergedProci)
                        );

                        if (pp.name() == fromProcString)
                        {
                            isConnected = true;
                            break;
                        }
                    }
                }

                if (isConnected)
                {
                    label meshFacei = pp.start();
                    forAll(pp, i)
                    {
                        addFaces.append(meshFacei++);
                    }
                }
            }
        }
        addFaces.shrink();

        return autoPtr<faceCoupleInfo>
        (
            new faceCoupleInfo
            (
                masterMesh,
                masterFaces,
                meshToAdd,
                addFaces,
                mergeDist,      // Absolute merging distance
                true,           // Matching faces identical?
                false,          // If perfect match are faces already ordered
                                // (e.g. processor patches)
                false           // are faces each on separate patch?
            )
        );
    }
}


autoPtr<mapPolyMesh> mergeSharedPoints
(
    const scalar mergeDist,
    polyMesh& mesh,
    labelListList& pointProcAddressing
)
{
    // Find out which sets of points get merged and create a map from
    // mesh point to unique point.
    Map<label> pointToMaster
    (
        fvMeshAdder::findSharedPoints
        (
            mesh,
            mergeDist
        )
    );

    Info<< "mergeSharedPoints : detected " << pointToMaster.size()
        << " points that are to be merged." << endl;

    if (returnReduce(pointToMaster.size(), sumOp<label>()) == 0)
    {
        return autoPtr<mapPolyMesh>(nullptr);
    }

    polyTopoChange meshMod(mesh);

    fvMeshAdder::mergePoints(mesh, pointToMaster, meshMod);

    // Change the mesh (no inflation). Note: parallel comms allowed.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false, true);

    // Update fields. No inflation, parallel sync.
    mesh.updateMesh(map);

    // pointProcAddressing give indices into the master mesh so adapt them
    // for changed point numbering.

    // Adapt constructMaps for merged points.
    forAll(pointProcAddressing, proci)
    {
        labelList& constructMap = pointProcAddressing[proci];

        forAll(constructMap, i)
        {
            label oldPointi = constructMap[i];

            // New label of point after changeMesh.
            label newPointi = map().reversePointMap()[oldPointi];

            if (newPointi < -1)
            {
                constructMap[i] = -newPointi-2;
            }
            else if (newPointi >= 0)
            {
                constructMap[i] = newPointi;
            }
            else
            {
                FatalErrorInFunction
                    << "Problem. oldPointi:" << oldPointi
                    << " newPointi:" << newPointi << abort(FatalError);
            }
        }
    }

    return map;
}


boundBox procBounds
(
    const argList& args,
    const PtrList<Time>& databases,
    const word& regionDir
)
{
    boundBox bb = boundBox::invertedBox;

    forAll(databases, proci)
    {
        fileName pointsInstance
        (
            databases[proci].findInstance
            (
                regionDir/polyMesh::meshSubDir,
                "points"
            )
        );

        if (pointsInstance != databases[proci].timeName())
        {
            FatalErrorInFunction
                << "Your time was specified as " << databases[proci].timeName()
                << " but there is no polyMesh/points in that time." << endl
                << "(there is a points file in " << pointsInstance
                << ")" << endl
                << "Please rerun with the correct time specified"
                << " (through the -constant, -time or -latestTime "
                << "(at your option)."
                << endl << exit(FatalError);
        }

        Info<< "Reading points from "
            << databases[proci].caseName()
            << " for time = " << databases[proci].timeName()
            << nl << endl;

        pointIOField points
        (
            IOobject
            (
                "points",
                databases[proci].findInstance
                (
                    regionDir/polyMesh::meshSubDir,
                    "points"
                ),
                regionDir/polyMesh::meshSubDir,
                databases[proci],
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        boundBox domainBb(points, false);

        bb.min() = min(bb.min(), domainBb.min());
        bb.max() = max(bb.max(), domainBb.max());
    }

    return bb;
}


void writeCellDistance
(
    Time& runTime,
    const fvMesh& masterMesh,
    const labelListList& cellProcAddressing

)
{
    // Write the decomposition as labelList for use with 'manual'
    // decomposition method.
    labelIOList cellDecomposition
    (
        IOobject
        (
            "cellDecomposition",
            masterMesh.facesInstance(),
            masterMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        masterMesh.nCells()
    );

    forAll(cellProcAddressing, proci)
    {
        const labelList& pCells = cellProcAddressing[proci];
        UIndirectList<label>(cellDecomposition, pCells) = proci;
    }

    cellDecomposition.write();

    Info<< nl << "Wrote decomposition to "
        << cellDecomposition.objectPath()
        << " for use in manual decomposition." << endl;


    // Write as volScalarField for postprocessing. Change time to 0
    // if was 'constant'
    {
        const scalar oldTime = runTime.value();
        const label oldIndex = runTime.timeIndex();
        if (runTime.timeName() == runTime.constant() && oldIndex == 0)
        {
            runTime.setTime(0, oldIndex+1);
        }

        volScalarField cellDist
        (
            IOobject
            (
                "cellDist",
                runTime.timeName(),
                masterMesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            masterMesh,
            dimensionedScalar("cellDist", dimless, 0),
            extrapolatedCalculatedFvPatchScalarField::typeName
        );

        forAll(cellDecomposition, celli)
        {
            cellDist[celli] = cellDecomposition[celli];
        }
        cellDist.correctBoundaryConditions();

        cellDist.write();

        Info<< nl << "Wrote decomposition as volScalarField to "
            << cellDist.name() << " for use in postprocessing."
            << endl;

        // Restore time
        runTime.setTime(oldTime, oldIndex);
    }
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "reconstruct a mesh using geometric information only"
    );

    // Enable -constant ... if someone really wants it
    // Enable -withZero to prevent accidentally trashing the initial fields
    timeSelector::addOptions(true, true);
    argList::noParallel();
    argList::addOption
    (
        "mergeTol",
        "scalar",
        "specify the merge distance relative to the bounding box size "
        "(default 1e-7)"
    );
    argList::addBoolOption
    (
        "fullMatch",
        "do (slower) geometric matching on all boundary faces"
    );
    argList::addBoolOption
    (
        "cellDist",
        "write cell distribution as a labelList - for use with 'manual' "
        "decomposition method or as a volScalarField for post-processing."
    );

    #include "addRegionOption.H"
    #include "addAllRegionsOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    Info<< "This is an experimental tool which tries to merge"
        << " individual processor" << nl
        << "meshes back into one master mesh. Use it if the original"
        << " master mesh has" << nl
        << "been deleted or if the processor meshes have been modified"
        << " (topology change)." << nl
        << "This tool will write the resulting mesh to a new time step"
        << " and construct" << nl
        << "xxxxProcAddressing files in the processor meshes so"
        << " reconstructPar can be" << nl
        << "used to regenerate the fields on the master mesh." << nl
        << nl
        << "Not well tested & use at your own risk!" << nl
        << endl;


    const wordList regionNames(selectRegionNames(args, runTime));
    if (regionNames.size() > 1)
    {
        Info<< "Operating on regions " << regionNames[0];
        for (label regioni = 1; regioni < regionNames.size() - 1; ++ regioni)
        {
            Info<< ", " << regionNames[regioni];
        }
        Info<< " and " << regionNames.last() << nl << endl;
    }
    else if (regionNames[0] != polyMesh::defaultRegion)
    {
        Info<< "Operating on region " << regionNames[0] << nl << endl;
    }


    scalar mergeTol = defaultMergeTol;
    args.optionReadIfPresent("mergeTol", mergeTol);

    scalar writeTol = Foam::pow(10.0, -scalar(IOstream::defaultPrecision()));

    Info<< "Merge tolerance : " << mergeTol << nl
        << "Write tolerance : " << writeTol << endl;

    if (runTime.writeFormat() == IOstream::ASCII && mergeTol < writeTol)
    {
        FatalErrorInFunction
            << "Your current settings specify ASCII writing with "
            << IOstream::defaultPrecision() << " digits precision." << endl
            << "Your merging tolerance (" << mergeTol << ") is finer than this."
            << endl
            << "Please change your writeFormat to binary"
            << " or increase the writePrecision" << endl
            << "or adjust the merge tolerance (-mergeTol)."
            << exit(FatalError);
    }


    const bool fullMatch = args.optionFound("fullMatch");

    if (fullMatch)
    {
        Info<< "Doing geometric matching on all boundary faces." << nl << endl;
    }
    else
    {
        Info<< "Doing geometric matching on correct procBoundaries only."
            << nl << "This assumes a correct decomposition." << endl;
    }

    bool writeCellDist = args.optionFound("cellDist");


    label nProcs = fileHandler().nProcs(args.path());

    Info<< "Found " << nProcs << " processor directories" << nl << endl;


    // Read all time databases
    PtrList<Time> databases(nProcs);

    forAll(databases, proci)
    {
        Info<< "Reading database "
            << args.caseName()/fileName(word("processor") + name(proci))
            << endl;

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

    // Loop over all times
    forAll(timeDirs, timeI)
    {
        // Set time for global database
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Set time for all databases
        forAll(databases, proci)
        {
            databases[proci].setTime(timeDirs[timeI], timeI);
        }

        forAll(regionNames, regioni)
        {
            const word& regionName = regionNames[regioni];
            const word regionDir =
                regionName == polyMesh::defaultRegion
              ? word::null
              : regionName;

            IOobject facesIO
            (
                "faces",
                databases[0].timeName(),
                regionDir/polyMesh::meshSubDir,
                databases[0],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            );


            // Problem: faceCompactIOList recognises both 'faceList' and
            //          'faceCompactList' so we should be lenient when doing
            //          typeHeaderOk
            if (!facesIO.typeHeaderOk<faceCompactIOList>(false))
            {
                Info<< "No mesh." << nl << endl;
                continue;
            }


            // Read point on individual processors to determine merge tolerance
            // (otherwise single cell domains might give problems)

            const boundBox bb = procBounds(args, databases, regionDir);
            const scalar mergeDist = mergeTol*bb.mag();

            Info<< "Overall mesh bounding box  : " << bb << nl
                << "Relative tolerance         : " << mergeTol << nl
                << "Absolute matching distance : " << mergeDist << nl
                << endl;


            // Addressing from processor to reconstructed case
            labelListList cellProcAddressing(nProcs);
            labelListList faceProcAddressing(nProcs);
            labelListList pointProcAddressing(nProcs);
            labelListList boundaryProcAddressing(nProcs);

            // Internal faces on the final reconstructed mesh
            label masterInternalFaces;

            // Owner addressing on the final reconstructed mesh
            labelList masterOwner;

            {
                // Construct empty mesh.
                // fvMesh** masterMesh = new fvMesh*[nProcs];
                PtrList<fvMesh> masterMesh(nProcs);

                for (label proci=0; proci<nProcs; proci++)
                {
                    masterMesh.set
                    (
                        proci,
                        new fvMesh
                        (
                            IOobject
                            (
                                regionName,
                                runTime.timeName(),
                                runTime,
                                IOobject::NO_READ
                            ),
                            xferCopy(pointField()),
                            xferCopy(faceList()),
                            xferCopy(cellList())
                        )
                    );

                    fvMesh meshToAdd
                    (
                        IOobject
                        (
                            regionName,
                            databases[proci].timeName(),
                            databases[proci]
                        )
                    );

                    // Initialize its addressing
                    cellProcAddressing[proci] = identity(meshToAdd.nCells());
                    faceProcAddressing[proci] = identity(meshToAdd.nFaces());
                    pointProcAddressing[proci] = identity(meshToAdd.nPoints());
                    boundaryProcAddressing[proci] =
                        identity(meshToAdd.boundaryMesh().size());

                    // Find geometrically shared points/faces.
                    autoPtr<faceCoupleInfo> couples = determineCoupledFaces
                    (
                        fullMatch,
                        proci,
                        proci,
                        masterMesh[proci],
                        proci,
                        proci,
                        meshToAdd,
                        mergeDist
                    );

                    // Add elements to mesh
                    autoPtr<mapAddedPolyMesh> map = fvMeshAdder::add
                    (
                        masterMesh[proci],
                        meshToAdd,
                        couples
                    );

                    // Added processor
                    renumber(map().addedCellMap(), cellProcAddressing[proci]);
                    renumber(map().addedFaceMap(), faceProcAddressing[proci]);
                    renumber(map().addedPointMap(), pointProcAddressing[proci]);
                    renumber
                    (
                        map().addedPatchMap(),
                        boundaryProcAddressing[proci]
                    );
                }
                for (label step=2; step<nProcs*2; step*=2)
                {
                    for (label proci=0; proci<nProcs; proci+=step)
                    {
                        label next = proci + step/2;
                        if(next >= nProcs)
                        {
                            continue;
                        }

                        Info<< "Merging mesh " << proci << " with " << next
                            << endl;

                        // Find geometrically shared points/faces.
                        autoPtr<faceCoupleInfo> couples = determineCoupledFaces
                        (
                            fullMatch,
                            proci,
                            next,
                            masterMesh[proci],
                            next,
                            proci+step,
                            masterMesh[next],
                            mergeDist
                        );

                        // Add elements to mesh
                        autoPtr<mapAddedPolyMesh> map = fvMeshAdder::add
                        (
                            masterMesh[proci],
                            masterMesh[next],
                            couples
                        );

                        // Processors that were already in masterMesh
                        for (label mergedI=proci; mergedI<next; mergedI++)
                        {
                            renumber
                            (
                                map().oldCellMap(),
                                cellProcAddressing[mergedI]
                            );

                            renumber
                            (
                                map().oldFaceMap(),
                                faceProcAddressing[mergedI]
                            );

                            renumber
                            (
                                map().oldPointMap(),
                                pointProcAddressing[mergedI]
                            );

                            // Note: boundary is special since can contain -1.
                            renumber
                            (
                                map().oldPatchMap(),
                                boundaryProcAddressing[mergedI]
                            );
                        }

                        // Added processor
                        for
                        (
                            label addedI=next;
                            addedI<min(proci+step, nProcs);
                            addedI++
                        )
                        {
                            renumber
                            (
                                map().addedCellMap(),
                                cellProcAddressing[addedI]
                            );

                            renumber
                            (
                                map().addedFaceMap(),
                                faceProcAddressing[addedI]
                            );

                            renumber
                            (
                                map().addedPointMap(),
                                pointProcAddressing[addedI]
                            );

                            renumber
                            (
                                map().addedPatchMap(),
                                boundaryProcAddressing[addedI]
                            );
                        }

                        masterMesh.set(next, nullptr);
                    }
                }

                for (label proci=0; proci<nProcs; proci++)
                {
                    Info<< "Reading mesh to add from "
                        << databases[proci].caseName()
                        << " for time = " << databases[proci].timeName()
                        << nl << nl << endl;
                }

                // See if any points on the mastermesh have become connected
                // because of connections through processor meshes.
                mergeSharedPoints
                (
                    mergeDist,
                    masterMesh[0],
                    pointProcAddressing
                );

                // Save some properties on the reconstructed mesh
                masterInternalFaces = masterMesh[0].nInternalFaces();
                masterOwner = masterMesh[0].faceOwner();


                Info<< "\nWriting merged mesh to "
                    << runTime.path()/runTime.timeName()
                    << nl << endl;

                if (!masterMesh[0].write())
                {
                    FatalErrorInFunction
                        << "Failed writing polyMesh."
                        << exit(FatalError);
                }

                if (writeCellDist)
                {
                    writeCellDistance
                    (
                        runTime,
                        masterMesh[0],
                        cellProcAddressing
                    );
                }
            }


            // Write the addressing

            Info<< "Reconstructing the addressing from the processor meshes"
                << " to the newly reconstructed mesh" << nl << endl;

            forAll(databases, proci)
            {
                Info<< "Reading processor " << proci << " mesh from "
                    << databases[proci].caseName() << endl;

                polyMesh procMesh
                (
                    IOobject
                    (
                        regionName,
                        databases[proci].timeName(),
                        databases[proci]
                    )
                );


                // From processor point to reconstructed mesh point

                Info<< "Writing pointProcAddressing to "
                    << databases[proci].caseName()
                      /procMesh.facesInstance()
                      /polyMesh::meshSubDir
                    << endl;

                labelIOList
                (
                    IOobject
                    (
                        "pointProcAddressing",
                        procMesh.facesInstance(),
                        polyMesh::meshSubDir,
                        procMesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false                       // Do not register
                    ),
                    pointProcAddressing[proci]
                ).write();


                // From processor face to reconstructed mesh face

                Info<< "Writing faceProcAddressing to "
                    << databases[proci].caseName()
                      /procMesh.facesInstance()
                      /polyMesh::meshSubDir
                    << endl;

                labelIOList faceProcAddr
                (
                    IOobject
                    (
                        "faceProcAddressing",
                        procMesh.facesInstance(),
                        polyMesh::meshSubDir,
                        procMesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false                       // Do not register
                    ),
                    faceProcAddressing[proci]
                );

                // Now add turning index to faceProcAddressing.
                // See reconstructPar for meaning of turning index.
                forAll(faceProcAddr, procFacei)
                {
                    const label masterFacei = faceProcAddr[procFacei];

                    if
                    (
                       !procMesh.isInternalFace(procFacei)
                     && masterFacei < masterInternalFaces
                    )
                    {
                        // proc face is now external but used to be internal
                        // face. Check if we have owner or neighbour.

                        label procOwn = procMesh.faceOwner()[procFacei];
                        label masterOwn = masterOwner[masterFacei];

                        if (cellProcAddressing[proci][procOwn] == masterOwn)
                        {
                            // No turning. Offset by 1.
                            faceProcAddr[procFacei]++;
                        }
                        else
                        {
                            // Turned face.
                            faceProcAddr[procFacei] =
                                -1 - faceProcAddr[procFacei];
                        }
                    }
                    else
                    {
                        // No turning. Offset by 1.
                        faceProcAddr[procFacei]++;
                    }
                }

                faceProcAddr.write();


                // From processor cell to reconstructed mesh cell

                Info<< "Writing cellProcAddressing to "
                    << databases[proci].caseName()
                      /procMesh.facesInstance()
                      /polyMesh::meshSubDir
                    << endl;

                labelIOList
                (
                    IOobject
                    (
                        "cellProcAddressing",
                        procMesh.facesInstance(),
                        polyMesh::meshSubDir,
                        procMesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false                       // Do not register
                    ),
                    cellProcAddressing[proci]
                ).write();



                // From processor patch to reconstructed mesh patch

                Info<< "Writing boundaryProcAddressing to "
                    << databases[proci].caseName()
                      /procMesh.facesInstance()
                      /polyMesh::meshSubDir
                    << endl;

                labelIOList
                (
                    IOobject
                    (
                        "boundaryProcAddressing",
                        procMesh.facesInstance(),
                        polyMesh::meshSubDir,
                        procMesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false                       // Do not register
                    ),
                    boundaryProcAddressing[proci]
                ).write();

                Info<< endl;
            }
        }
    }


    Info<< "End.\n" << endl;

    return 0;
}


// ************************************************************************* //
