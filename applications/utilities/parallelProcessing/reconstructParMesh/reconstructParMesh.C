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
#include "zeroGradientFvPatchFields.H"

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
// procBoundaries for procI. Assumes that masterMesh contains already merged
// all the processors < procI.
autoPtr<faceCoupleInfo> determineCoupledFaces
(
    const bool fullMatch,
    const label procI,
    const polyMesh& masterMesh,
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
                mergeDist,      // absolute merging distance
                true            // matching faces identical
            )
        );
    }
    else
    {
        // Pick up all patches on masterMesh ending in "toDDD" where DDD is
        // the processor number procI.

        const polyBoundaryMesh& masterPatches = masterMesh.boundaryMesh();

        const string toProcString("to" + name(procI));

        DynamicList<label> masterFaces
        (
            masterMesh.nFaces()
          - masterMesh.nInternalFaces()
        );

        forAll(masterPatches, patchI)
        {
            const polyPatch& pp = masterPatches[patchI];

            if
            (
                isA<processorPolyPatch>(pp)
             && (
                    pp.name().rfind(toProcString)
                 == (pp.name().size()-toProcString.size())
                )
            )
            {
                label meshFaceI = pp.start();
                forAll(pp, i)
                {
                    masterFaces.append(meshFaceI++);
                }
            }
        }
        masterFaces.shrink();


        // Pick up all patches on meshToAdd ending in "procBoundaryDDDtoYYY"
        // where DDD is the processor number procI and YYY is < procI.

        const polyBoundaryMesh& addPatches = meshToAdd.boundaryMesh();

        DynamicList<label> addFaces
        (
            meshToAdd.nFaces()
          - meshToAdd.nInternalFaces()
        );

        forAll(addPatches, patchI)
        {
            const polyPatch& pp = addPatches[patchI];

            if (isA<processorPolyPatch>(pp))
            {
                bool isConnected = false;

                for (label mergedProcI = 0; mergedProcI < procI; mergedProcI++)
                {
                    const string fromProcString
                    (
                        "procBoundary"
                      + name(procI)
                      + "to"
                      + name(mergedProcI)
                    );

                    if (pp.name() == fromProcString)
                    {
                        isConnected = true;
                        break;
                    }
                }

                if (isConnected)
                {
                    label meshFaceI = pp.start();
                    forAll(pp, i)
                    {
                        addFaces.append(meshFaceI++);
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
                mergeDist,      // absolute merging distance
                true,           // matching faces identical?
                false,          // if perfectmatch are faces already ordered
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
        return autoPtr<mapPolyMesh>(NULL);
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
    forAll(pointProcAddressing, procI)
    {
        labelList& constructMap = pointProcAddressing[procI];

        forAll(constructMap, i)
        {
            label oldPointI = constructMap[i];

            // New label of point after changeMesh.
            label newPointI = map().reversePointMap()[oldPointI];

            if (newPointI < -1)
            {
                constructMap[i] = -newPointI-2;
            }
            else if (newPointI >= 0)
            {
                constructMap[i] = newPointI;
            }
            else
            {
                FatalErrorIn("fvMeshDistribute::mergeSharedPoints()")
                    << "Problem. oldPointI:" << oldPointI
                    << " newPointI:" << newPointI << abort(FatalError);
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

    forAll(databases, procI)
    {
        fileName pointsInstance
        (
            databases[procI].findInstance
            (
                regionDir/polyMesh::meshSubDir,
                "points"
            )
        );

        if (pointsInstance != databases[procI].timeName())
        {
            FatalErrorIn(args.executable())
                << "Your time was specified as " << databases[procI].timeName()
                << " but there is no polyMesh/points in that time." << endl
                << "(there is a points file in " << pointsInstance
                << ")" << endl
                << "Please rerun with the correct time specified"
                << " (through the -constant, -time or -latestTime "
                << "(at your option)."
                << endl << exit(FatalError);
        }

        Info<< "Reading points from "
            << databases[procI].caseName()
            << " for time = " << databases[procI].timeName()
            << nl << endl;

        pointIOField points
        (
            IOobject
            (
                "points",
                databases[procI].findInstance
                (
                    regionDir/polyMesh::meshSubDir,
                    "points"
                ),
                regionDir/polyMesh::meshSubDir,
                databases[procI],
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

    forAll(cellProcAddressing, procI)
    {
        const labelList& pCells = cellProcAddressing[procI];
        UIndirectList<label>(cellDecomposition, pCells) = procI;
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
            zeroGradientFvPatchScalarField::typeName
        );

        forAll(cellDecomposition, cellI)
        {
            cellDist[cellI] = cellDecomposition[cellI];
        }

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

    #include "addTimeOptions.H"
    #include "addRegionOption.H"
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


    word regionName = polyMesh::defaultRegion;
    word regionDir = word::null;

    if
    (
        args.optionReadIfPresent("region", regionName)
     && regionName != polyMesh::defaultRegion
    )
    {
        regionDir = regionName;
        Info<< "Operating on region " << regionName << nl << endl;
    }

    scalar mergeTol = defaultMergeTol;
    args.optionReadIfPresent("mergeTol", mergeTol);

    scalar writeTol = Foam::pow(10.0, -scalar(IOstream::defaultPrecision()));

    Info<< "Merge tolerance : " << mergeTol << nl
        << "Write tolerance : " << writeTol << endl;

    if (runTime.writeFormat() == IOstream::ASCII && mergeTol < writeTol)
    {
        FatalErrorIn(args.executable())
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


    int nProcs = 0;

    while
    (
        isDir
        (
            args.rootPath()
          / args.caseName()
          / fileName(word("processor") + name(nProcs))
        )
    )
    {
        nProcs++;
    }

    Info<< "Found " << nProcs << " processor directories" << nl << endl;


    // Read all time databases
    PtrList<Time> databases(nProcs);

    forAll(databases, procI)
    {
        Info<< "Reading database "
            << args.caseName()/fileName(word("processor") + name(procI))
            << endl;

        databases.set
        (
            procI,
            new Time
            (
                Time::controlDictName,
                args.rootPath(),
                args.caseName()/fileName(word("processor") + name(procI))
            )
        );
    }


    // use the times list from the master processor
    // and select a subset based on the command-line options
    instantList Times = timeSelector::select
    (
        databases[0].times(),
        args
    );

    // set startTime and endTime depending on -time and -latestTime options
    #include "checkTimeOptions.H"

    if (Times.empty())
    {
        FatalErrorIn(args.executable())
            << "No times selected"
            << exit(FatalError);
    }


    // Loop over all times
    for (label timeI = startTime; timeI < endTime; timeI++)
    {
        // Set time for global database
        runTime.setTime(Times[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl << endl;

        // Set time for all databases
        forAll(databases, procI)
        {
            databases[procI].setTime(Times[timeI], timeI);
        }

        const fileName meshPath =
            databases[0].path()
           /databases[0].timeName()
           /regionDir
           /polyMesh::meshSubDir;

        if (!isFile(meshPath/"faces"))
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
            Info<< "Constructing empty mesh to add to." << nl << endl;
            fvMesh masterMesh
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
            );

            for (label procI = 0; procI < nProcs; procI++)
            {
                Info<< "Reading mesh to add from "
                    << databases[procI].caseName()
                    << " for time = " << databases[procI].timeName()
                    << nl << endl;

                fvMesh meshToAdd
                (
                    IOobject
                    (
                        regionName,
                        databases[procI].timeName(),
                        databases[procI]
                    )
                );

                // Initialize its addressing
                cellProcAddressing[procI] = identity(meshToAdd.nCells());
                faceProcAddressing[procI] = identity(meshToAdd.nFaces());
                pointProcAddressing[procI] = identity(meshToAdd.nPoints());
                boundaryProcAddressing[procI] =
                    identity(meshToAdd.boundaryMesh().size());


                // Find geometrically shared points/faces.
                autoPtr<faceCoupleInfo> couples = determineCoupledFaces
                (
                    fullMatch,
                    procI,
                    masterMesh,
                    meshToAdd,
                    mergeDist
                );


                // Add elements to mesh
                Info<< "Adding to master mesh" << nl << endl;

                autoPtr<mapAddedPolyMesh> map = fvMeshAdder::add
                (
                    masterMesh,
                    meshToAdd,
                    couples
                );

                // Update all addressing so xxProcAddressing points to correct
                // item in masterMesh.

                // Processors that were already in masterMesh
                for (label mergedI = 0; mergedI < procI; mergedI++)
                {
                    renumber(map().oldCellMap(), cellProcAddressing[mergedI]);
                    renumber(map().oldFaceMap(), faceProcAddressing[mergedI]);
                    renumber(map().oldPointMap(), pointProcAddressing[mergedI]);
                    // Note: boundary is special since can contain -1.
                    renumber
                    (
                        map().oldPatchMap(),
                        boundaryProcAddressing[mergedI]
                    );
                }

                // Added processor
                renumber(map().addedCellMap(), cellProcAddressing[procI]);
                renumber(map().addedFaceMap(), faceProcAddressing[procI]);
                renumber(map().addedPointMap(), pointProcAddressing[procI]);
                renumber(map().addedPatchMap(), boundaryProcAddressing[procI]);

                Info<< endl;
            }

            // See if any points on the mastermesh have become connected
            // because of connections through processor meshes.
            mergeSharedPoints(mergeDist, masterMesh, pointProcAddressing);

            // Save some properties on the reconstructed mesh
            masterInternalFaces = masterMesh.nInternalFaces();
            masterOwner = masterMesh.faceOwner();


            Info<< "\nWriting merged mesh to "
                << runTime.path()/runTime.timeName()
                << nl << endl;

            if (!masterMesh.write())
            {
                FatalErrorIn(args.executable())
                    << "Failed writing polyMesh."
                    << exit(FatalError);
            }

            if (writeCellDist)
            {
                writeCellDistance(runTime, masterMesh, cellProcAddressing);
            }
        }


        // Write the addressing

        Info<< "Reconstructing the addressing from the processor meshes"
            << " to the newly reconstructed mesh" << nl << endl;

        forAll(databases, procI)
        {
            Info<< "Reading processor " << procI << " mesh from "
                << databases[procI].caseName() << endl;

            polyMesh procMesh
            (
                IOobject
                (
                    regionName,
                    databases[procI].timeName(),
                    databases[procI]
                )
            );


            // From processor point to reconstructed mesh point

            Info<< "Writing pointProcAddressing to "
                << databases[procI].caseName()
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
                    false                       // do not register
                ),
                pointProcAddressing[procI]
            ).write();


            // From processor face to reconstructed mesh face

            Info<< "Writing faceProcAddressing to "
                << databases[procI].caseName()
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
                    false                       // do not register
                ),
                faceProcAddressing[procI]
            );

            // Now add turning index to faceProcAddressing.
            // See reconstructPar for meaning of turning index.
            forAll(faceProcAddr, procFaceI)
            {
                label masterFaceI = faceProcAddr[procFaceI];

                if
                (
                   !procMesh.isInternalFace(procFaceI)
                 && masterFaceI < masterInternalFaces
                )
                {
                    // proc face is now external but used to be internal face.
                    // Check if we have owner or neighbour.

                    label procOwn = procMesh.faceOwner()[procFaceI];
                    label masterOwn = masterOwner[masterFaceI];

                    if (cellProcAddressing[procI][procOwn] == masterOwn)
                    {
                        // No turning. Offset by 1.
                        faceProcAddr[procFaceI]++;
                    }
                    else
                    {
                        // Turned face.
                        faceProcAddr[procFaceI] =
                            -1 - faceProcAddr[procFaceI];
                    }
                }
                else
                {
                    // No turning. Offset by 1.
                    faceProcAddr[procFaceI]++;
                }
            }

            faceProcAddr.write();


            // From processor cell to reconstructed mesh cell

            Info<< "Writing cellProcAddressing to "
                << databases[procI].caseName()
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
                    false                       // do not register
                ),
                cellProcAddressing[procI]
            ).write();



            // From processor patch to reconstructed mesh patch

            Info<< "Writing boundaryProcAddressing to "
                << databases[procI].caseName()
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
                    false                       // do not register
                ),
                boundaryProcAddressing[procI]
            ).write();

            Info<< endl;
        }
    }


    Info<< "End.\n" << endl;

    return 0;
}


// ************************************************************************* //
