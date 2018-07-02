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

\*---------------------------------------------------------------------------*/

#include "backgroundMeshDecomposition.H"
#include "meshSearch.H"
#include "conformationSurfaces.H"
#include "zeroGradientFvPatchFields.H"
#include "Time.H"
#include "Random.H"
#include "pointConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(backgroundMeshDecomposition, 0);

}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::autoPtr<Foam::mapDistribute> Foam::backgroundMeshDecomposition::buildMap
(
    const List<label>& toProc
)
{
    // Determine send map
    // ~~~~~~~~~~~~~~~~~~

    // 1. Count
    labelList nSend(Pstream::nProcs(), 0);

    forAll(toProc, i)
    {
        label proci = toProc[i];

        nSend[proci]++;
    }


    // 2. Size sendMap
    labelListList sendMap(Pstream::nProcs());

    forAll(nSend, proci)
    {
        sendMap[proci].setSize(nSend[proci]);

        nSend[proci] = 0;
    }

    // 3. Fill sendMap
    forAll(toProc, i)
    {
        label proci = toProc[i];

        sendMap[proci][nSend[proci]++] = i;
    }

    // 4. Send over how many I need to receive
    labelList recvSizes;
    Pstream::exchangeSizes(sendMap, recvSizes);


    // Determine receive map
    // ~~~~~~~~~~~~~~~~~~~~~

    labelListList constructMap(Pstream::nProcs());

    // Local transfers first
    constructMap[Pstream::myProcNo()] = identity
    (
        sendMap[Pstream::myProcNo()].size()
    );

    label constructSize = constructMap[Pstream::myProcNo()].size();

    forAll(constructMap, proci)
    {
        if (proci != Pstream::myProcNo())
        {
            label nRecv = recvSizes[proci];

            constructMap[proci].setSize(nRecv);

            for (label i = 0; i < nRecv; i++)
            {
                constructMap[proci][i] = constructSize++;
            }
        }
    }

    return autoPtr<mapDistribute>
    (
        new mapDistribute
        (
            constructSize,
            sendMap.xfer(),
            constructMap.xfer()
        )
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::backgroundMeshDecomposition::initialRefinement()
{
    volScalarField cellWeights
    (
        IOobject
        (
            "cellWeights",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("one", dimless, 1.0),
        zeroGradientFvPatchScalarField::typeName
    );

    const conformationSurfaces& geometry = geometryToConformTo_;

    decompositionMethod& decomposer = decomposerPtr_();

    volScalarField::Internal& icellWeights = cellWeights;

    // For each cell in the mesh has it been determined if it is fully
    // inside, outside, or overlaps the surface
    List<volumeType> volumeStatus
    (
        mesh_.nCells(),
        volumeType::UNKNOWN
    );

    // Surface refinement
    {
        while (true)
        {
            // Determine/update the status of each cell
            forAll(volumeStatus, celli)
            {
                if (volumeStatus[celli] == volumeType::UNKNOWN)
                {
                    treeBoundBox cellBb
                    (
                        mesh_.cells()[celli].points
                        (
                            mesh_.faces(),
                            mesh_.points()
                        )
                    );

                    if (geometry.overlaps(cellBb))
                    {
                        volumeStatus[celli] = volumeType::MIXED;
                    }
                    else if (geometry.inside(cellBb.midpoint()))
                    {
                        volumeStatus[celli] = volumeType::INSIDE;
                    }
                    else
                    {
                        volumeStatus[celli] = volumeType::OUTSIDE;
                    }
                }
            }

            {
                labelList refCells = selectRefinementCells
                (
                    volumeStatus,
                    cellWeights
                );

                // Maintain 2:1 ratio
                labelList newCellsToRefine
                (
                    meshCutter_.consistentRefinement
                    (
                        refCells,
                        true                  // extend set
                    )
                );

                forAll(newCellsToRefine, nCTRI)
                {
                    label celli = newCellsToRefine[nCTRI];

                    if (volumeStatus[celli] == volumeType::MIXED)
                    {
                        volumeStatus[celli] = volumeType::UNKNOWN;
                    }

                    icellWeights[celli] = max
                    (
                        1.0,
                        icellWeights[celli]/8.0
                    );
                }

                if (returnReduce(newCellsToRefine.size(), sumOp<label>()) == 0)
                {
                    break;
                }

                // Mesh changing engine.
                polyTopoChange meshMod(mesh_);

                // Play refinement commands into mesh changer.
                meshCutter_.setRefinement(newCellsToRefine, meshMod);

                // Create mesh, return map from old to new mesh.
                autoPtr<mapPolyMesh> map = meshMod.changeMesh
                (
                    mesh_,
                    false,  // inflate
                    true,   // syncParallel
                    true,   // orderCells (to reduce cell transfers)
                    false   // orderPoints
                );

                // Update fields
                mesh_.updateMesh(map);

                // Update numbering of cells/vertices.
                meshCutter_.updateMesh(map);

                {
                    // Map volumeStatus

                    const labelList& cellMap = map().cellMap();

                    List<volumeType> newVolumeStatus(cellMap.size());

                    forAll(cellMap, newCelli)
                    {
                        label oldCelli = cellMap[newCelli];

                        if (oldCelli == -1)
                        {
                            newVolumeStatus[newCelli] = volumeType::UNKNOWN;
                        }
                        else
                        {
                            newVolumeStatus[newCelli] = volumeStatus[oldCelli];
                        }
                    }

                    volumeStatus.transfer(newVolumeStatus);
                }

                Info<< "    Background mesh refined from "
                    << returnReduce(map().nOldCells(), sumOp<label>())
                    << " to " << mesh_.globalData().nTotalCells()
                    << " cells." << endl;
            }

            // Determine/update the status of each cell
            forAll(volumeStatus, celli)
            {
                if (volumeStatus[celli] == volumeType::UNKNOWN)
                {
                    treeBoundBox cellBb
                    (
                        mesh_.cells()[celli].points
                        (
                            mesh_.faces(),
                            mesh_.points()
                        )
                    );

                    if (geometry.overlaps(cellBb))
                    {
                        volumeStatus[celli] = volumeType::MIXED;
                    }
                    else if (geometry.inside(cellBb.midpoint()))
                    {
                        volumeStatus[celli] = volumeType::INSIDE;
                    }
                    else
                    {
                        volumeStatus[celli] = volumeType::OUTSIDE;
                    }
                }
            }

            // Hard code switch for this stage for testing
            bool removeOutsideCells = false;

            if (removeOutsideCells)
            {
                DynamicList<label> cellsToRemove;

                forAll(volumeStatus, celli)
                {
                    if (volumeStatus[celli] == volumeType::OUTSIDE)
                    {
                        cellsToRemove.append(celli);
                    }
                }

                removeCells cellRemover(mesh_);

                // Mesh changing engine.
                polyTopoChange meshMod(mesh_);

                labelList exposedFaces = cellRemover.getExposedFaces
                (
                    cellsToRemove
                );

                // Play refinement commands into mesh changer.
                cellRemover.setRefinement
                (
                    cellsToRemove,
                    exposedFaces,
                    labelList(exposedFaces.size(), 0),  // patchID dummy
                    meshMod
                );

                // Create mesh, return map from old to new mesh.
                autoPtr<mapPolyMesh> map = meshMod.changeMesh
                (
                    mesh_,
                    false,  // inflate
                    true,   // syncParallel
                    true,   // orderCells (to reduce cell transfers)
                    false   // orderPoints
                );

                // Update fields
                mesh_.updateMesh(map);

                // Update numbering of cells/vertices.
                meshCutter_.updateMesh(map);
                cellRemover.updateMesh(map);

                {
                    // Map volumeStatus

                    const labelList& cellMap = map().cellMap();

                    List<volumeType> newVolumeStatus(cellMap.size());

                    forAll(cellMap, newCelli)
                    {
                        label oldCelli = cellMap[newCelli];

                        if (oldCelli == -1)
                        {
                            newVolumeStatus[newCelli] = volumeType::UNKNOWN;
                        }
                        else
                        {
                            newVolumeStatus[newCelli] = volumeStatus[oldCelli];
                        }
                    }

                    volumeStatus.transfer(newVolumeStatus);
                }

                Info<< "Removed "
                    << returnReduce(map().nOldCells(), sumOp<label>())
                     - mesh_.globalData().nTotalCells()
                    << " cells." << endl;
            }

            if (debug)
            {
                // const_cast<Time&>(mesh_.time())++;
                // Info<< "Time " << mesh_.time().timeName() << endl;
                meshCutter_.write();
                mesh_.write();
                cellWeights.write();
            }

            labelList newDecomp = decomposer.decompose
            (
                mesh_,
                mesh_.cellCentres(),
                icellWeights
            );

            fvMeshDistribute distributor(mesh_, mergeDist_);

            autoPtr<mapDistributePolyMesh> mapDist = distributor.distribute
            (
                newDecomp
            );

            meshCutter_.distribute(mapDist);

            mapDist().distributeCellData(volumeStatus);

            if (debug)
            {
                printMeshData(mesh_);

                // const_cast<Time&>(mesh_.time())++;
                // Info<< "Time " << mesh_.time().timeName() << endl;
                meshCutter_.write();
                mesh_.write();
                cellWeights.write();
            }
        }
    }

    if (debug)
    {
        // const_cast<Time&>(mesh_.time())++;
        // Info<< "Time " << mesh_.time().timeName() << endl;
        cellWeights.write();
        mesh_.write();
    }

    buildPatchAndTree();
}


void Foam::backgroundMeshDecomposition::printMeshData
(
    const polyMesh& mesh
) const
{
    // Collect all data on master

    globalIndex globalCells(mesh.nCells());
    // labelListList patchNeiProcNo(Pstream::nProcs());
    // labelListList patchSize(Pstream::nProcs());
    // const labelList& pPatches = mesh.globalData().processorPatches();
    // patchNeiProcNo[Pstream::myProcNo()].setSize(pPatches.size());
    // patchSize[Pstream::myProcNo()].setSize(pPatches.size());
    // forAll(pPatches, i)
    // {
    //     const processorPolyPatch& ppp = refCast<const processorPolyPatch>
    //     (
    //         mesh.boundaryMesh()[pPatches[i]]
    //     );
    //     patchNeiProcNo[Pstream::myProcNo()][i] = ppp.neighbProcNo();
    //     patchSize[Pstream::myProcNo()][i] = ppp.size();
    // }
    // Pstream::gatherList(patchNeiProcNo);
    // Pstream::gatherList(patchSize);


    // // Print stats

    // globalIndex globalBoundaryFaces(mesh.nFaces()-mesh.nInternalFaces());

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        Info<< "Processor " << proci << " "
            << "Number of cells = " << globalCells.localSize(proci)
            << endl;

        // label nProcFaces = 0;

        // const labelList& nei = patchNeiProcNo[proci];

        // forAll(patchNeiProcNo[proci], i)
        // {
        //     Info<< "    Number of faces shared with processor "
        //         << patchNeiProcNo[proci][i] << " = " << patchSize[proci][i]
        //         << endl;

        //     nProcFaces += patchSize[proci][i];
        // }

        // Info<< "    Number of processor patches = " << nei.size() << nl
        //     << "    Number of processor faces = " << nProcFaces << nl
        //     << "    Number of boundary faces = "
        //     << globalBoundaryFaces.localSize(proci) << endl;
    }
}


bool Foam::backgroundMeshDecomposition::refineCell
(
    label celli,
    volumeType volType,
    scalar& weightEstimate
) const
{
    // Sample the box to find an estimate of the min size, and a volume
    // estimate when overlapping == true.

//    const conformationSurfaces& geometry = geometryToConformTo_;

    treeBoundBox cellBb
    (
        mesh_.cells()[celli].points
        (
            mesh_.faces(),
            mesh_.points()
        )
    );

    weightEstimate = 1.0;

    if (volType == volumeType::MIXED)
    {
//        // Assess the cell size at the nearest point on the surface for the
//        // MIXED cells, if the cell is large with respect to the cell size,
//        // then refine it.
//
//        pointField samplePoints
//        (
//            volRes_*volRes_*volRes_,
//            Zero
//        );
//
//        // scalar sampleVol = cellBb.volume()/samplePoints.size();
//
//        vector delta = cellBb.span()/volRes_;
//
//        label pI = 0;
//
//        for (label i = 0; i < volRes_; i++)
//        {
//            for (label j = 0; j < volRes_; j++)
//            {
//                for (label k = 0; k < volRes_; k++)
//                {
//                    samplePoints[pI++] =
//                        cellBb.min()
//                      + vector
//                        (
//                            delta.x()*(i + 0.5),
//                            delta.y()*(j + 0.5),
//                            delta.z()*(k + 0.5)
//                        );
//                }
//            }
//        }
//
//        List<pointIndexHit> hitInfo;
//        labelList hitSurfaces;
//
//        geometry.findSurfaceNearest
//        (
//            samplePoints,
//            scalarField(samplePoints.size(), sqr(great)),
//            hitInfo,
//            hitSurfaces
//        );
//
//        // weightEstimate = 0.0;
//
//        scalar minCellSize = great;
//
//        forAll(samplePoints, i)
//        {
//            scalar s = cellShapeControl_.cellSize
//            (
//                hitInfo[i].hitPoint()
//            );
//
//            // Info<< "cellBb.midpoint() " << cellBb.midpoint() << nl
//            //     << samplePoints[i] << nl
//            //     << hitInfo[i] << nl
//            //     << "cellBb.span() " << cellBb.span() << nl
//            //     << "cellBb.mag() " << cellBb.mag() << nl
//            //     << s << endl;
//
//            if (s < minCellSize)
//            {
//                minCellSize = max(s, minCellSizeLimit_);
//            }
//
//            // Estimate the number of points in the cell by the surface size,
//            // this is likely to be too small, so reduce.
//            // weightEstimate += sampleVol/pow3(s);
//        }
//
//        if (sqr(spanScale_)*sqr(minCellSize) < magSqr(cellBb.span()))
//        {
//            return true;
//        }
    }
    else if (volType == volumeType::INSIDE)
    {
        // scalar s =
        //    foamyHexMesh_.cellShapeControl_.cellSize(cellBb.midpoint());

        // Estimate the number of points in the cell by the size at the cell
        // midpoint
        // weightEstimate = cellBb.volume()/pow3(s);

        return false;
    }
    // else
    // {
    //     weightEstimate = 1.0;

    //     return false;
    // }

    return false;
}


Foam::labelList Foam::backgroundMeshDecomposition::selectRefinementCells
(
    List<volumeType>& volumeStatus,
    volScalarField& cellWeights
) const
{
    volScalarField::Internal& icellWeights = cellWeights;

    labelHashSet cellsToRefine;

    // Determine/update the status of each cell
    forAll(volumeStatus, celli)
    {
        if (volumeStatus[celli] == volumeType::MIXED)
        {
            if (meshCutter_.cellLevel()[celli] < minLevels_)
            {
                cellsToRefine.insert(celli);
            }
        }

        if (volumeStatus[celli] != volumeType::OUTSIDE)
        {
            if
            (
                refineCell
                (
                    celli,
                    volumeStatus[celli],
                    icellWeights[celli]
                )
            )
            {
                cellsToRefine.insert(celli);
            }
        }
    }

    return cellsToRefine.toc();
}


void Foam::backgroundMeshDecomposition::buildPatchAndTree()
{
    primitivePatch tmpBoundaryFaces
    (
        SubList<face>
        (
            mesh_.faces(),
            mesh_.nFaces() - mesh_.nInternalFaces(),
            mesh_.nInternalFaces()
        ),
        mesh_.points()
    );

    boundaryFacesPtr_.reset
    (
        new bPatch
        (
            tmpBoundaryFaces.localFaces(),
            tmpBoundaryFaces.localPoints()
        )
    );

    // Overall bb
    treeBoundBox overallBb(boundaryFacesPtr_().localPoints());

    bFTreePtr_.reset
    (
        new indexedOctree<treeDataBPatch>
        (
            treeDataBPatch
            (
                false,
                boundaryFacesPtr_(),
                indexedOctree<treeDataBPatch>::perturbTol()
            ),
            overallBb.extend(1e-4),
            10, // maxLevel
            10, // leafSize
            3.0 // duplicity
        )
    );

    // Give the bounds of every processor to every other processor
    allBackgroundMeshBounds_[Pstream::myProcNo()] = overallBb;

    Pstream::gatherList(allBackgroundMeshBounds_);
    Pstream::scatterList(allBackgroundMeshBounds_);

    point bbMin(great, great, great);
    point bbMax(-great, -great, -great);

    forAll(allBackgroundMeshBounds_, proci)
    {
        bbMin = min(bbMin, allBackgroundMeshBounds_[proci].min());
        bbMax = max(bbMax, allBackgroundMeshBounds_[proci].max());
    }

    globalBackgroundBounds_ = treeBoundBox(bbMin, bbMax);

    if (false)
    {
        OFstream fStr
        (
            mesh_.time().path()
           /"backgroundMeshDecomposition_proc_"
          + name(Pstream::myProcNo())
          + "_boundaryFaces.obj"
        );

        const faceList& faces = boundaryFacesPtr_().localFaces();
        const List<point>& points = boundaryFacesPtr_().localPoints();

        Map<label> foamToObj(points.size());

        label vertI = 0;

        forAll(faces, i)
        {
            const face& f = faces[i];

            forAll(f, fPI)
            {
                if (foamToObj.insert(f[fPI], vertI))
                {
                    meshTools::writeOBJ(fStr, points[f[fPI]]);
                    vertI++;
                }
            }

            fStr<< 'f';

            forAll(f, fPI)
            {
                fStr<< ' ' << foamToObj[f[fPI]] + 1;
            }

            fStr<< nl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::backgroundMeshDecomposition::backgroundMeshDecomposition
(
    const Time& runTime,
    Random& rndGen,
    const conformationSurfaces& geometryToConformTo,
    const dictionary& coeffsDict
)
:
    runTime_(runTime),
    geometryToConformTo_(geometryToConformTo),
    rndGen_(rndGen),
    mesh_
    (
        IOobject
        (
            "backgroundMeshDecomposition",
            runTime_.timeName(),
            runTime_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            false
        )
    ),
    meshCutter_
    (
        mesh_,
        labelList(mesh_.nCells(), 0),
        labelList(mesh_.nPoints(), 0)
    ),
    boundaryFacesPtr_(),
    bFTreePtr_(),
    allBackgroundMeshBounds_(Pstream::nProcs()),
    globalBackgroundBounds_(),
    decomposeDict_
    (
        IOobject
        (
            "decomposeParDict",
            runTime_.system(),
            runTime_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    decomposerPtr_(decompositionMethod::New(decomposeDict_)),
    mergeDist_(1e-6*mesh_.bounds().mag()),
    spanScale_(readScalar(coeffsDict.lookup("spanScale"))),
    minCellSizeLimit_
    (
        coeffsDict.lookupOrDefault<scalar>("minCellSizeLimit", 0.0)
    ),
    minLevels_(readLabel(coeffsDict.lookup("minLevels"))),
    volRes_(readLabel(coeffsDict.lookup("sampleResolution"))),
    maxCellWeightCoeff_(readScalar(coeffsDict.lookup("maxCellWeightCoeff")))
{
    if (!Pstream::parRun())
    {
        FatalErrorInFunction
            << "This cannot be used when not running in parallel."
            << exit(FatalError);
    }

    if (!decomposerPtr_().parallelAware())
    {
        FatalErrorInFunction
            << "You have selected decomposition method "
            << decomposerPtr_().typeName
            << " which is not parallel aware." << endl
            << exit(FatalError);
    }

    Info<< nl << "Building initial background mesh decomposition" << endl;

    initialRefinement();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::backgroundMeshDecomposition::~backgroundMeshDecomposition()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mapDistributePolyMesh>
Foam::backgroundMeshDecomposition::distribute
(
    volScalarField& cellWeights
)
{
    if (debug)
    {
        // const_cast<Time&>(mesh_.time())++;
        // Info<< "Time " << mesh_.time().timeName() << endl;
        cellWeights.write();
        mesh_.write();
    }

    volScalarField::Internal& icellWeights = cellWeights;

    while (true)
    {
        // Refine large cells if necessary

        label nOccupiedCells = 0;

        forAll(icellWeights, cI)
        {
            if (icellWeights[cI] > 1 - small)
            {
                nOccupiedCells++;
            }
        }

        // Only look at occupied cells, as there is a possibility of runaway
        // refinement if the number of cells grows too fast.  Also, clip the
        // minimum cellWeightLimit at maxCellWeightCoeff_

        scalar cellWeightLimit = max
        (
            maxCellWeightCoeff_
           *sum(cellWeights).value()
           /returnReduce(nOccupiedCells, sumOp<label>()),
            maxCellWeightCoeff_
        );

        if (debug)
        {
            Info<< "    cellWeightLimit " << cellWeightLimit << endl;

            Pout<< "    sum(cellWeights) " << sum(cellWeights.primitiveField())
                << " max(cellWeights) " << max(cellWeights.primitiveField())
                << endl;
        }

        labelHashSet cellsToRefine;

        forAll(icellWeights, cWI)
        {
            if (icellWeights[cWI] > cellWeightLimit)
            {
                cellsToRefine.insert(cWI);
            }
        }

        if (returnReduce(cellsToRefine.size(), sumOp<label>()) == 0)
        {
            break;
        }

        // Maintain 2:1 ratio
        labelList newCellsToRefine
        (
            meshCutter_.consistentRefinement
            (
                cellsToRefine.toc(),
                true                  // extend set
            )
        );

        if (debug && !cellsToRefine.empty())
        {
            Pout<< "    cellWeights too large in " << cellsToRefine.size()
                << " cells" << endl;
        }

        forAll(newCellsToRefine, nCTRI)
        {
            label celli = newCellsToRefine[nCTRI];

            icellWeights[celli] /= 8.0;
        }

        // Mesh changing engine.
        polyTopoChange meshMod(mesh_);

        // Play refinement commands into mesh changer.
        meshCutter_.setRefinement(newCellsToRefine, meshMod);

        // Create mesh, return map from old to new mesh.
        autoPtr<mapPolyMesh> map = meshMod.changeMesh
        (
            mesh_,
            false,  // inflate
            true,   // syncParallel
            true,   // orderCells (to reduce cell motion)
            false   // orderPoints
        );

        // Update fields
        mesh_.updateMesh(map);

        // Update numbering of cells/vertices.
        meshCutter_.updateMesh(map);

        Info<< "    Background mesh refined from "
            << returnReduce(map().nOldCells(), sumOp<label>())
            << " to " << mesh_.globalData().nTotalCells()
            << " cells." << endl;

        if (debug)
        {
            // const_cast<Time&>(mesh_.time())++;
            // Info<< "Time " << mesh_.time().timeName() << endl;
            cellWeights.write();
            mesh_.write();
        }
    }

    if (debug)
    {
        printMeshData(mesh_);

        Pout<< "    Pre distribute sum(cellWeights) "
            << sum(icellWeights)
            << " max(cellWeights) "
            << max(icellWeights)
            << endl;
    }

    labelList newDecomp = decomposerPtr_().decompose
    (
        mesh_,
        mesh_.cellCentres(),
        icellWeights
    );

    Info<< "    Redistributing background mesh cells" << endl;

    fvMeshDistribute distributor(mesh_, mergeDist_);

    autoPtr<mapDistributePolyMesh> mapDist = distributor.distribute(newDecomp);

    meshCutter_.distribute(mapDist);

    if (debug)
    {
        printMeshData(mesh_);

        Pout<< "    Post distribute sum(cellWeights) "
            << sum(icellWeights)
            << " max(cellWeights) "
            << max(icellWeights)
            << endl;

        // const_cast<Time&>(mesh_.time())++;
        // Info<< "Time " << mesh_.time().timeName() << endl;
        mesh_.write();
        cellWeights.write();
    }

    buildPatchAndTree();

    return mapDist;
}


bool Foam::backgroundMeshDecomposition::positionOnThisProcessor
(
    const point& pt
) const
{
//    return bFTreePtr_().findAnyOverlap(pt, 0.0);

    return (bFTreePtr_().getVolumeType(pt) == volumeType::INSIDE);
}


Foam::boolList Foam::backgroundMeshDecomposition::positionOnThisProcessor
(
    const List<point>& pts
) const
{
    boolList posProc(pts.size(), true);

    forAll(pts, pI)
    {
        posProc[pI] = positionOnThisProcessor(pts[pI]);
    }

    return posProc;
}


bool Foam::backgroundMeshDecomposition::overlapsThisProcessor
(
    const treeBoundBox& box
) const
{
//    return !procBounds().contains(box);
    return !bFTreePtr_().findBox(box).empty();
}


bool Foam::backgroundMeshDecomposition::overlapsThisProcessor
(
    const point& centre,
    const scalar radiusSqr
) const
{
    // return bFTreePtr_().findAnyOverlap(centre, radiusSqr);

    return bFTreePtr_().findNearest(centre, radiusSqr).hit();
}


Foam::pointIndexHit Foam::backgroundMeshDecomposition::findLine
(
    const point& start,
    const point& end
) const
{
    return bFTreePtr_().findLine(start, end);
}


Foam::pointIndexHit Foam::backgroundMeshDecomposition::findLineAny
(
    const point& start,
    const point& end
) const
{
    return bFTreePtr_().findLineAny(start, end);
}


Foam::labelList Foam::backgroundMeshDecomposition::processorNearestPosition
(
    const List<point>& pts
) const
{
    DynamicList<label> toCandidateProc;
    DynamicList<point> testPoints;
    labelList ptBlockStart(pts.size(), -1);
    labelList ptBlockSize(pts.size(), -1);

    label nTotalCandidates = 0;

    forAll(pts, pI)
    {
        const point& pt = pts[pI];

        label nCandidates = 0;

        forAll(allBackgroundMeshBounds_, proci)
        {
            // Candidate points may lie just outside a processor box, increase
            // test range by using overlaps rather than contains
            if (allBackgroundMeshBounds_[proci].overlaps(pt, sqr(small*100)))
            {
                toCandidateProc.append(proci);
                testPoints.append(pt);

                nCandidates++;
            }
        }

        ptBlockStart[pI] = nTotalCandidates;
        ptBlockSize[pI] = nCandidates;

        nTotalCandidates += nCandidates;
    }

    // Needed for reverseDistribute
    label preDistributionToCandidateProcSize = toCandidateProc.size();

    autoPtr<mapDistribute> map(buildMap(toCandidateProc));

    map().distribute(testPoints);

    List<scalar> distanceSqrToCandidate(testPoints.size(), sqr(great));

    // Test candidate points on candidate processors
    forAll(testPoints, tPI)
    {
        pointIndexHit info = bFTreePtr_().findNearest
        (
            testPoints[tPI],
            sqr(great)
        );

        if (info.hit())
        {
            distanceSqrToCandidate[tPI] = magSqr
            (
                testPoints[tPI] - info.hitPoint()
            );
        }
    }

    map().reverseDistribute
    (
        preDistributionToCandidateProcSize,
        distanceSqrToCandidate
    );

    labelList ptNearestProc(pts.size(), -1);

    forAll(pts, pI)
    {
        // Extract the sub list of results for this point

        SubList<scalar> ptNearestProcResults
        (
            distanceSqrToCandidate,
            ptBlockSize[pI],
            ptBlockStart[pI]
        );

        scalar nearestProcDistSqr = great;

        forAll(ptNearestProcResults, pPRI)
        {
            if (ptNearestProcResults[pPRI] < nearestProcDistSqr)
            {
                nearestProcDistSqr = ptNearestProcResults[pPRI];

                ptNearestProc[pI] = toCandidateProc[ptBlockStart[pI] + pPRI];
            }
        }

        if (debug)
        {
            Pout<< pts[pI] << " nearestProcDistSqr " << nearestProcDistSqr
                << " ptNearestProc[pI] " << ptNearestProc[pI] << endl;
        }

        if (ptNearestProc[pI] < 0)
        {
            FatalErrorInFunction
                << "The position " << pts[pI]
                << " did not find a nearest point on the background mesh."
                << exit(FatalError);
        }
    }

    return ptNearestProc;
}



Foam::List<Foam::List<Foam::pointIndexHit>>
Foam::backgroundMeshDecomposition::intersectsProcessors
(
    const List<point>& starts,
    const List<point>& ends,
    bool includeOwnProcessor
) const
{
    DynamicList<label> toCandidateProc;
    DynamicList<point> testStarts;
    DynamicList<point> testEnds;
    labelList segmentBlockStart(starts.size(), -1);
    labelList segmentBlockSize(starts.size(), -1);

    label nTotalCandidates = 0;

    forAll(starts, sI)
    {
        const point& s = starts[sI];
        const point& e = ends[sI];

        // Dummy point for treeBoundBox::intersects
        point p(Zero);

        label nCandidates = 0;

        forAll(allBackgroundMeshBounds_, proci)
        {
            // It is assumed that the sphere in question overlaps the source
            // processor, so don't test it, unless includeOwnProcessor is true
            if
            (
                (includeOwnProcessor || proci != Pstream::myProcNo())
              && allBackgroundMeshBounds_[proci].intersects(s, e, p)
            )
            {
                toCandidateProc.append(proci);
                testStarts.append(s);
                testEnds.append(e);

                nCandidates++;
            }
        }

        segmentBlockStart[sI] = nTotalCandidates;
        segmentBlockSize[sI] = nCandidates;

        nTotalCandidates += nCandidates;
    }

    // Needed for reverseDistribute
    label preDistributionToCandidateProcSize = toCandidateProc.size();

    autoPtr<mapDistribute> map(buildMap(toCandidateProc));

    map().distribute(testStarts);
    map().distribute(testEnds);

    List<pointIndexHit> segmentIntersectsCandidate(testStarts.size());

    // Test candidate segments on candidate processors
    forAll(testStarts, sI)
    {
        const point& s = testStarts[sI];
        const point& e = testEnds[sI];

        // If the sphere finds some elements of the patch, then it overlaps
        segmentIntersectsCandidate[sI] = bFTreePtr_().findLine(s, e);
    }

    map().reverseDistribute
    (
        preDistributionToCandidateProcSize,
        segmentIntersectsCandidate
    );

    List<List<pointIndexHit>> segmentHitProcs(starts.size());

    // Working storage for assessing processors
    DynamicList<pointIndexHit> tmpProcHits;

    forAll(starts, sI)
    {
        tmpProcHits.clear();

        // Extract the sub list of results for this point

        SubList<pointIndexHit> segmentProcResults
        (
            segmentIntersectsCandidate,
            segmentBlockSize[sI],
            segmentBlockStart[sI]
        );

        forAll(segmentProcResults, sPRI)
        {
            if (segmentProcResults[sPRI].hit())
            {
                tmpProcHits.append(segmentProcResults[sPRI]);

                tmpProcHits.last().setIndex
                (
                    toCandidateProc[segmentBlockStart[sI] + sPRI]
                );
            }
        }

        segmentHitProcs[sI] = tmpProcHits;
    }

    return segmentHitProcs;
}


bool Foam::backgroundMeshDecomposition::overlapsOtherProcessors
(
    const point& centre,
    const scalar& radiusSqr
) const
{
    forAll(allBackgroundMeshBounds_, proci)
    {
        if (bFTreePtr_().findNearest(centre, radiusSqr).hit())
        {
            return true;
        }
    }

    return false;
}


Foam::labelList Foam::backgroundMeshDecomposition::overlapProcessors
(
    const point& centre,
    const scalar radiusSqr
) const
{
    DynamicList<label> toProc(Pstream::nProcs());

    forAll(allBackgroundMeshBounds_, proci)
    {
        // Test against the bounding box of the processor
        if
        (
            proci != Pstream::myProcNo()
         && allBackgroundMeshBounds_[proci].overlaps(centre, radiusSqr)
        )
        {
            // Expensive test
//            if (bFTreePtr_().findNearest(centre, radiusSqr).hit())
            {
                toProc.append(proci);
            }
        }
    }

    return toProc;
}


//Foam::labelListList Foam::backgroundMeshDecomposition::overlapsProcessors
//(
//    const List<point>& centres,
//    const List<scalar>& radiusSqrs,
//    const Delaunay& T,
//    bool includeOwnProcessor
//) const
//{
//    DynamicList<label> toCandidateProc;
//    DynamicList<point> testCentres;
//    DynamicList<scalar> testRadiusSqrs;
//    labelList sphereBlockStart(centres.size(), -1);
//    labelList sphereBlockSize(centres.size(), -1);
//
//    label nTotalCandidates = 0;
//
//    forAll(centres, sI)
//    {
//        const point& c = centres[sI];
//        scalar rSqr = radiusSqrs[sI];
//
//        label nCandidates = 0;
//
//        forAll(allBackgroundMeshBounds_, proci)
//        {
//            // It is assumed that the sphere in question overlaps the source
//            // processor, so don't test it, unless includeOwnProcessor is true
//            if
//            (
//                (includeOwnProcessor || proci != Pstream::myProcNo())
//             && allBackgroundMeshBounds_[proci].overlaps(c, rSqr)
//            )
//            {
//                if (bFTreePtr_().findNearest(c, rSqr).hit())
//                {
//                    toCandidateProc.append(proci);
//                    testCentres.append(c);
//                    testRadiusSqrs.append(rSqr);
//
//                    nCandidates++;
//                }
//            }
//        }
//
//        sphereBlockStart[sI] = nTotalCandidates;
//        sphereBlockSize[sI] = nCandidates;
//
//        nTotalCandidates += nCandidates;
//    }
//
//    // Needed for reverseDistribute
////    label preDistributionToCandidateProcSize = toCandidateProc.size();
////
////    autoPtr<mapDistribute> map(buildMap(toCandidateProc));
////
////    map().distribute(testCentres);
////    map().distribute(testRadiusSqrs);
//
//    // TODO: This is faster, but results in more vertices being referred
//    boolList sphereOverlapsCandidate(testCentres.size(), true);
////    boolList sphereOverlapsCandidate(testCentres.size(), false);
////
////    // Test candidate spheres on candidate processors
////    forAll(testCentres, sI)
////    {
////        const point& c = testCentres[sI];
////        const scalar rSqr = testRadiusSqrs[sI];
////
////        const bool flagOverlap = bFTreePtr_().findNearest(c, rSqr).hit();
////
////        if (flagOverlap)
////        {
////            // if (vertexOctree.findAnyOverlap(c, rSqr))
//////            if (vertexOctree.findNearest(c, rSqr*1.001).hit())
//////            {
//////                sphereOverlapsCandidate[sI] = true;
//////            }
////
//////            Vertex_handle nearestVertex = T.nearest_vertex
//////            (
//////                toPoint<Point>(c)
//////            );
//////
//////            const scalar distSqr = magSqr
//////            (
//////                topoint(nearestVertex->point()) - c
//////            );
//////
//////            if (distSqr <= rSqr)
//////            {
//////                // If the sphere finds a nearest element of the patch,
//////                // then it overlaps
////                sphereOverlapsCandidate[sI] = true;
//////            }
////        }
////    }
//
////    map().reverseDistribute
////    (
////        preDistributionToCandidateProcSize,
////        sphereOverlapsCandidate
////    );
//
//    labelListList sphereProcs(centres.size());
//
//    // Working storage for assessing processors
//    DynamicList<label> tmpProcs;
//
//    forAll(centres, sI)
//    {
//        tmpProcs.clear();
//
//        // Extract the sub list of results for this point
//
//        SubList<bool> sphereProcResults
//        (
//            sphereOverlapsCandidate,
//            sphereBlockSize[sI],
//            sphereBlockStart[sI]
//        );
//
//        forAll(sphereProcResults, sPRI)
//        {
//            if (sphereProcResults[sPRI])
//            {
//                tmpProcs.append(toCandidateProc[sphereBlockStart[sI] + sPRI]);
//            }
//        }
//
//        sphereProcs[sI] = tmpProcs;
//    }
//
//    return sphereProcs;
//}


//Foam::labelListList Foam::backgroundMeshDecomposition::overlapProcessors
//(
//    const point& cc,
//    const scalar rSqr
//) const
//{
//    DynamicList<label> toCandidateProc;
//    label sphereBlockStart(-1);
//    label sphereBlockSize(-1);
//
//    label nCandidates = 0;
//
//    forAll(allBackgroundMeshBounds_, proci)
//    {
//        // It is assumed that the sphere in question overlaps the source
//        // processor, so don't test it, unless includeOwnProcessor is true
//        if
//        (
//            (includeOwnProcessor || proci != Pstream::myProcNo())
//         && allBackgroundMeshBounds_[proci].overlaps(cc, rSqr)
//        )
//        {
//            toCandidateProc.append(proci);
//
//            nCandidates++;
//        }
//    }
//
//    sphereBlockSize = nCandidates;
//    nTotalCandidates += nCandidates;
//
//    // Needed for reverseDistribute
//    label preDistributionToCandidateProcSize = toCandidateProc.size();
//
//    autoPtr<mapDistribute> map(buildMap(toCandidateProc));
//
//    map().distribute(testCentres);
//    map().distribute(testRadiusSqrs);
//
//    // TODO: This is faster, but results in more vertices being referred
////    boolList sphereOverlapsCandidate(testCentres.size(), true);
//    boolList sphereOverlapsCandidate(testCentres.size(), false);
//
//    // Test candidate spheres on candidate processors
//    forAll(testCentres, sI)
//    {
//        const point& c = testCentres[sI];
//        const scalar rSqr = testRadiusSqrs[sI];
//
//        const bool flagOverlap = bFTreePtr_().findNearest(c, rSqr).hit();
//
//        if (flagOverlap)
//        {
//            // if (vertexOctree.findAnyOverlap(c, rSqr))
////            if (vertexOctree.findNearest(c, rSqr*1.001).hit())
////            {
////                sphereOverlapsCandidate[sI] = true;
////            }
//
////            Vertex_handle nearestVertex = T.nearest_vertex
////            (
////                toPoint<Point>(c)
////            );
////
////            const scalar distSqr = magSqr
////            (
////                topoint(nearestVertex->point()) - c
////            );
////
////            if (distSqr <= rSqr)
////            {
////                // If the sphere finds a nearest element of the patch, then
////                // it overlaps
//                sphereOverlapsCandidate[sI] = true;
////            }
//        }
//    }
//
//    map().reverseDistribute
//    (
//        preDistributionToCandidateProcSize,
//        sphereOverlapsCandidate
//    );
//
//    labelListList sphereProcs(centres.size());
//
//    // Working storage for assessing processors
//    DynamicList<label> tmpProcs;
//
//    forAll(centres, sI)
//    {
//        tmpProcs.clear();
//
//        // Extract the sub list of results for this point
//
//        SubList<bool> sphereProcResults
//        (
//            sphereOverlapsCandidate,
//            sphereBlockSize[sI],
//            sphereBlockStart[sI]
//        );
//
//        forAll(sphereProcResults, sPRI)
//        {
//            if (sphereProcResults[sPRI])
//            {
//                tmpProcs.append(toCandidateProc[sphereBlockStart[sI] + sPRI]);
//            }
//        }
//
//        sphereProcs[sI] = tmpProcs;
//    }
//
//    return sphereProcs;
//}


// ************************************************************************* //
