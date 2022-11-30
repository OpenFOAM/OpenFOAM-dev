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

\*---------------------------------------------------------------------------*/

#include "domainDecomposition.H"
#include "decompositionMethod.H"
#include "IOobjectList.H"
#include "cyclicFvPatch.H"
#include "processorCyclicFvPatch.H"
#include "nonConformalCyclicFvPatch.H"
#include "nonConformalProcessorCyclicFvPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::domainDecomposition::mark
(
    const labelList& zoneElems,
    const label zoneI,
    labelList& elementToZone
)
{
    forAll(zoneElems, i)
    {
        label pointi = zoneElems[i];

        if (elementToZone[pointi] == -1)
        {
            // First occurrence
            elementToZone[pointi] = zoneI;
        }
        else if (elementToZone[pointi] >= 0)
        {
            // Multiple zones
            elementToZone[pointi] = -2;
        }
    }
}


void Foam::domainDecomposition::addInterProcFace
(
    const label facei,
    const label ownerProc,
    const label nbrProc,
    const label subPatchID,
    List<Map<label>>& nbrToInterPatch,
    List<DynamicList<DynamicList<label>>>& interPatchFaces,
    List<labelListList>& subPatchIDs,
    List<labelListList>& subPatchStarts
) const
{
    // Create new interproc patches, if they do not already exist
    label toNbrProcPatchi = -1, toOwnerProcPatchi = -1;
    if (!nbrToInterPatch[ownerProc].found(nbrProc))
    {
        toNbrProcPatchi = nbrToInterPatch[ownerProc].size();
        nbrToInterPatch[ownerProc].insert(nbrProc, toNbrProcPatchi);
        interPatchFaces[ownerProc].append(DynamicList<label>());
        subPatchIDs[ownerProc].append(labelList(1, subPatchID));
        subPatchStarts[ownerProc].append(labelList(1, label(0)));

        if (facei != -1 && completeMesh().isInternalFace(facei))
        {
            toOwnerProcPatchi = nbrToInterPatch[nbrProc].size();
            nbrToInterPatch[nbrProc].insert(ownerProc, toOwnerProcPatchi);
            interPatchFaces[nbrProc].append(DynamicList<label>());
            subPatchIDs[nbrProc].append(labelList(1, subPatchID));
            subPatchStarts[nbrProc].append(labelList(1, label(0)));
        }
    }
    else
    {
        toNbrProcPatchi = nbrToInterPatch[ownerProc][nbrProc];

        if (facei != -1 && completeMesh().isInternalFace(facei))
        {
            toOwnerProcPatchi = nbrToInterPatch[nbrProc][ownerProc];
        }
    }

    // If the sub patch has changed then add new sub-patch entries
    if (subPatchIDs[ownerProc][toNbrProcPatchi].last() != subPatchID)
    {
        subPatchIDs[ownerProc][toNbrProcPatchi].append(subPatchID);
        subPatchStarts[ownerProc][toNbrProcPatchi].append
        (
            interPatchFaces[ownerProc][toNbrProcPatchi].size()
        );

        if (facei != -1 && completeMesh().isInternalFace(facei))
        {
            subPatchIDs[nbrProc][toOwnerProcPatchi].append(subPatchID);
            subPatchStarts[nbrProc][toOwnerProcPatchi].append
            (
                interPatchFaces[nbrProc][toOwnerProcPatchi].size()
            );
        }
    }

    // Add face to the inter-proc patches. Note use of turning index.
    if (facei != -1)
    {
        interPatchFaces[ownerProc][toNbrProcPatchi].append(facei + 1);

        if (completeMesh().isInternalFace(facei))
        {
            interPatchFaces[nbrProc][toOwnerProcPatchi].append(- facei - 1);
        }
    }
}


Foam::labelList Foam::domainDecomposition::distributeCells()
{
    Info<< "\nCalculating distribution of cells" << endl;

    cpuTime decompositionTime;

    const dictionary decomposeParDict =
        decompositionMethod::decomposeParDict(runTimes_.completeTime());

    scalarField cellWeights;
    if (decomposeParDict.found("weightField"))
    {
        const word weightName = decomposeParDict.lookup("weightField");

        volScalarField weights
        (
            IOobject
            (
                weightName,
                completeMesh().time().name(),
                completeMesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            completeMesh()
        );
        cellWeights = weights.primitiveField();
    }

    const labelList result =
        decompositionMethod::NewDecomposer(decomposeParDict)->decompose
        (
            completeMesh(),
            cellWeights
        );

    Info<< "\nFinished decomposition in "
        << decompositionTime.elapsedCpuTime()
        << " s" << endl;

    return result;
}


template<class BinaryOp>
inline void Foam::domainDecomposition::processInterCyclics
(
    const labelList& cellProc,
    const polyBoundaryMesh& patches,
    List<DynamicList<DynamicList<label>>>& interPatchFaces,
    List<Map<label>>& procNbrToInterPatch,
    List<labelListList>& subPatchIDs,
    List<labelListList>& subPatchStarts,
    bool owner,
    BinaryOp bop
) const
{
    // Processor boundaries from split cyclics
    forAll(patches, patchi)
    {
        if (isA<nonConformalCyclicPolyPatch>(patches[patchi]))
        {
            const nonConformalCyclicPolyPatch& nccPp =
                refCast<const nonConformalCyclicPolyPatch>(patches[patchi]);

            if (nccPp.owner() != owner) continue;

            const polyPatch& origPp = nccPp.origPatch();
            const polyPatch& nbrOrigPp = nccPp.nbrPatch().origPatch();

            const labelUList& origPatchFaceCells = origPp.faceCells();
            const labelUList& nbrOrigPatchFaceCells = nbrOrigPp.faceCells();

            // Add all possible interfaces between processors that contain the
            // original patches
            forAll(origPatchFaceCells, origFacei)
            {
                const label ownerProc =
                    cellProc[origPatchFaceCells[origFacei]];

                forAll(nbrOrigPatchFaceCells, nbrOrigFacei)
                {
                    const label nbrProc =
                        cellProc[nbrOrigPatchFaceCells[nbrOrigFacei]];

                    if (bop(ownerProc, nbrProc))
                    {
                        addInterProcFace
                        (
                            -1,
                            ownerProc,
                            nbrProc,
                            patchi,
                            procNbrToInterPatch,
                            interPatchFaces,
                            subPatchIDs,
                            subPatchStarts
                        );
                    }
                }
            }
        }
        else if (isA<cyclicPolyPatch>(patches[patchi]))
        {
            const cyclicPolyPatch& cPp =
                refCast<const cyclicPolyPatch>(patches[patchi]);

            if (cPp.owner() != owner) continue;

            const labelUList& patchFaceCells = cPp.faceCells();
            const labelUList& nbrPatchFaceCells = cPp.nbrPatch().faceCells();

            // Add faces with different owner and neighbour processors
            forAll(patchFaceCells, facei)
            {
                const label ownerProc = cellProc[patchFaceCells[facei]];
                const label nbrProc = cellProc[nbrPatchFaceCells[facei]];

                if (bop(ownerProc, nbrProc))
                {
                    addInterProcFace
                    (
                        cPp.start() + facei,
                        ownerProc,
                        nbrProc,
                        patchi,
                        procNbrToInterPatch,
                        interPatchFaces,
                        subPatchIDs,
                        subPatchStarts
                    );
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::domainDecomposition::decompose()
{
    // Decide which cell goes to which processor
    cellProc_ = distributeCells();

    // Distribute the cells according to the given processor label

    // calculate the addressing information for the original mesh
    Info<< "\nCalculating original mesh data" << endl;

    // set references to the original mesh
    const polyBoundaryMesh& patches = completeMesh().boundaryMesh();
    const faceList& fcs = completeMesh().faces();
    const labelList& owner = completeMesh().faceOwner();
    const labelList& neighbour = completeMesh().faceNeighbour();

    // loop through the list of processor labels for the cell and add the
    // cell shape to the list of cells for the appropriate processor

    Info<< "\nDistributing cells to processors" << endl;

    // Cells per processor
    procCellAddressing_ = invertOneToMany(nProcs(), cellProc_);

    Info<< "\nDistributing faces to processors" << endl;

    // Loop through all internal faces and decide which processor they belong to
    // First visit all internal faces. If cells at both sides belong to the
    // same processor, the face is an internal face. If they are different,
    // it belongs to both processors.
    List<DynamicList<label>> dynProcFaceAddressing(nProcs());

    // Internal faces
    forAll(neighbour, facei)
    {
        if (cellProc_[owner[facei]] == cellProc_[neighbour[facei]])
        {
            // Face internal to processor. Notice no turning index.
            dynProcFaceAddressing[cellProc_[owner[facei]]].append(facei+1);
        }
    }

    // for all processors, set the size of start index and patch size
    // lists to the number of patches in the mesh
    labelListList procPatchSize(nProcs());
    labelListList procPatchStartIndex(nProcs());
    forAll(procPatchSize, proci)
    {
        procPatchSize[proci].setSize(patches.size());
        procPatchStartIndex[proci].setSize(patches.size());
    }

    // Patch faces
    forAll(patches, patchi)
    {
        // Reset size and start index for all processors
        forAll(procPatchSize, proci)
        {
            procPatchSize[proci][patchi] = 0;
            procPatchStartIndex[proci][patchi] =
                dynProcFaceAddressing[proci].size();
        }

        const label patchStart = patches[patchi].start();

        if (!isA<cyclicPolyPatch>(patches[patchi]))
        {
            // Normal patch. Add faces to processor where the cell
            // next to the face lives.

            const labelUList& patchFaceCells = patches[patchi].faceCells();

            forAll(patchFaceCells, facei)
            {
                const label curProc = cellProc_[patchFaceCells[facei]];

                // add the face without turning index
                dynProcFaceAddressing[curProc].append(patchStart+facei+1);

                // increment the number of faces for this patch
                procPatchSize[curProc][patchi]++;
            }
        }
        else
        {
            // Cyclic patch. Add if the opposite sides are on the same
            // processor.

            const cyclicPolyPatch& pp =
                refCast<const cyclicPolyPatch>(patches[patchi]);

            const labelUList& patchFaceCells = pp.faceCells();
            const labelUList& nbrPatchFaceCells = pp.nbrPatch().faceCells();

            forAll(patchFaceCells, facei)
            {
                const label curProc = cellProc_[patchFaceCells[facei]];
                const label nbrProc = cellProc_[nbrPatchFaceCells[facei]];

                if (curProc == nbrProc)
                {
                    // add the face without turning index
                    dynProcFaceAddressing[curProc].append(patchStart+facei+1);

                    // increment the number of faces for this patch
                    procPatchSize[curProc][patchi]++;
                }
            }
        }
    }

    // Done internal bits of the new mesh and the ordinary patches.

    // Per processor, from neighbour processor to the inter-processor patch
    // that communicates with that neighbour
    List<Map<label>> procNbrToInterPatch(nProcs());

    // Per processor the faces per inter-processor patch
    List<DynamicList<DynamicList<label>>> interPatchFaces(nProcs());

    // Sub-patch ID-s and starts for inter-processor patches
    List<labelListList> subPatchIDs(nProcs());
    List<labelListList> subPatchStarts(nProcs());

    // Processor boundaries from internal faces
    forAll(neighbour, facei)
    {
        const label ownerProc = cellProc_[owner[facei]];
        const label nbrProc = cellProc_[neighbour[facei]];

        if (ownerProc != nbrProc)
        {
            // inter - processor patch face found.
            addInterProcFace
            (
                facei,
                ownerProc,
                nbrProc,
                -1,
                procNbrToInterPatch,
                interPatchFaces,
                subPatchIDs,
                subPatchStarts
            );
        }
    }

    // Special handling needed for the case that multiple processor cyclic
    // patches are created on each local processor domain, e.g. if a 3x3 case
    // is decomposed using the decomposition:
    //
    //              | 1 | 0 | 2 |
    //  cyclic left | 2 | 0 | 1 | cyclic right
    //              | 2 | 0 | 1 |
    //
    // - processors 1 and 2 will both have pieces of both cyclic left- and
    //   right sub-patches present
    // - the interface patch faces are stored in a single list, where each
    //   sub-patch is referenced into the list using a patch start index and
    //   size
    // - if the patches are in order (in the boundary file) of left, right
    //   - processor 1 will send: left, right
    //   - processor 1 will need to receive in reverse order: right, left
    //   - similarly for processor 2
    // - the sub-patches are therefore generated in 4 passes of the patch lists
    //   1. add faces from owner patch where local proc i < nbr proc i
    //   2. add faces from nbr patch where local proc i < nbr proc i
    //   3. add faces from nbr patch where local proc i > nbr proc i
    //   4. add faces from owner patch where local proc i > nbr proc i

    processInterCyclics
    (
        cellProc_,
        patches,
        interPatchFaces,
        procNbrToInterPatch,
        subPatchIDs,
        subPatchStarts,
        true,
        lessOp<label>()
    );

    processInterCyclics
    (
        cellProc_,
        patches,
        interPatchFaces,
        procNbrToInterPatch,
        subPatchIDs,
        subPatchStarts,
        false,
        lessOp<label>()
    );

    processInterCyclics
    (
        cellProc_,
        patches,
        interPatchFaces,
        procNbrToInterPatch,
        subPatchIDs,
        subPatchStarts,
        false,
        greaterOp<label>()
    );

    processInterCyclics
    (
        cellProc_,
        patches,
        interPatchFaces,
        procNbrToInterPatch,
        subPatchIDs,
        subPatchStarts,
        true,
        greaterOp<label>()
    );

    // Sort inter-proc patch by neighbour
    labelListList procNeighbourProcessors(nProcs());
    labelListList procProcessorPatchSize(nProcs());
    labelListList procProcessorPatchStartIndex(nProcs());
    List<labelListList> procProcessorPatchSubPatchIDs(nProcs());
    List<labelListList> procProcessorPatchSubPatchStarts(nProcs());
    labelList order;
    forAll(procNbrToInterPatch, proci)
    {
        label nInterfaces = procNbrToInterPatch[proci].size();

        procNeighbourProcessors[proci].setSize(nInterfaces);
        procProcessorPatchSize[proci].setSize(nInterfaces);
        procProcessorPatchStartIndex[proci].setSize(nInterfaces);
        procProcessorPatchSubPatchIDs[proci].setSize(nInterfaces);
        procProcessorPatchSubPatchStarts[proci].setSize(nInterfaces);

        // Get sorted neighbour processors
        const Map<label>& curNbrToInterPatch = procNbrToInterPatch[proci];
        labelList nbrs = curNbrToInterPatch.toc();

        sortedOrder(nbrs, order);

        DynamicList<DynamicList<label>>& curInterPatchFaces =
            interPatchFaces[proci];

        forAll(nbrs, i)
        {
            const label nbrProc = nbrs[i];
            const label interPatch = curNbrToInterPatch[nbrProc];

            procNeighbourProcessors[proci][i] = nbrProc;
            procProcessorPatchSize[proci][i] =
                curInterPatchFaces[interPatch].size();
            procProcessorPatchStartIndex[proci][i] =
                dynProcFaceAddressing[proci].size();

            // Add size as last element to substarts and transfer
            subPatchStarts[proci][interPatch].append
            (
                curInterPatchFaces[interPatch].size()
            );
            procProcessorPatchSubPatchIDs[proci][i].transfer
            (
                subPatchIDs[proci][interPatch]
            );
            procProcessorPatchSubPatchStarts[proci][i].transfer
            (
                subPatchStarts[proci][interPatch]
            );

            // And add all the face labels for interPatch
            DynamicList<label>& interPatchFaces =
                curInterPatchFaces[interPatch];

            forAll(interPatchFaces, j)
            {
                dynProcFaceAddressing[proci].append(interPatchFaces[j]);
            }

            interPatchFaces.clearStorage();
        }

        curInterPatchFaces.clearStorage();
        dynProcFaceAddressing[proci].shrink();
    }

    // Transfer face addressing to non-dynamic member data
    forAll(dynProcFaceAddressing, proci)
    {
        procFaceAddressing_[proci].transfer(dynProcFaceAddressing[proci]);
    }

    if (debug)
    {
        forAll(procPatchStartIndex, proci)
        {
            Info<< "Processor:" << proci << endl;
            Info<< "    total faces:" << procFaceAddressing_[proci].size()
                << endl;

            const labelList& curProcPatchStartIndex =
                procPatchStartIndex[proci];

            forAll(curProcPatchStartIndex, patchi)
            {
                Info<< "    patch:" << patchi
                    << "\tstart:" << curProcPatchStartIndex[patchi]
                    << "\tsize:" << procPatchSize[proci][patchi]
                    << endl;
            }
        }
        Info<< endl;

        forAll(procNeighbourProcessors, proci)
        {
            Info<< "Processor " << proci << endl;
            forAll(procNeighbourProcessors[proci], i)
            {
                Info<< "    nbr:" << procNeighbourProcessors[proci][i] << endl;
                Info<< "    size:" << procProcessorPatchSize[proci][i] << endl;
                Info<< "    start:" << procProcessorPatchStartIndex[proci][i]
                    << endl;
            }
        }
        Info<< endl;
    }

    if (debug > 1)
    {
        forAll(procFaceAddressing_, proci)
        {
            Info<< "Processor:" << proci << endl;
            Info<< "    faces:" << procFaceAddressing_[proci] << endl;
        }
    }

    Info<< "\nDistributing points to processors" << endl;

    // For every processor, loop through the list of faces for the processor.
    // For every face, loop through the list of points and mark the point as
    // used for the processor. Collect the list of used points for the
    // processor.

    forAll(procPointAddressing_, proci)
    {
        boolList pointLabels(completeMesh().nPoints(), false);

        // Get reference to list of used faces
        const labelList& procFaceLabels = procFaceAddressing_[proci];

        forAll(procFaceLabels, facei)
        {
            // Because of the turning index, some labels may be negative
            const labelList& facePoints = fcs[mag(procFaceLabels[facei]) - 1];

            forAll(facePoints, pointi)
            {
                // Mark the point as used
                pointLabels[facePoints[pointi]] = true;
            }
        }

        // Collect the used points
        labelList& procPointLabels = procPointAddressing_[proci];

        procPointLabels.setSize(pointLabels.size());

        label nUsedPoints = 0;

        forAll(pointLabels, pointi)
        {
            if (pointLabels[pointi])
            {
                procPointLabels[nUsedPoints] = pointi;

                nUsedPoints++;
            }
        }

        // Reset the size of used points
        procPointLabels.setSize(nUsedPoints);
    }

    Info<< "\nConstructing processor meshes" << endl;

    // Mark point/faces/cells that are in zones.
    // -1   : not in zone
    // -2   : in multiple zones
    // >= 0 : in single given zone
    // This will give direct lookup of elements that are in a single zone
    // and we'll only have to revert back to searching through all zones
    // for the duplicate elements

    // Point zones
    labelList pointToZone(completeMesh().points().size(), -1);
    forAll(completeMesh().pointZones(), zoneI)
    {
        mark(completeMesh().pointZones()[zoneI], zoneI, pointToZone);
    }

    // Face zones
    labelList faceToZone(completeMesh().faces().size(), -1);
    forAll(completeMesh().faceZones(), zoneI)
    {
        mark(completeMesh().faceZones()[zoneI], zoneI, faceToZone);
    }

    // Cell zones
    labelList cellToZone(completeMesh().nCells(), -1);
    forAll(completeMesh().cellZones(), zoneI)
    {
        mark(completeMesh().cellZones()[zoneI], zoneI, cellToZone);
    }

    // Initialise information for reporting
    label maxProcCells = 0;
    label totProcFaces = 0;
    label maxProcPatches = 0;
    label totProcPatches = 0;
    label maxProcFaces = 0;

    // Generate the meshes
    for (label proci = 0; proci < nProcs(); proci++)
    {
        // Create processor points
        const labelList& curPointLabels = procPointAddressing_[proci];
        const pointField& meshPoints = completeMesh().points();
        labelList pointLookup(completeMesh().nPoints(), -1);
        pointField procPoints(curPointLabels.size());
        forAll(curPointLabels, pointi)
        {
            procPoints[pointi] = meshPoints[curPointLabels[pointi]];
            pointLookup[curPointLabels[pointi]] = pointi;
        }

        // Create processor faces
        const labelList& curFaceLabels = procFaceAddressing_[proci];
        const faceList& meshFaces = completeMesh().faces();
        labelList faceLookup(completeMesh().nFaces(), -1);
        faceList procFaces(curFaceLabels.size());
        forAll(curFaceLabels, facei)
        {
            // Mark the original face as used
            // Remember to decrement the index by one (turning index)
            label curF = mag(curFaceLabels[facei]) - 1;

            faceLookup[curF] = facei;

            // get the original face
            labelList origFaceLabels;

            if (curFaceLabels[facei] >= 0)
            {
                // face not turned
                origFaceLabels = meshFaces[curF];
            }
            else
            {
                origFaceLabels = meshFaces[curF].reverseFace();
            }

            // translate face labels into local point list
            face& procFaceLabels = procFaces[facei];

            procFaceLabels.setSize(origFaceLabels.size());

            forAll(origFaceLabels, pointi)
            {
                procFaceLabels[pointi] = pointLookup[origFaceLabels[pointi]];
            }
        }

        // Create processor cells
        const labelList& curCellLabels = procCellAddressing_[proci];
        const cellList& meshCells = completeMesh().cells();
        cellList procCells(curCellLabels.size());
        forAll(curCellLabels, celli)
        {
            const labelList& origCellLabels = meshCells[curCellLabels[celli]];
            cell& curCell = procCells[celli];
            curCell.setSize(origCellLabels.size());
            forAll(origCellLabels, cellFacei)
            {
                curCell[cellFacei] = faceLookup[origCellLabels[cellFacei]];
            }
        }

        // Create a processor mesh without a boundary
        procMeshes_.set
        (
            proci,
            new fvMesh
            (
                IOobject
                (
                    regionName_,
                    completeMesh().facesInstance(),
                    runTimes_.procTimes()[proci]
                ),
                move(procPoints),
                move(procFaces),
                move(procCells)
            )
        );
        fvMesh& procMesh = procMeshes_[proci];
        procMesh.setPointsInstance(completeMesh().pointsInstance());

        // Create processor boundary patch information
        const labelList& curPatchSizes = procPatchSize[proci];
        const labelList& curPatchStarts = procPatchStartIndex[proci];
        const labelList& curNeighbourProcessors =
            procNeighbourProcessors[proci];
        const labelList& curProcessorPatchSizes =
            procProcessorPatchSize[proci];
        const labelList& curProcessorPatchStarts =
            procProcessorPatchStartIndex[proci];
        const labelListList& curSubPatchIDs =
            procProcessorPatchSubPatchIDs[proci];
        const labelListList& curSubStarts =
            procProcessorPatchSubPatchStarts[proci];
        const polyPatchList& meshPatches = completeMesh().boundaryMesh();

        // Count the number of inter-proc patches
        label nInterProcPatches = 0;
        forAll(curSubPatchIDs, procPatchi)
        {
            nInterProcPatches += curSubPatchIDs[procPatchi].size();
        }
        List<polyPatch*> procPatches
        (
            curPatchSizes.size() + nInterProcPatches,
            nullptr
        );

        label nPatches = 0;

        // Copy existing non-proc patches
        forAll(curPatchSizes, patchi)
        {
            const polyPatch& meshPatch = meshPatches[patchi];

            procPatches[nPatches] =
                meshPatch.clone
                (
                    procMesh.boundaryMesh(),
                    nPatches,
                    curPatchSizes[patchi],
                    curPatchStarts[patchi]
                ).ptr();

            nPatches++;
        }

        // Create new inter-proc patches
        forAll(curProcessorPatchSizes, procPatchi)
        {
            const labelList& subPatchID = curSubPatchIDs[procPatchi];
            const labelList& subStarts = curSubStarts[procPatchi];

            label curStart = curProcessorPatchStarts[procPatchi];

            forAll(subPatchID, i)
            {
                const label size =
                    i < subPatchID.size() - 1
                  ? subStarts[i+1] - subStarts[i]
                  : curProcessorPatchSizes[procPatchi] - subStarts[i];

                if (subPatchID[i] == -1)
                {
                    // Processor patch from internal faces
                    procPatches[nPatches] =
                        new processorPolyPatch
                        (
                            size,
                            curStart,
                            nPatches,
                            procMesh.boundaryMesh(),
                            proci,
                            curNeighbourProcessors[procPatchi]
                        );
                }
                else if
                (
                    isA<nonConformalCyclicPolyPatch>
                    (
                        completeMesh().boundaryMesh()[subPatchID[i]]
                    )
                )
                {
                    // Non-conformal processor cyclic patch from split
                    // non-conformal cyclic patch
                    const nonConformalCyclicPolyPatch& nccPp =
                        refCast<const nonConformalCyclicPolyPatch>
                        (
                            completeMesh().boundaryMesh()[subPatchID[i]]
                        );

                    procPatches[nPatches] =
                        new nonConformalProcessorCyclicPolyPatch
                        (
                            size,
                            curStart,
                            nPatches,
                            procMesh.boundaryMesh(),
                            proci,
                            curNeighbourProcessors[procPatchi],
                            nccPp.name(),
                            nccPp.origPatch().name()
                        );
                }
                else if
                (
                    isA<cyclicPolyPatch>
                    (
                        completeMesh().boundaryMesh()[subPatchID[i]]
                    )
                )
                {
                    // Processor cyclic patch from split cyclic faces
                    const cyclicPolyPatch& cPp =
                        refCast<const cyclicPolyPatch>
                        (
                            completeMesh().boundaryMesh()[subPatchID[i]]
                        );

                    procPatches[nPatches] =
                        new processorCyclicPolyPatch
                        (
                            size,
                            curStart,
                            nPatches,
                            procMesh.boundaryMesh(),
                            proci,
                            curNeighbourProcessors[procPatchi],
                            cPp.name()
                        );
                }
                else
                {
                    FatalErrorInFunction
                        << "Sub patch ID set for non-cyclic patch type"
                        << exit(FatalError);
                }

                curStart += size;

                nPatches++;
            }
        }

        // Add patches to the mesh
        procMesh.addFvPatches(procPatches);

        // Create point zones
        {
            const meshPointZones& pz = completeMesh().pointZones();

            // Go through all the zoned points and find out if they
            // belong to a zone.  If so, add it to the zone as
            // necessary
            List<DynamicList<label>> zonePoints(pz.size());

            // Estimate size
            forAll(zonePoints, zoneI)
            {
                zonePoints[zoneI].setCapacity(pz[zoneI].size()/nProcs());
            }

            // Use the pointToZone map to find out the single zone (if any),
            // use slow search only for shared points.
            forAll(curPointLabels, pointi)
            {
                label curPoint = curPointLabels[pointi];

                label zoneI = pointToZone[curPoint];

                if (zoneI >= 0)
                {
                    // Single zone.
                    zonePoints[zoneI].append(pointi);
                }
                else if (zoneI == -2)
                {
                    // Multiple zones. Lookup.
                    forAll(pz, zoneI)
                    {
                        label index = pz[zoneI].whichPoint(curPoint);

                        if (index != -1)
                        {
                            zonePoints[zoneI].append(pointi);
                        }
                    }
                }
            }

            procMesh.pointZones().clearAddressing();
            procMesh.pointZones().setSize(zonePoints.size());
            forAll(zonePoints, zoneI)
            {
                procMesh.pointZones().set
                (
                    zoneI,
                    pz[zoneI].clone
                    (
                        procMesh.pointZones(),
                        zoneI,
                        zonePoints[zoneI].shrink()
                    )
                );
            }

            if (pz.size())
            {
                // Force writing on all processors
                procMesh.pointZones().writeOpt() = IOobject::AUTO_WRITE;
            }
        }

        // Create face zones
        {
            const meshFaceZones& fz = completeMesh().faceZones();

            // Go through all the zoned face and find out if they
            // belong to a zone.  If so, add it to the zone as
            // necessary
            List<DynamicList<label>> zoneFaces(fz.size());
            List<DynamicList<bool>> zoneFaceFlips(fz.size());

            // Estimate size
            forAll(zoneFaces, zoneI)
            {
                label procSize = fz[zoneI].size()/nProcs();

                zoneFaces[zoneI].setCapacity(procSize);
                zoneFaceFlips[zoneI].setCapacity(procSize);
            }

            // Go through all the zoned faces and find out if they
            // belong to a zone.  If so, add it to the zone as
            // necessary
            forAll(curFaceLabels, facei)
            {
                // Remember to decrement the index by one (turning index)
                //
                label curF = mag(curFaceLabels[facei]) - 1;

                label zoneI = faceToZone[curF];

                if (zoneI >= 0)
                {
                    // Single zone. Add the face
                    zoneFaces[zoneI].append(facei);

                    label index = fz[zoneI].whichFace(curF);

                    bool flip = fz[zoneI].flipMap()[index];

                    if (curFaceLabels[facei] < 0)
                    {
                        flip = !flip;
                    }

                    zoneFaceFlips[zoneI].append(flip);
                }
                else if (zoneI == -2)
                {
                    // Multiple zones. Lookup.
                    forAll(fz, zoneI)
                    {
                        label index = fz[zoneI].whichFace(curF);

                        if (index != -1)
                        {
                            zoneFaces[zoneI].append(facei);

                            bool flip = fz[zoneI].flipMap()[index];

                            if (curFaceLabels[facei] < 0)
                            {
                                flip = !flip;
                            }

                            zoneFaceFlips[zoneI].append(flip);
                        }
                    }
                }
            }

            procMesh.faceZones().clearAddressing();
            procMesh.faceZones().setSize(zoneFaces.size());
            forAll(zoneFaces, zoneI)
            {
                procMesh.faceZones().set
                (
                    zoneI,
                    fz[zoneI].clone
                    (
                        zoneFaces[zoneI].shrink(),          // addressing
                        zoneFaceFlips[zoneI].shrink(),      // flipmap
                        zoneI,
                        procMesh.faceZones()
                    )
                );
            }

            if (fz.size())
            {
                // Force writing on all processors
                procMesh.faceZones().writeOpt() = IOobject::AUTO_WRITE;
            }
        }

        // Create cell zones
        {
            const meshCellZones& cz = completeMesh().cellZones();

            // Go through all the zoned cells and find out if they
            // belong to a zone.  If so, add it to the zone as
            // necessary
            List<DynamicList<label>> zoneCells(cz.size());

            // Estimate size
            forAll(zoneCells, zoneI)
            {
                zoneCells[zoneI].setCapacity(cz[zoneI].size()/nProcs());
            }

            forAll(curCellLabels, celli)
            {
                label curCelli = curCellLabels[celli];

                label zoneI = cellToZone[curCelli];

                if (zoneI >= 0)
                {
                    // Single zone.
                    zoneCells[zoneI].append(celli);
                }
                else if (zoneI == -2)
                {
                    // Multiple zones. Lookup.
                    forAll(cz, zoneI)
                    {
                        label index = cz[zoneI].whichCell(curCelli);

                        if (index != -1)
                        {
                            zoneCells[zoneI].append(celli);
                        }
                    }
                }
            }

            procMesh.cellZones().clearAddressing();
            procMesh.cellZones().setSize(zoneCells.size());
            forAll(zoneCells, zoneI)
            {
                procMesh.cellZones().set
                (
                    zoneI,
                    cz[zoneI].clone
                    (
                        zoneCells[zoneI].shrink(),
                        zoneI,
                        procMesh.cellZones()
                    )
                );
            }

            if (cz.size())
            {
                // Force writing on all processors
                procMesh.cellZones().writeOpt() = IOobject::AUTO_WRITE;
            }
        }

        // Report processor and update global statistics
        {
            Info<< endl
                << "Processor " << proci << nl
                << "    Number of cells = " << procMesh.nCells()
                << endl;

            maxProcCells = max(maxProcCells, procMesh.nCells());

            label nBoundaryFaces = 0;
            label nProcPatches = 0;
            label nProcFaces = 0;

            forAll(procMesh.boundaryMesh(), patchi)
            {
                if (isA<processorPolyPatch>(procMesh.boundaryMesh()[patchi]))
                {
                    const processorPolyPatch& ppp =
                    refCast<const processorPolyPatch>
                    (
                        procMesh.boundaryMesh()[patchi]
                    );

                    Info<< "    Number of faces shared with processor "
                        << ppp.neighbProcNo() << " = " << ppp.size() << endl;

                    nProcPatches++;
                    nProcFaces += ppp.size();
                }
                else
                {
                    nBoundaryFaces += procMesh.boundaryMesh()[patchi].size();
                }
            }

            Info<< "    Number of processor patches = " << nProcPatches << nl
                << "    Number of processor faces = " << nProcFaces << nl
                << "    Number of boundary faces = " << nBoundaryFaces << endl;

            totProcFaces += nProcFaces;
            totProcPatches += nProcPatches;
            maxProcPatches = max(maxProcPatches, nProcPatches);
            maxProcFaces = max(maxProcFaces, nProcFaces);
        }
    }

    // Determine the average number of processor elements
    scalar avgProcCells = scalar(completeMesh().nCells())/nProcs();
    scalar avgProcPatches = scalar(totProcPatches)/nProcs();
    scalar avgProcFaces = scalar(totProcFaces)/nProcs();

    // Prevent division by zero in the case of all faces on one processor
    if (totProcPatches == 0)
    {
        avgProcPatches = 1;
    }
    if (totProcFaces == 0)
    {
        avgProcFaces = 1;
    }

    Info<< nl
        << "Number of processor faces = " << totProcFaces/2 << nl
        << "Max number of cells = " << maxProcCells
        << " (" << 100.0*(maxProcCells-avgProcCells)/avgProcCells
        << "% above average " << avgProcCells << ")" << nl
        << "Max number of processor patches = " << maxProcPatches
        << " (" << 100.0*(maxProcPatches-avgProcPatches)/avgProcPatches
        << "% above average " << avgProcPatches << ")" << nl
        << "Max number of faces between processors = " << maxProcFaces
        << " (" << 100.0*(maxProcFaces-avgProcFaces)/avgProcFaces
        << "% above average " << avgProcFaces << ")" << nl
        << endl;

    // Unconform any non-conformal parts of the processor meshes
    unconform();
}


// ************************************************************************* //
