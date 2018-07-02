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

#include "InteractionLists.H"
#include "globalIndexAndTransform.H"
#include "indexedOctree.H"
#include "treeDataFace.H"
#include "treeDataCell.H"
#include "volFields.H"
#include "meshTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ParticleType>
void Foam::InteractionLists<ParticleType>::buildInteractionLists()
{
    Info<< "Building InteractionLists with interaction distance "
        << maxDistance_ << endl;

    const vector interactionVec = maxDistance_*vector::one;

    treeBoundBox procBb(treeBoundBox(mesh_.points()));

    treeBoundBox extendedProcBb
    (
        procBb.min() - interactionVec,
        procBb.max() + interactionVec
    );

    treeBoundBoxList allExtendedProcBbs(Pstream::nProcs());

    allExtendedProcBbs[Pstream::myProcNo()] = extendedProcBb;

    Pstream::gatherList(allExtendedProcBbs);

    Pstream::scatterList(allExtendedProcBbs);

    List<treeBoundBox> extendedProcBbsInRange;
    List<label> extendedProcBbsTransformIndex;
    List<label> extendedProcBbsOrigProc;

    findExtendedProcBbsInRange
    (
        procBb,
        allExtendedProcBbs,
        mesh_.globalData().globalTransforms(),
        extendedProcBbsInRange,
        extendedProcBbsTransformIndex,
        extendedProcBbsOrigProc
    );

    treeBoundBoxList cellBbs(mesh_.nCells());

    forAll(cellBbs, celli)
    {
        cellBbs[celli] = treeBoundBox
        (
            mesh_.cells()[celli].points
            (
                mesh_.faces(),
                mesh_.points()
            )
        );
    }

    const globalIndexAndTransform& globalTransforms =
        mesh_.globalData().globalTransforms();

    // Recording which cells are in range of an extended boundBox, as
    // only these cells will need to be tested to determine which
    // referred cells that they interact with.
    PackedBoolList cellInRangeOfCoupledPatch(mesh_.nCells(), false);

    // IAndT: index (=local cell index) and transform (from
    // globalIndexAndTransform)
    DynamicList<labelPair> cellIAndTToExchange;

    DynamicList<treeBoundBox> cellBbsToExchange;

    DynamicList<label> procToDistributeCellTo;

    forAll(extendedProcBbsInRange, ePBIRI)
    {
        const treeBoundBox& otherExtendedProcBb =
            extendedProcBbsInRange[ePBIRI];

        label transformIndex = extendedProcBbsTransformIndex[ePBIRI];

        label origProc = extendedProcBbsOrigProc[ePBIRI];

        forAll(cellBbs, celli)
        {
            const treeBoundBox& cellBb = cellBbs[celli];

            if (cellBb.overlaps(otherExtendedProcBb))
            {
                // This cell is in range of the Bb of the other
                // processor Bb, and so needs to be referred to it

                cellInRangeOfCoupledPatch[celli] = true;

                cellIAndTToExchange.append
                (
                    globalTransforms.encode(celli, transformIndex)
                );

                cellBbsToExchange.append(cellBb);

                procToDistributeCellTo.append(origProc);
            }
        }
    }

    buildMap(cellMapPtr_, procToDistributeCellTo);

    // Needed for reverseDistribute
    label preDistributionCellMapSize = procToDistributeCellTo.size();

    cellMap().distribute(cellBbsToExchange);

    cellMap().distribute(cellIAndTToExchange);

    // Determine labelList specifying only cells that are in range of
    // a coupled boundary to build an octree limited to these cells.
    DynamicList<label> coupledPatchRangeCells;

    forAll(cellInRangeOfCoupledPatch, celli)
    {
        if (cellInRangeOfCoupledPatch[celli])
        {
            coupledPatchRangeCells.append(celli);
        }
    }

    treeBoundBox procBbRndExt
    (
        treeBoundBox(mesh_.points()).extend(1e-4)
    );

    indexedOctree<treeDataCell> coupledPatchRangeTree
    (
        treeDataCell
        (
            true,                   // Cache cell bb
            mesh_,
            coupledPatchRangeCells, // Subset of mesh
            polyMesh::CELL_TETS     // Consistent with tracking
        ),
        procBbRndExt,
        8,              // maxLevel,
        10,             // leafSize,
        100.0           // duplicity
    );

    ril_.setSize(cellBbsToExchange.size());

    // This needs to be a boolList, not PackedBoolList if
    // reverseDistribute is called.
    boolList cellBbRequiredByAnyCell(cellBbsToExchange.size(), false);

    Info<< "    Building referred interaction lists" << endl;

    forAll(cellBbsToExchange, bbI)
    {
        const labelPair& ciat = cellIAndTToExchange[bbI];

        const vectorTensorTransform& transform = globalTransforms.transform
        (
            globalTransforms.transformIndex(ciat)
        );

        treeBoundBox tempTransformedBb
        (
            transform.invTransformPosition(cellBbsToExchange[bbI].points())
        );

        treeBoundBox extendedBb
        (
            tempTransformedBb.min() - interactionVec,
            tempTransformedBb.max() + interactionVec
        );

        // Find all elements intersecting box.
        labelList interactingElems
        (
            coupledPatchRangeTree.findBox(extendedBb)
        );

        if (!interactingElems.empty())
        {
            cellBbRequiredByAnyCell[bbI] = true;
        }

        ril_[bbI].setSize(interactingElems.size(), -1);

        forAll(interactingElems, i)
        {
            label elemI = interactingElems[i];

            // Here, a more detailed geometric test could be applied,
            // i.e. a more accurate bounding volume like a OBB or
            // convex hull or an exact geometrical test.

            label c = coupledPatchRangeTree.shapes().cellLabels()[elemI];

            ril_[bbI][i] = c;
        }
    }

    // Perform subset of ril_, to remove any referred cells that do
    // not interact.  They will not be sent from originating
    // processors.  This assumes that the disappearance of the cell
    // from the sending list of the source processor, simply removes
    // the referred cell from the ril_, all of the subsequent indices
    // shuffle down one, but the order and structure is preserved,
    // i.e. it, is as if the cell had never been referred in the first
    // place.

    inplaceSubset(cellBbRequiredByAnyCell, ril_);

    // Send information about which cells are actually required back
    // to originating processors.

    // At this point, cellBbsToExchange does not need to be maintained
    // or distributed as it is not longer needed.

    cellBbsToExchange.setSize(0);

    cellMap().reverseDistribute
    (
        preDistributionCellMapSize,
        cellBbRequiredByAnyCell
    );

    cellMap().reverseDistribute
    (
        preDistributionCellMapSize,
        cellIAndTToExchange
    );

    // Perform ordering of cells to send, this invalidates the
    // previous value of preDistributionCellMapSize.

    preDistributionCellMapSize = -1;

    // Move all of the used cells to refer to the start of the list
    // and truncate it

    inplaceSubset(cellBbRequiredByAnyCell, cellIAndTToExchange);

    inplaceSubset(cellBbRequiredByAnyCell, procToDistributeCellTo);

    preDistributionCellMapSize = procToDistributeCellTo.size();

    // Rebuild mapDistribute with only required referred cells
    buildMap(cellMapPtr_, procToDistributeCellTo);

    // Store cellIndexAndTransformToDistribute
    cellIndexAndTransformToDistribute_.transfer(cellIAndTToExchange);

    // Determine inverse addressing for referred cells

    rilInverse_.setSize(mesh_.nCells());

    // Temporary Dynamic lists for accumulation
    List<DynamicList<label>> rilInverseTemp(rilInverse_.size());

    // Loop over all referred cells
    forAll(ril_, refCelli)
    {
        const labelList& realCells = ril_[refCelli];

        // Loop over all real cells in that the referred cell is to
        // supply interactions to and record the index of this
        // referred cell in the real cells entry in rilInverse

        forAll(realCells, realCelli)
        {
            rilInverseTemp[realCells[realCelli]].append(refCelli);
        }
    }

    forAll(rilInverse_, celli)
    {
        rilInverse_[celli].transfer(rilInverseTemp[celli]);
    }

    // Determine which wall faces to refer

    // The referring of wall patch data relies on patches having the same
    // index on each processor.
    mesh_.boundaryMesh().checkParallelSync(true);

    // Determine the index of all of the wall faces on this processor
    DynamicList<label> localWallFaces;

    forAll(mesh_.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchi];

        if (isA<wallPolyPatch>(patch))
        {
            localWallFaces.append(identity(patch.size()) + patch.start());
        }
    }

    treeBoundBoxList wallFaceBbs(localWallFaces.size());

    forAll(wallFaceBbs, i)
    {
        wallFaceBbs[i] = treeBoundBox
        (
            mesh_.faces()[localWallFaces[i]].points(mesh_.points())
        );
    }

    // IAndT: index and transform
    DynamicList<labelPair> wallFaceIAndTToExchange;

    DynamicList<treeBoundBox> wallFaceBbsToExchange;

    DynamicList<label> procToDistributeWallFaceTo;

    forAll(extendedProcBbsInRange, ePBIRI)
    {
        const treeBoundBox& otherExtendedProcBb =
            extendedProcBbsInRange[ePBIRI];

        label transformIndex = extendedProcBbsTransformIndex[ePBIRI];

        label origProc = extendedProcBbsOrigProc[ePBIRI];

        forAll(wallFaceBbs, i)
        {
            const treeBoundBox& wallFaceBb = wallFaceBbs[i];

            if (wallFaceBb.overlaps(otherExtendedProcBb))
            {
                // This wall face is in range of the Bb of the other
                // processor Bb, and so needs to be referred to it

                label wallFacei = localWallFaces[i];

                wallFaceIAndTToExchange.append
                (
                    globalTransforms.encode(wallFacei, transformIndex)
                );

                wallFaceBbsToExchange.append(wallFaceBb);

                procToDistributeWallFaceTo.append(origProc);
            }
        }
    }

    buildMap(wallFaceMapPtr_, procToDistributeWallFaceTo);

    // Needed for reverseDistribute
    label preDistributionWallFaceMapSize = procToDistributeWallFaceTo.size();

    wallFaceMap().distribute(wallFaceBbsToExchange);

    wallFaceMap().distribute(wallFaceIAndTToExchange);

    indexedOctree<treeDataCell> allCellsTree
    (
        treeDataCell(true, mesh_, polyMesh::CELL_TETS),
        procBbRndExt,
        8,              // maxLevel,
        10,             // leafSize,
        100.0           // duplicity
    );

    rwfil_.setSize(wallFaceBbsToExchange.size());

    // This needs to be a boolList, not PackedBoolList if
    // reverseDistribute is called.
    boolList wallFaceBbRequiredByAnyCell(wallFaceBbsToExchange.size(), false);

    forAll(wallFaceBbsToExchange, bbI)
    {
        const labelPair& wfiat = wallFaceIAndTToExchange[bbI];

        const vectorTensorTransform& transform = globalTransforms.transform
        (
            globalTransforms.transformIndex(wfiat)
        );

        treeBoundBox tempTransformedBb
        (
            transform.invTransformPosition(wallFaceBbsToExchange[bbI].points())
        );

        treeBoundBox extendedBb
        (
            tempTransformedBb.min() - interactionVec,
            tempTransformedBb.max() + interactionVec
        );

        // Find all elements intersecting box.
        labelList interactingElems
        (
            coupledPatchRangeTree.findBox(extendedBb)
        );

        if (!interactingElems.empty())
        {
            wallFaceBbRequiredByAnyCell[bbI] = true;
        }

        rwfil_[bbI].setSize(interactingElems.size(), -1);

        forAll(interactingElems, i)
        {
            label elemI = interactingElems[i];

            // Here, a more detailed geometric test could be applied,
            // i.e. a more accurate bounding volume like a OBB or
            // convex hull or an exact geometrical test.

            label c = coupledPatchRangeTree.shapes().cellLabels()[elemI];

            rwfil_[bbI][i] = c;
        }
    }

    // Perform subset of rwfil_, to remove any referred wallFaces that do
    // not interact.  They will not be sent from originating
    // processors.  This assumes that the disappearance of the wallFace
    // from the sending list of the source processor, simply removes
    // the referred wallFace from the rwfil_, all of the subsequent indices
    // shuffle down one, but the order and structure is preserved,
    // i.e. it, is as if the wallFace had never been referred in the first
    // place.

    inplaceSubset(wallFaceBbRequiredByAnyCell, rwfil_);

    // Send information about which wallFaces are actually required
    // back to originating processors.

    // At this point, wallFaceBbsToExchange does not need to be
    // maintained or distributed as it is not longer needed.

    wallFaceBbsToExchange.setSize(0);

    wallFaceMap().reverseDistribute
    (
        preDistributionWallFaceMapSize,
        wallFaceBbRequiredByAnyCell
    );

    wallFaceMap().reverseDistribute
    (
        preDistributionWallFaceMapSize,
        wallFaceIAndTToExchange
    );

    // Perform ordering of wallFaces to send, this invalidates the
    // previous value of preDistributionWallFaceMapSize.

    preDistributionWallFaceMapSize = -1;

    // Move all of the used wallFaces to refer to the start of the
    // list and truncate it

    inplaceSubset(wallFaceBbRequiredByAnyCell, wallFaceIAndTToExchange);

    inplaceSubset(wallFaceBbRequiredByAnyCell, procToDistributeWallFaceTo);

    preDistributionWallFaceMapSize = procToDistributeWallFaceTo.size();

    // Rebuild mapDistribute with only required referred wallFaces
    buildMap(wallFaceMapPtr_, procToDistributeWallFaceTo);

    // Store wallFaceIndexAndTransformToDistribute
    wallFaceIndexAndTransformToDistribute_.transfer(wallFaceIAndTToExchange);

    // Determine inverse addressing for referred wallFaces

    rwfilInverse_.setSize(mesh_.nCells());

    // Temporary Dynamic lists for accumulation
    List<DynamicList<label>> rwfilInverseTemp(rwfilInverse_.size());

    // Loop over all referred wall faces
    forAll(rwfil_, refWallFacei)
    {
        const labelList& realCells = rwfil_[refWallFacei];

        // Loop over all real cells in that the referred wall face is
        // to supply interactions to and record the index of this
        // referred wall face in the real cells entry in rwfilInverse

        forAll(realCells, realCelli)
        {
            rwfilInverseTemp[realCells[realCelli]].append(refWallFacei);
        }
    }

    forAll(rwfilInverse_, celli)
    {
        rwfilInverse_[celli].transfer(rwfilInverseTemp[celli]);
    }

    // Refer wall faces to the appropriate processor
    referredWallFaces_.setSize(wallFaceIndexAndTransformToDistribute_.size());

    forAll(referredWallFaces_, rWFI)
    {
        const labelPair& wfiat = wallFaceIndexAndTransformToDistribute_[rWFI];

        label wallFaceIndex = globalTransforms.index(wfiat);

        const vectorTensorTransform& transform = globalTransforms.transform
        (
            globalTransforms.transformIndex(wfiat)
        );

        const face& f = mesh_.faces()[wallFaceIndex];

        label patchi = mesh_.boundaryMesh().patchID()
        [
            wallFaceIndex - mesh_.nInternalFaces()
        ];

        referredWallFaces_[rWFI] = referredWallFace
        (
            face(identity(f.size())),
            transform.invTransformPosition(f.points(mesh_.points())),
            patchi
        );
    }

    wallFaceMap().distribute(referredWallFaces_);

    if (writeCloud_)
    {
        writeReferredWallFaces();
    }

    // Direct interaction list and direct wall faces

    Info<< "    Building direct interaction lists" << endl;

    indexedOctree<treeDataFace> wallFacesTree
    (
        treeDataFace(true, mesh_, localWallFaces),
        procBbRndExt,
        8,              // maxLevel,
        10,             // leafSize,
        100.0
    );

    dil_.setSize(mesh_.nCells());

    dwfil_.setSize(mesh_.nCells());

    forAll(cellBbs, celli)
    {
        const treeBoundBox& cellBb = cellBbs[celli];

        treeBoundBox extendedBb
        (
            cellBb.min() - interactionVec,
            cellBb.max() + interactionVec
        );

        // Find all cells intersecting extendedBb
        labelList interactingElems(allCellsTree.findBox(extendedBb));

        // Reserve space to avoid multiple resizing
        DynamicList<label> cellDIL(interactingElems.size());

        forAll(interactingElems, i)
        {
            label elemI = interactingElems[i];

            label c = allCellsTree.shapes().cellLabels()[elemI];

            // Here, a more detailed geometric test could be applied,
            // i.e. a more accurate bounding volume like a OBB or
            // convex hull, or an exact geometrical test.

            // The higher index cell is added to the lower index
            // cell's DIL.  A cell is not added to its own DIL.
            if (c > celli)
            {
                cellDIL.append(c);
            }
        }

        dil_[celli].transfer(cellDIL);

        // Find all wall faces intersecting extendedBb
        interactingElems = wallFacesTree.findBox(extendedBb);

        dwfil_[celli].setSize(interactingElems.size(), -1);

        forAll(interactingElems, i)
        {
            label elemI = interactingElems[i];

            label f = wallFacesTree.shapes().faceLabels()[elemI];

            dwfil_[celli][i] = f;
        }
    }
}


template<class ParticleType>
void Foam::InteractionLists<ParticleType>::findExtendedProcBbsInRange
(
    const treeBoundBox& procBb,
    const List<treeBoundBox>& allExtendedProcBbs,
    const globalIndexAndTransform& globalTransforms,
    List<treeBoundBox>& extendedProcBbsInRange,
    List<label>& extendedProcBbsTransformIndex,
    List<label>& extendedProcBbsOrigProc
)
{
    extendedProcBbsInRange.setSize(0);
    extendedProcBbsTransformIndex.setSize(0);
    extendedProcBbsOrigProc.setSize(0);

    DynamicList<treeBoundBox> tmpExtendedProcBbsInRange;
    DynamicList<label> tmpExtendedProcBbsTransformIndex;
    DynamicList<label> tmpExtendedProcBbsOrigProc;

    label nTrans = globalTransforms.nIndependentTransforms();

    forAll(allExtendedProcBbs, proci)
    {
        List<label> permutationIndices(nTrans, 0);

        if (nTrans == 0 && proci != Pstream::myProcNo())
        {
            treeBoundBox extendedReferredProcBb = allExtendedProcBbs[proci];

            if (procBb.overlaps(extendedReferredProcBb))
            {
                tmpExtendedProcBbsInRange.append(extendedReferredProcBb);

                // Dummy index, there are no transforms, so there will
                // be no resultant transform when this is decoded.
                tmpExtendedProcBbsTransformIndex.append(0);

                tmpExtendedProcBbsOrigProc.append(proci);
            }
        }
        else if (nTrans == 3)
        {
            label& i = permutationIndices[0];
            label& j = permutationIndices[1];
            label& k = permutationIndices[2];

            for (i = -1; i <= 1; i++)
            {
                for (j = -1; j <= 1; j++)
                {
                    for (k = -1; k <= 1; k++)
                    {
                        if
                        (
                            i == 0
                         && j == 0
                         && k == 0
                         && proci == Pstream::myProcNo()
                        )
                        {
                            // Skip this processor's extended boundBox
                            // when it has no transformation
                            continue;
                        }

                        label transI = globalTransforms.encodeTransformIndex
                        (
                            permutationIndices
                        );

                        const vectorTensorTransform& transform =
                            globalTransforms.transform(transI);

                        treeBoundBox extendedReferredProcBb
                        (
                            transform.transformPosition
                            (
                                allExtendedProcBbs[proci].points()
                            )
                        );

                        if (procBb.overlaps(extendedReferredProcBb))
                        {
                            tmpExtendedProcBbsInRange.append
                            (
                                extendedReferredProcBb
                            );

                            tmpExtendedProcBbsTransformIndex.append(transI);

                            tmpExtendedProcBbsOrigProc.append(proci);
                        }
                    }
                }
            }
        }
        else if (nTrans == 2)
        {
            label& i = permutationIndices[0];
            label& j = permutationIndices[1];

            for (i = -1; i <= 1; i++)
            {
                for (j = -1; j <= 1; j++)
                {
                    if (i == 0 && j == 0 && proci == Pstream::myProcNo())
                    {
                        // Skip this processor's extended boundBox
                        // when it has no transformation
                        continue;
                    }

                    label transI = globalTransforms.encodeTransformIndex
                    (
                        permutationIndices
                    );

                    const vectorTensorTransform& transform =
                        globalTransforms.transform(transI);

                    treeBoundBox extendedReferredProcBb
                    (
                        transform.transformPosition
                        (
                            allExtendedProcBbs[proci].points()
                        )
                    );

                    if (procBb.overlaps(extendedReferredProcBb))
                    {
                        tmpExtendedProcBbsInRange.append
                        (
                            extendedReferredProcBb
                        );

                        tmpExtendedProcBbsTransformIndex.append(transI);

                        tmpExtendedProcBbsOrigProc.append(proci);
                    }
                }
            }
        }
        else if (nTrans == 1)
        {
            label& i = permutationIndices[0];

            for (i = -1; i <= 1; i++)
            {
                if (i == 0 && proci == Pstream::myProcNo())
                {
                    // Skip this processor's extended boundBox when it
                    // has no transformation
                    continue;
                }

                label transI = globalTransforms.encodeTransformIndex
                (
                    permutationIndices
                );

                const vectorTensorTransform& transform =
                    globalTransforms.transform(transI);

                treeBoundBox extendedReferredProcBb
                (
                    transform.transformPosition
                    (
                        allExtendedProcBbs[proci].points()
                    )
                );

                if (procBb.overlaps(extendedReferredProcBb))
                {
                    tmpExtendedProcBbsInRange.append
                    (
                        extendedReferredProcBb
                    );

                    tmpExtendedProcBbsTransformIndex.append(transI);

                    tmpExtendedProcBbsOrigProc.append(proci);
                }
            }
        }
    }

    extendedProcBbsInRange = tmpExtendedProcBbsInRange.xfer();
    extendedProcBbsTransformIndex = tmpExtendedProcBbsTransformIndex.xfer();
    extendedProcBbsOrigProc = tmpExtendedProcBbsOrigProc.xfer();
}


template<class ParticleType>
void Foam::InteractionLists<ParticleType>::buildMap
(
    autoPtr<mapDistribute>& mapPtr,
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

    mapPtr.reset
    (
        new mapDistribute
        (
            constructSize,
            sendMap.xfer(),
            constructMap.xfer()
        )
    );
}


template<class ParticleType>
void Foam::InteractionLists<ParticleType>::prepareParticlesToRefer
(
    const List<DynamicList<ParticleType*>>& cellOccupancy
)
{
    const globalIndexAndTransform& globalTransforms =
        mesh_.globalData().globalTransforms();

    referredParticles_.setSize(cellIndexAndTransformToDistribute_.size());

    // Clear all existing referred particles

    forAll(referredParticles_, i)
    {
        referredParticles_[i].clear();
    }

    // Clear all particles that may have been populated into the cloud
    cloud_.clear();

    forAll(cellIndexAndTransformToDistribute_, i)
    {
        const labelPair ciat = cellIndexAndTransformToDistribute_[i];

        label cellIndex = globalTransforms.index(ciat);

        List<ParticleType*> realParticles = cellOccupancy[cellIndex];

        IDLList<ParticleType>& particlesToRefer = referredParticles_[i];

        forAll(realParticles, rM)
        {
            const ParticleType& particle = *realParticles[rM];

            particlesToRefer.append(particle.clone().ptr());

            prepareParticleToBeReferred(particlesToRefer.last(), ciat);
        }
    }
}


template<class ParticleType>
void Foam::InteractionLists<ParticleType>::prepareParticleToBeReferred
(
    ParticleType* particle,
    labelPair ciat
)
{
    const globalIndexAndTransform& globalTransforms =
        mesh_.globalData().globalTransforms();

    const vectorTensorTransform& transform = globalTransforms.transform
    (
        globalTransforms.transformIndex(ciat)
    );

    particle->prepareForInteractionListReferral(transform);
}


template<class ParticleType>
void Foam::InteractionLists<ParticleType>::fillReferredParticleCloud()
{
    if (writeCloud_)
    {
        forAll(referredParticles_, refCelli)
        {
            const IDLList<ParticleType>& refCell =
                referredParticles_[refCelli];

            forAllConstIter(typename IDLList<ParticleType>, refCell, iter)
            {
                cloud_.addParticle
                (
                    static_cast<ParticleType*>(iter().clone().ptr())
                );
            }
        }
    }
}


template<class ParticleType>
void Foam::InteractionLists<ParticleType>::prepareWallDataToRefer()
{
    const globalIndexAndTransform& globalTransforms =
        mesh_.globalData().globalTransforms();

    referredWallData_.setSize
    (
        wallFaceIndexAndTransformToDistribute_.size()
    );

    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);

    forAll(referredWallData_, rWVI)
    {
        const labelPair& wfiat = wallFaceIndexAndTransformToDistribute_[rWVI];

        label wallFaceIndex = globalTransforms.index(wfiat);

        const vectorTensorTransform& transform = globalTransforms.transform
        (
            globalTransforms.transformIndex(wfiat)
        );

        label patchi = mesh_.boundaryMesh().patchID()
        [
            wallFaceIndex - mesh_.nInternalFaces()
        ];

        label patchFacei =
            wallFaceIndex
          - mesh_.boundaryMesh()[patchi].start();

        // Need to transform velocity when tensor transforms are
        // supported
        referredWallData_[rWVI] = U.boundaryField()[patchi][patchFacei];

        if (transform.hasR())
        {
            referredWallData_[rWVI] =
                transform.R().T() & referredWallData_[rWVI];
        }
    }
}


template<class ParticleType>
void Foam::InteractionLists<ParticleType>::writeReferredWallFaces() const
{
    if (referredWallFaces_.empty())
    {
        return;
    }

    fileName objDir = mesh_.time().timePath()/cloud::prefix;

    mkDir(objDir);

    fileName objFileName = "referredWallFaces.obj";

    OFstream str(objDir/objFileName);

    Info<< "    Writing "
        << mesh_.time().timeName()/cloud::prefix/objFileName
        << endl;

    label offset = 1;

    forAll(referredWallFaces_, rWFI)
    {
        const referredWallFace& rwf = referredWallFaces_[rWFI];

        forAll(rwf, fPtI)
        {
            meshTools::writeOBJ(str, rwf.points()[rwf[fPtI]]);
        }

        str<< 'f';

        forAll(rwf, fPtI)
        {
            str<< ' ' << fPtI + offset;
        }

        str<< nl;

        offset += rwf.size();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::InteractionLists<ParticleType>::InteractionLists(const polyMesh& mesh)
:
    mesh_(mesh),
    cloud_(mesh_, "nullptr_Cloud", IDLList<ParticleType>()),
    writeCloud_(false),
    cellMapPtr_(),
    wallFaceMapPtr_(),
    maxDistance_(0.0),
    dil_(),
    dwfil_(),
    ril_(),
    rilInverse_(),
    cellIndexAndTransformToDistribute_(),
    wallFaceIndexAndTransformToDistribute_(),
    referredWallFaces_(),
    UName_("unknown_U"),
    referredWallData_(),
    referredParticles_()
{}


template<class ParticleType>
Foam::InteractionLists<ParticleType>::InteractionLists
(
    const polyMesh& mesh,
    scalar maxDistance,
    Switch writeCloud,
    const word& UName
)
:
    mesh_(mesh),
    cloud_(mesh_, "referredParticleCloud", IDLList<ParticleType>()),
    writeCloud_(writeCloud),
    cellMapPtr_(),
    wallFaceMapPtr_(),
    maxDistance_(maxDistance),
    dil_(),
    dwfil_(),
    ril_(),
    rilInverse_(),
    cellIndexAndTransformToDistribute_(),
    wallFaceIndexAndTransformToDistribute_(),
    referredWallFaces_(),
    UName_(UName),
    referredWallData_(),
    referredParticles_()
{
    buildInteractionLists();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::InteractionLists<ParticleType>::~InteractionLists()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ParticleType>
void Foam::InteractionLists<ParticleType>::sendReferredData
(
    const List<DynamicList<ParticleType*>>& cellOccupancy,
    PstreamBuffers& pBufs
)
{
    if (mesh_.changing())
    {
        WarningInFunction
            << "Mesh changing, rebuilding InteractionLists form scratch."
            << endl;

        buildInteractionLists();
    }

    prepareWallDataToRefer();

    prepareParticlesToRefer(cellOccupancy);

    for (label domain = 0; domain < Pstream::nProcs(); domain++)
    {
        const labelList& subMap = cellMap().subMap()[domain];

        if (subMap.size())
        {
            UOPstream toDomain(domain, pBufs);

            UIndirectList<IDLList<ParticleType>> subMappedParticles
            (
                referredParticles_,
                subMap
            );

            forAll(subMappedParticles, i)
            {
                toDomain << subMappedParticles[i];
            }
        }
    }

    // Using the mapDistribute to start sending and receiving the
    // buffer but not block, i.e. it is calling
    //     pBufs.finishedSends(false);
    wallFaceMap().send(pBufs, referredWallData_);
}


template<class ParticleType>
void Foam::InteractionLists<ParticleType>::receiveReferredData
(
    PstreamBuffers& pBufs,
    const label startOfRequests
)
{
    Pstream::waitRequests(startOfRequests);

    referredParticles_.setSize(cellMap().constructSize());

    for (label domain = 0; domain < Pstream::nProcs(); domain++)
    {
        const labelList& constructMap = cellMap().constructMap()[domain];

        if (constructMap.size())
        {
            UIPstream str(domain, pBufs);

            forAll(constructMap, i)
            {
                referredParticles_[constructMap[i]] = IDLList<ParticleType>
                (
                    str,
                    typename ParticleType::iNew(mesh_)
                );
            }
        }
    }

    forAll(referredParticles_, refCelli)
    {
        IDLList<ParticleType>& refCell = referredParticles_[refCelli];
        forAllIter(typename IDLList<ParticleType>, refCell, iter)
        {
            iter().correctAfterInteractionListReferral(ril_[refCelli][0]);
        }
    }

    fillReferredParticleCloud();

    wallFaceMap().receive(pBufs, referredWallData_);
}


// ************************************************************************* //
