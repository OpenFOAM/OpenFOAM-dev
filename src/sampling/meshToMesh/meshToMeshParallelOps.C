/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "meshToMesh.H"
#include "OFstream.H"
#include "Time.H"
#include "globalIndex.H"
#include "mergePoints.H"
#include "processorPolyPatch.H"
#include "SubField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::meshToMesh::calcDistribution
(
    const polyMesh& src,
    const polyMesh& tgt
) const
{
    label proci = 0;

    if (Pstream::parRun())
    {
        List<label> cellsPresentOnProc(Pstream::nProcs(), 0);
        if ((src.nCells() > 0) || (tgt.nCells() > 0))
        {
            cellsPresentOnProc[Pstream::myProcNo()] = 1;
        }
        else
        {
            cellsPresentOnProc[Pstream::myProcNo()] = 0;
        }

        Pstream::gatherList(cellsPresentOnProc);
        Pstream::scatterList(cellsPresentOnProc);

        label nHaveCells = sum(cellsPresentOnProc);

        if (nHaveCells > 1)
        {
            proci = -1;
            if (debug)
            {
                InfoInFunction
                    << "Meshes split across multiple processors" << endl;
            }
        }
        else if (nHaveCells == 1)
        {
            proci = findIndex(cellsPresentOnProc, 1);
            if (debug)
            {
                InfoInFunction
                    << "Meshes local to processor" << proci << endl;
            }
        }
    }

    return proci;
}


Foam::label Foam::meshToMesh::calcOverlappingProcs
(
    const List<boundBox>& procBb,
    const boundBox& bb,
    boolList& overlaps
) const
{
    overlaps = false;

    label nOverlaps = 0;

    forAll(procBb, proci)
    {
        const boundBox& bbp = procBb[proci];

        if (bbp.overlaps(bb))
        {
            overlaps[proci] = true;
            nOverlaps++;
        }
    }

    return nOverlaps;
}


Foam::autoPtr<Foam::mapDistribute> Foam::meshToMesh::calcProcMap
(
    const polyMesh& src,
    const polyMesh& tgt
) const
{
    // get decomposition of cells on src mesh
    List<boundBox> procBb(Pstream::nProcs());

    if (src.nCells() > 0)
    {
        // bounding box for my mesh - do not parallel reduce
        procBb[Pstream::myProcNo()] = boundBox(src.points(), false);

        // slightly increase size of bounding boxes to allow for cases where
        // bounding boxes are perfectly aligned
        procBb[Pstream::myProcNo()].inflate(0.01);
    }
    else
    {
        procBb[Pstream::myProcNo()] = boundBox();
    }


    Pstream::gatherList(procBb);
    Pstream::scatterList(procBb);


    if (debug)
    {
        InfoInFunction
            << "Determining extent of src mesh per processor:" << nl
            << "\tproc\tbb" << endl;
        forAll(procBb, proci)
        {
            Info<< '\t' << proci << '\t' << procBb[proci] << endl;
        }
    }


    // determine which cells of tgt mesh overlaps src mesh per proc
    const cellList& cells = tgt.cells();
    const faceList& faces = tgt.faces();
    const pointField& points = tgt.points();

    labelListList sendMap;

    {
        // per processor indices into all segments to send
        List<DynamicList<label>> dynSendMap(Pstream::nProcs());
        label iniSize = floor(tgt.nCells()/Pstream::nProcs());

        forAll(dynSendMap, proci)
        {
            dynSendMap[proci].setCapacity(iniSize);
        }

        // work array - whether src processor bb overlaps the tgt cell bounds
        boolList procBbOverlaps(Pstream::nProcs());
        forAll(cells, celli)
        {
            const cell& c = cells[celli];

            // determine bounding box of tgt cell
            boundBox cellBb(point::max, point::min);
            forAll(c, facei)
            {
                const face& f = faces[c[facei]];
                forAll(f, fp)
                {
                    cellBb.min() = min(cellBb.min(), points[f[fp]]);
                    cellBb.max() = max(cellBb.max(), points[f[fp]]);
                }
            }

            // find the overlapping tgt cells on each src processor
            (void)calcOverlappingProcs(procBb, cellBb, procBbOverlaps);

            forAll(procBbOverlaps, proci)
            {
                if (procBbOverlaps[proci])
                {
                    dynSendMap[proci].append(celli);
                }
            }
        }

        // convert dynamicList to labelList
        sendMap.setSize(Pstream::nProcs());
        forAll(sendMap, proci)
        {
            sendMap[proci].transfer(dynSendMap[proci]);
        }
    }

    // debug printing
    if (debug)
    {
        Pout<< "Of my " << cells.size() << " target cells I need to send to:"
            << nl << "\tproc\tcells" << endl;
        forAll(sendMap, proci)
        {
            Pout<< '\t' << proci << '\t' << sendMap[proci].size() << endl;
        }
    }


    // send over how many tgt cells I need to receive from each processor
    labelListList sendSizes(Pstream::nProcs());
    sendSizes[Pstream::myProcNo()].setSize(Pstream::nProcs());
    forAll(sendMap, proci)
    {
        sendSizes[Pstream::myProcNo()][proci] = sendMap[proci].size();
    }
    Pstream::gatherList(sendSizes);
    Pstream::scatterList(sendSizes);


    // determine order of receiving
    labelListList constructMap(Pstream::nProcs());

    label segmentI = 0;
    forAll(constructMap, proci)
    {
        // what I need to receive is what other processor is sending to me
        label nRecv = sendSizes[proci][Pstream::myProcNo()];
        constructMap[proci].setSize(nRecv);

        for (label i = 0; i < nRecv; i++)
        {
            constructMap[proci][i] = segmentI++;
        }
    }

    autoPtr<mapDistribute> mapPtr
    (
        new mapDistribute
        (
            segmentI,       // size after construction
            sendMap.xfer(),
            constructMap.xfer()
        )
    );

    return mapPtr;
}


void Foam::meshToMesh::distributeCells
(
    const mapDistribute& map,
    const polyMesh& tgtMesh,
    const globalIndex& globalI,
    List<pointField>& points,
    List<label>& nInternalFaces,
    List<faceList>& faces,
    List<labelList>& faceOwner,
    List<labelList>& faceNeighbour,
    List<labelList>& cellIDs,
    List<labelList>& nbrProcIDs,
    List<labelList>& procLocalFaceIDs
) const
{
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    points.setSize(Pstream::nProcs());
    nInternalFaces.setSize(Pstream::nProcs(), 0);
    faces.setSize(Pstream::nProcs());
    faceOwner.setSize(Pstream::nProcs());
    faceNeighbour.setSize(Pstream::nProcs());
    cellIDs.setSize(Pstream::nProcs());

    nbrProcIDs.setSize(Pstream::nProcs());;
    procLocalFaceIDs.setSize(Pstream::nProcs());;


    for (label domain = 0; domain < Pstream::nProcs(); domain++)
    {
        const labelList& sendElems = map.subMap()[domain];

        if (sendElems.size())
        {
            // reverse cell map
            labelList reverseCellMap(tgtMesh.nCells(), -1);
            forAll(sendElems, subCelli)
            {
                reverseCellMap[sendElems[subCelli]] = subCelli;
            }

            DynamicList<face> subFaces(tgtMesh.nFaces());
            DynamicList<label> subFaceOwner(tgtMesh.nFaces());
            DynamicList<label> subFaceNeighbour(tgtMesh.nFaces());

            DynamicList<label> subNbrProcIDs(tgtMesh.nFaces());
            DynamicList<label> subProcLocalFaceIDs(tgtMesh.nFaces());

            label nInternal = 0;

            // internal faces
            forAll(tgtMesh.faceNeighbour(), facei)
            {
                label own = tgtMesh.faceOwner()[facei];
                label nbr = tgtMesh.faceNeighbour()[facei];
                label subOwn = reverseCellMap[own];
                label subNbr = reverseCellMap[nbr];

                if (subOwn != -1 && subNbr != -1)
                {
                    nInternal++;

                    if (subOwn < subNbr)
                    {
                        subFaces.append(tgtMesh.faces()[facei]);
                        subFaceOwner.append(subOwn);
                        subFaceNeighbour.append(subNbr);
                        subNbrProcIDs.append(-1);
                        subProcLocalFaceIDs.append(-1);
                    }
                    else
                    {
                        subFaces.append(tgtMesh.faces()[facei].reverseFace());
                        subFaceOwner.append(subNbr);
                        subFaceNeighbour.append(subOwn);
                        subNbrProcIDs.append(-1);
                        subProcLocalFaceIDs.append(-1);
                    }
                }
            }

            // boundary faces for new region
            forAll(tgtMesh.faceNeighbour(), facei)
            {
                label own = tgtMesh.faceOwner()[facei];
                label nbr = tgtMesh.faceNeighbour()[facei];
                label subOwn = reverseCellMap[own];
                label subNbr = reverseCellMap[nbr];

                if (subOwn != -1 && subNbr == -1)
                {
                    subFaces.append(tgtMesh.faces()[facei]);
                    subFaceOwner.append(subOwn);
                    subFaceNeighbour.append(subNbr);
                    subNbrProcIDs.append(-1);
                    subProcLocalFaceIDs.append(-1);
                }
                else if (subOwn == -1 && subNbr != -1)
                {
                    subFaces.append(tgtMesh.faces()[facei].reverseFace());
                    subFaceOwner.append(subNbr);
                    subFaceNeighbour.append(subOwn);
                    subNbrProcIDs.append(-1);
                    subProcLocalFaceIDs.append(-1);
                }
            }

            // boundary faces of existing region
            forAll(tgtMesh.boundaryMesh(), patchi)
            {
                const polyPatch& pp = tgtMesh.boundaryMesh()[patchi];

                label nbrProci = -1;

                // store info for faces on processor patches
                if (isA<processorPolyPatch>(pp))
                {
                    const processorPolyPatch& ppp =
                        dynamic_cast<const processorPolyPatch&>(pp);

                    nbrProci = ppp.neighbProcNo();
                }

                forAll(pp, i)
                {
                    label facei = pp.start() + i;
                    label own = tgtMesh.faceOwner()[facei];

                    if (reverseCellMap[own] != -1)
                    {
                        subFaces.append(tgtMesh.faces()[facei]);
                        subFaceOwner.append(reverseCellMap[own]);
                        subFaceNeighbour.append(-1);
                        subNbrProcIDs.append(nbrProci);
                        subProcLocalFaceIDs.append(i);
                    }
                }
            }

            // reverse point map
            labelList reversePointMap(tgtMesh.nPoints(), -1);
            DynamicList<point> subPoints(tgtMesh.nPoints());
            forAll(subFaces, subFacei)
            {
                face& f = subFaces[subFacei];
                forAll(f, fp)
                {
                    label pointi = f[fp];
                    if (reversePointMap[pointi] == -1)
                    {
                        reversePointMap[pointi] = subPoints.size();
                        subPoints.append(tgtMesh.points()[pointi]);
                    }

                    f[fp] = reversePointMap[pointi];
                }
            }

            // tgt cells into global numbering
            labelList globalElems(sendElems.size());
            forAll(sendElems, i)
            {
                if (debug)
                {
                    Pout<< "tgtProc:" << Pstream::myProcNo()
                        << " sending tgt cell " << sendElems[i]
                        << "[" << globalI.toGlobal(sendElems[i]) << "]"
                        << " to srcProc " << domain << endl;
                }

                globalElems[i] = globalI.toGlobal(sendElems[i]);
            }

            // pass data
            if (domain == Pstream::myProcNo())
            {
                // allocate my own data
                points[Pstream::myProcNo()] = subPoints;
                nInternalFaces[Pstream::myProcNo()] = nInternal;
                faces[Pstream::myProcNo()] = subFaces;
                faceOwner[Pstream::myProcNo()] = subFaceOwner;
                faceNeighbour[Pstream::myProcNo()] = subFaceNeighbour;
                cellIDs[Pstream::myProcNo()] = globalElems;
                nbrProcIDs[Pstream::myProcNo()] = subNbrProcIDs;
                procLocalFaceIDs[Pstream::myProcNo()] = subProcLocalFaceIDs;
            }
            else
            {
                // send data to other processor domains
                UOPstream toDomain(domain, pBufs);

                toDomain
                    << subPoints
                    << nInternal
                    << subFaces
                    << subFaceOwner
                    << subFaceNeighbour
                    << globalElems
                    << subNbrProcIDs
                    << subProcLocalFaceIDs;
            }
        }
    }

    // Start receiving
    pBufs.finishedSends();

    // Consume
    for (label domain = 0; domain < Pstream::nProcs(); domain++)
    {
        const labelList& recvElems = map.constructMap()[domain];

        if (domain != Pstream::myProcNo() && recvElems.size())
        {
            UIPstream str(domain, pBufs);

            str >> points[domain]
                >> nInternalFaces[domain]
                >> faces[domain]
                >> faceOwner[domain]
                >> faceNeighbour[domain]
                >> cellIDs[domain]
                >> nbrProcIDs[domain]
                >> procLocalFaceIDs[domain];
        }

        if (debug)
        {
            Pout<< "Target mesh send sizes[" << domain << "]"
                << ": points="<< points[domain].size()
                << ", faces=" << faces[domain].size()
                << ", nInternalFaces=" << nInternalFaces[domain]
                << ", faceOwn=" << faceOwner[domain].size()
                << ", faceNbr=" << faceNeighbour[domain].size()
                << ", cellIDs=" << cellIDs[domain].size() << endl;
        }
    }
}


void Foam::meshToMesh::distributeAndMergeCells
(
    const mapDistribute& map,
    const polyMesh& tgt,
    const globalIndex& globalI,
    pointField& tgtPoints,
    faceList& tgtFaces,
    labelList& tgtFaceOwners,
    labelList& tgtFaceNeighbours,
    labelList& tgtCellIDs
) const
{
    // Exchange per-processor data
    List<pointField> allPoints;
    List<label> allNInternalFaces;
    List<faceList> allFaces;
    List<labelList> allFaceOwners;
    List<labelList> allFaceNeighbours;
    List<labelList> allTgtCellIDs;

    // Per target mesh face the neighbouring proc and index in
    // processor patch (all -1 for normal boundary face)
    List<labelList> allNbrProcIDs;
    List<labelList> allProcLocalFaceIDs;

    distributeCells
    (
        map,
        tgt,
        globalI,
        allPoints,
        allNInternalFaces,
        allFaces,
        allFaceOwners,
        allFaceNeighbours,
        allTgtCellIDs,
        allNbrProcIDs,
        allProcLocalFaceIDs
    );

    // Convert lists into format that can be used to generate a valid polyMesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // Points and cells are collected into single flat lists:
    // - i.e. proc0, proc1 ... procN
    //
    // Faces need to be sorted after collection to that internal faces are
    // contiguous, followed by all boundary faces
    //
    // Processor patch faces between included cells on neighbouring processors
    // are converted into internal faces
    //
    // Face list structure:
    // - Per processor:
    //   - internal faces
    //   - processor faces that have been converted into internal faces
    // - Followed by all boundary faces
    //   - from 'normal' boundary faces
    //   - from singularly-sided processor patch faces


    // Number of internal+coupled faces
    labelList allNIntCoupledFaces(allNInternalFaces);

    // Starting offset for points
    label nPoints = 0;
    labelList pointOffset(Pstream::nProcs(), 0);
    forAll(allPoints, proci)
    {
        pointOffset[proci] = nPoints;
        nPoints += allPoints[proci].size();
    }

    // Starting offset for cells
    label nCells = 0;
    labelList cellOffset(Pstream::nProcs(), 0);
    forAll(allTgtCellIDs, proci)
    {
        cellOffset[proci] = nCells;
        nCells += allTgtCellIDs[proci].size();
    }

    // Count any coupled faces
    typedef FixedList<label, 3> label3;
    typedef HashTable<label, label3, label3::Hash<>> procCoupleInfo;
    procCoupleInfo procFaceToGlobalCell;

    forAll(allNbrProcIDs, proci)
    {
        const labelList& nbrProci = allNbrProcIDs[proci];
        const labelList& localFacei = allProcLocalFaceIDs[proci];

        forAll(nbrProci, i)
        {
            if (nbrProci[i] != -1 && localFacei[i] != -1)
            {
                label3 key;
                key[0] = min(proci, nbrProci[i]);
                key[1] = max(proci, nbrProci[i]);
                key[2] = localFacei[i];

                procCoupleInfo::const_iterator fnd =
                    procFaceToGlobalCell.find(key);

                if (fnd == procFaceToGlobalCell.end())
                {
                    procFaceToGlobalCell.insert(key, -1);
                }
                else
                {
                    if (debug)
                    {
                        Pout<< "Additional internal face between procs:"
                            << key[0] << " and " << key[1]
                            << " across local face " << key[2] << endl;
                    }

                    allNIntCoupledFaces[proci]++;
                }
            }
        }
    }


    // Starting offset for internal faces
    label nIntFaces = 0;
    label nFacesTotal = 0;
    labelList internalFaceOffset(Pstream::nProcs(), 0);
    forAll(allNIntCoupledFaces, proci)
    {
        label nCoupledFaces =
            allNIntCoupledFaces[proci] - allNInternalFaces[proci];

        internalFaceOffset[proci] = nIntFaces;
        nIntFaces += allNIntCoupledFaces[proci];
        nFacesTotal += allFaceOwners[proci].size() - nCoupledFaces;
    }

    tgtPoints.setSize(nPoints);
    tgtFaces.setSize(nFacesTotal);
    tgtFaceOwners.setSize(nFacesTotal);
    tgtFaceNeighbours.setSize(nFacesTotal);
    tgtCellIDs.setSize(nCells);

    // Insert points
    forAll(allPoints, proci)
    {
        const pointField& pts = allPoints[proci];
        SubList<point>(tgtPoints, pts.size(), pointOffset[proci]) = pts;
    }

    // Insert cellIDs
    forAll(allTgtCellIDs, proci)
    {
        const labelList& cellIDs = allTgtCellIDs[proci];
        SubList<label>(tgtCellIDs, cellIDs.size(), cellOffset[proci]) = cellIDs;
    }


    // Insert internal faces (from internal faces)
    forAll(allFaces, proci)
    {
        const faceList& fcs = allFaces[proci];
        const labelList& faceOs = allFaceOwners[proci];
        const labelList& faceNs = allFaceNeighbours[proci];

        SubList<face> slice
        (
            tgtFaces,
            allNInternalFaces[proci],
            internalFaceOffset[proci]
        );
        slice = SubList<face>(fcs, allNInternalFaces[proci]);
        forAll(slice, i)
        {
            add(slice[i], pointOffset[proci]);
        }

        SubField<label> ownSlice
        (
            tgtFaceOwners,
            allNInternalFaces[proci],
            internalFaceOffset[proci]
        );
        ownSlice = SubField<label>(faceOs, allNInternalFaces[proci]);
        add(ownSlice, cellOffset[proci]);

        SubField<label> nbrSlice
        (
            tgtFaceNeighbours,
            allNInternalFaces[proci],
            internalFaceOffset[proci]
        );
        nbrSlice = SubField<label>(faceNs, allNInternalFaces[proci]);
        add(nbrSlice, cellOffset[proci]);

        internalFaceOffset[proci] += allNInternalFaces[proci];
    }


    // Insert internal faces (from coupled face-pairs)
    forAll(allNbrProcIDs, proci)
    {
        const labelList& nbrProci = allNbrProcIDs[proci];
        const labelList& localFacei = allProcLocalFaceIDs[proci];
        const labelList& faceOs = allFaceOwners[proci];
        const faceList& fcs = allFaces[proci];

        forAll(nbrProci, i)
        {
            if (nbrProci[i] != -1 && localFacei[i] != -1)
            {
                label3 key;
                key[0] = min(proci, nbrProci[i]);
                key[1] = max(proci, nbrProci[i]);
                key[2] = localFacei[i];

                procCoupleInfo::iterator fnd = procFaceToGlobalCell.find(key);

                if (fnd != procFaceToGlobalCell.end())
                {
                    label tgtFacei = fnd();
                    if (tgtFacei == -1)
                    {
                        // on first visit store the new cell on this side
                        fnd() = cellOffset[proci] + faceOs[i];
                    }
                    else
                    {
                        // get owner and neighbour in new cell numbering
                        label newOwn = cellOffset[proci] + faceOs[i];
                        label newNbr = fnd();
                        label tgtFacei = internalFaceOffset[proci]++;

                        if (debug)
                        {
                            Pout<< "    proc " << proci
                                << "\tinserting face:" << tgtFacei
                                << " connection between owner " << newOwn
                                << " and neighbour " << newNbr
                                << endl;
                        }

                        if (newOwn < newNbr)
                        {
                            // we have correct orientation
                            tgtFaces[tgtFacei] = fcs[i];
                            tgtFaceOwners[tgtFacei] = newOwn;
                            tgtFaceNeighbours[tgtFacei] = newNbr;
                        }
                        else
                        {
                            // reverse orientation
                            tgtFaces[tgtFacei] = fcs[i].reverseFace();
                            tgtFaceOwners[tgtFacei] = newNbr;
                            tgtFaceNeighbours[tgtFacei] = newOwn;
                        }

                        add(tgtFaces[tgtFacei], pointOffset[proci]);

                        // mark with unique value
                        fnd() = -2;
                    }
                }
            }
        }
    }


    forAll(allNbrProcIDs, proci)
    {
        const labelList& nbrProci = allNbrProcIDs[proci];
        const labelList& localFacei = allProcLocalFaceIDs[proci];
        const labelList& faceOs = allFaceOwners[proci];
        const labelList& faceNs = allFaceNeighbours[proci];
        const faceList& fcs = allFaces[proci];

        forAll(nbrProci, i)
        {
            // coupled boundary face
            if (nbrProci[i] != -1 && localFacei[i] != -1)
            {
                label3 key;
                key[0] = min(proci, nbrProci[i]);
                key[1] = max(proci, nbrProci[i]);
                key[2] = localFacei[i];

                label tgtFacei = procFaceToGlobalCell[key];

                if (tgtFacei == -1)
                {
                    FatalErrorInFunction
                        << "Unvisited " << key
                        << abort(FatalError);
                }
                else if (tgtFacei != -2)
                {
                    label newOwn = cellOffset[proci] + faceOs[i];
                    label tgtFacei = nIntFaces++;

                    if (debug)
                    {
                        Pout<< "    proc " << proci
                            << "\tinserting boundary face:" << tgtFacei
                            << " from coupled face " << key
                            << endl;
                    }

                    tgtFaces[tgtFacei] = fcs[i];
                    add(tgtFaces[tgtFacei], pointOffset[proci]);

                    tgtFaceOwners[tgtFacei] = newOwn;
                    tgtFaceNeighbours[tgtFacei] = -1;
                }
            }
            // normal boundary face
            else
            {
                label own = faceOs[i];
                label nbr = faceNs[i];
                if ((own != -1) && (nbr == -1))
                {
                    label newOwn = cellOffset[proci] + faceOs[i];
                    label tgtFacei = nIntFaces++;

                    tgtFaces[tgtFacei] = fcs[i];
                    add(tgtFaces[tgtFacei], pointOffset[proci]);

                    tgtFaceOwners[tgtFacei] = newOwn;
                    tgtFaceNeighbours[tgtFacei] = -1;
                }
            }
        }
    }


    if (debug)
    {
        // only merging points in debug mode

        labelList oldToNew;
        pointField newTgtPoints;
        bool hasMerged = mergePoints
        (
            tgtPoints,
            small,
            false,
            oldToNew,
            newTgtPoints
        );

        if (hasMerged)
        {
            if (debug)
            {
                Pout<< "Merged from " << tgtPoints.size()
                    << " down to " << newTgtPoints.size() << " points" << endl;
            }

            tgtPoints.transfer(newTgtPoints);
            forAll(tgtFaces, i)
            {
                inplaceRenumber(oldToNew, tgtFaces[i]);
            }
        }
    }
}


// ************************************************************************* //
