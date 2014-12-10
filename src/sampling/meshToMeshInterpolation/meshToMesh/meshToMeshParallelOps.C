/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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
    label procI = 0;

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
            procI = -1;
            if (debug)
            {
                Info<< "meshToMesh::calcDistribution: "
                    << "Meshes split across multiple processors" << endl;
            }
        }
        else if (nHaveCells == 1)
        {
            procI = findIndex(cellsPresentOnProc, 1);
            if (debug)
            {
                Info<< "meshToMesh::calcDistribution: "
                    << "Meshes local to processor" << procI << endl;
            }
        }
    }

    return procI;
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

    forAll(procBb, procI)
    {
        const boundBox& bbp = procBb[procI];

        if (bbp.overlaps(bb))
        {
            overlaps[procI] = true;
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
        // bounding boxes are perfectly alligned
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
        Info<< "Determining extent of src mesh per processor:" << nl
            << "\tproc\tbb" << endl;
        forAll(procBb, procI)
        {
            Info<< '\t' << procI << '\t' << procBb[procI] << endl;
        }
    }


    // determine which cells of tgt mesh overlaps src mesh per proc
    const cellList& cells = tgt.cells();
    const faceList& faces = tgt.faces();
    const pointField& points = tgt.points();

    labelListList sendMap;

    {
        // per processor indices into all segments to send
        List<DynamicList<label> > dynSendMap(Pstream::nProcs());
        label iniSize = floor(tgt.nCells()/Pstream::nProcs());

        forAll(dynSendMap, procI)
        {
            dynSendMap[procI].setCapacity(iniSize);
        }

        // work array - whether src processor bb overlaps the tgt cell bounds
        boolList procBbOverlaps(Pstream::nProcs());
        forAll(cells, cellI)
        {
            const cell& c = cells[cellI];

            // determine bounding box of tgt cell
            boundBox cellBb(point::max, point::min);
            forAll(c, faceI)
            {
                const face& f = faces[c[faceI]];
                forAll(f, fp)
                {
                    cellBb.min() = min(cellBb.min(), points[f[fp]]);
                    cellBb.max() = max(cellBb.max(), points[f[fp]]);
                }
            }

            // find the overlapping tgt cells on each src processor
            (void)calcOverlappingProcs(procBb, cellBb, procBbOverlaps);

            forAll(procBbOverlaps, procI)
            {
                if (procBbOverlaps[procI])
                {
                    dynSendMap[procI].append(cellI);
                }
            }
        }

        // convert dynamicList to labelList
        sendMap.setSize(Pstream::nProcs());
        forAll(sendMap, procI)
        {
            sendMap[procI].transfer(dynSendMap[procI]);
        }
    }

    // debug printing
    if (debug)
    {
        Pout<< "Of my " << cells.size() << " target cells I need to send to:"
            << nl << "\tproc\tcells" << endl;
        forAll(sendMap, procI)
        {
            Pout<< '\t' << procI << '\t' << sendMap[procI].size() << endl;
        }
    }


    // send over how many tgt cells I need to receive from each processor
    labelListList sendSizes(Pstream::nProcs());
    sendSizes[Pstream::myProcNo()].setSize(Pstream::nProcs());
    forAll(sendMap, procI)
    {
        sendSizes[Pstream::myProcNo()][procI] = sendMap[procI].size();
    }
    Pstream::gatherList(sendSizes);
    Pstream::scatterList(sendSizes);


    // determine order of receiving
    labelListList constructMap(Pstream::nProcs());

    label segmentI = 0;
    forAll(constructMap, procI)
    {
        // what I need to receive is what other processor is sending to me
        label nRecv = sendSizes[procI][Pstream::myProcNo()];
        constructMap[procI].setSize(nRecv);

        for (label i = 0; i < nRecv; i++)
        {
            constructMap[procI][i] = segmentI++;
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
    PstreamBuffers pBufs(Pstream::nonBlocking);

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
            forAll(sendElems, subCellI)
            {
                reverseCellMap[sendElems[subCellI]] = subCellI;
            }

            DynamicList<face> subFaces(tgtMesh.nFaces());
            DynamicList<label> subFaceOwner(tgtMesh.nFaces());
            DynamicList<label> subFaceNeighbour(tgtMesh.nFaces());

            DynamicList<label> subNbrProcIDs(tgtMesh.nFaces());
            DynamicList<label> subProcLocalFaceIDs(tgtMesh.nFaces());

            label nInternal = 0;

            // internal faces
            forAll(tgtMesh.faceNeighbour(), faceI)
            {
                label own = tgtMesh.faceOwner()[faceI];
                label nbr = tgtMesh.faceNeighbour()[faceI];
                label subOwn = reverseCellMap[own];
                label subNbr = reverseCellMap[nbr];

                if (subOwn != -1 && subNbr != -1)
                {
                    nInternal++;

                    if (subOwn < subNbr)
                    {
                        subFaces.append(tgtMesh.faces()[faceI]);
                        subFaceOwner.append(subOwn);
                        subFaceNeighbour.append(subNbr);
                        subNbrProcIDs.append(-1);
                        subProcLocalFaceIDs.append(-1);
                    }
                    else
                    {
                        subFaces.append(tgtMesh.faces()[faceI].reverseFace());
                        subFaceOwner.append(subNbr);
                        subFaceNeighbour.append(subOwn);
                        subNbrProcIDs.append(-1);
                        subProcLocalFaceIDs.append(-1);
                    }
                }
            }

            // boundary faces for new region
            forAll(tgtMesh.faceNeighbour(), faceI)
            {
                label own = tgtMesh.faceOwner()[faceI];
                label nbr = tgtMesh.faceNeighbour()[faceI];
                label subOwn = reverseCellMap[own];
                label subNbr = reverseCellMap[nbr];

                if (subOwn != -1 && subNbr == -1)
                {
                    subFaces.append(tgtMesh.faces()[faceI]);
                    subFaceOwner.append(subOwn);
                    subFaceNeighbour.append(subNbr);
                    subNbrProcIDs.append(-1);
                    subProcLocalFaceIDs.append(-1);
                }
                else if (subOwn == -1 && subNbr != -1)
                {
                    subFaces.append(tgtMesh.faces()[faceI].reverseFace());
                    subFaceOwner.append(subNbr);
                    subFaceNeighbour.append(subOwn);
                    subNbrProcIDs.append(-1);
                    subProcLocalFaceIDs.append(-1);
                }
            }

            // boundary faces of existing region
            forAll(tgtMesh.boundaryMesh(), patchI)
            {
                const polyPatch& pp = tgtMesh.boundaryMesh()[patchI];

                label nbrProcI = -1;

                // store info for faces on processor patches
                if (isA<processorPolyPatch>(pp))
                {
                    const processorPolyPatch& ppp =
                        dynamic_cast<const processorPolyPatch&>(pp);

                    nbrProcI = ppp.neighbProcNo();
                }

                forAll(pp, i)
                {
                    label faceI = pp.start() + i;
                    label own = tgtMesh.faceOwner()[faceI];

                    if (reverseCellMap[own] != -1)
                    {
                        subFaces.append(tgtMesh.faces()[faceI]);
                        subFaceOwner.append(reverseCellMap[own]);
                        subFaceNeighbour.append(-1);
                        subNbrProcIDs.append(nbrProcI);
                        subProcLocalFaceIDs.append(i);
                    }
                }
            }

            // reverse point map
            labelList reversePointMap(tgtMesh.nPoints(), -1);
            DynamicList<point> subPoints(tgtMesh.nPoints());
            forAll(subFaces, subFaceI)
            {
                face& f = subFaces[subFaceI];
                forAll(f, fp)
                {
                    label pointI = f[fp];
                    if (reversePointMap[pointI] == -1)
                    {
                        reversePointMap[pointI] = subPoints.size();
                        subPoints.append(tgtMesh.points()[pointI]);
                    }

                    f[fp] = reversePointMap[pointI];
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
    forAll(allPoints, procI)
    {
        pointOffset[procI] = nPoints;
        nPoints += allPoints[procI].size();
    }

    // Starting offset for cells
    label nCells = 0;
    labelList cellOffset(Pstream::nProcs(), 0);
    forAll(allTgtCellIDs, procI)
    {
        cellOffset[procI] = nCells;
        nCells += allTgtCellIDs[procI].size();
    }

    // Count any coupled faces
    typedef FixedList<label, 3> label3;
    typedef HashTable<label, label3, label3::Hash<> > procCoupleInfo;
    procCoupleInfo procFaceToGlobalCell;

    forAll(allNbrProcIDs, procI)
    {
        const labelList& nbrProcI = allNbrProcIDs[procI];
        const labelList& localFaceI = allProcLocalFaceIDs[procI];

        forAll(nbrProcI, i)
        {
            if (nbrProcI[i] != -1 && localFaceI[i] != -1)
            {
                label3 key;
                key[0] = min(procI, nbrProcI[i]);
                key[1] = max(procI, nbrProcI[i]);
                key[2] = localFaceI[i];

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

                    allNIntCoupledFaces[procI]++;
                }
            }
        }
    }


    // Starting offset for internal faces
    label nIntFaces = 0;
    label nFacesTotal = 0;
    labelList internalFaceOffset(Pstream::nProcs(), 0);
    forAll(allNIntCoupledFaces, procI)
    {
        label nCoupledFaces =
            allNIntCoupledFaces[procI] - allNInternalFaces[procI];

        internalFaceOffset[procI] = nIntFaces;
        nIntFaces += allNIntCoupledFaces[procI];
        nFacesTotal += allFaceOwners[procI].size() - nCoupledFaces;
    }

    tgtPoints.setSize(nPoints);
    tgtFaces.setSize(nFacesTotal);
    tgtFaceOwners.setSize(nFacesTotal);
    tgtFaceNeighbours.setSize(nFacesTotal);
    tgtCellIDs.setSize(nCells);

    // Insert points
    forAll(allPoints, procI)
    {
        const pointField& pts = allPoints[procI];
        SubList<point>(tgtPoints, pts.size(), pointOffset[procI]).assign(pts);
    }

    // Insert cellIDs
    forAll(allTgtCellIDs, procI)
    {
        const labelList& cellIDs = allTgtCellIDs[procI];
        SubList<label>(tgtCellIDs, cellIDs.size(), cellOffset[procI]).assign
        (
            cellIDs
        );
    }


    // Insert internal faces (from internal faces)
    forAll(allFaces, procI)
    {
        const faceList& fcs = allFaces[procI];
        const labelList& faceOs = allFaceOwners[procI];
        const labelList& faceNs = allFaceNeighbours[procI];

        SubList<face> slice
        (
            tgtFaces,
            allNInternalFaces[procI],
            internalFaceOffset[procI]
        );
        slice.assign(SubList<face>(fcs, allNInternalFaces[procI]));
        forAll(slice, i)
        {
            add(slice[i], pointOffset[procI]);
        }

        SubField<label> ownSlice
        (
            tgtFaceOwners,
            allNInternalFaces[procI],
            internalFaceOffset[procI]
        );
        ownSlice.assign(SubField<label>(faceOs, allNInternalFaces[procI]));
        add(ownSlice, cellOffset[procI]);

        SubField<label> nbrSlice
        (
            tgtFaceNeighbours,
            allNInternalFaces[procI],
            internalFaceOffset[procI]
        );
        nbrSlice.assign(SubField<label>(faceNs, allNInternalFaces[procI]));
        add(nbrSlice, cellOffset[procI]);

        internalFaceOffset[procI] += allNInternalFaces[procI];
    }


    // Insert internal faces (from coupled face-pairs)
    forAll(allNbrProcIDs, procI)
    {
        const labelList& nbrProcI = allNbrProcIDs[procI];
        const labelList& localFaceI = allProcLocalFaceIDs[procI];
        const labelList& faceOs = allFaceOwners[procI];
        const faceList& fcs = allFaces[procI];

        forAll(nbrProcI, i)
        {
            if (nbrProcI[i] != -1 && localFaceI[i] != -1)
            {
                label3 key;
                key[0] = min(procI, nbrProcI[i]);
                key[1] = max(procI, nbrProcI[i]);
                key[2] = localFaceI[i];

                procCoupleInfo::iterator fnd = procFaceToGlobalCell.find(key);

                if (fnd != procFaceToGlobalCell.end())
                {
                    label tgtFaceI = fnd();
                    if (tgtFaceI == -1)
                    {
                        // on first visit store the new cell on this side
                        fnd() = cellOffset[procI] + faceOs[i];
                    }
                    else
                    {
                        // get owner and neighbour in new cell numbering
                        label newOwn = cellOffset[procI] + faceOs[i];
                        label newNbr = fnd();
                        label tgtFaceI = internalFaceOffset[procI]++;

                        if (debug)
                        {
                            Pout<< "    proc " << procI
                                << "\tinserting face:" << tgtFaceI
                                << " connection between owner " << newOwn
                                << " and neighbour " << newNbr
                                << endl;
                        }

                        if (newOwn < newNbr)
                        {
                            // we have correct orientation
                            tgtFaces[tgtFaceI] = fcs[i];
                            tgtFaceOwners[tgtFaceI] = newOwn;
                            tgtFaceNeighbours[tgtFaceI] = newNbr;
                        }
                        else
                        {
                            // reverse orientation
                            tgtFaces[tgtFaceI] = fcs[i].reverseFace();
                            tgtFaceOwners[tgtFaceI] = newNbr;
                            tgtFaceNeighbours[tgtFaceI] = newOwn;
                        }

                        add(tgtFaces[tgtFaceI], pointOffset[procI]);

                        // mark with unique value
                        fnd() = -2;
                    }
                }
            }
        }
    }


    forAll(allNbrProcIDs, procI)
    {
        const labelList& nbrProcI = allNbrProcIDs[procI];
        const labelList& localFaceI = allProcLocalFaceIDs[procI];
        const labelList& faceOs = allFaceOwners[procI];
        const labelList& faceNs = allFaceNeighbours[procI];
        const faceList& fcs = allFaces[procI];

        forAll(nbrProcI, i)
        {
            // coupled boundary face
            if (nbrProcI[i] != -1 && localFaceI[i] != -1)
            {
                label3 key;
                key[0] = min(procI, nbrProcI[i]);
                key[1] = max(procI, nbrProcI[i]);
                key[2] = localFaceI[i];

                label tgtFaceI = procFaceToGlobalCell[key];

                if (tgtFaceI == -1)
                {
                    FatalErrorIn
                    (
                        "void Foam::meshToMesh::"
                        "distributeAndMergeCells"
                        "("
                            "const mapDistribute&, "
                            "const polyMesh&, "
                            "const globalIndex&, "
                            "pointField&, "
                            "faceList&, "
                            "labelList&, "
                            "labelList&, "
                            "labelList&"
                        ") const"
                    )
                        << "Unvisited " << key
                        << abort(FatalError);
                }
                else if (tgtFaceI != -2)
                {
                    label newOwn = cellOffset[procI] + faceOs[i];
                    label tgtFaceI = nIntFaces++;

                    if (debug)
                    {
                        Pout<< "    proc " << procI
                            << "\tinserting boundary face:" << tgtFaceI
                            << " from coupled face " << key
                            << endl;
                    }

                    tgtFaces[tgtFaceI] = fcs[i];
                    add(tgtFaces[tgtFaceI], pointOffset[procI]);

                    tgtFaceOwners[tgtFaceI] = newOwn;
                    tgtFaceNeighbours[tgtFaceI] = -1;
                }
            }
            // normal boundary face
            else
            {
                label own = faceOs[i];
                label nbr = faceNs[i];
                if ((own != -1) && (nbr == -1))
                {
                    label newOwn = cellOffset[procI] + faceOs[i];
                    label tgtFaceI = nIntFaces++;

                    tgtFaces[tgtFaceI] = fcs[i];
                    add(tgtFaces[tgtFaceI], pointOffset[procI]);

                    tgtFaceOwners[tgtFaceI] = newOwn;
                    tgtFaceNeighbours[tgtFaceI] = -1;
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
            SMALL,
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
