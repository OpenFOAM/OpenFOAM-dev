/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2023 OpenFOAM Foundation
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

#include "cellsToCells.H"
#include "processorPolyPatch.H"
#include "SubField.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

inline void offset(label& lst, const label o)
{
    lst += o;
}

template<class ListType>
void offset(ListType& lst, const label o)
{
    forAll(lst, i)
    {
        offset(lst[i], o);
    }
}

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelListList Foam::cellsToCells::tgtMeshSendCells
(
    const polyMesh& srcMesh,
    const polyMesh& tgtMesh
) const
{
    // Calculate and communicate the bound boxes for the source meshes
    List<boundBox> procBb(Pstream::nProcs(), boundBox());
    if (srcMesh.nCells() > 0)
    {
        // Bounding box for this mesh. Do not reduce.
        procBb[Pstream::myProcNo()] = boundBox(srcMesh.points(), false);

        // Slightly increase size to allow for cases where boxes are aligned
        procBb[Pstream::myProcNo()].inflate(0.01);
    }
    Pstream::gatherList(procBb);
    Pstream::scatterList(procBb);

    // per processor indices into all segments to send
    List<DynamicList<label>> resultDyn
    (
        Pstream::nProcs(),
        DynamicList<label>(tgtMesh.nCells()/Pstream::nProcs())
    );

    // Send a target cell to a process if it overlaps the source bound box
    // for that process
    forAll(tgtMesh.cells(), tgtCelli)
    {
        const boundBox cellBb =
            tgtMesh.cells()[tgtCelli].bb(tgtMesh.points(), tgtMesh.faces());

        forAll(procBb, proci)
        {
            if (procBb[proci].overlaps(cellBb))
            {
                resultDyn[proci].append(tgtCelli);
            }
        }
    }

    // Transfer to permanent storage and return
    labelListList result(Pstream::nProcs());
    forAll(result, proci)
    {
        result[proci].transfer(resultDyn[proci]);
    }

    return result;
}


Foam::List<Foam::remote> Foam::cellsToCells::distributeMesh
(
    const distributionMap& map,
    const polyMesh& mesh,
    autoPtr<polyMesh>& localMeshPtr
)
{
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // Exchange raw mesh topology/geometry
    List<pointField> allPoints(Pstream::nProcs());
    labelList allNInternalFaces(Pstream::nProcs(), 0);
    List<faceList> allFaces(Pstream::nProcs());
    List<labelList> allFaceOwners(Pstream::nProcs());
    List<labelList> allFaceNeighbours(Pstream::nProcs());
    List<List<remote>> allProcCells(Pstream::nProcs());
    List<List<remote>> allProcProcessorPatchFaces(Pstream::nProcs());
    for (label domain = 0; domain < Pstream::nProcs(); domain++)
    {
        const labelList& sendElems = map.subMap()[domain];

        if (sendElems.size())
        {
            // reverse cell map
            labelList reverseCellMap(mesh.nCells(), -1);
            forAll(sendElems, subCelli)
            {
                reverseCellMap[sendElems[subCelli]] = subCelli;
            }

            DynamicList<face> subFaces(mesh.nFaces());
            DynamicList<label> subFaceOwner(mesh.nFaces());
            DynamicList<label> subFaceNeighbour(mesh.nFaces());
            DynamicList<remote> subProcProcessorPatchFaces(mesh.nFaces());

            label nInternal = 0;

            // internal faces
            forAll(mesh.faceNeighbour(), facei)
            {
                const label own = mesh.faceOwner()[facei];
                const label nbr = mesh.faceNeighbour()[facei];
                const label subOwn = reverseCellMap[own];
                const label subNbr = reverseCellMap[nbr];

                if (subOwn != -1 && subNbr != -1)
                {
                    nInternal++;

                    if (subOwn < subNbr)
                    {
                        subFaces.append(mesh.faces()[facei]);
                        subFaceOwner.append(subOwn);
                        subFaceNeighbour.append(subNbr);
                        subProcProcessorPatchFaces.append(remote());
                    }
                    else
                    {
                        subFaces.append(mesh.faces()[facei].reverseFace());
                        subFaceOwner.append(subNbr);
                        subFaceNeighbour.append(subOwn);
                        subProcProcessorPatchFaces.append(remote());
                    }
                }
            }

            // boundary faces for new region
            forAll(mesh.faceNeighbour(), facei)
            {
                const label own = mesh.faceOwner()[facei];
                const label nbr = mesh.faceNeighbour()[facei];
                const label subOwn = reverseCellMap[own];
                const label subNbr = reverseCellMap[nbr];

                if (subOwn != -1 && subNbr == -1)
                {
                    subFaces.append(mesh.faces()[facei]);
                    subFaceOwner.append(subOwn);
                    subFaceNeighbour.append(subNbr);
                    subProcProcessorPatchFaces.append(remote());
                }
                else if (subOwn == -1 && subNbr != -1)
                {
                    subFaces.append(mesh.faces()[facei].reverseFace());
                    subFaceOwner.append(subNbr);
                    subFaceNeighbour.append(subOwn);
                    subProcProcessorPatchFaces.append(remote());
                }
            }

            // boundary faces of existing region
            forAll(mesh.boundaryMesh(), patchi)
            {
                const polyPatch& pp = mesh.boundaryMesh()[patchi];

                const label nbrProci =
                    isType<processorPolyPatch>(pp)
                  ? refCast<const processorPolyPatch>(pp).neighbProcNo()
                  : -1;

                forAll(pp, i)
                {
                    const label facei = pp.start() + i;
                    const label own = mesh.faceOwner()[facei];

                    if (reverseCellMap[own] != -1)
                    {
                        subFaces.append(mesh.faces()[facei]);
                        subFaceOwner.append(reverseCellMap[own]);
                        subFaceNeighbour.append(-1);
                        subProcProcessorPatchFaces.append
                        (
                            nbrProci != -1 ? remote(nbrProci, i) : remote()
                        );
                    }
                }
            }

            // reverse point map
            labelList reversePointMap(mesh.nPoints(), -1);
            DynamicList<point> subPoints(mesh.nPoints());
            forAll(subFaces, subFacei)
            {
                face& f = subFaces[subFacei];
                forAll(f, fp)
                {
                    label pointi = f[fp];
                    if (reversePointMap[pointi] == -1)
                    {
                        reversePointMap[pointi] = subPoints.size();
                        subPoints.append(mesh.points()[pointi]);
                    }

                    f[fp] = reversePointMap[pointi];
                }
            }

            // cell indices
            List<remote> subProcCells(sendElems.size());
            forAll(sendElems, i)
            {
                subProcCells[i] = remote(Pstream::myProcNo(), sendElems[i]);
            }

            // pass data
            if (domain == Pstream::myProcNo())
            {
                // allocate my own data
                allPoints[Pstream::myProcNo()] = subPoints;
                allNInternalFaces[Pstream::myProcNo()] = nInternal;
                allFaces[Pstream::myProcNo()] = subFaces;
                allFaceOwners[Pstream::myProcNo()] = subFaceOwner;
                allFaceNeighbours[Pstream::myProcNo()] = subFaceNeighbour;
                allProcCells[Pstream::myProcNo()] = subProcCells;
                allProcProcessorPatchFaces[Pstream::myProcNo()] =
                    subProcProcessorPatchFaces;
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
                    << subProcCells
                    << subProcProcessorPatchFaces;
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

            str >> allPoints[domain]
                >> allNInternalFaces[domain]
                >> allFaces[domain]
                >> allFaceOwners[domain]
                >> allFaceNeighbours[domain]
                >> allProcCells[domain]
                >> allProcProcessorPatchFaces[domain];
        }
    }

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

    // Count any coupled faces
    typedef FixedList<label, 3> label3;
    typedef HashTable<label, label3, label3::Hash<>> procCoupleInfo;
    procCoupleInfo procFaceToGlobalCell;
    forAll(allProcProcessorPatchFaces, proci)
    {
        const List<remote>& procProcessorPatchFaces =
            allProcProcessorPatchFaces[proci];

        forAll(procProcessorPatchFaces, i)
        {
            if (procProcessorPatchFaces[i] != remote())
            {
                const label3 key
                ({
                    min(proci, procProcessorPatchFaces[i].proci),
                    max(proci, procProcessorPatchFaces[i].proci),
                    procProcessorPatchFaces[i].elementi
                });

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

    // Starting offset for cells
    label nCells = 0;
    labelList cellOffset(Pstream::nProcs(), 0);
    forAll(allProcCells, proci)
    {
        cellOffset[proci] = nCells;
        nCells += allProcCells[proci].size();
    }

    // Allocate
    List<remote> localProcCells(nCells);
    pointField localPoints(nPoints);
    faceList localFaces(nFacesTotal);
    labelList localFaceOwners(nFacesTotal);
    labelList localFaceNeighbours(nFacesTotal);

    // Insert proc-cells
    forAll(allProcCells, proci)
    {
        const List<remote>& procCells = allProcCells[proci];
        SubList<remote>(localProcCells, procCells.size(), cellOffset[proci]) =
            procCells;
    }

    // Insert points
    forAll(allPoints, proci)
    {
        const pointField& pts = allPoints[proci];
        SubList<point>(localPoints, pts.size(), pointOffset[proci]) = pts;
    }

    // Insert internal faces (from internal faces)
    forAll(allFaces, proci)
    {
        const faceList& fcs = allFaces[proci];
        const labelList& faceOs = allFaceOwners[proci];
        const labelList& faceNs = allFaceNeighbours[proci];

        SubList<face> slice
        (
            localFaces,
            allNInternalFaces[proci],
            internalFaceOffset[proci]
        );
        slice = SubList<face>(fcs, allNInternalFaces[proci]);
        offset(slice, pointOffset[proci]);

        SubField<label> ownSlice
        (
            localFaceOwners,
            allNInternalFaces[proci],
            internalFaceOffset[proci]
        );
        ownSlice = SubField<label>(faceOs, allNInternalFaces[proci]);
        offset(ownSlice, cellOffset[proci]);

        SubField<label> nbrSlice
        (
            localFaceNeighbours,
            allNInternalFaces[proci],
            internalFaceOffset[proci]
        );
        nbrSlice = SubField<label>(faceNs, allNInternalFaces[proci]);
        offset(nbrSlice, cellOffset[proci]);

        internalFaceOffset[proci] += allNInternalFaces[proci];
    }

    // Insert internal faces (from coupled face-pairs)
    forAll(allProcProcessorPatchFaces, proci)
    {
        const List<remote>& procProcessorPatchFaces =
            allProcProcessorPatchFaces[proci];
        const labelList& faceOs = allFaceOwners[proci];
        const faceList& fcs = allFaces[proci];

        forAll(procProcessorPatchFaces, i)
        {
            if (procProcessorPatchFaces[i] != remote())
            {
                const label3 key
                ({
                    min(proci, procProcessorPatchFaces[i].proci),
                    max(proci, procProcessorPatchFaces[i].proci),
                    procProcessorPatchFaces[i].elementi
                });

                procCoupleInfo::iterator fnd = procFaceToGlobalCell.find(key);

                if (fnd != procFaceToGlobalCell.end())
                {
                    if (fnd() == -1)
                    {
                        // on first visit store the new cell on this side
                        fnd() = cellOffset[proci] + faceOs[i];
                    }
                    else
                    {
                        // get owner and neighbour in new cell numbering
                        const label newOwn = cellOffset[proci] + faceOs[i];
                        const label newNbr = fnd();
                        const label localFacei = internalFaceOffset[proci]++;

                        if (debug)
                        {
                            Pout<< "    proc " << proci
                                << "\tinserting face:" << localFacei
                                << " connection between owner " << newOwn
                                << " and neighbour " << newNbr
                                << endl;
                        }

                        if (newOwn < newNbr)
                        {
                            // we have correct orientation
                            localFaces[localFacei] = fcs[i];
                            localFaceOwners[localFacei] = newOwn;
                            localFaceNeighbours[localFacei] = newNbr;
                        }
                        else
                        {
                            // reverse orientation
                            localFaces[localFacei] = fcs[i].reverseFace();
                            localFaceOwners[localFacei] = newNbr;
                            localFaceNeighbours[localFacei] = newOwn;
                        }

                        offset(localFaces[localFacei], pointOffset[proci]);

                        // mark with unique value
                        fnd() = -2;
                    }
                }
            }
        }
    }

    forAll(allProcProcessorPatchFaces, proci)
    {
        const List<remote>& procProcessorPatchFaces =
            allProcProcessorPatchFaces[proci];
        const labelList& faceOs = allFaceOwners[proci];
        const labelList& faceNs = allFaceNeighbours[proci];
        const faceList& fcs = allFaces[proci];

        forAll(procProcessorPatchFaces, i)
        {
            // coupled boundary face
            if (procProcessorPatchFaces[i] != remote())
            {
                const label3 key
                ({
                    min(proci, procProcessorPatchFaces[i].proci),
                    max(proci, procProcessorPatchFaces[i].proci),
                    procProcessorPatchFaces[i].elementi
                });

                if (procFaceToGlobalCell[key] == -1)
                {
                    FatalErrorInFunction
                        << "Unvisited " << key
                        << abort(FatalError);
                }
                else if (procFaceToGlobalCell[key] != -2)
                {
                    const label newOwn = cellOffset[proci] + faceOs[i];
                    const label localFacei = nIntFaces++;

                    if (debug)
                    {
                        Pout<< "    proc " << proci
                            << "\tinserting boundary face:" << localFacei
                            << " from coupled face " << key
                            << endl;
                    }

                    localFaces[localFacei] = fcs[i];
                    offset(localFaces[localFacei], pointOffset[proci]);

                    localFaceOwners[localFacei] = newOwn;
                    localFaceNeighbours[localFacei] = -1;
                }
            }

            // normal boundary face
            else
            {
                const label own = faceOs[i];
                const label nbr = faceNs[i];

                if ((own != -1) && (nbr == -1))
                {
                    const label newOwn = cellOffset[proci] + faceOs[i];
                    const label localFacei = nIntFaces++;

                    localFaces[localFacei] = fcs[i];
                    offset(localFaces[localFacei], pointOffset[proci]);

                    localFaceOwners[localFacei] = newOwn;
                    localFaceNeighbours[localFacei] = -1;
                }
            }
        }
    }

    // Create the local mesh
    localMeshPtr.reset
    (
        new polyMesh
        (
            IOobject
            (
                "local" + mesh.name().capitalise(),
                mesh.time().name(),
                mesh.time(),
                IOobject::NO_READ
            ),
            move(localPoints),
            move(localFaces),
            move(localFaceOwners),
            move(localFaceNeighbours),
            false
        )
    );

    // Add a dummy patch to the target mesh
    List<polyPatch*> patches(1);
    patches[0] = new polyPatch
    (
        "defaultFaces",
        localMeshPtr().nFaces() - localMeshPtr().nInternalFaces(),
        localMeshPtr().nInternalFaces(),
        0,
        localMeshPtr().boundaryMesh(),
        word::null
    );
    localMeshPtr().addPatches(patches);

    // Force calculation of tet-base points used for point-in-cell
    (void) localMeshPtr().tetBasePtIs();

    return localProcCells;
}


void Foam::cellsToCells::trimLocalTgt()
{
    // Determine which local target cells are actually used
    boolList oldLocalTgtCellIsUsed(localTgtProcCellsPtr_().size(), false);
    forAll(srcLocalTgtCells_, srcCelli)
    {
        forAll(srcLocalTgtCells_[srcCelli], i)
        {
            oldLocalTgtCellIsUsed[srcLocalTgtCells_[srcCelli][i]] = true;
        }
    }

    // Trim the target map
    labelList oldToNewLocalTgtCell, newToOldLocalTgtCell;
    patchToPatchTools::trimDistributionMap
    (
        oldLocalTgtCellIsUsed,
        tgtMapPtr_(),
        oldToNewLocalTgtCell,
        newToOldLocalTgtCell
    );

    if (debug)
    {
        Pout<< "Trim from " << oldToNewLocalTgtCell.size() << " to "
            << newToOldLocalTgtCell.size() << " cells" << endl;
    }

    // Renumber the source addressing
    forAll(srcLocalTgtCells_, srcCelli)
    {
        forAll(srcLocalTgtCells_[srcCelli], i)
        {
            srcLocalTgtCells_[srcCelli][i] =
                oldToNewLocalTgtCell[srcLocalTgtCells_[srcCelli][i]];
        }
    }

    // Trim the target addressing
    tgtLocalSrcCells_ =
        labelListList(tgtLocalSrcCells_, newToOldLocalTgtCell);
    localTgtProcCellsPtr_() =
        List<remote>(localTgtProcCellsPtr_(), newToOldLocalTgtCell);

    // Trim the target weights
    tgtWeights_ = scalarListList(tgtWeights_, newToOldLocalTgtCell);

    // Trim the stored local target mesh
    const polyMesh& oldLocalTgtMesh = localTgtMeshPtr_();
    const labelList& oldLocalTgtFaceOwner =
        oldLocalTgtMesh.faceOwner();
    labelList oldLocalTgtFaceNeighbour(oldLocalTgtFaceOwner.size(), -1);
    SubList<label>(oldLocalTgtFaceNeighbour, oldLocalTgtMesh.nInternalFaces()) =
        oldLocalTgtMesh.faceNeighbour();

    // ...
    labelList newToOldLocalTgtFace(identityMap(oldLocalTgtMesh.nFaces()));
    labelList oldToNewLocalTgtFace;
    {
        label i0 = 0;
        label i1 = newToOldLocalTgtFace.size();
        label iEnd = newToOldLocalTgtFace.size();

        while (i0 < i1)
        {
            label& oldLocalTgtFacei0 = newToOldLocalTgtFace[i0];
            label& oldLocalTgtFacei1 = newToOldLocalTgtFace[i1 - 1];
            label& oldLocalTgtFaceiEnd = newToOldLocalTgtFace[iEnd - 1];

            const label newLocalTgtOwni0 =
                oldLocalTgtFaceOwner[oldLocalTgtFacei0] != -1
              ? oldToNewLocalTgtCell
                [oldLocalTgtFaceOwner[oldLocalTgtFacei0]]
              : -1;
            const label newLocalTgtOwni1 =
                oldLocalTgtFaceOwner[oldLocalTgtFacei1] != -1
              ? oldToNewLocalTgtCell
                [oldLocalTgtFaceOwner[oldLocalTgtFacei1]]
              : -1;

            const label newLocalTgtNbri0 =
                oldLocalTgtFaceNeighbour[oldLocalTgtFacei0] != -1
              ? oldToNewLocalTgtCell
                [oldLocalTgtFaceNeighbour[oldLocalTgtFacei0]]
              : -1;
            const label newLocalTgtNbri1 =
                oldLocalTgtFaceNeighbour[oldLocalTgtFacei1] != -1
              ? oldToNewLocalTgtCell
                [oldLocalTgtFaceNeighbour[oldLocalTgtFacei1]]
              : -1;

            const bool used0 =
                newLocalTgtOwni0 != -1 || newLocalTgtNbri0 != -1;
            const bool used1 =
                newLocalTgtOwni1 != -1 || newLocalTgtNbri1 != -1;

            const bool internal0 =
                newLocalTgtOwni0 != -1 && newLocalTgtNbri0 != -1;
            const bool internal1 =
                newLocalTgtOwni1 != -1 && newLocalTgtNbri1 != -1;

            // If face 0 is not used, move it to the end to remove it
            if (!used0)
            {
                Swap(oldLocalTgtFacei0, oldLocalTgtFaceiEnd);
                if (i1 == iEnd) i1 --;
                iEnd --;
            }

            // If face 1 is not used, move it to the end to remove it
            else if (!used1)
            {
                Swap(oldLocalTgtFacei1, oldLocalTgtFaceiEnd);
                if (i1 == iEnd) i1 --;
                iEnd --;
            }

            // Both are internal faces. Face 0 is fine, but face 1 might be out
            // of order. So move to the next face 0.
            else if (internal0 && internal1)
            {
                i0 ++;
            }

            // Both are boundary faces. Face 0 might be out of order, but face
            // 1 is fine. So move to the next face 1.
            else if (!internal0 && !internal1)
            {
                i1 --;
            }

            // Face 0 is an internal face and face 1 is a boundary face. Both
            // are fine. So move to the next two faces.
            else if (internal0 && !internal1)
            {
                i0 ++;
                i1 --;
            }

            // Face 0 is a boundary face and face 1 is an internal face. Both
            // are out of order. So swap and then move to the next two faces.
            else if (!internal0 && internal1)
            {
                Swap(oldLocalTgtFacei0, oldLocalTgtFacei1);
                i0 ++;
                i1 --;
            }
        }

        newToOldLocalTgtFace.resize(iEnd);

        oldToNewLocalTgtFace =
            invert(oldLocalTgtMesh.nFaces(), newToOldLocalTgtFace);
    }

    if (debug)
    {
        Pout<< "Trim from " << oldToNewLocalTgtFace.size() << " to "
            << newToOldLocalTgtFace.size() << " faces" << endl;
    }

    // Create trimmed mesh primitives
    pointField newLocalTgtPoints(oldLocalTgtMesh.points());
    faceList newLocalTgtFaces(oldLocalTgtMesh.faces(), newToOldLocalTgtFace);
    labelList newLocalTgtFaceOwner
    (
        oldLocalTgtFaceOwner,
        static_cast<const labelUList&>(newToOldLocalTgtFace)
    );
    inplaceRenumber(oldToNewLocalTgtCell, newLocalTgtFaceOwner);
    labelList newLocalTgtFaceNeighbour
    (
        oldLocalTgtFaceNeighbour,
        static_cast<const labelUList&>(newToOldLocalTgtFace)
    );
    inplaceRenumber(oldToNewLocalTgtCell, newLocalTgtFaceNeighbour);

    // Check the trimmed mesh structure and flip any faces that have a
    // neighbour but not an owner
    {
        label newLocalTgtNInternalFaces = 0;
        bool internal0 = true;

        forAll(newLocalTgtFaces, newLocalTgtFacei)
        {
            face& newLocalTgtF = newLocalTgtFaces[newLocalTgtFacei];
            label& newLocalTgtOwni = newLocalTgtFaceOwner[newLocalTgtFacei];
            label& newLocalTgtNbri = newLocalTgtFaceNeighbour[newLocalTgtFacei];

            // Check that internal and boundary faces are in order
            const bool internal =
                newLocalTgtOwni != -1 && newLocalTgtNbri != -1;
            if (internal0 && !internal)
            {
                newLocalTgtNInternalFaces = newLocalTgtFacei;
                internal0 = false;
            }
            if (!internal0 && internal)
            {
                FatalErrorInFunction
                    << "Trimmed mesh has boundary faces before internal faces"
                    << exit(FatalError);
            }

            // Flip any face that has a neighbour but not an owner
            const bool flip =
                newLocalTgtOwni == -1 && newLocalTgtNbri != -1;
            if (flip)
            {
                newLocalTgtF = newLocalTgtF.reverseFace();
                Swap(newLocalTgtOwni, newLocalTgtNbri);
            }
        }

        newLocalTgtFaceNeighbour.resize(newLocalTgtNInternalFaces);
    }

    // Create the local mesh
    localTgtMeshPtr_.reset
    (
        new polyMesh
        (
            IOobject
            (
                "trimmed" + oldLocalTgtMesh.name().capitalise(),
                oldLocalTgtMesh.time().name(),
                oldLocalTgtMesh.time(),
                IOobject::NO_READ
            ),
            move(newLocalTgtPoints),
            move(newLocalTgtFaces),
            move(newLocalTgtFaceOwner),
            move(newLocalTgtFaceNeighbour),
            false
        )
    );

    // Add a dummy patch to the target mesh
    List<polyPatch*> patches(1);
    patches[0] = new polyPatch
    (
        "defaultFaces",
        localTgtMeshPtr_().nFaces() - localTgtMeshPtr_().nInternalFaces(),
        localTgtMeshPtr_().nInternalFaces(),
        0,
        localTgtMeshPtr_().boundaryMesh(),
        word::null
    );
    localTgtMeshPtr_().addPatches(patches);

    // Force calculation of tet-base points used for point-in-cell
    (void) localTgtMeshPtr_().tetBasePtIs();
}


// ************************************************************************* //
