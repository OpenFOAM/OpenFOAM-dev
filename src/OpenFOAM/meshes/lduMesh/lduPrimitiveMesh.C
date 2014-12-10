/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "lduPrimitiveMesh.H"
#include "processorLduInterface.H"
#include "EdgeMap.H"
#include "labelPair.H"
#include "processorGAMGInterface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lduPrimitiveMesh, 0);

    //- Less operator for pairs of \<processor\>\<index\>
    class procLess
    {
        const labelPairList& lst_;

    public:

        procLess(const labelPairList& lst)
        :
            lst_(lst)
        {}

        bool operator()(const label a, const label b)
        {
            return lst_[a].first() < lst_[b].first();
        }
    };

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::lduPrimitiveMesh::checkUpperTriangular
(
    const label size,
    const labelUList& l,
    const labelUList& u
)
{
    forAll(l, faceI)
    {
        if (u[faceI] < l[faceI])
        {
            FatalErrorIn
            (
                "checkUpperTriangular"
                "(const label, const labelUList&, const labelUList&)"
            )   << "Reversed face. Problem at face " << faceI
                << " l:" << l[faceI] << " u:" << u[faceI]
                << abort(FatalError);
        }
        if (l[faceI] < 0 || u[faceI] < 0 || u[faceI] >= size)
        {
            FatalErrorIn
            (
                "checkUpperTriangular"
                "(const label, const labelUList&, const labelUList&)"
            )   << "Illegal cell label. Problem at face " << faceI
                << " l:" << l[faceI] << " u:" << u[faceI]
                << abort(FatalError);
        }
    }

    for (label faceI=1; faceI < l.size(); faceI++)
    {
        if (l[faceI-1] > l[faceI])
        {
            FatalErrorIn
            (
                "checkUpperTriangular"
                "(const label, const labelUList&, const labelUList&)"
            )   << "Lower not in incremental cell order."
                << " Problem at face " << faceI
                << " l:" << l[faceI] << " u:" << u[faceI]
                << " previous l:" << l[faceI-1]
                << abort(FatalError);
        }
        else if (l[faceI-1] == l[faceI])
        {
            // Same cell.
            if (u[faceI-1] > u[faceI])
            {
                FatalErrorIn
                (
                    "checkUpperTriangular"
                    "(const label, const labelUList&, const labelUList&)"
                )   << "Upper not in incremental cell order."
                    << " Problem at face " << faceI
                    << " l:" << l[faceI] << " u:" << u[faceI]
                    << " previous u:" << u[faceI-1]
                    << abort(FatalError);
            }
        }
    }
}


Foam::label Foam::lduPrimitiveMesh::totalSize
(
    const PtrList<lduPrimitiveMesh>& meshes
)
{
    label size = 0;

    forAll(meshes, i)
    {
        size += meshes[i].lduAddr().size();
    }
    return size;
}


Foam::labelList Foam::lduPrimitiveMesh::upperTriOrder
(
    const label nCells,
    const labelUList& lower,
    const labelUList& upper
)
{
    labelList nNbrs(nCells, 0);

    // Count number of upper neighbours
    forAll(lower, faceI)
    {
        if (upper[faceI] < lower[faceI])
        {
            FatalErrorIn("lduPrimitiveMesh::upperTriOrder(..)")
                << "Problem at face:" << faceI
                << " lower:" << lower[faceI]
                << " upper:" << upper[faceI]
                << exit(FatalError);
        }
        nNbrs[lower[faceI]]++;
    }

    // Construct cell-upper cell addressing
    labelList offsets(nCells+1);
    offsets[0] = 0;
    forAll(nNbrs, cellI)
    {
        offsets[cellI+1] = offsets[cellI]+nNbrs[cellI];
    }

    nNbrs = offsets;

    labelList cellToFaces(offsets.last());
    forAll(upper, faceI)
    {
        label cellI = lower[faceI];
        cellToFaces[nNbrs[cellI]++] = faceI;
    }

    // Sort

    labelList oldToNew(lower.size());

    labelList order;
    labelList nbr;

    label newFaceI = 0;

    for (label cellI = 0; cellI < nCells; cellI++)
    {
        label startOfCell = offsets[cellI];
        label nNbr = offsets[cellI+1] - startOfCell;

        nbr.setSize(nNbr);
        order.setSize(nNbr);
        forAll(order, i)
        {
            nbr[i] = upper[cellToFaces[offsets[cellI]+i]];
        }
        sortedOrder(nbr, order);

        forAll(order, i)
        {
            label index = order[i];
            oldToNew[cellToFaces[startOfCell + index]] = newFaceI++;
        }
    }

    return oldToNew;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lduPrimitiveMesh::lduPrimitiveMesh
(
    const label nCells,
    labelList& l,
    labelList& u,
    const label comm,
    bool reUse
)
:
    lduAddressing(nCells),
    lowerAddr_(l, reUse),
    upperAddr_(u, reUse),
    comm_(comm)
{}


void Foam::lduPrimitiveMesh::addInterfaces
(
    lduInterfacePtrsList& interfaces,
    const lduSchedule& ps
)
{
    interfaces_ = interfaces;
    patchSchedule_ = ps;

    // Create interfaces
    primitiveInterfaces_.setSize(interfaces_.size());
    forAll(interfaces_, i)
    {
        if (interfaces_.set(i))
        {
            primitiveInterfaces_.set(i, &interfaces_[i]);
        }
    }
}


Foam::lduPrimitiveMesh::lduPrimitiveMesh
(
    const label nCells,
    labelList& l,
    labelList& u,
    PtrList<const lduInterface>& primitiveInterfaces,
    const lduSchedule& ps,
    const label comm
)
:
    lduAddressing(nCells),
    lowerAddr_(l, true),
    upperAddr_(u, true),
    primitiveInterfaces_(0),
    patchSchedule_(ps),
    comm_(comm)
{
    primitiveInterfaces_.transfer(primitiveInterfaces);

    // Create interfaces
    interfaces_.setSize(primitiveInterfaces_.size());
    forAll(primitiveInterfaces_, i)
    {
        if (primitiveInterfaces_.set(i))
        {
            interfaces_.set(i, &primitiveInterfaces_[i]);
        }
    }
}


Foam::lduPrimitiveMesh::lduPrimitiveMesh
(
    const label comm,
    const labelList& procAgglomMap,

    const labelList& procIDs,
    const lduMesh& myMesh,
    const PtrList<lduPrimitiveMesh>& otherMeshes,

    labelList& cellOffsets,
    labelList& faceOffsets,
    labelListList& faceMap,
    labelListList& boundaryMap,
    labelListListList& boundaryFaceMap
)
:
    lduAddressing(myMesh.lduAddr().size() + totalSize(otherMeshes)),
    lowerAddr_(0),
    upperAddr_(0),
    interfaces_(0),
    patchSchedule_(0),
    comm_(comm)
{
    const label currentComm = myMesh.comm();

    forAll(otherMeshes, i)
    {
        if (otherMeshes[i].comm() != currentComm)
        {
            WarningIn
            (
                "lduPrimitiveMesh::lduPrimitiveMesh(..)"
            )   << "Communicator " << otherMeshes[i].comm()
                << " at index " << i
                << " differs from that of predecessor "
                << currentComm
                << endl;    //exit(FatalError);
        }
    }

    const label nMeshes = otherMeshes.size()+1;

    const label myAgglom = procAgglomMap[UPstream::myProcNo(currentComm)];

    if (lduPrimitiveMesh::debug)
    {
        Pout<< "I am " << UPstream::myProcNo(currentComm)
            << " agglomerating into " << myAgglom
            << " as are " << findIndices(procAgglomMap, myAgglom)
            << endl;
    }


    forAll(procIDs, i)
    {
        if (procAgglomMap[procIDs[i]] != procAgglomMap[procIDs[0]])
        {
            FatalErrorIn
            (
                "lduPrimitiveMesh::lduPrimitiveMesh(..)"
            )   << "Processor " << procIDs[i]
                << " agglomerates to " << procAgglomMap[procIDs[i]]
                << " whereas other processors " << procIDs
                << " agglomerate to "
                << UIndirectList<label>(procAgglomMap, procIDs)
                << exit(FatalError);
        }
    }


    // Cells get added in order.
    cellOffsets.setSize(nMeshes+1);
    cellOffsets[0] = 0;
    for (label procMeshI = 0; procMeshI < nMeshes; procMeshI++)
    {
        const lduMesh& procMesh = mesh(myMesh, otherMeshes, procMeshI);

        cellOffsets[procMeshI+1] =
            cellOffsets[procMeshI]
          + procMesh.lduAddr().size();
    }


    // Faces initially get added in order, sorted later
    labelList internalFaceOffsets(nMeshes+1);
    internalFaceOffsets[0] = 0;
    for (label procMeshI = 0; procMeshI < nMeshes; procMeshI++)
    {
        const lduMesh& procMesh = mesh(myMesh, otherMeshes, procMeshI);

        internalFaceOffsets[procMeshI+1] =
            internalFaceOffsets[procMeshI]
          + procMesh.lduAddr().lowerAddr().size();
    }

    // Count how faces get added. Interfaces inbetween get merged.

    // Merged interfaces: map from two coarse processors back to
    // - procMeshes
    // - interface in procMesh
    // (estimate size from number of patches of mesh0)
    EdgeMap<labelPairList> mergedMap(2*myMesh.interfaces().size());

    // Unmerged interfaces: map from two coarse processors back to
    // - procMeshes
    // - interface in procMesh
    EdgeMap<labelPairList> unmergedMap(mergedMap.size());

    boundaryMap.setSize(nMeshes);
    boundaryFaceMap.setSize(nMeshes);


    label nOtherInterfaces = 0;
    labelList nCoupledFaces(nMeshes, 0);

    for (label procMeshI = 0; procMeshI < nMeshes; procMeshI++)
    {
        const lduInterfacePtrsList interfaces =
            mesh(myMesh, otherMeshes, procMeshI).interfaces();

        // Inialise all boundaries as merged
        boundaryMap[procMeshI].setSize(interfaces.size(), -1);
        boundaryFaceMap[procMeshI].setSize(interfaces.size());

        // Get sorted order of processors
        forAll(interfaces, intI)
        {
            if (interfaces.set(intI))
            {
                const lduInterface& ldui = interfaces[intI];

                if (isA<processorLduInterface>(ldui))
                {
                    const processorLduInterface& pldui =
                        refCast<const processorLduInterface>(ldui);

                    label agglom0 = procAgglomMap[pldui.myProcNo()];
                    label agglom1 = procAgglomMap[pldui.neighbProcNo()];

                    const edge procEdge(agglom0, agglom1);

                    if (agglom0 != myAgglom && agglom1 != myAgglom)
                    {
                        FatalErrorIn("lduPrimitiveMesh::lduPrimitiveMesh(..)")
                            << "At mesh from processor " << procIDs[procMeshI]
                            << " have interface " << intI
                            << " with myProcNo:" << pldui.myProcNo()
                            << " with neighbProcNo:" << pldui.neighbProcNo()
                            << exit(FatalError);
                    }
                    else if (agglom0 == myAgglom && agglom1 == myAgglom)
                    {
                        // Merged interface
                        if (debug)
                        {
                            Pout<< "merged interface: myProcNo:"
                                << pldui.myProcNo()
                                << " nbr:" << pldui.neighbProcNo()
                                << " size:" << ldui.faceCells().size()
                                << endl;
                        }

                        label nbrProcMeshI = findIndex
                        (
                            procIDs,
                            pldui.neighbProcNo()
                        );

                        if (procMeshI < nbrProcMeshI)
                        {
                            // I am 'master' since get lowest numbered cells
                            nCoupledFaces[procMeshI] += ldui.faceCells().size();
                        }

                        EdgeMap<labelPairList>::iterator iter =
                            mergedMap.find(procEdge);

                        if (iter != mergedMap.end())
                        {
                            iter().append(labelPair(procMeshI, intI));
                        }
                        else
                        {
                            mergedMap.insert
                            (
                                procEdge,
                                labelPairList(1, labelPair(procMeshI, intI))
                            );
                        }
                    }
                    else
                    {
                        if (debug)
                        {
                            Pout<< "external interface: myProcNo:"
                                << pldui.myProcNo()
                                << " nbr:" << pldui.neighbProcNo()
                                << " size:" << ldui.faceCells().size()
                                << endl;
                        }

                        EdgeMap<labelPairList>::iterator iter =
                            unmergedMap.find(procEdge);

                        if (iter != unmergedMap.end())
                        {
                            iter().append(labelPair(procMeshI, intI));
                        }
                        else
                        {
                            unmergedMap.insert
                            (
                                procEdge,
                                labelPairList(1, labelPair(procMeshI, intI))
                            );
                        }
                    }
                }
                else
                {
                    // Still external (non proc) interface
                    FatalErrorIn("lduPrimitiveMesh::lduPrimitiveMesh(..)")
                        << "At mesh from processor " << procIDs[procMeshI]
                        << " have interface " << intI
                        << " of unhandled type " << interfaces[intI].type()
                        << exit(FatalError);

                    nOtherInterfaces++;
                }
            }
        }
    }



    if (debug)
    {
        Pout<< "Remaining interfaces:" << endl;
        forAllConstIter(EdgeMap<labelPairList>, unmergedMap, iter)
        {
            Pout<< "    agglom procEdge:" << iter.key() << endl;
            const labelPairList& elems = iter();
            forAll(elems, i)
            {
                label procMeshI = elems[i][0];
                label interfaceI = elems[i][1];
                const lduInterfacePtrsList interfaces =
                    mesh(myMesh, otherMeshes, procMeshI).interfaces();

                const processorLduInterface& pldui =
                    refCast<const processorLduInterface>
                    (
                        interfaces[interfaceI]
                    );

                Pout<< "        proc:" << procIDs[procMeshI]
                    << " interfaceI:" << interfaceI
                    << " between:" << pldui.myProcNo()
                    << " and:" << pldui.neighbProcNo()
                    << endl;
            }
            Pout<< endl;
        }
    }
    if (debug)
    {
        Pout<< "Merged interfaces:" << endl;
        forAllConstIter(EdgeMap<labelPairList>, mergedMap, iter)
        {
            Pout<< "    agglom procEdge:" << iter.key() << endl;
            const labelPairList& elems = iter();

            forAll(elems, i)
            {
                label procMeshI = elems[i][0];
                label interfaceI = elems[i][1];
                const lduInterfacePtrsList interfaces =
                    mesh(myMesh, otherMeshes, procMeshI).interfaces();
                const processorLduInterface& pldui =
                    refCast<const processorLduInterface>
                    (
                        interfaces[interfaceI]
                    );

                Pout<< "        proc:" << procIDs[procMeshI]
                    << " interfaceI:" << interfaceI
                    << " between:" << pldui.myProcNo()
                    << " and:" << pldui.neighbProcNo()
                    << endl;
            }
            Pout<< endl;
        }
    }


    // Adapt faceOffsets for internal interfaces
    faceOffsets.setSize(nMeshes+1);
    faceOffsets[0] = 0;
    faceMap.setSize(nMeshes);
    for (label procMeshI = 0; procMeshI < nMeshes; procMeshI++)
    {
        const lduMesh& procMesh = mesh(myMesh, otherMeshes, procMeshI);
        label nInternal = procMesh.lduAddr().lowerAddr().size();

        faceOffsets[procMeshI+1] =
            faceOffsets[procMeshI]
          + nInternal
          + nCoupledFaces[procMeshI];

        labelList& map = faceMap[procMeshI];
        map.setSize(nInternal);
        forAll(map, i)
        {
            map[i] = faceOffsets[procMeshI] + i;
        }
    }


    // Combine upper and lower
    lowerAddr_.setSize(faceOffsets.last(), -1);
    upperAddr_.setSize(lowerAddr_.size(), -1);


    // Old internal faces and resolved coupled interfaces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    for (label procMeshI = 0; procMeshI < nMeshes; procMeshI++)
    {
        const lduMesh& procMesh = mesh(myMesh, otherMeshes, procMeshI);

        const labelUList& l = procMesh.lduAddr().lowerAddr();
        const labelUList& u = procMesh.lduAddr().upperAddr();

        // Add internal faces
        label allFaceI = faceOffsets[procMeshI];

        forAll(l, faceI)
        {
            lowerAddr_[allFaceI] = cellOffsets[procMeshI]+l[faceI];
            upperAddr_[allFaceI] = cellOffsets[procMeshI]+u[faceI];
            allFaceI++;
        }


        // Add merged interfaces
        const lduInterfacePtrsList interfaces = procMesh.interfaces();

        forAll(interfaces, intI)
        {
            if (interfaces.set(intI))
            {
                const lduInterface& ldui = interfaces[intI];

                if (isA<processorLduInterface>(ldui))
                {
                    const processorLduInterface& pldui =
                        refCast<const processorLduInterface>(ldui);

                    // Look up corresponding interfaces
                    label myP = pldui.myProcNo();
                    label nbrP = pldui.neighbProcNo();
                    label nbrProcMeshI = findIndex(procIDs, nbrP);

                    if (procMeshI < nbrProcMeshI)
                    {
                        // I am 'master' since my cell numbers will be lower
                        // since cells get added in procMeshI order.

                        label agglom0 = procAgglomMap[myP];
                        label agglom1 = procAgglomMap[nbrP];

                        EdgeMap<labelPairList>::const_iterator fnd =
                            mergedMap.find(edge(agglom0, agglom1));

                        if (fnd != mergedMap.end())
                        {
                            const labelPairList& elems = fnd();

                            // Find nbrP in elems
                            label nbrIntI = -1;
                            forAll(elems, i)
                            {
                                label procI = elems[i][0];
                                label interfaceI = elems[i][1];
                                const lduInterfacePtrsList interfaces =
                                    mesh
                                    (
                                        myMesh,
                                        otherMeshes,
                                        procI
                                    ).interfaces();
                                const processorLduInterface& pldui =
                                    refCast<const processorLduInterface>
                                    (
                                        interfaces[interfaceI]
                                    );

                                if
                                (
                                    elems[i][0] == nbrProcMeshI
                                 && pldui.neighbProcNo() == procIDs[procMeshI]
                                )
                                {
                                    nbrIntI = elems[i][1];
                                    break;
                                }
                            }


                            if (nbrIntI == -1)
                            {
                                FatalErrorIn
                                (
                                    "lduPrimitiveMesh::lduPrimitiveMesh(..)"
                                )   << "elems:" << elems << abort(FatalError);
                            }


                            const lduInterfacePtrsList nbrInterfaces = mesh
                            (
                                myMesh,
                                otherMeshes,
                                nbrProcMeshI
                            ).interfaces();


                            const labelUList& faceCells =
                                ldui.faceCells();
                            const labelUList& nbrFaceCells =
                                nbrInterfaces[nbrIntI].faceCells();

                            if (faceCells.size() != nbrFaceCells.size())
                            {
                                FatalErrorIn
                                (
                                    "lduPrimitiveMesh::lduPrimitiveMesh(..)"
                                )   << "faceCells:" << faceCells
                                    << " nbrFaceCells:" << nbrFaceCells
                                    << abort(FatalError);
                            }


                            labelList& bfMap =
                                boundaryFaceMap[procMeshI][intI];
                            labelList& nbrBfMap =
                                boundaryFaceMap[nbrProcMeshI][nbrIntI];

                            bfMap.setSize(faceCells.size());
                            nbrBfMap.setSize(faceCells.size());

                            forAll(faceCells, pfI)
                            {
                                lowerAddr_[allFaceI] =
                                    cellOffsets[procMeshI]+faceCells[pfI];
                                bfMap[pfI] = allFaceI;
                                upperAddr_[allFaceI] =
                                    cellOffsets[nbrProcMeshI]+nbrFaceCells[pfI];
                                nbrBfMap[pfI] = (-allFaceI-1);
                                allFaceI++;
                            }
                        }
                    }
                }
            }
        }
    }


    // Sort upper-tri order
    {
        labelList oldToNew
        (
            upperTriOrder
            (
                cellOffsets.last(), //nCells
                lowerAddr_,
                upperAddr_
            )
        );

        forAll(faceMap, procMeshI)
        {
            labelList& map = faceMap[procMeshI];
            forAll(map, i)
            {
                if (map[i] >= 0)
                {
                    map[i] = oldToNew[map[i]];
                }
                else
                {
                    label allFaceI = -map[i]-1;
                    map[i] = -oldToNew[allFaceI]-1;
                }
            }
        }


        inplaceReorder(oldToNew, lowerAddr_);
        inplaceReorder(oldToNew, upperAddr_);

        forAll(boundaryFaceMap, procI)
        {
            const labelList& bMap = boundaryMap[procI];
            forAll(bMap, intI)
            {
                if (bMap[intI] == -1)
                {
                    // Merged interface
                    labelList& bfMap = boundaryFaceMap[procI][intI];

                    forAll(bfMap, i)
                    {
                        if (bfMap[i] >= 0)
                        {
                            bfMap[i] = oldToNew[bfMap[i]];
                        }
                        else
                        {
                            label allFaceI = -bfMap[i]-1;
                            bfMap[i] = (-oldToNew[allFaceI]-1);
                        }
                    }
                }
            }
        }
    }


    // Kept interfaces
    // ~~~~~~~~~~~~~~~

    interfaces_.setSize(unmergedMap.size() + nOtherInterfaces);
    primitiveInterfaces_.setSize(interfaces_.size());

    label allInterfaceI = 0;

    forAllConstIter(EdgeMap<labelPairList>, unmergedMap, iter)
    {
        const labelPairList& elems = iter();

        // Sort processors in increasing order so both sides walk through in
        // same order.
        labelPairList procPairs(elems.size());
        forAll(elems, i)
        {
            const labelPair& elem = elems[i];
            label procMeshI = elem[0];
            label interfaceI = elem[1];
            const lduInterfacePtrsList interfaces = mesh
            (
                myMesh,
                otherMeshes,
                procMeshI
            ).interfaces();

            const processorLduInterface& pldui =
                refCast<const processorLduInterface>
                (
                    interfaces[interfaceI]
                );
            label myProcNo = pldui.myProcNo();
            label nbrProcNo = pldui.neighbProcNo();
            procPairs[i] = labelPair
            (
                min(myProcNo, nbrProcNo),
                max(myProcNo, nbrProcNo)
            );
        }

        labelList order;
        sortedOrder(procPairs, order);

        // Count
        label n = 0;

        forAll(order, i)
        {
            const labelPair& elem = elems[order[i]];
            label procMeshI = elem[0];
            label interfaceI = elem[1];
            const lduInterfacePtrsList interfaces = mesh
            (
                myMesh,
                otherMeshes,
                procMeshI
            ).interfaces();

            n += interfaces[interfaceI].faceCells().size();
        }

        // Size
        labelField allFaceCells(n);
        labelField allFaceRestrictAddressing(n);
        n = 0;

        // Fill
        forAll(order, i)
        {
            const labelPair& elem = elems[order[i]];
            label procMeshI = elem[0];
            label interfaceI = elem[1];
            const lduInterfacePtrsList interfaces = mesh
            (
                myMesh,
                otherMeshes,
                procMeshI
            ).interfaces();

            boundaryMap[procMeshI][interfaceI] = allInterfaceI;
            labelList& bfMap = boundaryFaceMap[procMeshI][interfaceI];

            const labelUList& l = interfaces[interfaceI].faceCells();
            bfMap.setSize(l.size());

            forAll(l, faceI)
            {
                allFaceCells[n] = cellOffsets[procMeshI]+l[faceI];
                allFaceRestrictAddressing[n] = n;
                bfMap[faceI] = n;
                n++;
            }
        }


        // Find out local and remote processor in new communicator

        label neighbProcNo = -1;

        // See what the two processors map onto

        if (iter.key()[0] == myAgglom)
        {
            if (iter.key()[1] == myAgglom)
            {
                FatalErrorIn
                (
                    "lduPrimitiveMesh::lduPrimitiveMesh(..)"
                )   << "problem procEdge:" << iter.key()
                    << exit(FatalError);
            }

            neighbProcNo = iter.key()[1];
        }
        else
        {
            if (iter.key()[1] != myAgglom)
            {
                FatalErrorIn
                (
                    "lduPrimitiveMesh::lduPrimitiveMesh(..)"
                )   << "problem procEdge:" << iter.key()
                    << exit(FatalError);
            }

            neighbProcNo = iter.key()[0];
        }

        primitiveInterfaces_.set
        (
            allInterfaceI,
            new processorGAMGInterface
            (
                allInterfaceI,
                interfaces_,
                allFaceCells,
                allFaceRestrictAddressing,
                comm_,
                myAgglom,
                neighbProcNo,
                tensorField(),          // forwardT
                Pstream::msgType()      // tag
            )
        );
        interfaces_.set(allInterfaceI, &primitiveInterfaces_[allInterfaceI]);

        if (debug)
        {
            Pout<< "Created " << interfaces_[allInterfaceI].type()
                << " interface at " << allInterfaceI
                << " comm:" << comm_
                << " myProcNo:" << myAgglom
                << " neighbProcNo:" << neighbProcNo
                << " nFaces:" << allFaceCells.size()
                << endl;
        }


        allInterfaceI++;
    }


    patchSchedule_ = nonBlockingSchedule<processorGAMGInterface>(interfaces_);

    if (debug)
    {
        checkUpperTriangular(cellOffsets.last(), lowerAddr_, upperAddr_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::lduMesh& Foam::lduPrimitiveMesh::mesh
(
    const lduMesh& myMesh,
    const PtrList<lduPrimitiveMesh>& otherMeshes,
    const label meshI
)
{
    return (meshI == 0 ? myMesh : otherMeshes[meshI-1]);
}


void Foam::lduPrimitiveMesh::gather
(
    const label comm,
    const lduMesh& mesh,
    const labelList& procIDs,
    PtrList<lduPrimitiveMesh>& otherMeshes
)
{
    // Force calculation of schedule (since does parallel comms)
    (void)mesh.lduAddr().patchSchedule();

    if (Pstream::myProcNo(comm) == procIDs[0])
    {
        otherMeshes.setSize(procIDs.size()-1);

        // Slave meshes
        for (label i = 1; i < procIDs.size(); i++)
        {
            //Pout<< "on master :"
            //    << " receiving from slave " << procIDs[i] << endl;

            IPstream fromSlave
            (
                Pstream::scheduled,
                procIDs[i],
                0,          // bufSize
                Pstream::msgType(),
                comm
            );

            label nCells = readLabel(fromSlave);
            labelList lowerAddr(fromSlave);
            labelList upperAddr(fromSlave);
            boolList validInterface(fromSlave);


            // Construct mesh without interfaces
            otherMeshes.set
            (
                i-1,
                new lduPrimitiveMesh
                (
                    nCells,
                    lowerAddr,
                    upperAddr,
                    comm,
                    true    // reuse
                )
            );

            // Construct GAMGInterfaces
            lduInterfacePtrsList newInterfaces(validInterface.size());
            forAll(validInterface, intI)
            {
                if (validInterface[intI])
                {
                    word coupleType(fromSlave);

                    newInterfaces.set
                    (
                        intI,
                        GAMGInterface::New
                        (
                            coupleType,
                            intI,
                            otherMeshes[i-1].rawInterfaces(),
                            fromSlave
                        ).ptr()
                    );
                }
            }

            otherMeshes[i-1].addInterfaces
            (
                newInterfaces,
                nonBlockingSchedule<processorGAMGInterface>
                (
                    newInterfaces
                )
            );
       }
    }
    else if (findIndex(procIDs, Pstream::myProcNo(comm)) != -1)
    {
        // Send to master

        const lduAddressing& addressing = mesh.lduAddr();
        lduInterfacePtrsList interfaces(mesh.interfaces());
        boolList validInterface(interfaces.size());
        forAll(interfaces, intI)
        {
            validInterface[intI] = interfaces.set(intI);
        }

        OPstream toMaster
        (
            Pstream::scheduled,
            procIDs[0],
            0,
            Pstream::msgType(),
            comm
        );

        toMaster
            << addressing.size()
            << addressing.lowerAddr()
            << addressing.upperAddr()
            << validInterface;

        forAll(interfaces, intI)
        {
            if (interfaces.set(intI))
            {
                const GAMGInterface& interface = refCast<const GAMGInterface>
                (
                    interfaces[intI]
                );

                toMaster << interface.type();
                interface.write(toMaster);
            }
        }
    }
}


// ************************************************************************* //
