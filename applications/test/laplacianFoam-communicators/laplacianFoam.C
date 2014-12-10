/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    laplacianFoam

Description
    Solves a simple Laplace equation, e.g. for thermal diffusion in a solid.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"
#include "globalIndex.H"
#include "lduPrimitiveMesh.H"
#include "processorGAMGInterface.H"
#include "GAMGInterfaceField.H"
#include "processorLduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void checkUpperTriangular
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
                << " l:" << l[faceI] << " u:" << u[faceI] << abort(FatalError);
        }
        if (l[faceI] < 0 || u[faceI] < 0 || u[faceI] >= size)
        {
            FatalErrorIn
            (
                "checkUpperTriangular"
                "(const label, const labelUList&, const labelUList&)"
            )   << "Illegal cell label. Problem at face " << faceI
                << " l:" << l[faceI] << " u:" << u[faceI] << abort(FatalError);
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
                << " previous l:" << l[faceI-1] << abort(FatalError);
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
                    << " previous u:" << u[faceI-1] << abort(FatalError);
            }
        }
    }
}


void print(const string& msg, const lduMesh& mesh)
{
    const lduAddressing& addressing = mesh.lduAddr();
    const lduInterfacePtrsList interfaces = mesh.interfaces();

    Pout<< "Mesh:" << msg.c_str() << nl
        << "    cells:" << addressing.size() << nl
        << "    faces:" << addressing.lowerAddr().size() << nl
        << "    patches:" << interfaces.size() << nl;


    const labelUList& l = addressing.lowerAddr();
    const labelUList& startAddr = addressing.losortStartAddr();
    const labelUList& addr = addressing.losortAddr();

    forAll(addressing, cellI)
    {
        Pout<< "    cell:" << cellI << nl;

        label start = startAddr[cellI];
        label end = startAddr[cellI+1];

        for (label index = start; index < end; index++)
        {
            Pout<< "        nbr:" << l[addr[index]] << nl;
        }
    }

    Pout<< "    Patches:" << nl;
    forAll(interfaces, i)
    {
        if (interfaces.set(i))
        {
            if (isA<processorLduInterface>(interfaces[i]))
            {
                const processorLduInterface& pldui =
                    refCast<const processorLduInterface>(interfaces[i]);
                Pout<< "        " << i
                    << " me:" << pldui.myProcNo()
                    << " nbr:" << pldui.neighbProcNo()
                    << " comm:" << pldui.comm()
                    << " tag:" << pldui.tag()
                    << nl;
            }

            {
                Pout<< "        " << i << " addressing:" << nl;
                const labelUList& faceCells = interfaces[i].faceCells();
                forAll(faceCells, i)
                {
                    Pout<< "\t\t" << i << '\t' << faceCells[i] << nl;
                }
            }
        }
    }
}


template<class ProcPatch>
lduSchedule nonBlockingSchedule
(
    const lduInterfacePtrsList& interfaces
)
{
    lduSchedule schedule(2*interfaces.size());
    label slotI = 0;

    forAll(interfaces, i)
    {
        if (interfaces.set(i) && !isA<ProcPatch>(interfaces[i]))
        {
            schedule[slotI].patch = i;
            schedule[slotI].init = true;
            slotI++;
            schedule[slotI].patch = i;
            schedule[slotI].init = false;
            slotI++;
        }
    }

    forAll(interfaces, i)
    {
        if (interfaces.set(i) && isA<ProcPatch>(interfaces[i]))
        {
            schedule[slotI].patch = i;
            schedule[slotI].init = true;
            slotI++;
        }
    }

    forAll(interfaces, i)
    {
        if (interfaces.set(i) && isA<ProcPatch>(interfaces[i]))
        {
            schedule[slotI].patch = i;
            schedule[slotI].init = false;
            slotI++;
        }
    }

    return schedule;
}


void sendReceive
(
    const label comm,
    const label tag,
    const globalIndex& offsets,
    const scalarField& field,

    scalarField& allField
)
{
    label nProcs = Pstream::nProcs(comm);

    if (Pstream::master(comm))
    {
        allField.setSize(offsets.size());

        // Assign master slot
        SubList<scalar>
        (
            allField,
            offsets.localSize(0),
            offsets.offset(0)
        ).assign(field);

        // Assign slave slots
        for (label procI = 1; procI < nProcs; procI++)
        {
            SubList<scalar> procSlot
            (
                allField,
                offsets.localSize(procI),
                offsets.offset(procI)
            );

            Pout<< "Receiving allField from " << procI
                << " at offset:" << offsets.offset(procI)
                << " size:" << offsets.size()
                << endl;

            IPstream::read
            (
                Pstream::nonBlocking,
                procI,
                reinterpret_cast<char*>(procSlot.begin()),
                procSlot.byteSize(),
                tag,
                comm
            );
        }
    }
    else
    {
        OPstream::write
        (
            Pstream::nonBlocking,
            0,                          // master
            reinterpret_cast<const char*>(field.begin()),
            field.byteSize(),
            tag,
            comm
        );
    }
}


void sendReceive
(
    const label comm,
    const label tag,
    const globalIndex& offsets,
    const FieldField<Field, scalar>& field,

    FieldField<Field, scalar>& allField
)
{
    PstreamBuffers pBufs(Pstream::nonBlocking, Pstream::msgType(), comm);

    if (!Pstream::master(comm))
    {
        UOPstream toMaster(Pstream::masterNo(), pBufs);

        Pout<< "To 0 sending " << field.size()
            << " fields." << endl;

        forAll(field, intI)
        {
            toMaster << field[intI];
        }
    }
    pBufs.finishedSends();
    if (Pstream::master(comm))
    {
        allField.setSize(offsets.size());
        forAll(allField, i)
        {
            allField.set(i, new scalarField(0));
        }

        // Insert master values
        forAll(field, intI)
        {
            allField[intI] = field[intI];
        }


        // Receive and insert slave values
        label nProcs = Pstream::nProcs(comm);

        for (label procI = 1; procI < nProcs; procI++)
        {
            UIPstream fromSlave(procI, pBufs);

            label nSlaveInts = offsets.localSize(procI);

            Pout<< "From " << procI << " receiving "
                << nSlaveInts << " fields." << endl;

            for (label intI = 0; intI < nSlaveInts; intI++)
            {
                label slotI = offsets.toGlobal(procI, intI);

                Pout<< "    int:" << intI << " goes into slot " << slotI
                    << endl;

                fromSlave >> allField[slotI];
            }
        }
    }
}



void collect
(
    const label comm,
    const globalIndex& cellOffsets,
    const globalIndex& faceOffsets,

    const scalarField& diagonal,
    const scalarField& upper,
    const scalarField& lower,

    scalarField& allDiagonal,
    scalarField& allUpper,
    scalarField& allLower
)
{
    label nOutstanding = Pstream::nRequests();
    int allDiagonalTag = Pstream::allocateTag("allDiagonal:" __FILE__);
    int allUpperTag = Pstream::allocateTag("allUpper:" __FILE__);
    int allLowerTag = Pstream::allocateTag("allLower:" __FILE__);


    sendReceive
    (
        comm,
        allDiagonalTag,
        cellOffsets,
        diagonal,
        allDiagonal
    );

    sendReceive
    (
        comm,
        allUpperTag,
        faceOffsets,
        upper,
        allUpper
    );

    sendReceive
    (
        comm,
        allLowerTag,
        faceOffsets,
        lower,
        allLower
    );

    Pstream::waitRequests(nOutstanding);

    Pstream::freeTag("allDiagonal:" __FILE__, allDiagonalTag);
    Pstream::freeTag("allUpper:" __FILE__, allUpperTag);
    Pstream::freeTag("allLower:" __FILE__, allLowerTag);
}


void setCommunicator(fvMesh& mesh, const label newComm)
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // The current state is consistent with the mesh so check where the new
    // communicator is and adjust accordingly.

    forAll(pbm, patchI)
    {
        if (isA<processorPolyPatch>(pbm[patchI]))
        {
            processorPolyPatch& ppp = const_cast<processorPolyPatch&>
            (
                refCast
                <
                    const processorPolyPatch
                >(pbm[patchI])
            );

            label thisRank = UPstream::procNo
            (
                newComm,
                ppp.comm(),
                ppp.myProcNo()
            );
            label nbrRank = UPstream::procNo
            (
                newComm,
                ppp.comm(),
                ppp.neighbProcNo()
            );

            //ppp.comm() = newComm;
            ppp.myProcNo() = thisRank;
            ppp.neighbProcNo() = nbrRank;
        }
    }
    mesh.polyMesh::comm() = newComm;
}


namespace Foam
{
    typedef UPtrList<const GAMGInterfaceField> GAMGInterfaceFieldPtrsList;
}


// Gather matrices from processors procIDs[1..] on procIDs[0]
void gatherMatrices
(
    const labelList& procIDs,
    const PtrList<lduMesh>& procMeshes,

    const lduMatrix& mat,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,

    PtrList<lduMatrix>& otherMats,
    PtrList<FieldField<Field, scalar> >& otherBouCoeffs,
    PtrList<FieldField<Field, scalar> >& otherIntCoeffs,
    PtrList<GAMGInterfaceFieldPtrsList>& otherInterfaces
)
{
    const label meshComm = mat.mesh().comm();

    //lduInterfacePtrsList interfaces(mesh.interfaces());

    if (Pstream::myProcNo(meshComm) == procIDs[0])
    {
        // Master.
        otherMats.setSize(procIDs.size()-1);
        otherBouCoeffs.setSize(procIDs.size()-1);
        otherIntCoeffs.setSize(procIDs.size()-1);
        otherInterfaces.setSize(procIDs.size()-1);

        for (label i = 1; i < procIDs.size(); i++)
        {
            const lduMesh& procMesh = procMeshes[i];
            const lduInterfacePtrsList procInterfaces = procMesh.interfaces();

            IPstream fromSlave
            (
                Pstream::scheduled,
                procIDs[i],
                0,          // bufSize
                Pstream::msgType(),
                meshComm
            );

            otherMats.set(i-1, new lduMatrix(procMesh, fromSlave));

            // Receive number of/valid interfaces
            boolList validTransforms(fromSlave);
            List<int> validRanks(fromSlave);

            // Size coefficients
            otherBouCoeffs.set
            (
                i-1,
                new FieldField<Field, scalar>(validTransforms.size())
            );
            otherIntCoeffs.set
            (
                i-1,
                new FieldField<Field, scalar>(validTransforms.size())
            );
            otherInterfaces.set
            (
                i-1,
                new GAMGInterfaceFieldPtrsList(validTransforms.size())
            );

            forAll(validTransforms, intI)
            {
                if (validTransforms[intI])
                {
                    const processorGAMGInterface& interface =
                        refCast<const processorGAMGInterface>
                        (
                            procInterfaces[intI]
                        );


                    otherBouCoeffs[i-1].set(intI, new scalarField(fromSlave));
                    otherIntCoeffs[i-1].set(intI, new scalarField(fromSlave));
                    otherInterfaces[i-1].set
                    (
                        intI,
                        GAMGInterfaceField::New
                        (
                            interface,  //procInterfaces[intI],
                            validTransforms[intI],
                            validRanks[intI]
                        ).ptr()
                    );
                }
            }
        }
    }
    else
    {
        // Send to master
        OPstream toMaster
        (
            Pstream::scheduled,
            procIDs[0],
            0,
            Pstream::msgType(),
            meshComm
        );


        // Count valid interfaces
        boolList validTransforms(interfaceBouCoeffs.size(), false);
        List<int> validRanks(interfaceBouCoeffs.size(), -1);
        forAll(interfaces, intI)
        {
            if (interfaces.set(intI))
            {
                const processorLduInterfaceField& interface =
                    refCast<const processorLduInterfaceField>
                    (
                        interfaces[intI]
                    );

                validTransforms[intI] = interface.doTransform();
                validRanks[intI] = interface.rank();
            }
        }

        toMaster << mat << validTransforms << validRanks;
        forAll(validTransforms, intI)
        {
            if (validTransforms[intI])
            {
                toMaster
                    << interfaceBouCoeffs[intI]
                    << interfaceIntCoeffs[intI];
            }
        }
    }
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    simpleControl simple(mesh);

    //const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;


    // Get a subset of processors
    labelList subProcs(3);
    subProcs[0] = 0;
    subProcs[1] = 1;
    subProcs[2] = 2;


    const UPstream::communicator newComm
    (
        UPstream::worldComm,
        subProcs,
        true
    );


    Pout<< "procIDs world  :" << UPstream::procID(UPstream::worldComm) << endl;
    Pout<< "procIDs newComm:" << UPstream::procID(newComm) << endl;


//// On ALL processors: get the interfaces:
//lduInterfacePtrsList interfaces(mesh.interfaces());
//PtrList<lduMesh> procMeshes;
//
//if (Pstream::myProcNo(newComm) != -1)
//{
//    print("InitialMesh", mesh);
//
//    labelList procIDs(3);
//    procIDs[0] = 0;
//    procIDs[1] = 1;
//    procIDs[2] = 2;
//
////XXXXXX
//    // Collect meshes from procs 0,1 (in newComm) on 1.
//    lduPrimitiveMesh::gather(mesh, procIDs, procMeshes);
//
//    if (Pstream::myProcNo(newComm) == procIDs[0])
//    {
//        // Print a bit
//        forAll(procMeshes, i)
//        {
//            const lduMesh& pMesh = procMeshes[i];
//            print("procMesh" + Foam::name(i), pMesh);
//
//            const lduAddressing& addr = pMesh.lduAddr();
//            checkUpperTriangular
//            (
//                addr.size(),
//                addr.lowerAddr(),
//                addr.upperAddr()
//            );
//        }
//
//
//        // Combine
//        labelList cellOffsets;
//        labelListList faceMap;
//        labelListList boundaryMap;
//        labelListListList boundaryFaceMap;
//        //autoPtr<lduPrimitiveMesh> allMeshPtr = combineMeshes
//        //(
//        //    newComm,
//        //    procIDs,
//        //    procMeshes,
//        //
//        //    cellOffsets,        // per mesh the starting cell
//        //    faceMap,            // per mesh the starting face
//        //    boundaryMap,        // per mesh,per interface the starting face
//        //    boundaryFaceMap
//        //);
//        //const lduPrimitiveMesh& allMesh = allMeshPtr();
//        const lduPrimitiveMesh allMesh
//        (
//            newComm,
//            procIDs,
//            procMeshes,
//
//            cellOffsets,
//            faceMap,
//            boundaryMap,
//            boundaryFaceMap
//        );
//
//
//        print("ALLMESH", allMesh);
//
//        forAll(procMeshes, procMeshI)
//        {
//            const lduMesh& pMesh = procMeshes[procMeshI];
//            //const lduAddressing& pAddressing = pMesh.lduAddr();
//
//            Pout<< "procMesh" << procMeshI << endl
//                << "    cells start at:" << cellOffsets[procMeshI] << endl
//                << "    faces to to:" << faceMap[procMeshI] << endl;
//
//            lduInterfacePtrsList interfaces = pMesh.interfaces();
//            forAll(interfaces, intI)
//            {
//                Pout<< "    patch:" << intI
//                    << " becomes patch:" << boundaryMap[procMeshI][intI]
//                    << endl;
//
//                Pout<< "    patch:" << intI
//                    << " faces become faces:"
//                    << boundaryFaceMap[procMeshI][intI]
//                    << endl;
//            }
//        }
//    }
//
//
//    // Construct test data
//    // ~~~~~~~~~~~~~~~~~~~
//
//    GAMGInterfaceFieldPtrsList interfaces(interfaces.size());
//    FieldField<Field, scalar> interfaceBouCoeffs(interfaces.size());
//    FieldField<Field, scalar> interfaceIntCoeffs(interfaces.size());
//
//    forAll(interfaces, intI)
//    {
//        if (interfaces.set(intI))
//        {
//            label size = interfaces[intI].size();
//
//            interfaces.set
//            (
//                intI,
//                GAMGInterfaceField::New
//                (
//                    mesh.interfaces()[intI],
//                    interfaces[intI]
//                )
//            );
//            interfaceBouCoeffs.set(intI, new scalarField(size, 111));
//            interfaceIntCoeffs.set(intI, new scalarField(size, 222));
//        }
//    }
//
//
//    PtrList<lduMatrix> otherMats;
//    PtrList<FieldField<Field, scalar> > otherBouCoeffs;
//    PtrList<FieldField<Field, scalar> > otherIntCoeffs;
//    PtrList<GAMGInterfaceFieldPtrsList> otherInterfaces;
//    gatherMatrices
//    (
//        procIDs,
//        procMeshes,
//
//        mat,
//        interfaceBouCoeffs,
//        interfaceIntCoeffs,
//        interfaces,
//
//        otherMats,
//        otherBouCoeffs,
//        otherIntCoeffs,
//        otherInterfaces
//    );
////XXXXXX
//}



    {
        Pout<< "World:" << UPstream::worldComm
            << " procID:" << 2
            << " subComm:" << newComm
            << " rank1:" << UPstream::procNo(newComm, UPstream::worldComm, 1)
            << " rank2:" << UPstream::procNo(newComm, UPstream::worldComm, 2)
            << endl;
    }


    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix Teqn
            (
                //fvm::ddt(T) - fvm::laplacian(DT, T)
                fvm::laplacian(DT, T)
            );


            {
                label oldWarn = UPstream::warnComm;
                UPstream::warnComm = newComm;

                label oldComm = mesh.comm();
                setCommunicator(mesh, newComm);
                Pout<< "** oldcomm:" << oldComm
                    << "  newComm:" << mesh.comm() << endl;

                if (Pstream::myProcNo(mesh.comm()) != -1)
                {
                    solve(Teqn);
                }

                setCommunicator(mesh, oldComm);
                Pout<< "** reset mesh to:" << mesh.comm() << endl;

                UPstream::warnComm = oldWarn;
            }
        }

        #include "write.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Pout<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
