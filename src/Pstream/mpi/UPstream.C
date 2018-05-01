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

Note
    The const_cast used in this file is a temporary hack for to work around
    bugs in OpenMPI for versions < 1.7.4

\*---------------------------------------------------------------------------*/

#include "UPstream.H"
#include "PstreamReduceOps.H"
#include "OSspecific.H"
#include "PstreamGlobals.H"
#include "SubList.H"
#include "allReduce.H"

#include <mpi.h>

#include <cstring>
#include <cstdlib>
#include <csignal>

#if defined(WM_SP)
    #define MPI_SCALAR MPI_FLOAT
#elif defined(WM_DP)
    #define MPI_SCALAR MPI_DOUBLE
#elif defined(WM_LP)
    #define MPI_SCALAR MPI_LONG_DOUBLE
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// NOTE:
// valid parallel options vary between implementations, but flag common ones.
// if they are not removed by MPI_Init(), the subsequent argument processing
// will notice that they are wrong
void Foam::UPstream::addValidParOptions(HashTable<string>& validParOptions)
{
    validParOptions.insert("np", "");
    validParOptions.insert("p4pg", "PI file");
    validParOptions.insert("p4wd", "directory");
    validParOptions.insert("p4amslave", "");
    validParOptions.insert("p4yourname", "hostname");
    validParOptions.insert("machinefile", "machine file");
}


bool Foam::UPstream::init(int& argc, char**& argv, const bool needsThread)
{
    // MPI_Init(&argc, &argv);
    int provided_thread_support;
    MPI_Init_thread
    (
        &argc,
        &argv,
        (
            needsThread
          ? MPI_THREAD_MULTIPLE
          : MPI_THREAD_SINGLE
        ),
        &provided_thread_support
    );

    // int numprocs;
    // MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    // int myRank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    MPI_Comm_split
    (
        MPI_COMM_WORLD,
        1,
        myGlobalRank,
        &PstreamGlobals::MPI_COMM_FOAM
    );

    int numprocs;
    MPI_Comm_size(PstreamGlobals::MPI_COMM_FOAM, &numprocs);
    int myRank;
    MPI_Comm_rank(PstreamGlobals::MPI_COMM_FOAM, &myRank);

    if (debug)
    {
        Pout<< "UPstream::init : initialised with numProcs:" << numprocs
            << " myRank:" << myRank << endl;
    }

    if (numprocs <= 1)
    {
        FatalErrorInFunction
            << "bool IPstream::init(int& argc, char**& argv) : "
               "attempt to run parallel on 1 processor"
            << Foam::abort(FatalError);
    }


    // Initialise parallel structure
    setParRun(numprocs, provided_thread_support == MPI_THREAD_MULTIPLE);

    #ifndef SGIMPI
    string bufferSizeName = getEnv("MPI_BUFFER_SIZE");

    if (bufferSizeName.size())
    {
        int bufferSize = atoi(bufferSizeName.c_str());

        if (bufferSize)
        {
            MPI_Buffer_attach(new char[bufferSize], bufferSize);
        }
    }
    else
    {
        FatalErrorInFunction
            << "UPstream::init(int& argc, char**& argv) : "
            << "environment variable MPI_BUFFER_SIZE not defined"
            << Foam::abort(FatalError);
    }
    #endif

    // int processorNameLen;
    // char processorName[MPI_MAX_PROCESSOR_NAME];
    //
    // MPI_Get_processor_name(processorName, &processorNameLen);
    // processorName[processorNameLen] = '\0';
    // Pout<< "Processor name:" << processorName << endl;

    return true;
}


void Foam::UPstream::exit(int errnum)
{
    if (debug)
    {
        Pout<< "UPstream::exit." << endl;
    }

    #ifndef SGIMPI
    int size;
    char* buff;
    MPI_Buffer_detach(&buff, &size);
    delete[] buff;
    #endif

    if (PstreamGlobals::outstandingRequests_.size())
    {
        label n = PstreamGlobals::outstandingRequests_.size();
        PstreamGlobals::outstandingRequests_.clear();

        WarningInFunction
            << "There are still " << n << " outstanding MPI_Requests." << endl
            << "This means that your code exited before doing a"
            << " UPstream::waitRequests()." << endl
            << "This should not happen for a normal code exit."
            << endl;
    }

    // Clean mpi communicators
    forAll(myProcNo_, communicator)
    {
        if (myProcNo_[communicator] != -1)
        {
            freePstreamCommunicator(communicator);
        }
    }

    if (errnum == 0)
    {
        MPI_Finalize();
        ::exit(errnum);
    }
    else
    {
        MPI_Abort(PstreamGlobals::MPI_COMM_FOAM, errnum);
    }
}


void Foam::UPstream::abort()
{
    MPI_Abort(PstreamGlobals::MPI_COMM_FOAM, 1);
}


void Foam::reduce
(
    scalar& Value,
    const sumOp<scalar>& bop,
    const int tag,
    const label communicator
)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }
    allReduce(Value, 1, MPI_SCALAR, MPI_SUM, bop, tag, communicator);
}


void Foam::reduce
(
    scalar& Value,
    const minOp<scalar>& bop,
    const int tag,
    const label communicator
)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }
    allReduce(Value, 1, MPI_SCALAR, MPI_MIN, bop, tag, communicator);
}


void Foam::reduce
(
    vector2D& Value,
    const sumOp<vector2D>& bop,
    const int tag,
    const label communicator
)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }
    allReduce(Value, 2, MPI_SCALAR, MPI_SUM, bop, tag, communicator);
}


void Foam::sumReduce
(
    scalar& Value,
    label& Count,
    const int tag,
    const label communicator
)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }
    vector2D twoScalars(Value, scalar(Count));
    reduce(twoScalars, sumOp<vector2D>(), tag, communicator);

    Value = twoScalars.x();
    Count = twoScalars.y();
}


void Foam::reduce
(
    scalar& Value,
    const sumOp<scalar>& bop,
    const int tag,
    const label communicator,
    label& requestID
)
{
#ifdef MPIX_COMM_TYPE_SHARED
    // Assume mpich2 with non-blocking collectives extensions. Once mpi3
    // is available this will change.
    MPI_Request request;
    scalar v = Value;
    MPIX_Ireduce
    (
        &v,
        &Value,
        1,
        MPI_SCALAR,
        MPI_SUM,
        0,              // root
        PstreamGlobals::MPICommunicators_[communicator],
        &request
    );

    requestID = PstreamGlobals::outstandingRequests_.size();
    PstreamGlobals::outstandingRequests_.append(request);

    if (UPstream::debug)
    {
        Pout<< "UPstream::allocateRequest for non-blocking reduce"
            << " : request:" << requestID
            << endl;
    }
#else
    // Non-blocking not yet implemented in mpi
    reduce(Value, bop, tag, communicator);
    requestID = -1;
#endif
}


void Foam::UPstream::allToAll
(
    const labelUList& sendData,
    labelUList& recvData,
    const label communicator
)
{
    label np = nProcs(communicator);

    if (sendData.size() != np || recvData.size() != np)
    {
        FatalErrorInFunction
            << "Size of sendData " << sendData.size()
            << " or size of recvData " << recvData.size()
            << " is not equal to the number of processors in the domain "
            << np
            << Foam::abort(FatalError);
    }

    if (!UPstream::parRun())
    {
        recvData.deepCopy(sendData);
    }
    else
    {
        if
        (
            MPI_Alltoall
            (
                const_cast<label*>(sendData.begin()),
                sizeof(label),
                MPI_BYTE,
                recvData.begin(),
                sizeof(label),
                MPI_BYTE,
                PstreamGlobals::MPICommunicators_[communicator]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Alltoall failed for " << sendData
                << " on communicator " << communicator
                << Foam::abort(FatalError);
        }
    }
}


void Foam::UPstream::allToAll
(
    const char* sendData,
    const UList<int>& sendSizes,
    const UList<int>& sendOffsets,

    char* recvData,
    const UList<int>& recvSizes,
    const UList<int>& recvOffsets,

    const label communicator
)
{
    label np = nProcs(communicator);

    if
    (
        sendSizes.size() != np
     || sendOffsets.size() != np
     || recvSizes.size() != np
     || recvOffsets.size() != np
    )
    {
        FatalErrorInFunction
            << "Size of sendSize " << sendSizes.size()
            << ", sendOffsets " << sendOffsets.size()
            << ", recvSizes " << recvSizes.size()
            << " or recvOffsets " << recvOffsets.size()
            << " is not equal to the number of processors in the domain "
            << np
            << Foam::abort(FatalError);
    }

    if (!UPstream::parRun())
    {
        if (recvSizes[0] != sendSizes[0])
        {
            FatalErrorInFunction
                << "Bytes to send " << sendSizes[0]
                << " does not equal bytes to receive " << recvSizes[0]
                << Foam::abort(FatalError);
        }
        memmove(recvData, &sendData[sendOffsets[0]], recvSizes[0]);
    }
    else
    {
        if
        (
            MPI_Alltoallv
            (
                const_cast<char*>(sendData),
                const_cast<int*>(sendSizes.begin()),
                const_cast<int*>(sendOffsets.begin()),
                MPI_BYTE,
                recvData,
                const_cast<int*>(recvSizes.begin()),
                const_cast<int*>(recvOffsets.begin()),
                MPI_BYTE,
                PstreamGlobals::MPICommunicators_[communicator]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Alltoallv failed for sendSizes " << sendSizes
                << " recvSizes " << recvSizes
                << " communicator " << communicator
                << Foam::abort(FatalError);
        }
    }
}


void Foam::UPstream::gather
(
    const char* sendData,
    int sendSize,

    char* recvData,
    const UList<int>& recvSizes,
    const UList<int>& recvOffsets,
    const label communicator
)
{
    label np = nProcs(communicator);

    if
    (
        UPstream::master(communicator)
     && (recvSizes.size() != np || recvOffsets.size() < np)
    )
    {
        // Note: allow recvOffsets to be e.g. 1 larger than np so we
        // can easily loop over the result

        FatalErrorInFunction
            << "Size of recvSizes " << recvSizes.size()
            << " or recvOffsets " << recvOffsets.size()
            << " is not equal to the number of processors in the domain "
            << np
            << Foam::abort(FatalError);
    }

    if (!UPstream::parRun())
    {
        memmove(recvData, sendData, sendSize);
    }
    else
    {
        if
        (
            MPI_Gatherv
            (
                const_cast<char*>(sendData),
                sendSize,
                MPI_BYTE,
                recvData,
                const_cast<int*>(recvSizes.begin()),
                const_cast<int*>(recvOffsets.begin()),
                MPI_BYTE,
                0,
                MPI_Comm(PstreamGlobals::MPICommunicators_[communicator])
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Gatherv failed for sendSize " << sendSize
                << " recvSizes " << recvSizes
                << " communicator " << communicator
                << Foam::abort(FatalError);
        }
    }
}


void Foam::UPstream::scatter
(
    const char* sendData,
    const UList<int>& sendSizes,
    const UList<int>& sendOffsets,

    char* recvData,
    int recvSize,
    const label communicator
)
{
    label np = nProcs(communicator);

    if
    (
        UPstream::master(communicator)
     && (sendSizes.size() != np || sendOffsets.size() != np)
    )
    {
        FatalErrorInFunction
            << "Size of sendSizes " << sendSizes.size()
            << " or sendOffsets " << sendOffsets.size()
            << " is not equal to the number of processors in the domain "
            << np
            << Foam::abort(FatalError);
    }

    if (!UPstream::parRun())
    {
        memmove(recvData, sendData, recvSize);
    }
    else
    {
        if
        (
            MPI_Scatterv
            (
                const_cast<char*>(sendData),
                const_cast<int*>(sendSizes.begin()),
                const_cast<int*>(sendOffsets.begin()),
                MPI_BYTE,
                recvData,
                recvSize,
                MPI_BYTE,
                0,
                MPI_Comm(PstreamGlobals::MPICommunicators_[communicator])
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Scatterv failed for sendSizes " << sendSizes
                << " sendOffsets " << sendOffsets
                << " communicator " << communicator
                << Foam::abort(FatalError);
        }
    }
}


void Foam::UPstream::allocatePstreamCommunicator
(
    const label parentIndex,
    const label index
)
{
    if (index == PstreamGlobals::MPIGroups_.size())
    {
        // Extend storage with dummy values
        MPI_Group newGroup = MPI_GROUP_NULL;
        PstreamGlobals::MPIGroups_.append(newGroup);
        MPI_Comm newComm = MPI_COMM_NULL;
        PstreamGlobals::MPICommunicators_.append(newComm);
    }
    else if (index > PstreamGlobals::MPIGroups_.size())
    {
        FatalErrorInFunction
            << "PstreamGlobals out of sync with UPstream data. Problem."
            << Foam::exit(FatalError);
    }


    if (parentIndex == -1)
    {
        // Allocate world communicator

        if (index != UPstream::worldComm)
        {
            FatalErrorInFunction
                << "world communicator should always be index "
                << UPstream::worldComm << Foam::exit(FatalError);
        }

        PstreamGlobals::MPICommunicators_[index] =
            PstreamGlobals::MPI_COMM_FOAM;
        MPI_Comm_group
        (
            PstreamGlobals::MPI_COMM_FOAM,
            &PstreamGlobals::MPIGroups_[index]
        );
        MPI_Comm_rank
        (
            PstreamGlobals::MPICommunicators_[index],
           &myProcNo_[index]
        );

        // Set the number of processes to the actual number
        int numProcs;
        MPI_Comm_size(PstreamGlobals::MPICommunicators_[index], &numProcs);

        // procIDs_[index] = identity(numProcs);
        procIDs_[index].setSize(numProcs);
        forAll(procIDs_[index], i)
        {
            procIDs_[index][i] = i;
        }
    }
    else
    {
        // Create new group
        MPI_Group_incl
        (
            PstreamGlobals::MPIGroups_[parentIndex],
            procIDs_[index].size(),
            procIDs_[index].begin(),
           &PstreamGlobals::MPIGroups_[index]
        );

        // Create new communicator
        MPI_Comm_create
        (
            PstreamGlobals::MPICommunicators_[parentIndex],
            PstreamGlobals::MPIGroups_[index],
           &PstreamGlobals::MPICommunicators_[index]
        );

        if (PstreamGlobals::MPICommunicators_[index] == MPI_COMM_NULL)
        {
            myProcNo_[index] = -1;
        }
        else
        {
            if
            (
                MPI_Comm_rank
                (
                    PstreamGlobals::MPICommunicators_[index],
                   &myProcNo_[index]
                )
            )
            {
                FatalErrorInFunction
                    << "Problem :"
                    << " when allocating communicator at " << index
                    << " from ranks " << procIDs_[index]
                    << " of parent " << parentIndex
                    << " cannot find my own rank"
                    << Foam::exit(FatalError);
            }
        }
    }
}


void Foam::UPstream::freePstreamCommunicator(const label communicator)
{
    if (communicator != UPstream::worldComm)
    {
        if (PstreamGlobals::MPICommunicators_[communicator] != MPI_COMM_NULL)
        {
            // Free communicator. Sets communicator to MPI_COMM_NULL
            MPI_Comm_free(&PstreamGlobals::MPICommunicators_[communicator]);
        }
        if (PstreamGlobals::MPIGroups_[communicator] != MPI_GROUP_NULL)
        {
            // Free greoup. Sets group to MPI_GROUP_NULL
            MPI_Group_free(&PstreamGlobals::MPIGroups_[communicator]);
        }
    }
}


Foam::label Foam::UPstream::nRequests()
{
    return PstreamGlobals::outstandingRequests_.size();
}


void Foam::UPstream::resetRequests(const label i)
{
    if (i < PstreamGlobals::outstandingRequests_.size())
    {
        PstreamGlobals::outstandingRequests_.setSize(i);
    }
}


void Foam::UPstream::waitRequests(const label start)
{
    if (debug)
    {
        Pout<< "UPstream::waitRequests : starting wait for "
            << PstreamGlobals::outstandingRequests_.size()-start
            << " outstanding requests starting at " << start << endl;
    }

    if (PstreamGlobals::outstandingRequests_.size())
    {
        SubList<MPI_Request> waitRequests
        (
            PstreamGlobals::outstandingRequests_,
            PstreamGlobals::outstandingRequests_.size() - start,
            start
        );

        if
        (
            MPI_Waitall
            (
                waitRequests.size(),
                waitRequests.begin(),
                MPI_STATUSES_IGNORE
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Waitall returned with error" << Foam::endl;
        }

        resetRequests(start);
    }

    if (debug)
    {
        Pout<< "UPstream::waitRequests : finished wait." << endl;
    }
}


void Foam::UPstream::waitRequest(const label i)
{
    if (debug)
    {
        Pout<< "UPstream::waitRequest : starting wait for request:" << i
            << endl;
    }

    if (i >= PstreamGlobals::outstandingRequests_.size())
    {
        FatalErrorInFunction
            << "There are " << PstreamGlobals::outstandingRequests_.size()
            << " outstanding send requests and you are asking for i=" << i
            << nl
            << "Maybe you are mixing blocking/non-blocking comms?"
            << Foam::abort(FatalError);
    }

    if
    (
        MPI_Wait
        (
           &PstreamGlobals::outstandingRequests_[i],
            MPI_STATUS_IGNORE
        )
    )
    {
        FatalErrorInFunction
            << "MPI_Wait returned with error" << Foam::endl;
    }

    if (debug)
    {
        Pout<< "UPstream::waitRequest : finished wait for request:" << i
            << endl;
    }
}


bool Foam::UPstream::finishedRequest(const label i)
{
    if (debug)
    {
        Pout<< "UPstream::finishedRequest : checking request:" << i
            << endl;
    }

    if (i >= PstreamGlobals::outstandingRequests_.size())
    {
        FatalErrorInFunction
            << "There are " << PstreamGlobals::outstandingRequests_.size()
            << " outstanding send requests and you are asking for i=" << i
            << nl
            << "Maybe you are mixing blocking/non-blocking comms?"
            << Foam::abort(FatalError);
    }

    int flag;
    MPI_Test
    (
       &PstreamGlobals::outstandingRequests_[i],
       &flag,
        MPI_STATUS_IGNORE
    );

    if (debug)
    {
        Pout<< "UPstream::finishedRequest : finished request:" << i
            << endl;
    }

    return flag != 0;
}


int Foam::UPstream::allocateTag(const char* s)
{
    int tag;
    if (PstreamGlobals::freedTags_.size())
    {
        tag = PstreamGlobals::freedTags_.remove();
    }
    else
    {
        tag = PstreamGlobals::nTags_++;
    }

    if (debug)
    {
        // if (UPstream::lateBlocking > 0)
        //{
        //    string& poutp = Pout.prefix();
        //    poutp[poutp.size()-2*(UPstream::lateBlocking+2)+tag] = 'X';
        //    Perr.prefix() = Pout.prefix();
        //}
        Pout<< "UPstream::allocateTag " << s
            << " : tag:" << tag
            << endl;
    }

    return tag;
}


int Foam::UPstream::allocateTag(const word& s)
{
    int tag;
    if (PstreamGlobals::freedTags_.size())
    {
        tag = PstreamGlobals::freedTags_.remove();
    }
    else
    {
        tag = PstreamGlobals::nTags_++;
    }

    if (debug)
    {
        // if (UPstream::lateBlocking > 0)
        //{
        //    string& poutp = Pout.prefix();
        //    poutp[poutp.size()-2*(UPstream::lateBlocking+2)+tag] = 'X';
        //    Perr.prefix() = Pout.prefix();
        //}
        Pout<< "UPstream::allocateTag " << s
            << " : tag:" << tag
            << endl;
    }

    return tag;
}


void Foam::UPstream::freeTag(const char* s, const int tag)
{
    if (debug)
    {
        // if (UPstream::lateBlocking > 0)
        //{
        //    string& poutp = Pout.prefix();
        //    poutp[poutp.size()-2*(UPstream::lateBlocking+2)+tag] = ' ';
        //    Perr.prefix() = Pout.prefix();
        //}
        Pout<< "UPstream::freeTag " << s << " tag:" << tag << endl;
    }
    PstreamGlobals::freedTags_.append(tag);
}


void Foam::UPstream::freeTag(const word& s, const int tag)
{
    if (debug)
    {
        // if (UPstream::lateBlocking > 0)
        //{
        //    string& poutp = Pout.prefix();
        //    poutp[poutp.size()-2*(UPstream::lateBlocking+2)+tag] = ' ';
        //    Perr.prefix() = Pout.prefix();
        //}
        Pout<< "UPstream::freeTag " << s << " tag:" << tag << endl;
    }
    PstreamGlobals::freedTags_.append(tag);
}


// ************************************************************************* //
