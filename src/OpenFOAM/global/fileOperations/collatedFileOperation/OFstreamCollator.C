/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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

#include "OFstreamCollator.H"
#include "OFstream.H"
#include "decomposedBlockData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(OFstreamCollator, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::OFstreamCollator::writeFile
(
    const label comm,
    const word& typeName,
    const fileName& fName,
    const string& masterData,
    const labelUList& recvSizes,
    const bool haveSlaveData,           // does master have slaveData
    const UList<char>& slaveData,       // on master: slave data
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool append
)
{
    if (debug)
    {
        Pout<< "OFstreamCollator : Writing " << masterData.size()
            << " bytes to " << fName
            << " using comm " << comm << endl;
    }

    autoPtr<OSstream> osPtr;
    if (UPstream::master(comm))
    {
        Foam::mkDir(fName.path());
        osPtr.reset
        (
            new OFstream
            (
                fName,
                fmt,
                ver,
                cmp,
                append
            )
        );

        // We don't have IOobject so cannot use IOobject::writeHeader
        OSstream& os = osPtr();
        decomposedBlockData::writeHeader
        (
            os,
            ver,
            fmt,
            typeName,
            "",
            fName,
            fName.name()
        );
    }


    UList<char> slice
    (
        const_cast<char*>(masterData.data()),
        label(masterData.size())
    );

    // Assuming threaded writing hides any slowness so we
    // can use scheduled communication to send the data to
    // the master processor in order. However can be unstable
    // for some mpi so default is non-blocking.

    List<std::streamoff> start;
    decomposedBlockData::writeBlocks
    (
        comm,
        osPtr,
        start,
        slice,
        recvSizes,
        haveSlaveData,
        slaveData,
        UPstream::commsTypes::nonBlocking,  //scheduled,
        false       // do not reduce return state
    );

    if (osPtr.valid() && !osPtr().good())
    {
        FatalIOErrorInFunction(osPtr())
            << "Failed writing to " << fName << exit(FatalIOError);
    }

    if (debug)
    {
        Pout<< "OFstreamCollator : Finished writing " << masterData.size()
            << " bytes";
        if (UPstream::master(comm))
        {
            off_t sum = 0;
            forAll(recvSizes, i)
            {
                sum += recvSizes[i];
            }
            Pout<< " (overall " << sum << ")";
        }
        Pout<< " to " << fName
            << " using comm " << comm << endl;
    }

    return true;
}


void* Foam::OFstreamCollator::writeAll(void *threadarg)
{
    OFstreamCollator& handler = *static_cast<OFstreamCollator*>(threadarg);

    // Consume stack
    while (true)
    {
        writeData* ptr = nullptr;

        lockMutex(handler.mutex_);
        if (handler.objects_.size())
        {
            ptr = handler.objects_.pop();
        }
        unlockMutex(handler.mutex_);

        if (!ptr)
        {
            break;
        }
        else
        {
            bool ok = writeFile
            (
                ptr->comm_,
                ptr->typeName_,
                ptr->pathName_,
                ptr->data_,
                ptr->sizes_,
                ptr->haveSlaveData_,
                ptr->slaveData_,

                ptr->format_,
                ptr->version_,
                ptr->compression_,
                ptr->append_
            );
            if (!ok)
            {
                FatalIOErrorInFunction(ptr->pathName_)
                    << "Failed writing " << ptr->pathName_
                    << exit(FatalIOError);
            }

            delete ptr;
        }
        //sleep(1);
    }

    if (debug)
    {
        Pout<< "OFstreamCollator : Exiting write thread " << endl;
    }

    lockMutex(handler.mutex_);
    handler.threadRunning_ = false;
    unlockMutex(handler.mutex_);

    return nullptr;
}


void Foam::OFstreamCollator::waitForBufferSpace(const off_t wantedSize) const
{
    while (true)
    {
        // Count files to be written
        off_t totalSize = 0;

        lockMutex(mutex_);
        forAllConstIter(FIFOStack<writeData*>, objects_, iter)
        {
            totalSize += iter()->size();
        }
        unlockMutex(mutex_);

        if (totalSize == 0 || (totalSize+wantedSize) <= maxBufferSize_)
        {
            break;
        }

        if (debug)
        {
            lockMutex(mutex_);
            Pout<< "OFstreamCollator : Waiting for buffer space."
                << " Currently in use:" << totalSize
                << " limit:" << maxBufferSize_
                << " files:" << objects_.size()
                << endl;
            unlockMutex(mutex_);
        }

        sleep(5);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::OFstreamCollator::OFstreamCollator(const off_t maxBufferSize)
:
    maxBufferSize_(maxBufferSize),
    mutex_
    (
        maxBufferSize_ > 0
      ? allocateMutex()
      : -1
    ),
    thread_
    (
        maxBufferSize_ > 0
      ? allocateThread()
      : -1
    ),
    threadRunning_(false),
    comm_
    (
        UPstream::allocateCommunicator
        (
            UPstream::worldComm,
            identity(UPstream::nProcs(UPstream::worldComm))
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::OFstreamCollator::~OFstreamCollator()
{
    if (threadRunning_)
    {
        if (debug)
        {
            Pout<< "~OFstreamCollator : Waiting for write thread" << endl;
        }

        joinThread(thread_);
    }
    if (thread_ != -1)
    {
        freeThread(thread_);
    }
    if (mutex_ != -1)
    {
        freeMutex(mutex_);
    }
    if (comm_ != -1)
    {
        UPstream::freeCommunicator(comm_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::OFstreamCollator::write
(
    const word& typeName,
    const fileName& fName,
    const string& data,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool append
)
{
    // Determine (on master) sizes to receive. Note: do NOT use thread
    // communicator
    labelList recvSizes;
    decomposedBlockData::gather(Pstream::worldComm, data.size(), recvSizes);
    off_t totalSize = 0;
    label maxLocalSize = 0;
    {
        for (label proci = 0; proci < recvSizes.size(); proci++)
        {
            totalSize += recvSizes[proci];
            maxLocalSize = max(maxLocalSize, recvSizes[proci]);
        }
        Pstream::scatter(totalSize, Pstream::msgType(), Pstream::worldComm);
        Pstream::scatter(maxLocalSize, Pstream::msgType(), Pstream::worldComm);
    }

    if (maxBufferSize_ == 0 || maxLocalSize > maxBufferSize_)
    {
        if (debug)
        {
            Pout<< "OFstreamCollator : non-thread gather and write of " << fName
                << " using worldComm" << endl;
        }
        // Direct collating and writing (so master blocks until all written!)
        const List<char> dummySlaveData;
        return writeFile
        (
            UPstream::worldComm,
            typeName,
            fName,
            data,
            recvSizes,
            false,              // no slave data provided yet
            dummySlaveData,
            fmt,
            ver,
            cmp,
            append
        );
    }
    else if (totalSize <= maxBufferSize_)
    {
        // Total size can be stored locally so receive all data now and only
        // do the writing in the thread

        if (debug)
        {
            Pout<< "OFstreamCollator : non-thread gather; thread write of "
                << fName << endl;
        }

        if (Pstream::master())
        {
            waitForBufferSpace(totalSize);
        }

        // Allocate local buffer for all collated data
        autoPtr<writeData> fileAndDataPtr
        (
            new writeData
            (
                comm_,      // Note: comm not actually used anymore
                typeName,
                fName,
                data,
                recvSizes,
                true,       // have slave data (collected below)
                fmt,
                ver,
                cmp,
                append
            )
        );
        writeData& fileAndData = fileAndDataPtr();

        // Gather the slave data and insert into fileAndData
        UList<char> slice(const_cast<char*>(data.data()), label(data.size()));
        List<int> slaveOffsets;
        decomposedBlockData::gatherSlaveData
        (
            Pstream::worldComm,         // Note: using simulation thread
            slice,
            recvSizes,

            1,                          // startProc,
            Pstream::nProcs()-1,        // n procs

            slaveOffsets,
            fileAndData.slaveData_
        );

        // Append to thread buffer
        lockMutex(mutex_);
        objects_.push(fileAndDataPtr.ptr());
        unlockMutex(mutex_);

        // Start thread if not running
        lockMutex(mutex_);
        if (!threadRunning_)
        {
            createThread(thread_, writeAll, this);
            if (debug)
            {
                Pout<< "OFstreamCollator : Started write thread "
                    << thread_ << endl;
            }
            threadRunning_ = true;
        }
        unlockMutex(mutex_);

        return true;
    }
    else
    {
        if (debug)
        {
            Pout<< "OFstreamCollator : thread gather and write of " << fName
                << " in thread " << thread_
                << " using communicator " << comm_ << endl;
        }

        if (!UPstream::haveThreads())
        {
            FatalErrorInFunction
                << "mpi does not seem to have thread support."
                << "Please increase the buffer size 'maxThreadFileBufferSize'"
                << " to at least " << totalSize
                << " to be able to do the collating before threading."
                << exit(FatalError);
        }

        if (Pstream::master())
        {
            waitForBufferSpace(data.size());
        }

        lockMutex(mutex_);
        // Push all file info on buffer. Note that no slave data provided
        // so it will trigger communication inside the thread
        objects_.push
        (
            new writeData
            (
                comm_,
                typeName,
                fName,
                data,
                recvSizes,
                false,          // Have no slave data; collect in thread
                fmt,
                ver,
                cmp,
                append
            )
        );
        unlockMutex(mutex_);

        lockMutex(mutex_);
        if (!threadRunning_)
        {
            createThread(thread_, writeAll, this);
            if (debug)
            {
                Pout<< "OFstreamCollator : Started write thread " << endl;
            }
            threadRunning_ = true;
        }
        unlockMutex(mutex_);

        return true;
    }
}


// ************************************************************************* //
