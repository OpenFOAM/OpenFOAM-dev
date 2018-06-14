/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
#include "masterUncollatedFileOperation.H"

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
    const PtrList<SubList<char>>& slaveData,    // optional slave data
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool append
)
{
    if (debug)
    {
        Pout<< "OFstreamCollator : Writing master " << masterData.size()
            << " bytes to " << fName
            << " using comm " << comm << endl;
        if (slaveData.size())
        {
            Pout<< "OFstreamCollator :  Slave data" << endl;
            forAll(slaveData, proci)
            {
                if (slaveData.set(proci))
                {
                    Pout<< "    " << proci
                        << " size:" << slaveData[proci].size()
                        << endl;
                }
            }
        }
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
        if (!append)
        {
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
        slaveData,
        (
            fileOperations::masterUncollatedFileOperation::
                maxMasterFileBufferSize == 0
          ? UPstream::commsTypes::scheduled
          : UPstream::commsTypes::nonBlocking
        ),
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
            // Use ostringstream to display long int (until writing these is
            // supported)
            std::ostringstream os;
            os << sum;
            Pout<< " (overall " << os.str() << ")";
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

        {
            std::lock_guard<std::mutex> guard(handler.mutex_);
            if (handler.objects_.size())
            {
                ptr = handler.objects_.pop();
            }
        }

        if (!ptr)
        {
            break;
        }
        else
        {
            // Convert storage to pointers
            PtrList<SubList<char>> slaveData;
            if (ptr->slaveData_.size())
            {
                slaveData.setSize(ptr->slaveData_.size());
                forAll(slaveData, proci)
                {
                    if (ptr->slaveData_.set(proci))
                    {
                        slaveData.set
                        (
                            proci,
                            new SubList<char>
                            (
                                ptr->slaveData_[proci],
                                ptr->sizes_[proci]
                            )
                        );
                    }
                }
            }

            bool ok = writeFile
            (
                ptr->comm_,
                ptr->typeName_,
                ptr->pathName_,
                ptr->data_,
                ptr->sizes_,
                slaveData,
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
        // sleep(1);
    }

    if (debug)
    {
        Pout<< "OFstreamCollator : Exiting write thread " << endl;
    }

    {
        std::lock_guard<std::mutex> guard(handler.mutex_);
        handler.threadRunning_ = false;
    }

    return nullptr;
}


void Foam::OFstreamCollator::waitForBufferSpace(const off_t wantedSize) const
{
    while (true)
    {
        // Count files to be written
        off_t totalSize = 0;

        {
            std::lock_guard<std::mutex> guard(mutex_);
            forAllConstIter(FIFOStack<writeData*>, objects_, iter)
            {
                totalSize += iter()->size();
            }
        }

        if
        (
            totalSize == 0
         || (wantedSize >= 0 && (totalSize+wantedSize) <= maxBufferSize_)
        )
        {
            break;
        }

        if (debug)
        {
            std::lock_guard<std::mutex> guard(mutex_);
            Pout<< "OFstreamCollator : Waiting for buffer space."
                << " Currently in use:" << totalSize
                << " limit:" << maxBufferSize_
                << " files:" << objects_.size()
                << endl;
        }

        sleep(5);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::OFstreamCollator::OFstreamCollator(const off_t maxBufferSize)
:
    maxBufferSize_(maxBufferSize),
    threadRunning_(false),
    localComm_(UPstream::worldComm),
    threadComm_
    (
        UPstream::allocateCommunicator
        (
            localComm_,
            identity(UPstream::nProcs(localComm_))
        )
    )
{}


Foam::OFstreamCollator::OFstreamCollator
(
    const off_t maxBufferSize,
    const label comm
)
:
    maxBufferSize_(maxBufferSize),
    threadRunning_(false),
    localComm_(comm),
    threadComm_
    (
        UPstream::allocateCommunicator
        (
            localComm_,
            identity(UPstream::nProcs(localComm_))
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::OFstreamCollator::~OFstreamCollator()
{
    if (thread_.valid())
    {
        if (debug)
        {
            Pout<< "~OFstreamCollator : Waiting for write thread" << endl;
        }
        thread_().join();
        thread_.clear();
    }

    if (threadComm_ != -1)
    {
        UPstream::freeCommunicator(threadComm_);
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
    const bool append,
    const bool useThread
)
{
    // Determine (on master) sizes to receive. Note: do NOT use thread
    // communicator
    labelList recvSizes;
    decomposedBlockData::gather(localComm_, label(data.size()), recvSizes);

    off_t totalSize = 0;
    label maxLocalSize = 0;
    {
        for (label proci = 0; proci < recvSizes.size(); proci++)
        {
            totalSize += recvSizes[proci];
            maxLocalSize = max(maxLocalSize, recvSizes[proci]);
        }
        Pstream::scatter(totalSize, Pstream::msgType(), localComm_);
        Pstream::scatter(maxLocalSize, Pstream::msgType(), localComm_);
    }

    if (!useThread || maxBufferSize_ == 0 || maxLocalSize > maxBufferSize_)
    {
        if (debug)
        {
            Pout<< "OFstreamCollator : non-thread gather and write of " << fName
                << " using local comm " << localComm_ << endl;
        }
        // Direct collating and writing (so master blocks until all written!)
        const PtrList<SubList<char>> dummySlaveData;
        return writeFile
        (
            localComm_,
            typeName,
            fName,
            data,
            recvSizes,
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

        if (Pstream::master(localComm_))
        {
            waitForBufferSpace(totalSize);
        }


        // Receive in chunks of labelMax (2^31-1) since this is the maximum
        // size that a List can be

        autoPtr<writeData> fileAndDataPtr
        (
            new writeData
            (
                threadComm_,        // Note: comm not actually used anymore
                typeName,
                fName,
                (
                    Pstream::master(localComm_)
                  ? data            // Only used on master
                  : string::null
                ),
                recvSizes,
                fmt,
                ver,
                cmp,
                append
            )
        );
        writeData& fileAndData = fileAndDataPtr();

        PtrList<List<char>>& slaveData = fileAndData.slaveData_;

        UList<char> slice(const_cast<char*>(data.data()), label(data.size()));

        slaveData.setSize(recvSizes.size());

        // Gather all data onto master. Is done in local communicator since
        // not in write thread. Note that we do not store in contiguous
        // buffer since that would limit to 2G chars.
        label startOfRequests = Pstream::nRequests();
        if (Pstream::master(localComm_))
        {
            for (label proci = 1; proci < slaveData.size(); proci++)
            {
                slaveData.set(proci, new List<char>(recvSizes[proci]));
                UIPstream::read
                (
                    UPstream::commsTypes::nonBlocking,
                    proci,
                    reinterpret_cast<char*>(slaveData[proci].begin()),
                    slaveData[proci].byteSize(),
                    Pstream::msgType(),
                    localComm_
                );
            }
        }
        else
        {
            if
            (
               !UOPstream::write
                (
                    UPstream::commsTypes::nonBlocking,
                    0,
                    reinterpret_cast<const char*>(slice.begin()),
                    slice.byteSize(),
                    Pstream::msgType(),
                    localComm_
                )
            )
            {
                FatalErrorInFunction
                    << "Cannot send outgoing message. "
                    << "to:" << 0 << " nBytes:"
                    << label(slice.byteSize())
                    << Foam::abort(FatalError);
            }
        }
        Pstream::waitRequests(startOfRequests);

        {
            std::lock_guard<std::mutex> guard(mutex_);

            // Append to thread buffer
            objects_.push(fileAndDataPtr.ptr());

            // Start thread if not running
            if (!threadRunning_)
            {
                if (thread_.valid())
                {
                    if (debug)
                    {
                        Pout<< "OFstreamCollator : Waiting for write thread"
                            << endl;
                    }
                    thread_().join();
                }

                if (debug)
                {
                    Pout<< "OFstreamCollator : Starting write thread"
                        << endl;
                }
                thread_.reset(new std::thread(writeAll, this));
                threadRunning_ = true;
            }
        }

        return true;
    }
    else
    {
        if (debug)
        {
            Pout<< "OFstreamCollator : thread gather and write of " << fName
                << " using communicator " << threadComm_ << endl;
        }

        if (!UPstream::haveThreads())
        {
            FatalErrorInFunction
                << "mpi does not seem to have thread support."
                << " Make sure to set buffer size 'maxThreadFileBufferSize'"
                << " to at least " << totalSize
                << " to be able to do the collating before threading."
                << exit(FatalError);
        }

        if (Pstream::master(localComm_))
        {
            waitForBufferSpace(data.size());
        }

        {
            std::lock_guard<std::mutex> guard(mutex_);

            // Push all file info on buffer. Note that no slave data provided
            // so it will trigger communication inside the thread
            objects_.push
            (
                new writeData
                (
                    threadComm_,
                    typeName,
                    fName,
                    data,
                    recvSizes,
                    fmt,
                    ver,
                    cmp,
                    append
                )
            );

            if (!threadRunning_)
            {
                if (thread_.valid())
                {
                    if (debug)
                    {
                        Pout<< "OFstreamCollator : Waiting for write thread"
                            << endl;
                    }
                    thread_().join();
                }

                if (debug)
                {
                    Pout<< "OFstreamCollator : Starting write thread" << endl;
                }
                thread_.reset(new std::thread(writeAll, this));
                threadRunning_ = true;
            }
        }

        return true;
    }
}


void Foam::OFstreamCollator::waitAll()
{
    // Wait for all buffer space to be available i.e. wait for all jobs
    // to finish
    if (Pstream::master(localComm_))
    {
        if (debug)
        {
            Pout<< "OFstreamCollator : waiting for thread to have consumed all"
                << endl;
        }
        waitForBufferSpace(-1);
    }
}


// ************************************************************************* //
