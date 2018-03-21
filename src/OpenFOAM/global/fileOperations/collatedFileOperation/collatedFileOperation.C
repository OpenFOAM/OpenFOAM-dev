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

#include "collatedFileOperation.H"
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "Time.H"
#include "threadedCollatedOFstream.H"
#include "decomposedBlockData.H"
#include "registerSwitch.H"
#include "masterOFstream.H"
#include "OFstream.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(collatedFileOperation, 0);
    addToRunTimeSelectionTable
    (
        fileOperation,
        collatedFileOperation,
        word
    );

    float collatedFileOperation::maxThreadFileBufferSize
    (
        debug::floatOptimisationSwitch("maxThreadFileBufferSize", 1e9)
    );
    registerOptSwitch
    (
        "maxThreadFileBufferSize",
        float,
        collatedFileOperation::maxThreadFileBufferSize
    );

    // Mark as needing threaded mpi
    addNamedToRunTimeSelectionTable
    (
        fileOperationInitialise,
        collatedFileOperationInitialise,
        word,
        collated
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::fileOperations::collatedFileOperation::ioRanks()
{
    labelList ioRanks;

    string ioRanksString(getEnv("FOAM_IORANKS"));
    if (!ioRanksString.empty())
    {
        IStringStream is(ioRanksString);
        is >> ioRanks;
    }

    return ioRanks;
}


bool Foam::fileOperations::collatedFileOperation::isMasterRank
(
    const label proci
)
const
{
    if (Pstream::parRun())
    {
        return Pstream::master(comm_);
    }
    else
    {
        // Use any IO ranks
        if (ioRanks_.size())
        {
            // Find myself in IO rank
            return findIndex(ioRanks_, proci) != -1;
        }
        else
        {
            // Assume all in single communicator
            return proci == 0;
        }
    }
}


bool Foam::fileOperations::collatedFileOperation::appendObject
(
    const regIOobject& io,
    const fileName& pathName,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    // Append to processors/ file

    label proci = detectProcessorPath(io.objectPath());

    if (debug)
    {
        Pout<< "collatedFileOperation::writeObject :"
            << " For local object : " << io.name()
            << " appending processor " << proci
            << " data to " << pathName << endl;
    }

    if (proci == -1)
    {
        FatalErrorInFunction
            << "Not a valid processor path " << pathName
            << exit(FatalError);
    }

    const bool isMaster = isMasterRank(proci);

    // Determine the local rank if the pathName is a per-rank one
    label localProci = proci;
    {
        fileName path, procDir, local;
        label groupStart, groupSize, nProcs;
        splitProcessorPath
        (
            pathName,
            path,
            procDir,
            local,
            groupStart,
            groupSize,
            nProcs
        );
        if (groupSize > 0 && groupStart != -1)
        {
            localProci = proci-groupStart;
        }
    }


    // Create string from all data to write
    string buf;
    {
        OStringStream os(fmt, ver);
        if (isMaster)
        {
            if (!io.writeHeader(os))
            {
                return false;
            }
        }

        // Write the data to the Ostream
        if (!io.writeData(os))
        {
            return false;
        }

        if (isMaster)
        {
            IOobject::writeEndDivider(os);
        }

        buf = os.str();
    }


    // Note: cannot do append + compression. This is a limitation
    // of ogzstream (or rather most compressed formats)

    OFstream os
    (
        pathName,
        IOstream::BINARY,
        ver,
        IOstream::UNCOMPRESSED, // no compression
        !isMaster
    );

    if (!os.good())
    {
        FatalIOErrorInFunction(os)
            << "Cannot open for appending"
            << exit(FatalIOError);
    }

    if (isMaster)
    {
        IOobject::writeBanner(os)
            << "FoamFile\n{\n"
            << "    version     " << os.version() << ";\n"
            << "    format      " << os.format() << ";\n"
            << "    class       " << decomposedBlockData::typeName
            << ";\n"
            << "    location    " << pathName << ";\n"
            << "    object      " << pathName.name() << ";\n"
            << "}" << nl;
        IOobject::writeDivider(os) << nl;
    }

    // Write data
    UList<char> slice
    (
        const_cast<char*>(buf.data()),
        label(buf.size())
    );
    os << nl << "// Processor" << localProci << nl << slice << nl;

    return os.good();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::collatedFileOperation::collatedFileOperation
(
    const bool verbose
)
:
    masterUncollatedFileOperation
    (
        (
            ioRanks().size()
          ? UPstream::allocateCommunicator
            (
                UPstream::worldComm,
                subRanks(Pstream::nProcs())
            )
          : UPstream::worldComm
        ),
        false
    ),
    myComm_(comm_),
    writer_(maxThreadFileBufferSize, comm_),
    nProcs_(Pstream::nProcs()),
    ioRanks_(ioRanks())
{
    if (verbose)
    {
        Info<< "I/O    : " << typeName
            << " (maxThreadFileBufferSize " << maxThreadFileBufferSize
            << ')' << endl;

        if (maxThreadFileBufferSize == 0)
        {
            Info<< "         Threading not activated "
                   "since maxThreadFileBufferSize = 0." << nl
                << "         Writing may run slowly for large file sizes."
                << endl;
        }
        else
        {
            Info<< "         Threading activated "
                   "since maxThreadFileBufferSize > 0." << nl
                << "         Requires large enough buffer to collect all data"
                    " or thread support " << nl
                << "         enabled in MPI. If thread support cannot be "
                   "enabled, deactivate" << nl
                << "         threading by setting maxThreadFileBufferSize "
                    "to 0 in" << nl
                << "         $FOAM_ETC/controlDict"
                << endl;
        }

        if (ioRanks_.size())
        {
            // Print a bit of information
            stringList ioRanks(Pstream::nProcs());
            if (Pstream::master(comm_))
            {
                ioRanks[Pstream::myProcNo()] = hostName()+"."+name(pid());
            }
            Pstream::gatherList(ioRanks);

            Info<< "         IO nodes:" << endl;
            forAll(ioRanks, proci)
            {
                if (!ioRanks[proci].empty())
                {
                    Info<< "             " << ioRanks[proci] << endl;
                }
            }
        }


        if
        (
            regIOobject::fileModificationChecking
         == regIOobject::inotifyMaster
        )
        {
            WarningInFunction
                << "Resetting fileModificationChecking to inotify" << endl;
        }

        if
        (
            regIOobject::fileModificationChecking
         == regIOobject::timeStampMaster
        )
        {
            WarningInFunction
                << "Resetting fileModificationChecking to timeStamp" << endl;
        }
    }
}


Foam::fileOperations::collatedFileOperation::collatedFileOperation
(
    const label comm,
    const labelList& ioRanks,
    const word& typeName,
    const bool verbose
)
:
    masterUncollatedFileOperation(comm, false),
    myComm_(-1),
    writer_(maxThreadFileBufferSize, comm),
    nProcs_(Pstream::nProcs()),
    ioRanks_(ioRanks)
{
    if (verbose)
    {
        Info<< "I/O    : " << typeName
            << " (maxThreadFileBufferSize " << maxThreadFileBufferSize
            << ')' << endl;

        if (maxThreadFileBufferSize == 0)
        {
            Info<< "         Threading not activated "
                   "since maxThreadFileBufferSize = 0." << nl
                << "         Writing may run slowly for large file sizes."
                << endl;
        }
        else
        {
            Info<< "         Threading activated "
                   "since maxThreadFileBufferSize > 0." << nl
                << "         Requires large enough buffer to collect all data"
                    " or thread support " << nl
                << "         enabled in MPI. If thread support cannot be "
                   "enabled, deactivate" << nl
                << "         threading by setting maxThreadFileBufferSize "
                    "to 0 in" << nl
                << "         $FOAM_ETC/controlDict"
                << endl;
        }

        if
        (
            regIOobject::fileModificationChecking
         == regIOobject::inotifyMaster
        )
        {
            WarningInFunction
                << "Resetting fileModificationChecking to inotify" << endl;
        }

        if
        (
            regIOobject::fileModificationChecking
         == regIOobject::timeStampMaster
        )
        {
            WarningInFunction
                << "Resetting fileModificationChecking to timeStamp" << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::collatedFileOperation::~collatedFileOperation()
{
    if (myComm_ != -1)
    {
        UPstream::freeCommunicator(myComm_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::fileOperations::collatedFileOperation::objectPath
(
    const IOobject& io,
    const word& typeName
) const
{
    // Replacement for objectPath
    if (io.time().processorCase())
    {
        return masterUncollatedFileOperation::localObjectPath
        (
            io,
            fileOperation::PROCOBJECT,
            "dummy",        // not used for processorsobject
            io.instance()
        );
    }
    else
    {
        return masterUncollatedFileOperation::localObjectPath
        (
            io,
            fileOperation::OBJECT,
            word::null,
            io.instance()
        );
    }
}


bool Foam::fileOperations::collatedFileOperation::writeObject
(
    const regIOobject& io,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool valid
) const
{
    const Time& tm = io.time();
    const fileName& inst = io.instance();

    if (inst.isAbsolute() || !tm.processorCase())
    {
        mkDir(io.path());
        fileName pathName(io.objectPath());

        if (debug)
        {
            Pout<< "collatedFileOperation::writeObject :"
                << " For object : " << io.name()
                << " falling back to master-only output to " << io.path()
                << endl;
        }

        masterOFstream os
        (
            pathName,
            fmt,
            ver,
            cmp,
            false,
            valid
        );

        // If any of these fail, return (leave error handling to Ostream class)
        if (!os.good())
        {
            return false;
        }
        if (!io.writeHeader(os))
        {
            return false;
        }
        // Write the data to the Ostream
        if (!io.writeData(os))
        {
            return false;
        }
        IOobject::writeEndDivider(os);

        return true;
    }
    else
    {
        // Construct the equivalent processors/ directory
        fileName path(processorsPath(io, inst, processorsDir(io)));

        mkDir(path);
        fileName pathName(path/io.name());

        if (io.global())
        {
            if (debug)
            {
                Pout<< "collatedFileOperation::writeObject :"
                    << " For global object : " << io.name()
                    << " falling back to master-only output to " << pathName
                    << endl;
            }

            masterOFstream os
            (
                pathName,
                fmt,
                ver,
                cmp,
                false,
                valid
            );

            // If any of these fail, return (leave error handling to Ostream
            // class)
            if (!os.good())
            {
                return false;
            }
            if (!io.writeHeader(os))
            {
                return false;
            }
            // Write the data to the Ostream
            if (!io.writeData(os))
            {
                return false;
            }
            IOobject::writeEndDivider(os);

            return true;
        }
        else if (!Pstream::parRun())
        {
            // Special path for e.g. decomposePar. Append to
            // processorsDDD/ file
            if (debug)
            {
                Pout<< "collatedFileOperation::writeObject :"
                    << " For object : " << io.name()
                    << " appending to " << pathName << endl;
            }

            return appendObject(io, pathName, fmt, ver, cmp);
        }
        else
        {
            if (debug)
            {
                Pout<< "collatedFileOperation::writeObject :"
                    << " For object : " << io.name()
                    << " starting collating output to " << pathName << endl;
            }

            threadedCollatedOFstream os(writer_, pathName, fmt, ver, cmp);

            // If any of these fail, return (leave error handling to Ostream
            // class)
            if (!os.good())
            {
                return false;
            }
            if (Pstream::master(comm_) && !io.writeHeader(os))
            {
                return false;
            }
            // Write the data to the Ostream
            if (!io.writeData(os))
            {
                return false;
            }
            if (Pstream::master(comm_))
            {
                IOobject::writeEndDivider(os);
            }

            return true;
        }
    }
}


Foam::word Foam::fileOperations::collatedFileOperation::processorsDir
(
    const fileName& fName
) const
{
    if (Pstream::parRun())
    {
        const List<int>& procs(UPstream::procID(comm_));

        word procDir(processorsBaseDir+Foam::name(Pstream::nProcs()));

        if (procs.size() != Pstream::nProcs())
        {
            procDir +=
              + "_"
              + Foam::name(procs[0])
              + "-"
              + Foam::name(procs.last());
        }
        return procDir;
    }
    else
    {
        word procDir(processorsBaseDir+Foam::name(nProcs_));

        if (ioRanks_.size())
        {
            // Detect current processor number
            label proci = detectProcessorPath(fName);

            if (proci != -1)
            {
                // Find lowest io rank
                label minProc = 0;
                label maxProc = nProcs_-1;
                forAll(ioRanks_, i)
                {
                    if (ioRanks_[i] >= nProcs_)
                    {
                        break;
                    }
                    else if (ioRanks_[i] <= proci)
                    {
                        minProc = ioRanks_[i];
                    }
                    else
                    {
                        maxProc = ioRanks_[i]-1;
                        break;
                    }
                }
                procDir +=
                  + "_"
                  + Foam::name(minProc)
                  + "-"
                  + Foam::name(maxProc);
            }
        }

        return procDir;
    }
}


Foam::word Foam::fileOperations::collatedFileOperation::processorsDir
(
    const IOobject& io
) const
{
    return processorsDir(io.objectPath());
}


void Foam::fileOperations::collatedFileOperation::setNProcs(const label nProcs)
{
    nProcs_ = nProcs;

    if (debug)
    {
        Pout<< "collatedFileOperation::setNProcs :"
            << " Setting number of processors to " << nProcs_ << endl;
    }
}


// ************************************************************************* //
