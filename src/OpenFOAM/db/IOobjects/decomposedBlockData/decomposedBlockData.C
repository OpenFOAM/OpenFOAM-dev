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

#include "decomposedBlockData.H"
#include "OPstream.H"
#include "IPstream.H"
#include "PstreamBuffers.H"
#include "OFstream.H"
#include "IFstream.H"
#include "IStringStream.H"
#include "dictionary.H"
#include "objectRegistry.H"
#include "SubList.H"
#include "labelPair.H"
#include "masterUncollatedFileOperation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(decomposedBlockData, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decomposedBlockData::decomposedBlockData
(
    const label comm,
    const IOobject& io,
    const UPstream::commsTypes commsType
)
:
    regIOobject(io),
    commsType_(commsType),
    comm_(comm)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningInFunction
            << "decomposedBlockData " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but decomposedBlockData does not support automatic rereading."
            << endl;
    }
    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        read();
    }
}


Foam::decomposedBlockData::decomposedBlockData
(
    const label comm,
    const IOobject& io,
    const UList<char>& list,
    const UPstream::commsTypes commsType
)
:
    regIOobject(io),
    commsType_(commsType),
    comm_(comm)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningInFunction
            << "decomposedBlockData " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but decomposedBlockData does not support automatic rereading."
            << endl;
    }

    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        read();
    }
    else
    {
        List<char>::operator=(list);
    }
}


Foam::decomposedBlockData::decomposedBlockData
(
    const label comm,
    const IOobject& io,
    const Xfer<List<char>>& list,
    const UPstream::commsTypes commsType
)
:
    regIOobject(io),
    commsType_(commsType),
    comm_(comm)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningInFunction
            << "decomposedBlockData " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but decomposedBlockData does not support automatic rereading."
            << endl;
    }

    List<char>::transfer(list());

    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        read();
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::decomposedBlockData::~decomposedBlockData()
{}


// * * * * * * * * * * * * * * * Members Functions * * * * * * * * * * * * * //

bool Foam::decomposedBlockData::readMasterHeader(IOobject& io, Istream& is)
{
    if (debug)
    {
        Pout<< "decomposedBlockData::readMasterHeader:"
            << " stream:" << is.name() << endl;
    }

    // Master-only reading of header
    is.fatalCheck("read(Istream&)");

    List<char> data(is);
    is.fatalCheck("read(Istream&) : reading entry");
    string buf(data.begin(), data.size());
    IStringStream str(is.name(), buf);

    return io.readHeader(str);
}


void Foam::decomposedBlockData::writeHeader
(
    Ostream& os,
    const IOstream::versionNumber version,
    const IOstream::streamFormat format,
    const word& type,
    const string& note,
    const fileName& location,
    const word& name
)
{
    IOobject::writeBanner(os)
        << "FoamFile\n{\n"
        << "    version     " << version << ";\n"
        << "    format      " << format << ";\n"
        << "    class       " << type << ";\n";
    if (note.size())
    {
        os  << "    note        " << note << ";\n";
    }

    if (location.size())
    {
        os  << "    location    " << location << ";\n";
    }

    os  << "    object      " << name << ";\n"
        << "}" << nl;

    IOobject::writeDivider(os) << nl;
}


Foam::autoPtr<Foam::ISstream> Foam::decomposedBlockData::readBlock
(
    const label blocki,
    Istream& is,
    IOobject& headerIO
)
{
    if (debug)
    {
        Pout<< "decomposedBlockData::readBlock:"
            << " stream:" << is.name() << " attempt to read block " << blocki
            << endl;
    }

    is.fatalCheck("read(Istream&)");

    List<char> data;
    autoPtr<ISstream> realIsPtr;

    if (blocki == 0)
    {
        is >> data;
        is.fatalCheck("read(Istream&) : reading entry");

        string buf(data.begin(), data.size());
        realIsPtr = new IStringStream(is.name(), buf);

        // Read header
        if (!headerIO.readHeader(realIsPtr()))
        {
            FatalIOErrorInFunction(realIsPtr())
                << "problem while reading header for object "
                << is.name() << exit(FatalIOError);
        }
    }
    else
    {
        // Read master for header
        is >> data;
        is.fatalCheck("read(Istream&) : reading entry");

        IOstream::versionNumber ver(IOstream::currentVersion);
        IOstream::streamFormat fmt;
        {
            string buf(data.begin(), data.size());
            IStringStream headerStream(is.name(), buf);

            // Read header
            if (!headerIO.readHeader(headerStream))
            {
                FatalIOErrorInFunction(headerStream)
                    << "problem while reading header for object "
                    << is.name() << exit(FatalIOError);
            }
            ver = headerStream.version();
            fmt = headerStream.format();
        }

        for (label i = 1; i < blocki+1; i++)
        {
            // Read data, override old data
            is >> data;
            is.fatalCheck("read(Istream&) : reading entry");
        }
        string buf(data.begin(), data.size());
        realIsPtr = new IStringStream(is.name(), buf);

        // Apply master stream settings to realIsPtr
        realIsPtr().format(fmt);
        realIsPtr().version(ver);
    }
    return realIsPtr;
}


bool Foam::decomposedBlockData::readBlocks
(
    const label comm,
    autoPtr<ISstream>& isPtr,
    List<char>& data,
    const UPstream::commsTypes commsType
)
{
    if (debug)
    {
        Pout<< "decomposedBlockData::readBlocks:"
            << " stream:" << (isPtr.valid() ? isPtr().name() : "invalid")
            << " commsType:" << Pstream::commsTypeNames[commsType]
            << " comm:" << comm << endl;
    }

    bool ok = false;

    if (commsType == UPstream::commsTypes::scheduled)
    {
        if (UPstream::master(comm))
        {
            Istream& is = isPtr();
            is.fatalCheck("read(Istream&)");

            // Read master data
            {
                is >> data;
                is.fatalCheck("read(Istream&) : reading entry");
            }

            // Read slave data
            for
            (
                label proci = 1;
                proci < UPstream::nProcs(comm);
                proci++
            )
            {
                List<char> elems(is);
                is.fatalCheck("read(Istream&) : reading entry");

                OPstream os
                (
                    UPstream::commsTypes::scheduled,
                    proci,
                    0,
                    UPstream::msgType(),
                    comm
                );
                os << elems;
            }

            ok = is.good();
        }
        else
        {
            IPstream is
            (
                UPstream::commsTypes::scheduled,
                UPstream::masterNo(),
                0,
                UPstream::msgType(),
                comm
            );
            is >> data;
        }
    }
    else
    {
        PstreamBuffers pBufs
        (
            UPstream::commsTypes::nonBlocking,
            UPstream::msgType(),
            comm
        );

        if (UPstream::master(comm))
        {
            Istream& is = isPtr();
            is.fatalCheck("read(Istream&)");

            // Read master data
            {
                is >> data;
                is.fatalCheck("read(Istream&) : reading entry");
            }

            // Read slave data
            for
            (
                label proci = 1;
                proci < UPstream::nProcs(comm);
                proci++
            )
            {
                List<char> elems(is);
                is.fatalCheck("read(Istream&) : reading entry");

                UOPstream os(proci, pBufs);
                os << elems;
            }
        }

        labelList recvSizes;
        pBufs.finishedSends(recvSizes);

        if (!UPstream::master(comm))
        {
            UIPstream is(UPstream::masterNo(), pBufs);
            is >> data;
        }
    }

    Pstream::scatter(ok, Pstream::msgType(), comm);

    return ok;
}


Foam::autoPtr<Foam::ISstream> Foam::decomposedBlockData::readBlocks
(
    const label comm,
    const fileName& fName,
    autoPtr<ISstream>& isPtr,
    IOobject& headerIO,
    const UPstream::commsTypes commsType
)
{
    if (debug)
    {
        Pout<< "decomposedBlockData::readBlocks:"
            << " stream:" << (isPtr.valid() ? isPtr().name() : "invalid")
            << " commsType:" << Pstream::commsTypeNames[commsType] << endl;
    }

    bool ok = false;

    List<char> data;
    autoPtr<ISstream> realIsPtr;

    if (commsType == UPstream::commsTypes::scheduled)
    {
        if (UPstream::master(comm))
        {
            Istream& is = isPtr();
            is.fatalCheck("read(Istream&)");

            // Read master data
            {
                is >> data;
                is.fatalCheck("read(Istream&) : reading entry");

                string buf(data.begin(), data.size());
                realIsPtr = new IStringStream(fName, buf);

                // Read header
                if (!headerIO.readHeader(realIsPtr()))
                {
                    FatalIOErrorInFunction(realIsPtr())
                        << "problem while reading header for object "
                        << is.name() << exit(FatalIOError);
                }
            }

            // Read slave data
            for
            (
                label proci = 1;
                proci < UPstream::nProcs(comm);
                proci++
            )
            {
                is >> data;
                is.fatalCheck("read(Istream&) : reading entry");

                OPstream os
                (
                    UPstream::commsTypes::scheduled,
                    proci,
                    0,
                    UPstream::msgType(),
                    comm
                );
                os << data;
            }

            ok = is.good();
        }
        else
        {
            IPstream is
            (
                UPstream::commsTypes::scheduled,
                UPstream::masterNo(),
                0,
                UPstream::msgType(),
                comm
            );
            is >> data;

            string buf(data.begin(), data.size());
            realIsPtr = new IStringStream(fName, buf);
        }
    }
    else
    {
        PstreamBuffers pBufs
        (
            UPstream::commsTypes::nonBlocking,
            UPstream::msgType(),
            comm
        );

        if (UPstream::master(comm))
        {
            Istream& is = isPtr();
            is.fatalCheck("read(Istream&)");

            // Read master data
            {
                is >> data;
                is.fatalCheck("read(Istream&) : reading entry");

                string buf(data.begin(), data.size());
                realIsPtr = new IStringStream(fName, buf);

                // Read header
                if (!headerIO.readHeader(realIsPtr()))
                {
                    FatalIOErrorInFunction(realIsPtr())
                        << "problem while reading header for object "
                        << is.name() << exit(FatalIOError);
                }
            }

            // Read slave data
            for
            (
                label proci = 1;
                proci < UPstream::nProcs(comm);
                proci++
            )
            {
                List<char> elems(is);
                is.fatalCheck("read(Istream&) : reading entry");

                UOPstream os(proci, pBufs);
                os << elems;
            }

            ok = is.good();
        }

        labelList recvSizes;
        pBufs.finishedSends(recvSizes);

        if (!UPstream::master(comm))
        {
            UIPstream is(UPstream::masterNo(), pBufs);
            is >> data;

            string buf(data.begin(), data.size());
            realIsPtr = new IStringStream(fName, buf);
        }
    }

    Pstream::scatter(ok, Pstream::msgType(), comm);

    // version
    string versionString(realIsPtr().version().str());
    Pstream::scatter(versionString,  Pstream::msgType(), comm);
    realIsPtr().version(IStringStream(versionString)());

    // stream
    {
        OStringStream os;
        os << realIsPtr().format();
        string formatString(os.str());
        Pstream::scatter(formatString,  Pstream::msgType(), comm);
        realIsPtr().format(formatString);
    }

    word name(headerIO.name());
    Pstream::scatter(name, Pstream::msgType(), comm);
    headerIO.rename(name);
    Pstream::scatter(headerIO.headerClassName(), Pstream::msgType(), comm);
    Pstream::scatter(headerIO.note(), Pstream::msgType(), comm);
    //Pstream::scatter(headerIO.instance(), Pstream::msgType(), comm);
    //Pstream::scatter(headerIO.local(), Pstream::msgType(), comm);

    return realIsPtr;
}


void Foam::decomposedBlockData::gather
(
    const label comm,
    const label data,
    labelList& datas
)
{
    const label nProcs = UPstream::nProcs(comm);
    datas.setSize(nProcs);

    char* data0Ptr = reinterpret_cast<char*>(datas.begin());

    List<int> recvOffsets;
    List<int> recvSizes;
    if (UPstream::master(comm))
    {
        recvOffsets.setSize(nProcs);
        forAll(recvOffsets, proci)
        {
            // Note: truncating long int to int since UPstream::gather limited
            // to ints
            recvOffsets[proci] =
                int(reinterpret_cast<char*>(&datas[proci]) - data0Ptr);
        }
        recvSizes.setSize(nProcs, sizeof(label));
    }

    UPstream::gather
    (
        reinterpret_cast<const char*>(&data),
        sizeof(label),
        data0Ptr,
        recvSizes,
        recvOffsets,
        comm
    );
}


void Foam::decomposedBlockData::gatherSlaveData
(
    const label comm,
    const UList<char>& data,
    const labelUList& recvSizes,

    const label startProc,
    const label nProcs,

    List<int>& sliceOffsets,
    List<char>& recvData
)
{
    // Calculate master data
    List<int> sliceSizes;
    if (UPstream::master(comm))
    {
        const label numProcs = UPstream::nProcs(comm);

        sliceSizes.setSize(numProcs, 0);
        sliceOffsets.setSize(numProcs+1, 0);

        int totalSize = 0;
        label proci = startProc;
        for (label i = 0; i < nProcs; i++)
        {
            sliceSizes[proci] = int(recvSizes[proci]);
            sliceOffsets[proci] = totalSize;
            totalSize += sliceSizes[proci];
            proci++;
        }
        sliceOffsets[proci] = totalSize;
        recvData.setSize(totalSize);
    }

    int nSend = 0;
    if
    (
       !UPstream::master(comm)
     && (UPstream::myProcNo(comm) >= startProc)
     && (UPstream::myProcNo(comm) < startProc+nProcs)
    )
    {
        // Note: UPstream::gather limited to int
        nSend = int(data.byteSize());
    }

    UPstream::gather
    (
        data.begin(),
        nSend,

        recvData.begin(),
        sliceSizes,
        sliceOffsets,
        comm
    );
}


Foam::label Foam::decomposedBlockData::calcNumProcs
(
    const label comm,
    const off_t maxBufferSize,
    const labelUList& recvSizes,
    const label startProci
)
{
    const label nProcs = UPstream::nProcs(comm);

    label nSendProcs = -1;
    if (UPstream::master(comm))
    {
        off_t totalSize = recvSizes[startProci];
        label proci = startProci+1;
        while (proci < nProcs && (totalSize+recvSizes[proci] < maxBufferSize))
        {
            totalSize += recvSizes[proci];
            proci++;
        }

        nSendProcs = proci-startProci;
    }

    // Scatter nSendProcs
    label n;
    UPstream::scatter
    (
        reinterpret_cast<const char*>(&nSendProcs),
        List<int>(nProcs, sizeof(nSendProcs)),
        List<int>(nProcs, 0),
        reinterpret_cast<char*>(&n),
        sizeof(n),
        comm
    );

    return n;
}


bool Foam::decomposedBlockData::writeBlocks
(
    const label comm,
    autoPtr<OSstream>& osPtr,
    List<std::streamoff>& start,
    const UList<char>& data,

    const labelUList& recvSizes,
    const PtrList<SubList<char>>& slaveData,

    const UPstream::commsTypes commsType,
    const bool syncReturnState
)
{
    if (debug)
    {
        Pout<< "decomposedBlockData::writeBlocks:"
            << " stream:" << (osPtr.valid() ? osPtr().name() : "invalid")
            << " data:" << data.size()
            << " (master only) slaveData:" << slaveData.size()
            << " commsType:" << Pstream::commsTypeNames[commsType] << endl;
    }

    const label nProcs = UPstream::nProcs(comm);

    bool ok = true;

    if (slaveData.size())
    {
        // Already have gathered the slave data. communicator only used to
        // check who is the master

        if (UPstream::master(comm))
        {
            OSstream& os = osPtr();

            start.setSize(nProcs);

            // Write master data
            {
                os << nl << "// Processor" << UPstream::masterNo() << nl;
                start[UPstream::masterNo()] = os.stdStream().tellp();
                os << data;
            }

            // Write slaves

            label slaveOffset = 0;

            for (label proci = 1; proci < nProcs; proci++)
            {
                os << nl << nl << "// Processor" << proci << nl;
                start[proci] = os.stdStream().tellp();

                os << slaveData[proci];
                slaveOffset += recvSizes[proci];
            }

            ok = os.good();
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        if (UPstream::master(comm))
        {
            start.setSize(nProcs);

            OSstream& os = osPtr();

            // Write master data
            {
                os << nl << "// Processor" << UPstream::masterNo() << nl;
                start[UPstream::masterNo()] = os.stdStream().tellp();
                os << data;
            }
            // Write slaves
            List<char> elems;
            for (label proci = 1; proci < nProcs; proci++)
            {
                elems.setSize(recvSizes[proci]);
                IPstream::read
                (
                    UPstream::commsTypes::scheduled,
                    proci,
                    elems.begin(),
                    elems.size(),
                    Pstream::msgType(),
                    comm
                );

                os << nl << nl << "// Processor" << proci << nl;
                start[proci] = os.stdStream().tellp();
                os << elems;
            }

            ok = os.good();
        }
        else
        {
            UOPstream::write
            (
                UPstream::commsTypes::scheduled,
                UPstream::masterNo(),
                data.begin(),
                data.byteSize(),
                Pstream::msgType(),
                comm
            );
        }
    }
    else
    {
        // Write master data
        if (UPstream::master(comm))
        {
            start.setSize(nProcs);

            OSstream& os = osPtr();

            os << nl << "// Processor" << UPstream::masterNo() << nl;
            start[UPstream::masterNo()] = os.stdStream().tellp();
            os << data;
        }


        // Find out how many processor can be received into
        // maxMasterFileBufferSize

        // Starting slave processor and number of processors
        label startProc = 1;
        label nSendProcs = nProcs-1;

        while (nSendProcs > 0)
        {
            nSendProcs = calcNumProcs
            (
                comm,
                off_t
                (
                    fileOperations::masterUncollatedFileOperation::
                    maxMasterFileBufferSize
                ),
                recvSizes,
                startProc
            );

            if (startProc == nProcs || nSendProcs == 0)
            {
                break;
            }


            // Gather data from (a slice of) the slaves
            List<int> sliceOffsets;
            List<char> recvData;
            gatherSlaveData
            (
                comm,
                data,
                recvSizes,

                startProc,    // startProc,
                nSendProcs,    // nProcs,

                sliceOffsets,
                recvData
            );

            if (UPstream::master(comm))
            {
                OSstream& os = osPtr();

                // Write slaves
                for
                (
                    label proci = startProc;
                    proci < startProc+nSendProcs;
                    proci++
                )
                {
                    os << nl << nl << "// Processor" << proci << nl;
                    start[proci] = os.stdStream().tellp();

                    os <<
                        SubList<char>
                        (
                            recvData,
                            sliceOffsets[proci+1]-sliceOffsets[proci],
                            sliceOffsets[proci]
                        );
                }
            }

            startProc += nSendProcs;
        }

        if (UPstream::master(comm))
        {
            ok = osPtr().good();
        }
    }

    if (syncReturnState)
    {
        //- Enable to get synchronised error checking. Is the one that keeps
        //  slaves as slow as the master (which does all the writing)
        Pstream::scatter(ok, Pstream::msgType(), comm);
    }

    return ok;
}


bool Foam::decomposedBlockData::read()
{
    autoPtr<ISstream> isPtr;
    fileName objPath(fileHandler().filePath(false, *this, word::null));
    if (UPstream::master(comm_))
    {
        isPtr.reset(new IFstream(objPath));
        IOobject::readHeader(isPtr());
    }

    List<char>& data = *this;
    return readBlocks(comm_, isPtr, data, commsType_);
}


bool Foam::decomposedBlockData::writeData(Ostream& os) const
{
    const List<char>& data = *this;

    string str
    (
        reinterpret_cast<const char*>(data.cbegin()),
        data.byteSize()
    );

    IOobject io(*this);
    if (Pstream::master(comm_))
    {
        IStringStream is(name(), str);
        io.readHeader(is);
    }

    // Scatter header information

    // version
    string versionString(os.version().str());
    Pstream::scatter(versionString, Pstream::msgType(), comm_);

    // stream
    string formatString;
    {
        OStringStream os;
        os << os.format();
        formatString  = os.str();
        Pstream::scatter(formatString, Pstream::msgType(), comm_);
    }

    //word masterName(name());
    //Pstream::scatter(masterName, Pstream::msgType(), comm_);

    Pstream::scatter(io.headerClassName(), Pstream::msgType(), comm_);
    Pstream::scatter(io.note(), Pstream::msgType(), comm_);
    //Pstream::scatter(io.instance(), Pstream::msgType(), comm);
    //Pstream::scatter(io.local(), Pstream::msgType(), comm);

    fileName masterLocation(instance()/db().dbDir()/local());
    Pstream::scatter(masterLocation, Pstream::msgType(), comm_);

    if (!Pstream::master(comm_))
    {
        writeHeader
        (
            os,
            IOstream::versionNumber(IStringStream(versionString)()),
            IOstream::formatEnum(formatString),
            io.headerClassName(),
            io.note(),
            masterLocation,
            name()
        );
    }

    os.writeQuoted(str, false);

    if (!Pstream::master(comm_))
    {
        IOobject::writeEndDivider(os);
    }

    return os.good();
}


bool Foam::decomposedBlockData::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool valid
) const
{
    autoPtr<OSstream> osPtr;
    if (UPstream::master(comm_))
    {
        // Note: always write binary. These are strings so readable
        //       anyway. They have already be tokenised on the sending side.
        osPtr.reset(new OFstream(objectPath(), IOstream::BINARY, ver, cmp));
        IOobject::writeHeader(osPtr());
    }

    labelList recvSizes;
    gather(comm_, label(this->byteSize()), recvSizes);

    List<std::streamoff> start;
    PtrList<SubList<char>> slaveData;  // dummy slave data
    return writeBlocks
    (
        comm_,
        osPtr,
        start,
        *this,
        recvSizes,
        slaveData,
        commsType_
    );
}


Foam::label Foam::decomposedBlockData::numBlocks(const fileName& fName)
{
    label nBlocks = 0;

    IFstream is(fName);
    is.fatalCheck("decomposedBlockData::numBlocks(const fileName&)");

    if (!is.good())
    {
        return nBlocks;
    }

    // Skip header
    token firstToken(is);

    if
    (
        is.good()
     && firstToken.isWord()
     && firstToken.wordToken() == "FoamFile"
    )
    {
        dictionary headerDict(is);
        is.version(headerDict.lookup("version"));
        is.format(headerDict.lookup("format"));
    }

    List<char> data;
    while (is.good())
    {
        token sizeToken(is);
        if (!sizeToken.isLabel())
        {
            return nBlocks;
        }
        is.putBack(sizeToken);

        is >> data;
        nBlocks++;
    }

    return nBlocks;
}


// ************************************************************************* //
