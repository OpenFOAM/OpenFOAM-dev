/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
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

#include "masterUncollatedFileOperation.H"
#include "Time.H"
#include "masterOFstream.H"
#include "decomposedBlockData.H"
#include "dummyISstream.H"
#include "SubList.H"
#include "PackedBoolList.H"
#include "gzstream.h"
#include "addToRunTimeSelectionTable.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(masterUncollatedFileOperation, 0);
    addToRunTimeSelectionTable
    (
        fileOperation,
        masterUncollatedFileOperation,
        word
    );

    float masterUncollatedFileOperation::maxMasterFileBufferSize
    (
        Foam::debug::floatOptimisationSwitch("maxMasterFileBufferSize", 1e9)
    );

    // Mark as not needing threaded mpi
    addNamedToRunTimeSelectionTable
    (
        fileOperationInitialise,
        masterUncollatedFileOperationInitialise,
        word,
        masterUncollated
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::fileOperations::masterUncollatedFileOperation::subRanks
(
    const label n
)
{
    string ioRanksString(getEnv("FOAM_IORANKS"));
    if (ioRanksString.empty())
    {
        return identity(n);
    }
    else
    {
        DynamicList<label> subRanks(n);

        IStringStream is(ioRanksString);
        labelList ioRanks(is);

        if (findIndex(ioRanks, 0) == -1)
        {
            FatalErrorInFunction
                << "Rank 0 (master) should be in the IO ranks. Currently "
                << ioRanks << exit(FatalError);
        }

        // The lowest numbered rank is the IO rank
        PackedBoolList isIOrank(n);
        isIOrank.set(ioRanks);

        for (label proci = Pstream::myProcNo(); proci >= 0; --proci)
        {
            if (isIOrank[proci])
            {
                // Found my master. Collect all processors with same master
                subRanks.append(proci);
                for
                (
                    label rank = proci+1;
                    rank < n && !isIOrank[rank];
                    ++rank
                )
                {
                    subRanks.append(rank);
                }
                break;
            }
        }
        return move(subRanks);
    }
}


Foam::word
Foam::fileOperations::masterUncollatedFileOperation::findInstancePath
(
    const instantList& timeDirs,
    const instant& t
)
{
    // Note:
    // - times will include constant (with value 0) as first element.
    //   For backwards compatibility make sure to find 0 in preference
    //   to constant.
    // - list is sorted so could use binary search

    forAllReverse(timeDirs, i)
    {
        if (t.equal(timeDirs[i].value()))
        {
            return timeDirs[i].name();
        }
    }

    return word::null;
}


Foam::fileName
Foam::fileOperations::masterUncollatedFileOperation::filePathInfo
(
    const bool checkGlobal,
    const bool isFile,
    const IOobject& io,
    pathType& searchType,
    word& procsDir,
    word& newInstancePath
) const
{
    procsDir = word::null;
    newInstancePath = word::null;

    if (io.instance().isAbsolute())
    {
        fileName objPath = io.instance()/io.name();

        if (isFileOrDir(isFile, objPath))
        {
            searchType = fileOperation::ABSOLUTE;
            return objPath;
        }
        else
        {
            searchType = fileOperation::NOTFOUND;
            return fileName::null;
        }
    }
    else
    {
        // 1. Check the writing fileName
        fileName writePath(objectPath(io, io.headerClassName()));

        if (isFileOrDir(isFile, writePath))
        {
            searchType = fileOperation::WRITEOBJECT;
            return writePath;
        }

        // 2. Check processors/
        if (io.time().processorCase())
        {
            tmpNrc<dirIndexList> pDirs(lookupProcessorsPath(io.objectPath()));
            forAll(pDirs(), i)
            {
                const fileName& pDir = pDirs()[i].first();
                fileName objPath =
                    processorsPath(io, io.instance(), pDir)
                   /io.name();
                if (objPath != writePath && isFileOrDir(isFile, objPath))
                {
                    searchType = pDirs()[i].second().first();
                    procsDir = pDir;
                    return objPath;
                }
            }
        }
        {
            // 3. Check local
            fileName localPath = io.objectPath();

            if
            (
                localPath != writePath
            &&  isFileOrDir(isFile, localPath)
            )
            {
                searchType = fileOperation::OBJECT;
                return localPath;
            }
        }



        // Any global checks
        if
        (
            checkGlobal
         && io.time().processorCase()
         && (
                io.instance() == io.time().system()
             || io.instance() == io.time().constant()
            )
        )
        {
            fileName parentPath =
                io.rootPath()/io.time().globalCaseName()
               /io.instance()/io.db().dbDir()/io.local()/io.name();

            if (isFileOrDir(isFile, parentPath))
            {
                searchType = fileOperation::PARENTOBJECT;
                return parentPath;
            }
        }

        // Check for approximately same time. E.g. if time = 1e-2 and
        // directory is 0.01 (due to different time formats)
        HashPtrTable<instantList>::const_iterator pathFnd
        (
            times_.find
            (
                io.time().path()
            )
        );
        if (pathFnd != times_.end())
        {
            newInstancePath = findInstancePath
            (
                *pathFnd(),
                instant(io.instance())
            );

            if (newInstancePath.size() && newInstancePath != io.instance())
            {
                // 1. Try processors equivalent
                tmpNrc<dirIndexList> pDirs
                (
                    lookupProcessorsPath(io.objectPath())
                );
                forAll(pDirs(), i)
                {
                    const fileName& pDir = pDirs()[i].first();

                    fileName fName
                    (
                        processorsPath(io, newInstancePath, pDir)
                       /io.name()
                    );
                    if (isFileOrDir(isFile, fName))
                    {
                        switch (pDirs()[i].second().first())
                        {
                            case fileOperation::PROCUNCOLLATED:
                            {
                                searchType =
                                    fileOperation::PROCUNCOLLATEDINSTANCE;
                            }
                            break;
                            case fileOperation::PROCBASEOBJECT:
                            {
                                searchType = fileOperation::PROCBASEINSTANCE;
                            }
                            break;
                            case fileOperation::PROCOBJECT:
                            {
                                searchType = fileOperation::PROCINSTANCE;
                            }
                            break;
                            default:
                            break;
                        }
                        procsDir = pDir;
                        return fName;
                    }
                }


                // 2. Check local
                fileName fName
                (
                   io.rootPath()/io.caseName()
                  /newInstancePath/io.db().dbDir()/io.local()/io.name()
                );
                if (isFileOrDir(isFile, fName))
                {
                    searchType = fileOperation::FINDINSTANCE;
                    return fName;
                }
            }
        }

        searchType = fileOperation::NOTFOUND;
        return fileName::null;
    }
}


Foam::fileName
Foam::fileOperations::masterUncollatedFileOperation::localObjectPath
(
    const IOobject& io,
    const pathType& searchType,
    const word& procDir,
    const word& instancePath
) const
{
    // Replacement for IOobject::objectPath()

    switch (searchType)
    {
        case fileOperation::ABSOLUTE:
        {
            return io.instance()/io.name();
        }
        break;

        case fileOperation::OBJECT:
        {
            return io.path()/io.name();
        }
        break;

        case fileOperation::WRITEOBJECT:
        {
            return objectPath(io, io.headerClassName());
        }
        break;

        case fileOperation::PROCUNCOLLATED:
        {
            // Uncollated type, e.g. processor1
            const word procName
            (
                "processor"
               +Foam::name(Pstream::myProcNo(Pstream::worldComm))
            );
            return
                processorsPath
                (
                    io,
                    io.instance(),
                    (
                        Pstream::parRun()
                      ? procName
                      : procDir
                    )
                )
               /io.name();
        }
        break;

        case fileOperation::PROCBASEOBJECT:
        {
            // Collated, e.g. processors4
            return
                processorsPath(io, io.instance(), procDir)
               /io.name();
        }
        break;

        case fileOperation::PROCOBJECT:
        {
            // Processors directory locally provided by the fileHandler itself
            return
                processorsPath(io, io.instance(), processorsDir(io))
               /io.name();
        }
        break;

        case fileOperation::PARENTOBJECT:
        {
            return
                io.rootPath()/io.time().globalCaseName()
               /io.instance()/io.db().dbDir()/io.local()/io.name();
        }
        break;

        case fileOperation::FINDINSTANCE:
        {
            return
                io.rootPath()/io.caseName()
               /instancePath/io.db().dbDir()/io.local()/io.name();
        }
        break;

        case fileOperation::PROCUNCOLLATEDINSTANCE:
        {
            // Uncollated type, e.g. processor1
            const word procName
            (
                "processor"
               +Foam::name(Pstream::myProcNo(Pstream::worldComm))
            );
            return
                processorsPath
                (
                    io,
                    instancePath,
                    (
                        Pstream::parRun()
                      ? procName
                      : procDir
                    )
                )
               /io.name();
        }
        break;

        case fileOperation::PROCBASEINSTANCE:
        {
            // Collated, e.g. processors4
            return
                processorsPath(io, instancePath, procDir)
               /io.name();
        }
        break;

        case fileOperation::PROCINSTANCE:
        {
            // Processors directory locally provided by the fileHandler itself
            return
                processorsPath(io, instancePath, processorsDir(io))
               /io.name();
        }
        break;

        case fileOperation::NOTFOUND:
        {
            return fileName::null;
        }
        break;

        default:
        {
            NotImplemented;
            return fileName::null;
        }
    }
}


bool Foam::fileOperations::masterUncollatedFileOperation::uniformFile
(
    const fileNameList& filePaths
)
{
    const fileName& object0 = filePaths[0];

    for (label i = 1; i < filePaths.size(); i++)
    {
        if (filePaths[i] != object0)
        {
            return false;
        }
    }
    return true;
}


void Foam::fileOperations::masterUncollatedFileOperation::readAndSend
(
    const fileName& filePath,
    const labelUList& procs,
    PstreamBuffers& pBufs
)
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << ": Opening " << filePath << endl;
    }

    IFstream is(filePath, IOstream::streamFormat::BINARY);

    if (!is.good())
    {
        FatalIOErrorInFunction(filePath) << "Cannot open file " << filePath
            << exit(FatalIOError);
    }

    if (isA<igzstream>(is.stdStream()))
    {
        if (debug)
        {
            Pout<< FUNCTION_NAME << ": Reading compressed" << endl;
        }

        std::ostringstream stringStr;
        stringStr << is.stdStream().rdbuf();
        string buf(stringStr.str());

        forAll(procs, i)
        {
            UOPstream os(procs[i], pBufs);
            os.write(&buf[0], buf.size());
        }
    }
    else
    {
        off_t count(Foam::fileSize(filePath));

        if (debug)
        {
            Pout<< FUNCTION_NAME << " : Reading " << count << " bytes " << endl;
        }

        List<char> buf(static_cast<label>(count));
        is.stdStream().read(buf.begin(), count);

        forAll(procs, i)
        {
            UOPstream os(procs[i], pBufs);
            os.write(buf.begin(), count);
        }
    }
}


Foam::autoPtr<Foam::ISstream>
Foam::fileOperations::masterUncollatedFileOperation::read
(
    IOobject& io,
    const label comm,
    const bool uniform,             // on comms master only
    const fileNameList& filePaths,  // on comms master only
    const boolList& read            // on comms master only
)
{
    autoPtr<ISstream> isPtr;

    // const bool uniform = uniformFile(filePaths);

    PstreamBuffers pBufs
    (
        Pstream::commsTypes::nonBlocking,
        Pstream::msgType(),
        comm
    );

    if (Pstream::master(comm))
    {
        if (uniform)
        {
            if (read[0])
            {
                if (filePaths[0].empty())
                {
                    FatalIOErrorInFunction(filePaths[0])
                        << "cannot find file " << io.objectPath()
                        << exit(FatalIOError);
                }

                DynamicList<label> validProcs(Pstream::nProcs(comm));
                for
                (
                    label proci = 0;
                    proci < Pstream::nProcs(comm);
                    proci++
                )
                {
                    if (read[proci])
                    {
                        validProcs.append(proci);
                    }
                }

                // Read on master and send to all processors (including
                // master for simplicity)
                if (debug)
                {
                    Pout<< "masterUncollatedFileOperation::readStream :"
                        << " For uniform file " << filePaths[0]
                        << " sending to " << validProcs
                        << " in comm:" << comm << endl;
                }
                readAndSend(filePaths[0], validProcs, pBufs);
            }
        }
        else
        {
            if (read[0])
            {
                if (filePaths[0].empty())
                {
                    FatalIOErrorInFunction(filePaths[0])
                        << "cannot find file " << io.objectPath()
                        << exit(FatalIOError);
                }

                autoPtr<IFstream> ifsPtr(new IFstream(filePaths[0]));

                // Read header
                if (!io.readHeader(ifsPtr()))
                {
                    FatalIOErrorInFunction(ifsPtr())
                        << "problem while reading header for object "
                        << io.name() << exit(FatalIOError);
                }

                // Open master (steal from ifsPtr)
                isPtr.reset(ifsPtr.ptr());
            }

            // Read slave files
            for
            (
                label proci = 1;
                proci < Pstream::nProcs(comm);
                proci++
            )
            {
                if (debug)
                {
                    Pout<< "masterUncollatedFileOperation::readStream :"
                        << " For processor " << proci
                        << " opening " << filePaths[proci] << endl;
                }

                const fileName& fPath = filePaths[proci];

                if (read[proci] && !fPath.empty())
                {
                    // Note: handle compression ourselves since size cannot
                    // be determined without actually uncompressing
                    readAndSend(fPath, labelList(1, proci), pBufs);
                }
            }
        }
    }

    labelList recvSizes;
    pBufs.finishedSends(recvSizes);

    // isPtr will be valid on master and will be the unbuffered
    // IFstream. Else the information is in the PstreamBuffers (and
    // the special case of a uniform file)

    if (read[Pstream::myProcNo(comm)])
    {
        // This processor needs to return something

        if (!isPtr.valid())
        {
            UIPstream is(Pstream::masterNo(), pBufs);
            string buf(recvSizes[Pstream::masterNo()], '\0');
            if (recvSizes[Pstream::masterNo()] > 0)
            {
                is.read(&buf[0], recvSizes[Pstream::masterNo()]);
            }

            if (debug)
            {
                Pout<< "masterUncollatedFileOperation::readStream :"
                    << " Done reading " << buf.size() << " bytes" << endl;
            }
            const fileName& fName = filePaths[Pstream::myProcNo(comm)];
            isPtr.reset(new IStringStream(fName, buf, IOstream::BINARY));

            if (!io.readHeader(isPtr()))
            {
                FatalIOErrorInFunction(isPtr())
                    << "problem while reading header for object "
                    << io.name() << exit(FatalIOError);
            }
        }
    }
    else
    {
        isPtr.reset(new dummyISstream());
    }

    return isPtr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::masterUncollatedFileOperation::
masterUncollatedFileOperation
(
    const bool verbose
)
:
    fileOperation
    (
        UPstream::allocateCommunicator
        (
            UPstream::worldComm,
            subRanks(Pstream::nProcs())
        )
    ),
    myComm_(comm_)
{
    if (verbose)
    {
        InfoHeader
            << "I/O    : " << typeName
            << " (maxMasterFileBufferSize " << maxMasterFileBufferSize << ')'
            << endl;
    }

    if (regIOobject::fileModificationChecking == regIOobject::timeStampMaster)
    {
        if (verbose)
        {
            WarningInFunction
                << "Resetting fileModificationChecking to timeStamp" << endl;
        }
        regIOobject::fileModificationChecking = regIOobject::timeStamp;
    }
    else if
    (
        regIOobject::fileModificationChecking
     == regIOobject::inotifyMaster
    )
    {
        if (verbose)
        {
            WarningInFunction
                << "Resetting fileModificationChecking to inotify"
                << endl;
        }
        regIOobject::fileModificationChecking = regIOobject::inotify;
    }
}


Foam::fileOperations::masterUncollatedFileOperation::
masterUncollatedFileOperation
(
    const label comm,
    const bool verbose
)
:
    fileOperation(comm),
    myComm_(-1)
{
    if (verbose)
    {
        InfoHeader
            << "I/O    : " << typeName
            << " (maxMasterFileBufferSize " << maxMasterFileBufferSize << ')'
            << endl;
    }

    if (regIOobject::fileModificationChecking == regIOobject::timeStampMaster)
    {
        if (verbose)
        {
            WarningInFunction
                << "Resetting fileModificationChecking to timeStamp" << endl;
        }
        regIOobject::fileModificationChecking = regIOobject::timeStamp;
    }
    else if
    (
        regIOobject::fileModificationChecking
     == regIOobject::inotifyMaster
    )
    {
        if (verbose)
        {
            WarningInFunction
                << "Resetting fileModificationChecking to inotify"
                << endl;
        }
        regIOobject::fileModificationChecking = regIOobject::inotify;
    }
}


Foam::fileOperations::masterUncollatedFileOperationInitialise::
masterUncollatedFileOperationInitialise(int& argc, char**& argv)
:
    unthreadedInitialise(argc, argv)
{
    // Filter out any of my arguments
    const string s("-ioRanks");

    int index = -1;
    for (int i=1; i<argc-1; i++)
    {
        if (argv[i] == s)
        {
            index = i;
            setEnv("FOAM_IORANKS", argv[i+1], true);
            break;
        }
    }

    if (index != -1)
    {
        for (int i=index+2; i<argc; i++)
        {
            argv[i-2] = argv[i];
        }
        argc -= 2;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::masterUncollatedFileOperation::
~masterUncollatedFileOperation()
{
    if (myComm_ != -1)
    {
        UPstream::freeCommunicator(myComm_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileOperations::masterUncollatedFileOperation::mkDir
(
    const fileName& dir,
    mode_t mode
) const
{
    return masterOp<mode_t, mkDirOp>
    (
        dir,
        mkDirOp(mode),
        Pstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::masterUncollatedFileOperation::chMod
(
    const fileName& fName,
    mode_t mode
) const
{
    return masterOp<mode_t, chModOp>
    (
        fName,
        chModOp(mode),
        Pstream::msgType(),
        comm_
    );
}


mode_t Foam::fileOperations::masterUncollatedFileOperation::mode
(
    const fileName& fName,
    const bool checkVariants,
    const bool followLink
) const
{
    return masterOp<mode_t, modeOp>
    (
        fName,
        modeOp(checkVariants, followLink),
        Pstream::msgType(),
        comm_
    );
}


Foam::fileType Foam::fileOperations::masterUncollatedFileOperation::type
(
    const fileName& fName,
    const bool checkVariants,
    const bool followLink
) const
{
    return fileType
    (
        masterOp<label, typeOp>
        (
            fName,
            typeOp(checkVariants, followLink),
            Pstream::msgType(),
            comm_
        )
    );
}


bool Foam::fileOperations::masterUncollatedFileOperation::exists
(
    const fileName& fName,
    const bool checkVariants,
    const bool followLink
) const
{
    return masterOp<bool, existsOp>
    (
        fName,
        existsOp(checkVariants, followLink),
        Pstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::masterUncollatedFileOperation::isDir
(
    const fileName& fName,
    const bool followLink
) const
{
    return masterOp<bool, isDirOp>
    (
        fName,
        isDirOp(followLink),
        Pstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::masterUncollatedFileOperation::isFile
(
    const fileName& fName,
    const bool checkVariants,
    const bool followLink
) const
{
    return masterOp<bool, isFileOp>
    (
        fName,
        isFileOp(checkVariants, followLink),
        Pstream::msgType(),
        comm_
    );
}


off_t Foam::fileOperations::masterUncollatedFileOperation::fileSize
(
    const fileName& fName,
    const bool checkVariants,
    const bool followLink
) const
{
    return masterOp<off_t, fileSizeOp>
    (
        fName,
        fileSizeOp(checkVariants, followLink),
        Pstream::msgType(),
        comm_
    );
}


time_t Foam::fileOperations::masterUncollatedFileOperation::lastModified
(
    const fileName& fName,
    const bool checkVariants,
    const bool followLink
) const
{
    return masterOp<time_t, lastModifiedOp>
    (
        fName,
        lastModifiedOp(checkVariants, followLink),
        Pstream::msgType(),
        comm_
    );
}


double Foam::fileOperations::masterUncollatedFileOperation::highResLastModified
(
    const fileName& fName,
    const bool checkVariants,
    const bool followLink
) const
{
    return masterOp<double, lastModifiedHROp>
    (
        fName,
        lastModifiedHROp(checkVariants, followLink),
        Pstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::masterUncollatedFileOperation::mvBak
(
    const fileName& fName,
    const std::string& ext
) const
{
    return masterOp<bool, mvBakOp>
    (
        fName,
        mvBakOp(ext),
        Pstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::masterUncollatedFileOperation::rm
(
    const fileName& fName
) const
{
    return masterOp<bool, rmOp>
    (
        fName,
        rmOp(),
        Pstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::masterUncollatedFileOperation::rmDir
(
    const fileName& dir
) const
{
    return masterOp<bool, rmDirOp>
    (
        dir,
        rmDirOp(),
        Pstream::msgType(),
        comm_
    );
}


Foam::fileNameList Foam::fileOperations::masterUncollatedFileOperation::readDir
(
    const fileName& dir,
    const fileType type,
    const bool filtergz,
    const bool followLink
) const
{
    return masterOp<fileNameList, readDirOp>
    (
        dir,
        readDirOp(type, filtergz, followLink),
        Pstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::masterUncollatedFileOperation::cp
(
    const fileName& src,
    const fileName& dst,
    const bool followLink
) const
{
    return masterOp<bool, cpOp>
    (
        src,
        dst,
        cpOp(followLink),
        Pstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::masterUncollatedFileOperation::ln
(
    const fileName& src,
    const fileName& dst
) const
{
    return masterOp<bool, lnOp>
    (
        src,
        dst,
        lnOp(),
        Pstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::masterUncollatedFileOperation::mv
(
    const fileName& src,
    const fileName& dst,
    const bool followLink
) const
{
    return masterOp<bool, mvOp>
    (
        src,
        dst,
        mvOp(followLink),
        Pstream::msgType(),
        comm_
    );
}


Foam::fileName Foam::fileOperations::masterUncollatedFileOperation::filePath
(
    const bool checkGlobal,
    const IOobject& io,
    const word& typeName
) const
{
    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::filePath :"
            << " objectPath:" << io.objectPath()
            << " checkGlobal:" << checkGlobal << endl;
    }

    // Now that we have an IOobject path use it to detect & cache
    // processor directory naming
    (void)lookupProcessorsPath(io.objectPath());

    // Trigger caching of times
    (void)findTimes(io.time().path(), io.time().constant());


    // Determine master filePath and scatter

    fileName objPath;
    pathType searchType = NOTFOUND;
    word procsDir;
    word newInstancePath;

    if (Pstream::master(comm_))
    {
        // All masters search locally. Note that global objects might
        // fail (except on master). This gets handled later on (in PARENTOBJECT)
        objPath = filePathInfo
        (
            checkGlobal,
            true,
            io,
            searchType,
            procsDir,
            newInstancePath
        );

        if (debug)
        {
            Pout<< "masterUncollatedFileOperation::filePath :"
                << " master objPath:" << objPath
                << " searchType:" << fileOperation::pathTypeNames_[searchType]
                << " procsDir:" << procsDir << " instance:" << newInstancePath
                << endl;
        }

    }

    // Scatter the information about where the master found the object
    // Note: use the worldComm to make sure all processors decide
    //       the same type. Only procsDir is allowed to differ; searchType
    //       and instance have to be same
    {
        label masterType(searchType);
        Pstream::scatter(masterType);
        searchType = pathType(masterType);
    }
    Pstream::scatter(newInstancePath);

    if
    (
        checkGlobal
     || searchType == fileOperation::PARENTOBJECT
     || searchType == fileOperation::PROCBASEOBJECT
     || searchType == fileOperation::PROCBASEINSTANCE
     || io.local() == "uniform"
    )
    {
            // Distribute master path. This makes sure it is seen as uniform
            // and only gets read from the master.
            Pstream::scatter(objPath);
            Pstream::scatter(procsDir);
    }
    else
    {
        Pstream::scatter(procsDir, Pstream::msgType(), comm_);

        // Use the master type to determine if additional information is
        // needed to construct the local equivalent
        switch (searchType)
        {
            case fileOperation::PARENTOBJECT:
            case fileOperation::PROCBASEOBJECT:
            case fileOperation::PROCBASEINSTANCE:
            {
                // Already handled above
            }
            break;

            case fileOperation::ABSOLUTE:
            case fileOperation::WRITEOBJECT:
            case fileOperation::PROCUNCOLLATED:
            case fileOperation::PROCOBJECT:
            case fileOperation::FINDINSTANCE:
            case fileOperation::PROCUNCOLLATEDINSTANCE:
            case fileOperation::PROCINSTANCE:
            {
                // Construct equivalent local path
                objPath = localObjectPath
                (
                    io,
                    searchType,
                    procsDir,
                    newInstancePath
                );
            }
            break;

            case fileOperation::OBJECT:
            case fileOperation::NOTFOUND:
            {
                // Retest all processors separately since some processors might
                // have the file and some not (e.g. lagrangian data)

                objPath = masterOp<fileName, fileOrNullOp>
                (
                    io.objectPath(),
                    fileOrNullOp(true),
                    Pstream::msgType(),
                    comm_
                );
            }
            break;
        }
    }

    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::filePath :"
            << " Returning from file searching:" << endl
            << "    objectPath:" << io.objectPath() << endl
            << "    filePath  :" << objPath << endl << endl;
    }
    return objPath;
}


Foam::fileName Foam::fileOperations::masterUncollatedFileOperation::dirPath
(
    const bool checkGlobal,
    const IOobject& io
) const
{
    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::dirPath :"
            << " objectPath:" << io.objectPath()
            << " checkGlobal:" << checkGlobal << endl;
    }

    // Now that we have an IOobject path use it to detect & cache
    // processor directory naming
    (void)lookupProcessorsPath(io.objectPath());

    // Determine master dirPath and scatter

    fileName objPath;
    pathType searchType = NOTFOUND;
    word procsDir;
    word newInstancePath;

    if (Pstream::master(comm_))
    {
        objPath = filePathInfo
        (
            checkGlobal,
            false,
            io,
            searchType,
            procsDir,
            newInstancePath
        );
    }

    {
        label masterType(searchType);
        Pstream::scatter(masterType);   //, Pstream::msgType(), comm_);
        searchType = pathType(masterType);
    }
    Pstream::scatter(newInstancePath);  //, Pstream::msgType(), comm_);

    if
    (
        checkGlobal
     || searchType == fileOperation::PARENTOBJECT
     || searchType == fileOperation::PROCBASEOBJECT
     || searchType == fileOperation::PROCBASEINSTANCE
     || io.local() == "uniform"
    )
    {
            // Distribute master path. This makes sure it is seen as uniform
            // and only gets read from the master.
            Pstream::scatter(objPath);
            Pstream::scatter(procsDir);
    }
    else
    {
        Pstream::scatter(procsDir, Pstream::msgType(), comm_);

        // Use the master type to determine if additional information is
        // needed to construct the local equivalent
        switch (searchType)
        {
            case fileOperation::PARENTOBJECT:
            case fileOperation::PROCBASEOBJECT:
            case fileOperation::PROCBASEINSTANCE:
            {
                // Already handled above
            }
            break;

            case fileOperation::ABSOLUTE:
            case fileOperation::WRITEOBJECT:
            case fileOperation::PROCUNCOLLATED:
            case fileOperation::PROCOBJECT:
            case fileOperation::FINDINSTANCE:
            case fileOperation::PROCUNCOLLATEDINSTANCE:
            case fileOperation::PROCINSTANCE:
            {
                // Construct equivalent local path
                objPath = localObjectPath
                (
                    io,
                    searchType,
                    procsDir,
                    newInstancePath
                );
            }
            break;

            case fileOperation::OBJECT:
            case fileOperation::NOTFOUND:
            {
                // Retest all processors separately since some processors might
                // have the file and some not (e.g. lagrangian data)
                objPath = masterOp<fileName, fileOrNullOp>
                (
                    io.objectPath(),
                    fileOrNullOp(false),
                    Pstream::msgType(),
                    comm_
                );
            }
            break;
        }
    }

    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::dirPath :"
            << " Returning from file searching:" << endl
            << "    objectPath:" << io.objectPath() << endl
            << "    filePath  :" << objPath << endl << endl;
    }
    return objPath;
}


bool Foam::fileOperations::masterUncollatedFileOperation::exists
(
    const dirIndexList& pDirs,
    IOobject& io
) const
{
    // Cut-down version of filePathInfo that does not look for
    // different instance or parent directory

    const bool isFile = !io.name().empty();

    // Generate output filename for object
    const fileName writePath(objectPath(io, word::null));

    // 1. Test writing name for either directory or a (valid) file
    if (isFileOrDir(isFile, writePath))
    {
        return true;
    }

    // 2. Check processors/
    if (io.time().processorCase())
    {
        forAll(pDirs, i)
        {
            const fileName& pDir = pDirs[i].first();
            fileName procPath =
                processorsPath(io, io.instance(), pDir)
               /io.name();
            if (procPath != writePath && isFileOrDir(isFile, procPath))
            {
                return true;
            }
        }
    }

    // 3. Check local
    fileName localPath = io.objectPath();

    if (localPath != writePath && isFileOrDir(isFile, localPath))
    {
        return true;
    }

    return false;
}


Foam::IOobject
Foam::fileOperations::masterUncollatedFileOperation::findInstance
(
    const IOobject& startIO,
    const scalar startValue,
    const word& stopInstance
) const
{
    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::findInstance :"
            << " Starting searching for name:" << startIO.name()
            << " local:" << startIO.local()
            << " from instance:" << startIO.instance()
            << endl;
    }


    const Time& time = startIO.time();

    IOobject io(startIO);

    // Note: - if name is empty, just check the directory itself
    //       - check both for isFile and headerOk since the latter does a
    //         filePath so searches for the file.
    //       - check for an object with local file scope (so no looking up in
    //         parent directory in case of parallel)


    tmpNrc<dirIndexList> pDirs(lookupProcessorsPath(io.objectPath()));

    word foundInstance;

    // if (Pstream::master(comm_))
    if (Pstream::master(UPstream::worldComm))
    {
        if (exists(pDirs, io))
        {
            foundInstance = io.instance();
        }
    }

    // Do parallel early exit to avoid calling time.times()
    // Pstream::scatter(foundInstance, Pstream::msgType(), comm_);
    Pstream::scatter(foundInstance, Pstream::msgType(), UPstream::worldComm);
    if (!foundInstance.empty())
    {
        io.instance() = foundInstance;
        if (debug)
        {
            Pout<< "masterUncollatedFileOperation::findInstance :"
                << " for name:" << io.name() << " local:" << io.local()
                << " found starting instance:" << io.instance() << endl;
        }
        return io;
    }


    // Search back through the time directories to find the time
    // closest to and lower than current time

    instantList ts = time.times();
    // if (Pstream::master(comm_))
    if (Pstream::master(UPstream::worldComm))
    {
        label instanceI;

        for (instanceI = ts.size()-1; instanceI >= 0; --instanceI)
        {
            if (ts[instanceI].value() <= startValue)
            {
                break;
            }
        }

        // continue searching from here
        for (; instanceI >= 0; --instanceI)
        {
            // Shortcut: if actual directory is the timeName we've
            // already tested it
            if (ts[instanceI].name() == time.timeName())
            {
                continue;
            }

            io.instance() = ts[instanceI].name();
            if (exists(pDirs, io))
            {
                foundInstance = io.instance();
                if (debug)
                {
                    Pout<< "masterUncollatedFileOperation::findInstance :"
                        << " for name:" << io.name() << " local:" << io.local()
                        << " found at:" << io.instance()
                        << endl;
                }
                break;
            }

            // Check if hit minimum instance
            if (ts[instanceI].name() == stopInstance)
            {
                if
                (
                    startIO.readOpt() == IOobject::MUST_READ
                 || startIO.readOpt() == IOobject::MUST_READ_IF_MODIFIED
                )
                {
                    if (io.name().empty())
                    {
                        FatalErrorInFunction
                            << "Cannot find directory "
                            << io.local() << " in times " << time.timeName()
                            << " down to " << stopInstance
                            << exit(FatalError);
                    }
                    else
                    {
                        FatalErrorInFunction
                            << "Cannot find file \"" << io.name()
                            << "\" in directory " << io.local()
                            << " in times " << time.timeName()
                            << " down to " << stopInstance
                            << exit(FatalError);
                    }
                }
                foundInstance = io.instance();
                if (debug)
                {
                    Pout<< "masterUncollatedFileOperation::findInstance :"
                        << " name:" << io.name() << " local:" << io.local()
                        << " found at stopinstance:" << io.instance() << endl;
                }
                break;
            }
        }


        if (foundInstance.empty())
        {
            // times() usually already includes the constant() so would
            // have been checked above. Re-test if
            // - times() is empty. Sometimes this can happen (e.g. decomposePar
            //   with collated)
            // - times()[0] is not constant
            if (!ts.size() || ts[0].name() != time.constant())
            {
                // Note. This needs to be a hard-coded constant, rather than the
                // constant function of the time, because the latter points to
                // the case constant directory in parallel cases

                io.instance() = time.constant();
                if (exists(pDirs, io))
                {
                    if (debug)
                    {
                        Pout<< "masterUncollatedFileOperation::findInstance :"
                            << " name:" << io.name()
                            << " local:" << io.local()
                            << " found at:" << io.instance() << endl;
                    }
                    foundInstance = io.instance();
                }
            }
        }

        if (foundInstance.empty())
        {
            if
            (
                startIO.readOpt() == IOobject::MUST_READ
             || startIO.readOpt() == IOobject::MUST_READ_IF_MODIFIED
            )
            {
                FatalErrorInFunction
                    << "Cannot find file \"" << io.name() << "\" in directory "
                    << io.local() << " in times " << startIO.instance()
                    << " down to " << time.constant()
                    << exit(FatalError);
            }
            else
            {
                foundInstance = time.constant();
            }
        }
    }

    // Pstream::scatter(foundInstance, Pstream::msgType(), comm_);
    Pstream::scatter(foundInstance, Pstream::msgType(), UPstream::worldComm);
    io.instance() = foundInstance;
    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::findInstance :"
            << " name:" << io.name() << " local:" << io.local()
            << " returning instance:" << io.instance() << endl;
    }
    return io;
}


Foam::fileNameList
Foam::fileOperations::masterUncollatedFileOperation::readObjects
(
    const objectRegistry& db,
    const fileName& instance,
    const fileName& local,
    word& newInstance
) const
{
    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::readObjects :"
            << " db:" << db.objectPath()
            << " local:" << local << " instance:" << instance << endl;
    }

    fileNameList objectNames;
    newInstance = word::null;

    // Note: readObjects uses WORLD to make sure order of objects is the
    //       same everywhere

    if (Pstream::master())  // comm_))
    {
        // Avoid fileOperation::readObjects from triggering parallel ops
        // (through call to filePath which triggers parallel )
        bool oldParRun = UPstream::parRun();
        UPstream::parRun() = false;

        //- Use non-time searching version
        objectNames = fileOperation::readObjects
        (
            db,
            instance,
            local,
            newInstance
        );

        if (newInstance.empty())
        {
            // Find similar time

            // Copy of Time::findInstancePath. We want to avoid the
            // parallel call to findTimes. Alternative is to have
            // version of findInstancePath that takes instantList ...
            const instantList timeDirs
            (
                fileOperation::findTimes
                (
                    db.time().path(),
                    db.time().constant()
                )
            );

            const instant t(instance);
            forAllReverse(timeDirs, i)
            {
                if (t.equal(timeDirs[i].value()))
                {
                    objectNames = fileOperation::readObjects
                    (
                        db,
                        timeDirs[i].name(),     // newly found time
                        local,
                        newInstance
                    );
                    break;
                }
            }
        }

        UPstream::parRun() = oldParRun;
    }

    Pstream::scatter(newInstance);  //, Pstream::msgType(), comm_);
    Pstream::scatter(objectNames);  //, Pstream::msgType(), comm_);

    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::readObjects :"
            << " newInstance:" << newInstance
            << " objectNames:" << objectNames << endl;
    }

    return objectNames;
}


bool Foam::fileOperations::masterUncollatedFileOperation::readHeader
(
    IOobject& io,
    const fileName& fName,
    const word& typeName
) const
{
    bool ok = false;

    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::readHeader :" << endl
            << "    objectPath:" << io.objectPath() << endl
            << "    fName     :" << fName << endl;
    }

    // Get filePaths on world master
    fileNameList filePaths(Pstream::nProcs(Pstream::worldComm));
    filePaths[Pstream::myProcNo(Pstream::worldComm)] = fName;
    Pstream::gatherList(filePaths, Pstream::msgType(), Pstream::worldComm);
    bool uniform = uniformFile(filePaths);
    Pstream::scatter(uniform, Pstream::msgType(), Pstream::worldComm);

    if (uniform)
    {
        if (Pstream::master(Pstream::worldComm))
        {
            if (!fName.empty())
            {
                IFstream is(fName);

                if (is.good())
                {
                    ok = io.readHeader(is);
                    if (io.headerClassName() == decomposedBlockData::typeName)
                    {
                        // Read the header inside the container (master data)
                        ok = decomposedBlockData::readMasterHeader(io, is);
                    }
                }
            }
        }
        Pstream::scatter(ok, Pstream::msgType(), Pstream::worldComm);
        Pstream::scatter
        (
            io.headerClassName(),
            Pstream::msgType(),
            Pstream::worldComm
        );
        Pstream::scatter(io.note(), Pstream::msgType(), Pstream::worldComm);
    }
    else
    {
        if (Pstream::nProcs(comm_) != Pstream::nProcs(Pstream::worldComm))
        {
            // Re-gather file paths on local master
            filePaths.setSize(Pstream::nProcs(comm_));
            filePaths[Pstream::myProcNo(comm_)] = fName;
            Pstream::gatherList(filePaths, Pstream::msgType(), comm_);
        }

        boolList result(Pstream::nProcs(comm_), false);
        wordList headerClassName(Pstream::nProcs(comm_));
        stringList note(Pstream::nProcs(comm_));
        if (Pstream::master(comm_))
        {
            forAll(filePaths, proci)
            {
                if (!filePaths[proci].empty())
                {
                    if (proci > 0 && filePaths[proci] == filePaths[proci-1])
                    {
                        result[proci] = result[proci-1];
                        headerClassName[proci] = headerClassName[proci-1];
                        note[proci] = note[proci-1];
                    }
                    else
                    {
                        IFstream is(filePaths[proci]);

                        if (is.good())
                        {
                            result[proci] = io.readHeader(is);
                            if
                            (
                                io.headerClassName()
                             == decomposedBlockData::typeName
                            )
                            {
                                // Read the header inside the container (master
                                // data)
                                result[proci] = decomposedBlockData::
                                readMasterHeader
                                (
                                    io,
                                    is
                                );
                            }
                            headerClassName[proci] = io.headerClassName();
                            note[proci] = io.note();
                        }
                    }
                }
            }
        }
        ok = scatterList(result, Pstream::msgType(), comm_);
        io.headerClassName() = scatterList
        (
            headerClassName,
            Pstream::msgType(),
            comm_
        );
        io.note() = scatterList(note, Pstream::msgType(), comm_);
    }

    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::readHeader :" << " ok:" << ok
            << " class:" << io.headerClassName() << endl;
    }
    return ok;
}


Foam::autoPtr<Foam::ISstream>
Foam::fileOperations::masterUncollatedFileOperation::readStream
(
    regIOobject& io,
    const fileName& fName,
    const word& typeName,
    const bool read
) const
{
    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::readStream :"
            << " object : " << io.name()
            << " global : " << io.global()
            << " fName : " << fName << " read:" << read << endl;
    }


    autoPtr<ISstream> isPtr;
    bool isCollated = false;
    IOobject headerIO(io);

    // Detect collated format. This could be done on the local communicator
    // but we do it on the master node only for now.
    if (UPstream::master()) // comm_))
    {
        if (!fName.empty())
        {
            // This can happen in lagrangian field reading some processors
            // have no file to read from. This will only happen when using
            // normal writing since then the fName for the valid processors is
            // processorDDD/<instance>/.. . In case of collocated writing
            // the fName is already rewritten to processors/.

            isPtr.reset(new IFstream(fName));

            if (isPtr().good())
            {
                // Read header data (on copy)
                headerIO.readHeader(isPtr());

                if (headerIO.headerClassName() == decomposedBlockData::typeName)
                {
                    isCollated = true;
                }
                else if (!Pstream::parRun())
                {
                    // Short circuit: non-collated format. No parallel bits.
                    // Copy header and return.
                    if (debug)
                    {
                        Pout<< "masterUncollatedFileOperation::readStream :"
                            << " For object : " << io.name()
                            << " doing straight IFstream input from "
                            << fName << endl;
                    }
                    io = headerIO;
                    return isPtr;
                }
            }

            if (!isCollated)
            {
                // Close file. Reopened below.
                isPtr.clear();
            }
        }
    }

    Pstream::scatter(isCollated);   //, Pstream::msgType(), comm_);

    if (isCollated)
    {
        if (debug)
        {
            Pout<< "masterUncollatedFileOperation::readStream :"
                << " For object : " << io.name()
                << " starting collating input from " << fName << endl;
        }


        // Analyse the file path (on (co)master) to see the processors type
        fileName path, procDir, local;
        label groupStart, groupSize, nProcs;
        splitProcessorPath
        (
            fName,
            path,
            procDir,
            local,
            groupStart,
            groupSize,
            nProcs
        );


        List<char> data;
        if (!Pstream::parRun())
        {
            // Analyse the objectpath to find out the processor we're trying
            // to access
            label proci = detectProcessorPath(io.objectPath());

            if (proci == -1)
            {
                FatalIOErrorInFunction(isPtr())
                    << "Could not detect processor number"
                    << " from objectPath:" << io.objectPath()
                    << exit(FatalIOError);
            }

            // Analyse the fileName for any processor subset. Note: this
            // should really be part of filePath() which should return
            // both file and index in file.
            if (groupStart != -1 && groupSize > 0)
            {
                proci = proci-groupStart;
            }

            if (debug)
            {
                Pout<< "masterUncollatedFileOperation::readStream :"
                    << " For object : " << io.name()
                    << " starting input from block " << proci
                    << " of " << isPtr().name() << endl;
            }

            return decomposedBlockData::readBlock(proci, isPtr(), io);
        }
        else
        {
            // Scatter master header info
            string versionString;
            string formatString;
            if (isPtr.valid())
            {
                versionString = isPtr().version().str();
                OStringStream os;
                os << isPtr().format();
                formatString = (os.str());
            }

            Pstream::scatter(versionString); //,  Pstream::msgType(), comm);
            Pstream::scatter(formatString); //,  Pstream::msgType(), comm);

            // Get size of file
            off_t sz = Foam::fileSize(fName);
            bool bigSize = sz > off_t(maxMasterFileBufferSize);
            Pstream::scatter(bigSize);

            // Are we reading from single-master file ('processors256') or
            // from multi-master files ('processors256_0-9')
            label readComm = -1;
            if (groupStart != -1 && groupSize > 0)
            {
                readComm = comm_;
                if (UPstream::master(comm_) && !isPtr.valid() && !fName.empty())
                {
                    // In multi-master mode also open the file on the other
                    // masters
                    isPtr.reset(new IFstream(fName));

                    if (isPtr().good())
                    {
                        // Read header data (on copy)
                        IOobject headerIO(io);
                        headerIO.readHeader(isPtr());
                    }
                }
            }
            else
            {
                // Single master so read on world
                readComm = Pstream::worldComm;
            }

            // Read my data
            return decomposedBlockData::readBlocks
            (
                readComm,
                fName,
                isPtr,
                io,
                (
                    bigSize
                  ? UPstream::commsTypes::scheduled
                  : UPstream::commsTypes::nonBlocking
                )
            );
        }
    }
    else
    {
        if (debug)
        {
            Pout<< "masterUncollatedFileOperation::readStream :"
                << " For object : " << io.name()
                << " starting separated input from " << fName << endl;
        }

        if (io.global())
        {
            // Use worldComm. Note: should not really need to gather filePaths
            // since we enforce sending from master anyway ...
            fileNameList filePaths(Pstream::nProcs());
            filePaths[Pstream::myProcNo()] = fName;
            Pstream::gatherList(filePaths);
            boolList procValid(Pstream::nProcs());
            procValid[Pstream::myProcNo()] = read;
            Pstream::gatherList(procValid);

            return this->read
            (
                io,
                Pstream::worldComm,
                true,
                filePaths,
                procValid
            );
        }
        else
        {
            // Use local communicator
            fileNameList filePaths(Pstream::nProcs(comm_));
            filePaths[Pstream::myProcNo(comm_)] = fName;
            Pstream::gatherList(filePaths, Pstream::msgType(), comm_);
            boolList procValid(Pstream::nProcs(comm_));
            procValid[Pstream::myProcNo(comm_)] = read;
            Pstream::gatherList(procValid, Pstream::msgType(), comm_);

            // Uniform in local comm
            bool uniform = uniformFile(filePaths);

            return this->read
            (
                io,
                comm_,
                uniform,
                filePaths,
                procValid
            );
        }
    }
}


bool Foam::fileOperations::masterUncollatedFileOperation::read
(
    regIOobject& io,
    const bool masterOnly,
    const IOstream::streamFormat defaultFormat,
    const word& typeName
) const
{
    bool ok = true;

    // Initialise format to the defaultFormat
    // but reset to ASCII if defaultFormat and file format are ASCII
    IOstream::streamFormat format = defaultFormat;

    if (io.global())
    {
        if (debug)
        {
            Pout<< "masterUncollatedFileOperation::read :"
                << " Reading global object " << io.name() << endl;
        }

        // Now that we have an IOobject path use it to detect & cache
        // processor directory naming
        (void)lookupProcessorsPath(io.objectPath());

        // Trigger caching of times
        (void)findTimes(io.time().path(), io.time().constant());

        bool ok = false;
        if (Pstream::master())  // comm_))
        {
            // Do master-only reading always.
            bool oldParRun = UPstream::parRun();
            UPstream::parRun() = false;

            // Open file and read header
            Istream& is = io.readStream(typeName);

            // Set format to ASCII if defaultFormat and file format are ASCII
            if (defaultFormat == IOstream::ASCII)
            {
                format = is.format();
            }

            // Read the data from the file
            ok = io.readData(is);

            // Close the file
            io.close();

            UPstream::parRun() = oldParRun;
        }

        Pstream::scatter(ok);
        Pstream::scatter(io.headerClassName());
        Pstream::scatter(io.note());

        if (defaultFormat == IOstream::ASCII)
        {
            std::underlying_type_t<IOstream::streamFormat> formatValue(format);
            Pstream::scatter(formatValue);
            format = IOstream::streamFormat(formatValue);
        }

        // scatter operation for regIOobjects

        // Get my communication order
        const List<Pstream::commsStruct>& comms =
        (
            (Pstream::nProcs() < Pstream::nProcsSimpleSum)
          ? Pstream::linearCommunication()
          : Pstream::treeCommunication()
        );
        const Pstream::commsStruct& myComm = comms[Pstream::myProcNo()];

        // Receive from up
        if (myComm.above() != -1)
        {
            IPstream fromAbove
            (
                Pstream::commsTypes::scheduled,
                myComm.above(),
                0,
                Pstream::msgType(),
                Pstream::worldComm,
                format
            );
            ok = io.readData(fromAbove);
        }

        // Send to my downstairs neighbours
        forAll(myComm.below(), belowI)
        {
            OPstream toBelow
            (
                Pstream::commsTypes::scheduled,
                myComm.below()[belowI],
                0,
                Pstream::msgType(),
                Pstream::worldComm,
                format
            );
            bool okWrite = io.writeData(toBelow);
            ok = ok && okWrite;
        }
    }
    else
    {
        if (debug)
        {
            Pout<< "masterUncollatedFileOperation::read :"
                << " Reading local object " << io.name() << endl;
        }

        ok = io.readData(io.readStream(typeName));
        io.close();
    }

    return ok;
}


bool Foam::fileOperations::masterUncollatedFileOperation::writeObject
(
    const regIOobject& io,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    fileName filePath(io.objectPath());

    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::writeObject :"
            << " io:" << filePath << " write:" << write << endl;
    }

    // Make sure to pick up any new times
    setTime(io.time());

    autoPtr<Ostream> osPtr
    (
        NewOFstream
        (
            filePath,
            fmt,
            ver,
            cmp,
            write
        )
    );
    Ostream& os = osPtr();

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


Foam::instantList Foam::fileOperations::masterUncollatedFileOperation::findTimes
(
    const fileName& directory,
    const word& constantName
) const
{
    HashPtrTable<instantList>::const_iterator iter = times_.find(directory);
    if (iter != times_.end())
    {
        if (debug)
        {
            Pout<< "masterUncollatedFileOperation::findTimes :"
                << " Found " << iter()->size() << " cached times" << endl;
        }
        return *iter();
    }
    else
    {
        instantList times;
        if (Pstream::master())  // comm_))
        {
            // Do master-only reading always.
            bool oldParRun = UPstream::parRun();
            UPstream::parRun() = false;
            times = fileOperation::findTimes(directory, constantName);
            UPstream::parRun() = oldParRun;
        }
        Pstream::scatter(times);    //, Pstream::msgType(), comm_);

        // Note: do we also cache if no times have been found since it might
        //       indicate a directory that is being filled later on ...

        instantList* tPtr = new instantList(move(times));

        times_.insert(directory, tPtr);

        if (debug)
        {
            Pout<< "masterUncollatedFileOperation::findTimes :"
                << " Caching times:" << *tPtr << nl
                << "    for directory:" << directory << endl;
        }
        return *tPtr;
    }
}


void Foam::fileOperations::masterUncollatedFileOperation::setTime
(
    const Time& tm
) const
{
    if (tm.subCycling())
    {
        return;
    }

    HashPtrTable<instantList>::const_iterator iter = times_.find(tm.path());
    if (iter != times_.end())
    {
        instantList& times = *iter();

        const instant timeNow(tm.value(), tm.timeName());

        if (times.size() > 0 && times[0].name() == tm.constant())
        {
            // Exclude constant
            SubList<instant> realTimes(times, times.size()-1, 1);
            if
            (
                findSortedIndex
                (
                    SubList<instant>(times, times.size()-1, 1),
                    timeNow
                )
             == -1
            )
            {
                if (debug)
                {
                    Pout<< "masterUncollatedFileOperation::setTime :"
                        << " Caching time " << tm.timeName()
                        << " for case:" << tm.path() << endl;
                }

                times.append(timeNow);
                SubList<instant> realTimes(times, times.size()-1, 1);
                Foam::stableSort(realTimes);
            }
        }
        else
        {
            if (findSortedIndex(times, timeNow) == -1)
            {
                if (debug)
                {
                    Pout<< "masterUncollatedFileOperation::setTime :"
                        << " Caching time " << tm.timeName()
                        << " for case:" << tm.path() << endl;
                }

                times.append(timeNow);
                Foam::stableSort(times);
            }
        }
    }
    fileOperation::setTime(tm);
}


Foam::autoPtr<Foam::ISstream>
Foam::fileOperations::masterUncollatedFileOperation::NewIFstream
(
    const fileName& filePath,
    IOstream::streamFormat format,
    IOstream::versionNumber version
) const
{
    if (Pstream::parRun())
    {
        // Insert logic of filePath. We assume that if a file is absolute
        // on the master it is absolute also on the slaves etc.

        fileNameList filePaths(Pstream::nProcs(Pstream::worldComm));
        filePaths[Pstream::myProcNo(Pstream::worldComm)] = filePath;
        Pstream::gatherList(filePaths, Pstream::msgType(), Pstream::worldComm);

        PstreamBuffers pBufs
        (
            Pstream::commsTypes::nonBlocking,
            Pstream::msgType(),
            Pstream::worldComm
        );

        if (Pstream::master(Pstream::worldComm))
        {
            const bool uniform = uniformFile(filePaths);

            if (uniform)
            {
                if (debug)
                {
                    Pout<< "masterUncollatedFileOperation::NewIFstream :"
                        << " Opening global file " << filePath << endl;
                }

                labelList procs(Pstream::nProcs(Pstream::worldComm)-1);
                for
                (
                    label proci = 1;
                    proci < Pstream::nProcs(Pstream::worldComm);
                    proci++
                )
                {
                    procs[proci-1] = proci;
                }

                readAndSend(filePath, procs, pBufs);
            }
            else
            {
                for
                (
                    label proci = 1;
                    proci < Pstream::nProcs(Pstream::worldComm);
                    proci++
                )
                {
                    readAndSend(filePaths[proci], labelList(1, proci), pBufs);
                }
            }
        }


        labelList recvSizes;
        pBufs.finishedSends(recvSizes);

        if (Pstream::master(Pstream::worldComm))
        {
            // Read myself
            return autoPtr<ISstream>
            (
                new IFstream(filePaths[Pstream::masterNo()], format, version)
            );
        }
        else
        {
            if (debug)
            {
                Pout<< "masterUncollatedFileOperation::NewIFstream :"
                    << " Reading " << filePath
                    << " from processor " << Pstream::masterNo() << endl;
            }

            UIPstream is(Pstream::masterNo(), pBufs);
            string buf(recvSizes[Pstream::masterNo()], '\0');
            is.read(&buf[0], recvSizes[Pstream::masterNo()]);

            if (debug)
            {
                Pout<< "masterUncollatedFileOperation::NewIFstream :"
                    << " Done reading " << buf.size() << " bytes" << endl;
            }

            // Note: IPstream is not an IStream so use a IStringStream to
            //       convert the buffer. Note that we construct with a string
            //       so it holds a copy of the buffer.
            return autoPtr<ISstream>
            (
                new IStringStream(filePath, buf, IOstream::BINARY)
            );
        }
    }
    else
    {
        // Read myself
        return autoPtr<ISstream>(new IFstream(filePath, format, version));
    }
}


Foam::autoPtr<Foam::Ostream>
Foam::fileOperations::masterUncollatedFileOperation::NewOFstream
(
    const fileName& filePath,
    IOstream::streamFormat format,
    IOstream::versionNumber version,
    IOstream::compressionType compression,
    const bool write
) const
{
    return autoPtr<Ostream>
    (
        new masterOFstream
        (
            filePath,
            format,
            version,
            compression,
            false,      // append
            write
        )
    );
}


void Foam::fileOperations::masterUncollatedFileOperation::flush() const
{
    fileOperation::flush();
    times_.clear();
}


Foam::label Foam::fileOperations::masterUncollatedFileOperation::addWatch
(
    const fileName& fName
) const
{
    label watchFd;
    if (Pstream::master())      // comm_))
    {
        watchFd = monitor().addWatch(fName);
    }
    Pstream::scatter(watchFd);  //, Pstream::msgType(), comm_);
    return watchFd;
}


bool Foam::fileOperations::masterUncollatedFileOperation::removeWatch
(
    const label watchIndex
) const
{
    bool ok;
    if (Pstream::master())  // comm_))
    {
        ok = monitor().removeWatch(watchIndex);
    }
    Pstream::scatter(ok);   //, Pstream::msgType(), comm_);
    return ok;
}


Foam::label Foam::fileOperations::masterUncollatedFileOperation::findWatch
(
    const labelList& watchIndices,
    const fileName& fName
) const
{
    label index = -1;

    if (Pstream::master())  // comm_))
    {
        forAll(watchIndices, i)
        {
            if (monitor().getFile(watchIndices[i]) == fName)
            {
                index = i;
                break;
            }
        }
    }
    Pstream::scatter(index);    //, Pstream::msgType(), comm_);
    return index;
}


void Foam::fileOperations::masterUncollatedFileOperation::addWatches
(
    regIOobject& rio,
    const fileNameList& files
) const
{
    const labelList& watchIndices = rio.watchIndices();

    DynamicList<label> newWatchIndices;
    labelHashSet removedWatches(watchIndices);

    forAll(files, i)
    {
        const fileName& f = files[i];
        label index = findWatch(watchIndices, f);

        if (index == -1)
        {
            newWatchIndices.append(addWatch(f));
        }
        else
        {
            // Existing watch
            newWatchIndices.append(watchIndices[index]);
            removedWatches.erase(index);
        }
    }

    // Remove any unused watches
    forAllConstIter(labelHashSet, removedWatches, iter)
    {
        removeWatch(watchIndices[iter.key()]);
    }

    rio.watchIndices() = newWatchIndices;
}


Foam::fileName Foam::fileOperations::masterUncollatedFileOperation::getFile
(
    const label watchIndex
) const
{
    fileName fName;
    if (Pstream::master())  // comm_))
    {
        fName = monitor().getFile(watchIndex);
    }
    Pstream::scatter(fName);    //, Pstream::msgType(), comm_);
    return fName;
}


void Foam::fileOperations::masterUncollatedFileOperation::updateStates
(
    const bool masterOnly,
    const bool syncPar
) const
{
    if (Pstream::master())  // comm_))
    {
        monitor().updateStates(true, false);
    }
}


Foam::fileMonitor::fileState
Foam::fileOperations::masterUncollatedFileOperation::getState
(
    const label watchFd
) const
{
    unsigned int state = fileMonitor::UNMODIFIED;
    if (Pstream::master())  // comm_))
    {
        state = monitor().getState(watchFd);
    }
    Pstream::scatter(state);    //, Pstream::msgType(), comm_);
    return fileMonitor::fileState(state);
}


void Foam::fileOperations::masterUncollatedFileOperation::setUnmodified
(
    const label watchFd
) const
{
    if (Pstream::master())  // comm_))
    {
        monitor().setUnmodified(watchFd);
    }
}


// ************************************************************************* //
