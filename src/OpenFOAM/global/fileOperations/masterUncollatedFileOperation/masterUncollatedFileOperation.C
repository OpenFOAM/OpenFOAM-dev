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

#include "masterUncollatedFileOperation.H"
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "Time.H"
#include "instant.H"
#include "IFstream.H"
#include "masterOFstream.H"
#include "decomposedBlockData.H"
#include "registerSwitch.H"
#include "dummyISstream.H"
#include "SubList.H"

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
    registerOptSwitch
    (
        "maxMasterFileBufferSize",
        float,
        masterUncollatedFileOperation::maxMasterFileBufferSize
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

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


Foam::fileName Foam::fileOperations::masterUncollatedFileOperation::filePathInfo
(
    const bool checkGlobal,
    const bool isFile,
    const IOobject& io,
    pathType& searchType,
    word& newInstancePath
) const
{
    newInstancePath = word::null;

    if (io.instance().isAbsolute())
    {
        fileName objectPath = io.instance()/io.name();

        if (isFileOrDir(isFile, objectPath))
        {
            searchType = fileOperation::ABSOLUTE;
            return objectPath;
        }
        else
        {
            searchType = fileOperation::NOTFOUND;
            return fileName::null;
        }
    }
    else
    {
        // 1. Check processors/
        if (io.time().processorCase())
        {
            fileName objectPath = processorsPath(io, io.instance())/io.name();
            if (isFileOrDir(isFile, objectPath))
            {
                searchType = fileOperation::PROCESSORSOBJECT;
                return objectPath;
            }
        }
        {
            // 2. Check local
            fileName localObjectPath = io.objectPath();

            if (isFileOrDir(isFile, localObjectPath))
            {
                searchType = fileOperation::OBJECT;
                return localObjectPath;
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
            fileName parentObjectPath =
                io.rootPath()/io.time().globalCaseName()
               /io.instance()/io.db().dbDir()/io.local()/io.name();

            if (isFileOrDir(isFile, parentObjectPath))
            {
                searchType = fileOperation::PARENTOBJECT;
                return parentObjectPath;
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

           if (newInstancePath.size())
           {
               // 1. Try processors equivalent

               fileName fName =
                   processorsPath(io, newInstancePath)
                  /io.name();
               if (isFileOrDir(isFile, fName))
               {
                   searchType = fileOperation::PROCESSORSFINDINSTANCE;
                   return fName;
               }

               fName =
                   io.rootPath()/io.caseName()
                  /newInstancePath/io.db().dbDir()/io.local()/io.name();

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
Foam::fileOperations::masterUncollatedFileOperation::processorsCasePath
(
    const IOobject& io
)
{
    return
        io.rootPath()
       /io.time().globalCaseName()
       /processorsDir;
}


Foam::fileName
Foam::fileOperations::masterUncollatedFileOperation::processorsPath
(
    const IOobject& io,
    const word& instance
)
{
    return
        processorsCasePath(io)
       /instance
       /io.db().dbDir()
       /io.local();
}


Foam::fileName
Foam::fileOperations::masterUncollatedFileOperation::processorsPath
(
    const fileName& dir
)
{
    // Check if directory is processorXXX
    word caseName(dir.name());

    std::string::size_type pos = caseName.find("processor");
    if (pos == 0)
    {
        return dir.path()/processorsDir;
    }
    else
    {
        return fileName::null;
    }
}


Foam::label
Foam::fileOperations::masterUncollatedFileOperation::splitProcessorPath
(
    const fileName& objectPath,
    fileName& path,
    fileName& local
)
{
    // Search for processor at start of line or /processor
    std::string::size_type pos = objectPath.find("processor");
    if (pos == string::npos)
    {
        return -1;
    }

    if (pos == 0)
    {
        path = "";
        local = objectPath.substr(pos+9);
    }
    else if (objectPath[pos-1] != '/')
    {
        return -1;
    }
    else
    {
        path = objectPath.substr(0, pos-1);
        local = objectPath.substr(pos+9);
    }

    pos = local.find('/');
    if (pos == string::npos)
    {
        // processorXXX without local
        label proci;
        if (Foam::read(local.c_str(), proci))
        {
            local.clear();
            return proci;
        }
        return -1;
    }
    string procName(local.substr(0, pos));
    label proci;
    if (Foam::read(procName.c_str(), proci))
    {
        local = local.substr(pos+1);
        return proci;
    }
    return -1;
}


Foam::fileName Foam::fileOperations::masterUncollatedFileOperation::objectPath
(
    const IOobject& io,
    const pathType& searchType,
    const word& instancePath
)
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

        case fileOperation::PROCESSORSOBJECT:
        {
            return processorsPath(io, io.instance())/io.name();
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

        case fileOperation::PROCESSORSFINDINSTANCE:
        {
            return processorsPath(io, instancePath)/io.name();
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
    const IOstream::compressionType cmp,
    const labelUList& procs,
    PstreamBuffers& pBufs
)
{
    if (cmp == IOstream::compressionType::COMPRESSED)
    {
        if (debug)
        {
            Pout<< "masterUncollatedFileOperation::readAndSend:"
                << " opening compressed " << filePath << endl;
        }

        IFstream is(filePath, IOstream::streamFormat::BINARY);

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
        IFstream is(filePath, IOstream::streamFormat::BINARY);

        if (debug)
        {
            Pout<< "masterUncollatedFileOperation::readStream:"
                << " From " << filePath <<  " reading " << label(count)
                << " bytes" << endl;
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::masterUncollatedFileOperation::
masterUncollatedFileOperation
(
    const bool verbose
)
{
    if (verbose)
    {
        Info<< "I/O    : " << typeName
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


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::masterUncollatedFileOperation::
~masterUncollatedFileOperation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileOperations::masterUncollatedFileOperation::mkDir
(
    const fileName& dir,
    mode_t mode
) const
{
    return masterOp<mode_t, mkDirOp>(dir, mkDirOp(mode));
}


bool Foam::fileOperations::masterUncollatedFileOperation::chMod
(
    const fileName& fName,
    mode_t mode
) const
{
    return masterOp<mode_t, chModOp>(fName, chModOp(mode));
}


mode_t Foam::fileOperations::masterUncollatedFileOperation::mode
(
    const fileName& fName,
    const bool followLink
) const
{
    return masterOp<mode_t, modeOp>(fName, modeOp(followLink));
}


Foam::fileName::Type Foam::fileOperations::masterUncollatedFileOperation::type
(
    const fileName& fName,
    const bool followLink
) const
{
    return fileName::Type(masterOp<label, typeOp>(fName, typeOp(followLink)));
}


bool Foam::fileOperations::masterUncollatedFileOperation::exists
(
    const fileName& fName,
    const bool checkGzip,
    const bool followLink
) const
{
    return masterOp<bool, existsOp>(fName, existsOp(checkGzip, followLink));
}


bool Foam::fileOperations::masterUncollatedFileOperation::isDir
(
    const fileName& fName,
    const bool followLink
) const
{
    return masterOp<bool, isDirOp>(fName, isDirOp(followLink));
}


bool Foam::fileOperations::masterUncollatedFileOperation::isFile
(
    const fileName& fName,
    const bool checkGzip,
    const bool followLink
) const
{
    return masterOp<bool, isFileOp>(fName, isFileOp(checkGzip, followLink));
}


off_t Foam::fileOperations::masterUncollatedFileOperation::fileSize
(
    const fileName& fName,
    const bool followLink
) const
{
    return masterOp<off_t, fileSizeOp>(fName, fileSizeOp(followLink));
}


time_t Foam::fileOperations::masterUncollatedFileOperation::lastModified
(
    const fileName& fName,
    const bool followLink
) const
{
    return masterOp<time_t, lastModifiedOp>
    (
        fName,
        lastModifiedOp(followLink)
    );
}


double Foam::fileOperations::masterUncollatedFileOperation::highResLastModified
(
    const fileName& fName,
    const bool followLink
) const
{
    return masterOp<double, lastModifiedHROp>
    (
        fName,
        lastModifiedHROp(followLink)
    );
}


bool Foam::fileOperations::masterUncollatedFileOperation::mvBak
(
    const fileName& fName,
    const std::string& ext
) const
{
    return masterOp<bool, mvBakOp>(fName, mvBakOp(ext));
}


bool Foam::fileOperations::masterUncollatedFileOperation::rm
(
    const fileName& fName
) const
{
    return masterOp<bool, rmOp>(fName, rmOp());
}


bool Foam::fileOperations::masterUncollatedFileOperation::rmDir
(
    const fileName& dir
) const
{
    return masterOp<bool, rmDirOp>(dir, rmDirOp());
}


Foam::fileNameList Foam::fileOperations::masterUncollatedFileOperation::readDir
(
    const fileName& dir,
    const fileName::Type type,
    const bool filtergz,
    const bool followLink
) const
{
    return masterOp<fileNameList, readDirOp>
    (
        dir,
        readDirOp(type, filtergz, followLink)
    );
}


bool Foam::fileOperations::masterUncollatedFileOperation::cp
(
    const fileName& src,
    const fileName& dst,
    const bool followLink
) const
{
    return masterOp<bool, cpOp>(src, dst, cpOp(followLink));
}


bool Foam::fileOperations::masterUncollatedFileOperation::ln
(
    const fileName& src,
    const fileName& dst
) const
{
    return masterOp<bool, lnOp>(src, dst, lnOp());
}


bool Foam::fileOperations::masterUncollatedFileOperation::mv
(
    const fileName& src,
    const fileName& dst,
    const bool followLink
) const
{
    return masterOp<bool, mvOp>(src, dst, mvOp(followLink));
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

    // Trigger caching of times
    (void)findTimes(io.time().path(), io.time().constant());

    // Determine master filePath and scatter

    fileName objPath;
    pathType searchType = NOTFOUND;
    word newInstancePath;

    if (Pstream::master())
    {
        objPath = filePathInfo
        (
            checkGlobal,
            true,
            io,
            searchType,
            newInstancePath
        );
    }

    {
        label masterType(searchType);
        Pstream::scatter(masterType);
        searchType = pathType(masterType);
    }

    Pstream::scatter(newInstancePath);


    // Use the master type to determine if additional information is
    // needed to construct the local equivalent
    switch (searchType)
    {
        case fileOperation::ABSOLUTE:
        case fileOperation::PROCESSORSOBJECT:
        case fileOperation::PARENTOBJECT:
        case fileOperation::FINDINSTANCE:
        case fileOperation::PROCESSORSFINDINSTANCE:
        {
            // Construct equivalent local path
            objPath = objectPath(io, searchType, newInstancePath);
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
                fileOrNullOp(true)
            );
        }
        break;
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

    // Determine master dirPath and scatter

    fileName objPath;
    pathType searchType = NOTFOUND;
    word newInstancePath;

    if (Pstream::master())
    {
        objPath = filePathInfo
        (
            checkGlobal,
            false,
            io,
            searchType,
            newInstancePath
        );
    }

    {
        label masterType(searchType);
        Pstream::scatter(masterType);
        searchType = pathType(masterType);
    }
    Pstream::scatter(newInstancePath);


    // Use the master type to determine if additional information is
    // needed to construct the local equivalent
    switch (searchType)
    {
        case fileOperation::ABSOLUTE:
        case fileOperation::PROCESSORSOBJECT:
        case fileOperation::PARENTOBJECT:
        case fileOperation::FINDINSTANCE:
        case fileOperation::PROCESSORSFINDINSTANCE:
        {
            // Construct equivalent local path
            objPath = objectPath(io, searchType, newInstancePath);
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
                fileOrNullOp(false)
            );
        }
        break;
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
            << " instance:" << instance << endl;
    }

    fileNameList objectNames;
    newInstance = word::null;

    if (Pstream::master())
    {
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
    }

    Pstream::scatter(newInstance);
    Pstream::scatter(objectNames);

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
        Pout<< "masterUncollatedFileOperation::readHeader:" << endl
            << "    objectPath:" << io.objectPath() << endl
            << "    fName     :" << fName << endl;
    }

    fileNameList filePaths(Pstream::nProcs());
    filePaths[Pstream::myProcNo()] = fName;
    Pstream::gatherList(filePaths);

    bool uniform = uniformFile(filePaths);
    Pstream::scatter(uniform);

    if (uniform)
    {
        if (Pstream::master())
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
        Pstream::scatter(ok);
        Pstream::scatter(io.headerClassName());
        Pstream::scatter(io.note());
    }
    else
    {
        boolList result(Pstream::nProcs(), false);
        wordList headerClassName(Pstream::nProcs());
        stringList note(Pstream::nProcs());
        if (Pstream::master())
        {
            forAll(filePaths, proci)
            {
                if (!filePaths[proci].empty())
                {
                    IFstream is(filePaths[proci]);

                    if (is.good())
                    {
                        result[proci] = io.readHeader(is);
                        headerClassName[proci] = io.headerClassName();
                        note[proci] = io.note();

                        if
                        (
                            io.headerClassName()
                         == decomposedBlockData::typeName
                        )
                        {
                            FatalErrorInFunction
                                << "Unexpected decomposedBlockData container"
                                << " for processor " << proci
                                << " file:" << filePaths[proci]
                                << ". A decomposedBlockData container should"
                                << " produce the same file name on all"
                                << " processors" << exit(FatalError);
                        }
                    }
                }
            }
        }
        ok = scatterList(result);
        io.headerClassName() = scatterList(headerClassName);
        io.note() = scatterList(note);
    }

    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::readHeader:" << " ok:" << ok
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
    const bool valid
) const
{
    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::readStream:"
            << " object : " << io.name()
            << " fName : " << fName << " valid:" << valid << endl;
    }


    autoPtr<ISstream> isPtr;
    bool isCollated = false;
    if (UPstream::master())
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
                IOobject headerIO(io);
                headerIO.readHeader(isPtr());

                if (headerIO.headerClassName() == decomposedBlockData::typeName)
                {
                    isCollated = true;
                }
            }

            if (!isCollated)
            {
                // Close file. Reopened below.
                isPtr.clear();
            }
        }
    }

    Pstream::scatter(isCollated);

    if (isCollated)
    {
        if (debug)
        {
            Pout<< "masterUncollatedFileOperation::readStream:"
                << " for object : " << io.name()
                << " starting collating input from " << fName << endl;
        }

        List<char> data;
        if (!Pstream::parRun())
        {
            // Analyse the objectpath to find out the processor we're trying
            // to access
            fileName path;
            fileName local;
            label proci = fileOperations::masterUncollatedFileOperation::
            splitProcessorPath
            (
                io.objectPath(),
                path,
                local
            );

            if (proci == -1)
            {
                FatalIOErrorInFunction(isPtr())
                    << "Could not detect processor number"
                    << " from objectPath:" << io.objectPath()
                    << exit(FatalIOError);
            }

            return decomposedBlockData::readBlock(proci, isPtr(), io);
        }
        else
        {
            // Get size of file (on master, scatter to slaves)
            off_t sz = fileSize(fName);

            // Read my data
            return decomposedBlockData::readBlocks
            (
                UPstream::worldComm,
                fName,
                isPtr,
                io,
                (
                    sz > off_t(maxMasterFileBufferSize)
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
            Pout<< "masterUncollatedFileOperation::readStream:"
                << " for object : " << io.name()
                << " starting separated input from " << fName << endl;
        }

        fileNameList filePaths(Pstream::nProcs());
        filePaths[Pstream::myProcNo()] = fName;
        Pstream::gatherList(filePaths);
        boolList procValid(Pstream::nProcs());
        procValid[Pstream::myProcNo()] = valid;
        Pstream::gatherList(procValid);

        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        if (Pstream::master())
        {
            //const bool uniform = uniformFile(filePaths);

            if (valid)
            {
                if (fName.empty())
                {
                    FatalErrorInFunction
                        << "cannot find file " << io.objectPath()
                        << exit(FatalError);
                }
                else
                {
                    autoPtr<IFstream> ifsPtr(new IFstream(fName));

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
            }

            // Read slave files
            for (label proci = 1; proci < Pstream::nProcs(); proci++)
            {
                if (debug)
                {
                    Pout<< "masterUncollatedFileOperation::readStream:"
                        << " For processor " << proci
                        << " opening " << filePaths[proci] << endl;
                }

                if (procValid[proci] && !filePaths[proci].empty())
                {
                    // Note: handle compression ourselves since size cannot
                    // be determined without actually uncompressing

                    if (Foam::exists(filePaths[proci]+".gz", false))
                    {
                        readAndSend
                        (
                            filePaths[proci],
                            IOstream::compressionType::COMPRESSED,
                            labelList(1, proci),
                            pBufs
                        );
                    }
                    else
                    {
                        readAndSend
                        (
                            filePaths[proci],
                            IOstream::compressionType::UNCOMPRESSED,
                            labelList(1, proci),
                            pBufs
                        );
                    }
                }
            }
        }

        labelList recvSizes;
        pBufs.finishedSends(recvSizes);

        // isPtr will be valid on master. Else the information is in the
        // PstreamBuffers

        if (Pstream::master())
        {
            if (!isPtr.valid())
            {
                return autoPtr<ISstream>(new dummyISstream());
            }
            else
            {
                return isPtr;
            }
        }
        else
        {
            if (valid)
            {
                UIPstream is(Pstream::masterNo(), pBufs);
                string buf(recvSizes[Pstream::masterNo()], '\0');
                if (recvSizes[Pstream::masterNo()] > 0)
                {
                    is.read(&buf[0], recvSizes[Pstream::masterNo()]);
                }

                if (debug)
                {
                    Pout<< "masterUncollatedFileOperation::readStream:"
                        << " Done reading " << buf.size() << " bytes" << endl;
                }
                isPtr.reset(new IStringStream(fName, buf));

                if (!io.readHeader(isPtr()))
                {
                    FatalIOErrorInFunction(isPtr())
                        << "problem while reading header for object "
                        << io.name() << exit(FatalIOError);
                }

                return isPtr;
            }
            else
            {
                return autoPtr<ISstream>(new dummyISstream());
            }
        }
    }
}


bool Foam::fileOperations::masterUncollatedFileOperation::read
(
    regIOobject& io,
    const bool masterOnly,
    const IOstream::streamFormat format,
    const word& typeName
) const
{
    bool ok = true;

    if (io.globalObject())
    {
        if (debug)
        {
            Pout<< "masterUncollatedFileOperation::read:"
                << "reading global object " << io.name() << endl;
        }

        bool ok = false;
        if (Pstream::master())
        {
            // Do master-only reading always.
            bool oldParRun = UPstream::parRun();
            UPstream::parRun() = false;

            ok = io.readData(io.readStream(typeName));
            io.close();

            UPstream::parRun() = oldParRun;
        }

        Pstream::scatter(ok);
        Pstream::scatter(io.headerClassName());
        Pstream::scatter(io.note());


        // scatter operation for regIOobjects

        // Get my communication order
        const List<Pstream::commsStruct>& comms =
        (
            (Pstream::nProcs() < Pstream::nProcsSimpleSum)
          ? Pstream::linearCommunication()
          : Pstream::treeCommunication()
        );
        const Pstream::commsStruct& myComm = comms[Pstream::myProcNo()];

        // Reveive from up
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
            Pout<< "masterUncollatedFileOperation::read:"
                << "reading local object " << io.name() << endl;
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
    const bool valid
) const
{
    fileName pathName(io.objectPath());

    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::writeObject:"
            << " io:" << pathName << " valid:" << valid << endl;
    }

    // Make sure to pick up any new times
    setTime(io.time());

    autoPtr<Ostream> osPtr
    (
        NewOFstream
        (
            pathName,
            fmt,
            ver,
            cmp,
            valid
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
    if (debug)
    {
        Pout<< "masterUncollatedFileOperation::findTimes:"
            << " Finding times in directory " << directory << endl;
    }

    HashPtrTable<instantList>::const_iterator iter = times_.find(directory);
    if (iter != times_.end())
    {
        if (debug)
        {
            Pout<< "masterUncollatedFileOperation::findTimes:"
                << " Found cached times:" << *iter() << endl;
        }
        return *iter();
    }
    else
    {
        instantList times;
        if (Pstream::master())
        {
            times = fileOperation::findTimes(directory, constantName);
        }
        Pstream::scatter(times);

        instantList* tPtr = new instantList(times.xfer());

        times_.insert(directory, tPtr);

        if (debug)
        {
            Pout<< "masterUncollatedFileOperation::findTimes:"
                << " Caching times:" << *tPtr << endl;
        }
        return *tPtr;
    }
}


void Foam::fileOperations::masterUncollatedFileOperation::setTime
(
    const Time& tm
) const
{
    HashPtrTable<instantList>::const_iterator iter = times_.find(tm.path());
    if (iter != times_.end())
    {
        if (debug)
        {
            Pout<< "masterUncollatedFileOperation::setTime:"
                << " Caching time " << tm.timeName()
                << " for case:" << tm.path() << endl;
        }

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
                times.append(timeNow);
                SubList<instant> realTimes(times, times.size()-1, 1);
                Foam::stableSort(realTimes);
            }
        }
        else
        {
            if (findSortedIndex(times, timeNow) == -1)
            {
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
    const fileName& filePath
) const
{
    if (Pstream::parRun())
    {
        // Insert logic of filePath. We assume that if a file is absolute
        // on the master it is absolute also on the slaves etc.

        fileNameList filePaths(Pstream::nProcs());
        filePaths[Pstream::myProcNo()] = filePath;
        Pstream::gatherList(filePaths);

        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        if (Pstream::master())
        {
            const bool uniform = uniformFile(filePaths);

            if (uniform)
            {
                if (debug)
                {
                    Pout<< "masterUncollatedFileOperation::NewIFstream:"
                        << " Opening global file " << filePath << endl;
                }

                IOstream::compressionType cmp
                (
                    Foam::exists(filePath+".gz", false)
                  ? IOstream::compressionType::COMPRESSED
                  : IOstream::compressionType::UNCOMPRESSED
                );

                labelList procs(Pstream::nProcs()-1);
                for (label proci = 1; proci < Pstream::nProcs(); proci++)
                {
                    procs[proci-1] = proci;
                }

                readAndSend(filePath, cmp, procs, pBufs);
            }
            else
            {
                for (label proci = 1; proci < Pstream::nProcs(); proci++)
                {
                    IOstream::compressionType cmp
                    (
                        Foam::exists(filePaths[proci]+".gz", false)
                      ? IOstream::compressionType::COMPRESSED
                      : IOstream::compressionType::UNCOMPRESSED
                    );

                    readAndSend
                    (
                        filePaths[proci],
                        cmp,
                        labelList(1, proci),
                        pBufs
                    );
                }
            }
        }


        labelList recvSizes;
        pBufs.finishedSends(recvSizes);

        if (Pstream::master())
        {
            // Read myself
            return autoPtr<ISstream>
            (
                new IFstream(filePaths[Pstream::masterNo()])
            );
        }
        else
        {
            if (debug)
            {
                Pout<< "masterUncollatedFileOperation::NewIFstream:"
                    << " Reading " << filePath
                    << " from processor " << Pstream::masterNo() << endl;
            }

            UIPstream is(Pstream::masterNo(), pBufs);
            string buf(recvSizes[Pstream::masterNo()], '\0');
            is.read(&buf[0], recvSizes[Pstream::masterNo()]);

            if (debug)
            {
                Pout<< "masterUncollatedFileOperation::NewIFstream:"
                    << " Done reading " << buf.size() << " bytes" << endl;
            }

            // Note: IPstream is not an IStream so use a IStringStream to
            //       convert the buffer. Note that we construct with a string
            //       so it holds a copy of the buffer.
            return autoPtr<ISstream>(new IStringStream(filePath, buf));
        }
    }
    else
    {
        // Read myself
        return autoPtr<ISstream>(new IFstream(filePath));
    }
}


Foam::autoPtr<Foam::Ostream>
Foam::fileOperations::masterUncollatedFileOperation::NewOFstream
(
    const fileName& pathName,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool valid
) const
{
    return autoPtr<Ostream>
    (
        new masterOFstream
        (
            pathName,
            fmt,
            ver,
            cmp,
            false,      // append
            valid
        )
    );
}


Foam::label Foam::fileOperations::masterUncollatedFileOperation::addWatch
(
    const fileName& fName
) const
{
    label watchFd;
    if (Pstream::master())
    {
        watchFd = monitor().addWatch(fName);
    }
    Pstream::scatter(watchFd);
    return watchFd;
}


bool Foam::fileOperations::masterUncollatedFileOperation::removeWatch
(
    const label watchIndex
) const
{
    bool ok;
    if (Pstream::master())
    {
        ok = monitor().removeWatch(watchIndex);
    }
    Pstream::scatter(ok);
    return ok;
}


Foam::label Foam::fileOperations::masterUncollatedFileOperation::findWatch
(
    const labelList& watchIndices,
    const fileName& fName
) const
{
    label index = -1;

    if (Pstream::master())
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
    Pstream::scatter(index);
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
    if (Pstream::master())
    {
        fName = monitor().getFile(watchIndex);
    }
    Pstream::scatter(fName);
    return fName;
}


void Foam::fileOperations::masterUncollatedFileOperation::updateStates
(
    const bool masterOnly,
    const bool syncPar
) const
{
    if (Pstream::master())
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
    if (Pstream::master())
    {
        state = monitor().getState(watchFd);
    }
    Pstream::scatter(state);
    return fileMonitor::fileState(state);
}


void Foam::fileOperations::masterUncollatedFileOperation::setUnmodified
(
    const label watchFd
) const
{
    if (Pstream::master())
    {
        monitor().setUnmodified(watchFd);
    }
}


// ************************************************************************* //
