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

#include "fileOperation.H"
#include "uncollatedFileOperation.H"
#include "regIOobject.H"
#include "argList.H"
#include "HashSet.H"
#include "objectRegistry.H"
#include "decomposedBlockData.H"
#include "polyMesh.H"
#include "registerSwitch.H"
#include "Time.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
    autoPtr<fileOperation> fileOperation::fileHandlerPtr_;

    defineTypeNameAndDebug(fileOperation, 0);
    defineRunTimeSelectionTable(fileOperation, word);

    template<>
    const char* Foam::NamedEnum
    <
        fileOperation::pathType,
        12
    >::names[] =
    {
        "notFound",
        "absolute",
        "objectPath",
        "writeObject",
        "uncollatedProc",
        "globalProc",
        "localProc",
        "parentObjectPath",
        "findInstance",
        "uncollatedProcInstance",
        "globalProcInstance",
        "localProcInstance"
    };
    const NamedEnum<fileOperation::pathType, 12> fileOperation::pathTypeNames_;

    word fileOperation::defaultFileHandler
    (
        debug::optimisationSwitches().lookupOrAddDefault
        (
            "fileHandler",
            // Foam::fileOperations::uncollatedFileOperation::typeName,
            word("uncollated"),
            false,
            false
        )
    );

    word fileOperation::processorsBaseDir = "processors";
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::fileMonitor& Foam::fileOperation::monitor() const
{
    if (!monitorPtr_.valid())
    {
        monitorPtr_.reset
        (
            new fileMonitor
            (
                regIOobject::fileModificationChecking == IOobject::inotify
             || regIOobject::fileModificationChecking == IOobject::inotifyMaster
            )
        );
    }
    return monitorPtr_();
}


Foam::instantList Foam::fileOperation::sortTimes
(
    const fileNameList& dirEntries,
    const word& constantName
)
{
    // Initialise instant list
    instantList Times(dirEntries.size() + 1);
    label nTimes = 0;

    // Check for "constant"
    bool haveConstant = false;
    forAll(dirEntries, i)
    {
        if (dirEntries[i] == constantName)
        {
            Times[nTimes].value() = 0;
            Times[nTimes].name() = dirEntries[i];
            nTimes++;
            haveConstant = true;
            break;
        }
    }

    // Read and parse all the entries in the directory
    forAll(dirEntries, i)
    {
        scalar timeValue;
        if (readScalar(dirEntries[i].c_str(), timeValue))
        {
            Times[nTimes].value() = timeValue;
            Times[nTimes].name() = dirEntries[i];
            nTimes++;
        }
    }

    // Reset the length of the times list
    Times.setSize(nTimes);

    if (haveConstant)
    {
        if (nTimes > 2)
        {
            std::sort(&Times[1], Times.end(), instant::less());
        }
    }
    else if (nTimes > 1)
    {
        std::sort(&Times[0], Times.end(), instant::less());
    }

    return Times;
}


void Foam::fileOperation::mergeTimes
(
    const instantList& extraTimes,
    const word& constantName,
    instantList& times
)
{
    if (extraTimes.size())
    {
        bool haveConstant =
        (
            times.size() > 0
         && times[0].name() == constantName
        );

        bool haveExtraConstant =
        (
            extraTimes.size() > 0
         && extraTimes[0].name() == constantName
        );

        // Combine times
        instantList combinedTimes(times.size()+extraTimes.size());
        label sz = 0;
        label extrai = 0;
        if (haveExtraConstant)
        {
            extrai = 1;
            if (!haveConstant)
            {
                combinedTimes[sz++] = extraTimes[0];    // constant
            }
        }
        forAll(times, i)
        {
            combinedTimes[sz++] = times[i];
        }
        for (; extrai < extraTimes.size(); extrai++)
        {
            combinedTimes[sz++] = extraTimes[extrai];
        }
        combinedTimes.setSize(sz);
        times.transfer(combinedTimes);

        // Sort
        if (times.size() > 1)
        {
            label starti = 0;
            if (times[0].name() == constantName)
            {
                starti = 1;
            }
            std::sort(&times[starti], times.end(), instant::less());

            // Filter out duplicates
            label newi = starti+1;
            for (label i = newi; i < times.size(); i++)
            {
                if (times[i].value() != times[i-1].value())
                {
                    if (newi != i)
                    {
                        times[newi] = times[i];
                    }
                    newi++;
                }
            }

            times.setSize(newi);
        }
    }
}


bool Foam::fileOperation::isFileOrDir(const bool isFile, const fileName& f)
{
    return
        (isFile && Foam::isFile(f))
     || (!isFile && Foam::isDir(f));
}


Foam::tmpNrc<Foam::fileOperation::dirIndexList>
Foam::fileOperation::lookupProcessorsPath(const fileName& fName) const
{
    // If path is local to a processor (e.g. contains 'processor2')
    // find the corresponding actual processor directory (e.g. 'processors4')
    // and index (2)

    fileName path;
    fileName pDir;
    fileName local;
    label gStart;
    label gSz;
    label numProcs;
    label proci =
        splitProcessorPath(fName, path, pDir, local, gStart, gSz, numProcs);

    if (proci != -1)
    {
        const fileName procPath(path/pDir);

        HashTable<dirIndexList>::const_iterator iter =
            procsDirs_.find(procPath);

        if (iter != procsDirs_.end())
        {
            return iter();
        }

        // Read all directories to see any beginning with processor
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        DynamicList<dirIndex> procDirs;

        // Note: use parallel synchronised reading so cache will be same
        //       order on all processors
        fileNameList dirNames(readDir(path, fileName::Type::DIRECTORY));

        // Extract info from processorsDDD or processorDDD:
        // - highest processor number
        // - directory+offset containing data for proci
        label maxProc = -1;

        forAll(dirNames, i)
        {
            const fileName& dirN = dirNames[i];

            // Analyse directory name
            fileName rp, rd, rl;
            label rStart, rSize, rNum;
            label readProci =
                splitProcessorPath(dirN, rp, rd, rl, rStart, rSize, rNum);
            maxProc = max(maxProc, readProci);

            if (proci == readProci)
            {
                // Found "processorDDD". No need for index.
                procDirs.append
                (
                    dirIndex
                    (
                        dirN,
                        Tuple2<pathType, label>(PROCUNCOLLATED, -1)
                    )
                );
            }
            else if (proci >= rStart && proci < rStart+rSize)
            {
                // "processorsDDD_start-end"
                // Found the file that contains the data for proci
                procDirs.append
                (
                    dirIndex
                    (
                        dirN,
                        Tuple2<pathType, label>(PROCOBJECT, proci-rStart)
                    )
                );
            }
            if (rNum != -1)
            {
                // Direct detection of processorsDDD
                maxProc = rNum-1;

                if (rStart == -1)
                {
                    // "processorsDDD"
                    procDirs.append
                    (
                        dirIndex
                        (
                            dirN,
                            Tuple2<pathType, label>(PROCBASEOBJECT, proci)
                        )
                    );
                }
            }
        }
        if (!Pstream::parRun())
        {
            // If (as a side effect) we found the number of decompositions
            // use it
            if (maxProc != -1)
            {
                const_cast<fileOperation&>(*this).setNProcs(maxProc+1);
            }
        }

        if (returnReduce(procDirs.size(), sumOp<label>()))
        {
            procsDirs_.insert(procPath, procDirs);

            if (debug)
            {
                Pout<< "fileOperation::lookupProcessorsPath : For:" << procPath
                    << " detected:" << procDirs << endl;
            }

            // Make sure to return a reference
            return procsDirs_[procPath];
        }
    }
    return tmpNrc<dirIndexList>(new dirIndexList(0, dirIndex()));
}


bool Foam::fileOperation::exists(IOobject& io) const
{
    // Generate output filename for object
    fileName objPath(objectPath(io, word::null));

    // Test for either directory or a (valid) file & IOobject
    bool ok;
    if (io.name().empty())
    {
        ok = isDir(objPath);
    }
    else
    {
        ok =
            isFile(objPath)
         && io.typeHeaderOk<IOList<label>>(false);// object with local scope
    }

    if (!ok)
    {
        // Re-test with searched for objectPath. This is for backwards
        // compatibility
        fileName originalPath(filePath(io.objectPath()));
        if (originalPath != objPath)
        {
            // Test for either directory or a (valid) file & IOobject
            if (io.name().empty())
            {
                ok = isDir(originalPath);
            }
            else
            {
                ok =
                    isFile(originalPath)
                 && io.typeHeaderOk<IOList<label>>(false);
            }
        }
    }

    return ok;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperation::fileOperation(label comm)
:
    comm_(comm)
{}


Foam::autoPtr<Foam::fileOperation> Foam::fileOperation::New
(
    const word& type,
    const bool verbose
)
{
    wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_->find(type);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown fileOperation type "
            << type << nl << nl
            << "Valid fileOperation types are" << endl
            << wordConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return autoPtr<fileOperation>(cstrIter()(verbose));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperation::~fileOperation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::fileOperation::objectPath
(
    const IOobject& io,
    const word& typeName
) const
{
    return io.objectPath();
}


bool Foam::fileOperation::writeObject
(
    const regIOobject& io,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool valid
) const
{
    if (valid)
    {
        fileName pathName(io.objectPath());

        mkDir(pathName.path());

        autoPtr<Ostream> osPtr
        (
            NewOFstream
            (
                pathName,
                fmt,
                ver,
                cmp
            )
        );

        if (!osPtr.valid())
        {
            return false;
        }

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
    }
    return true;
}


Foam::fileName Foam::fileOperation::filePath(const fileName& fName) const
{
    if (debug)
    {
        Pout<< "fileOperation::filePath :" << " fName:" << fName << endl;
    }

    fileName path;
    fileName pDir;
    fileName local;
    label gStart;
    label gSz;
    label numProcs;
    label proci =
        splitProcessorPath(fName, path, pDir, local, gStart, gSz, numProcs);

    if (numProcs != -1)
    {
        WarningInFunction << "Filename is already adapted:" << fName << endl;
    }

    // Give preference to processors variant
    if (proci != -1)
    {
        // Get all processor directories
        tmpNrc<dirIndexList> procDirs(lookupProcessorsPath(fName));
        forAll(procDirs(), i)
        {
            const fileName& procDir = procDirs()[i].first();

            fileName collatedName(path/procDir/local);
            if (exists(collatedName))
            {
                if (debug)
                {
                    Pout<< "fileOperation::filePath : " << collatedName << endl;
                }
                return collatedName;
            }
        }
    }

    if (exists(fName))
    {
        if (debug)
        {
            Pout<< "fileOperation::filePath : " << fName << endl;
        }
        return fName;
    }
    else
    {
        if (debug)
        {
            Pout<< "fileOperation::filePath : Not found" << endl;
        }
        return fileName::null;
    }
}


Foam::label Foam::fileOperation::addWatch(const fileName& fName) const
{
    return monitor().addWatch(fName);
}


bool Foam::fileOperation::removeWatch(const label watchIndex) const
{
    return monitor().removeWatch(watchIndex);
}


Foam::label Foam::fileOperation::findWatch
(
    const labelList& watchIndices,
    const fileName& fName
) const
{
    forAll(watchIndices, i)
    {
        if (getFile(watchIndices[i]) == fName)
        {
            return i;
        }
    }
    return -1;
}


void Foam::fileOperation::addWatches
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


Foam::fileName Foam::fileOperation::getFile(const label watchIndex) const
{
    return monitor().getFile(watchIndex);
}


void Foam::fileOperation::updateStates
(
    const bool masterOnly,
    const bool syncPar
) const
{
    monitor().updateStates(masterOnly, Pstream::parRun());
}


Foam::fileMonitor::fileState Foam::fileOperation::getState
(
    const label watchFd
) const
{
    return monitor().getState(watchFd);
}


void Foam::fileOperation::setUnmodified(const label watchFd) const
{
    monitor().setUnmodified(watchFd);
}


Foam::instantList Foam::fileOperation::findTimes
(
    const fileName& directory,
    const word& constantName
) const
{
    if (debug)
    {
        Pout<< "fileOperation::findTimes : Finding times in directory "
            << directory << endl;
    }

    // Read directory entries into a list
    fileNameList dirEntries
    (
        Foam::readDir
        (
            directory,
            fileName::DIRECTORY
        )
    );

    instantList times = sortTimes(dirEntries, constantName);


    // Get all processor directories
    tmpNrc<dirIndexList> procDirs(lookupProcessorsPath(directory));
    forAll(procDirs(), i)
    {
        const fileName& procDir = procDirs()[i].first();
        fileName collDir(processorsPath(directory, procDir));
        if (!collDir.empty() && collDir != directory)
        {
            fileNameList extraEntries
            (
                Foam::readDir
                (
                    collDir,
                    fileName::DIRECTORY
                )
            );
            mergeTimes
            (
                sortTimes(extraEntries, constantName),
                constantName,
                times
            );
        }
    }

    if (debug)
    {
        Pout<< "fileOperation::findTimes : Found times:" << times << endl;
    }
    return times;
}


Foam::IOobject Foam::fileOperation::findInstance
(
    const IOobject& startIO,
    const scalar startValue,
    const word& stopInstance
) const
{
    const Time& time = startIO.time();

    IOobject io(startIO);

    // Note: - if name is empty, just check the directory itself
    //       - check both for isFile and headerOk since the latter does a
    //         filePath so searches for the file.
    //       - check for an object with local file scope (so no looking up in
    //         parent directory in case of parallel)

    if (exists(io))
    {
        if (debug)
        {
            InfoInFunction
                << "Found exact match for \"" << io.name()
                << "\" in " << io.instance()/io.local()
                << endl;
        }

        return io;
    }

    // Search back through the time directories to find the time
    // closest to and lower than current time

    instantList ts = time.times();
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
        // Shortcut: if actual directory is the timeName we've already tested it
        if (ts[instanceI].name() == startIO.instance())
        {
            continue;
        }

        io.instance() = ts[instanceI].name();
        if (exists(io))
        {
            if (debug)
            {
                InfoInFunction
                    << "Found exact match for \"" << io.name()
                    << "\" in " << io.instance()/io.local()
                    << endl;
            }

            return io;
        }

        // Check if hit minimum instance
        if (ts[instanceI].name() == stopInstance)
        {
            if (debug)
            {
                InfoInFunction
                    << "Hit stopInstance " << stopInstance << endl;
            }

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
                        << io.local() << " in times " << startIO.instance()
                        << " down to " << stopInstance
                        << exit(FatalError);
                }
                else
                {
                    FatalErrorInFunction
                        << "Cannot find file \"" << io.name()
                        << "\" in directory " << io.local()
                        << " in times " << startIO.instance()
                        << " down to " << stopInstance
                        << exit(FatalError);
                }
            }

            return io;
        }
    }

    // times() usually already includes the constant() so would have been
    // checked above. Re-test if
    // - times() is empty. Sometimes this can happen (e.g. decomposePar with
    //   collated)
    // - times()[0] is not constant
    if (!ts.size() || ts[0].name() != time.constant())
    {
        // Note. This needs to be a hard-coded constant, rather than the
        // constant function of the time, because the latter points to
        // the case constant directory in parallel cases

        io.instance() = time.constant();
        if (exists(io))
        {
            if (debug)
            {
                InfoInFunction
                    << "Found constant match for \"" << io.name()
                    << "\" in " << io.instance()/io.local()
                    << endl;
            }
            return io;
        }
    }


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

    return io;
}


Foam::fileNameList Foam::fileOperation::readObjects
(
    const objectRegistry& db,
    const fileName& instance,
    const fileName& local,
    word& newInstance
) const
{
    if (debug)
    {
        Pout<< "fileOperation::readObjects :"
            << " db:" << db.objectPath()
            << " instance:" << instance << endl;
    }

    fileName path(db.path(instance, db.dbDir()/local));

    newInstance = word::null;
    fileNameList objectNames;

    if (Foam::isDir(path))
    {
        newInstance = instance;
        objectNames = Foam::readDir(path, fileName::FILE);
    }
    else
    {
        // Get processors equivalent of path
        fileName procsPath(filePath(path));

        if (!procsPath.empty())
        {
            newInstance = instance;
            objectNames = Foam::readDir(procsPath, fileName::FILE);
        }
    }
    return objectNames;
}


void Foam::fileOperation::setNProcs(const label nProcs)
{}


Foam::label Foam::fileOperation::nProcs
(
    const fileName& dir,
    const fileName& local
) const
{
    label nProcs = 0;
    if (Pstream::master(comm_))
    {
        fileNameList dirNames(Foam::readDir(dir, fileName::Type::DIRECTORY));

        // Detect any processorsDDD or processorDDD
        label maxProc = -1;
        forAll(dirNames, i)
        {
            const fileName& dirN = dirNames[i];

            fileName path, pDir, local;
            label start, size, n;
            maxProc = max
            (
                maxProc,
                splitProcessorPath(dirN, path, pDir, local, start, size, n)
            );
            if (n != -1)
            {
                // Direct detection of processorsDDD
                maxProc = n-1;
                break;
            }
        }
        nProcs = maxProc+1;


        if (nProcs == 0 && Foam::isDir(dir/processorsBaseDir))
        {
            fileName pointsFile
            (
                dir
               /processorsBaseDir
               /"constant"
               /local
               /polyMesh::meshSubDir
               /"points"
            );

            if (Foam::isFile(pointsFile))
            {
                nProcs = decomposedBlockData::numBlocks(pointsFile);
            }
            else
            {
                WarningInFunction << "Cannot read file " << pointsFile
                    << " to determine the number of decompositions."
                    << " Returning 1" << endl;
            }
        }
    }
    Pstream::scatter(nProcs, Pstream::msgType(), comm_);
    return nProcs;
}


void Foam::fileOperation::flush() const
{
    if (debug)
    {
        Pout<< "fileOperation::flush : clearing processor directories cache"
            << endl;
    }
    procsDirs_.clear();
}


Foam::fileName Foam::fileOperation::processorsCasePath
(
    const IOobject& io,
    const word& procsDir
) const
{
    return io.rootPath()/io.time().globalCaseName()/procsDir;
}


Foam::fileName Foam::fileOperation::processorsPath
(
    const IOobject& io,
    const word& instance,
    const word& procsDir
) const
{
    return
        processorsCasePath(io, procsDir)
       /instance
       /io.db().dbDir()
       /io.local();
}


Foam::fileName Foam::fileOperation::processorsPath
(
    const fileName& dir,
    const word& procsDir
) const
{
    // Check if directory is processorDDD
    word caseName(dir.name());

    std::string::size_type pos = caseName.find("processor");
    if (pos == 0)
    {
        if (caseName.size() <= 9 || caseName[9] == 's')
        {
            WarningInFunction << "Directory " << dir
                << " does not end in old-style processorDDD" << endl;
        }

        return dir.path()/procsDir;
    }
    else
    {
        return fileName::null;
    }
}


Foam::label Foam::fileOperation::splitProcessorPath
(
    const fileName& objectPath,
    fileName& path,
    fileName& procDir,
    fileName& local,

    label& groupStart,
    label& groupSize,

    label& nProcs
)
{
    path.clear();
    procDir.clear();
    local.clear();

    // Potentially detected start of number of processors in local group
    groupStart = -1;
    groupSize = 0;

    // Potentially detected number of processors
    nProcs = -1;

    // Search for processor at start of line or /processor
    std::string::size_type pos = objectPath.find("processor");
    if (pos == string::npos)
    {
        return -1;
    }

    // "processorDDD"
    // "processorsNNN"
    // "processorsNNN_AA-BB"


    if (pos > 0 && objectPath[pos-1] != '/')
    {
        // Directory not starting with "processor" e.g. "somenamewithprocessor"
        return -1;
    }

    procDir = objectPath;

    // Strip leading directory
    if (pos > 0)
    {
        path = objectPath.substr(0, pos-1);
        procDir = objectPath.substr(pos);
    }

    // Strip trailing local directory
    pos = procDir.find('/');
    if (pos != string::npos)
    {
        local = procDir.substr(pos+1);
        procDir = procDir.substr(0, pos);
    }

    // Now procDir is e.g.
    // - processor0
    // - processors0
    // - processorBananas

    // Look for number after "processor"

    fileName f(procDir.substr(9));

    if (f.size() && f[0] == 's')
    {
        // "processsorsNNN"

        f = f.substr(1);

        // Detect "processorsNNN_AA-BB"
        {
            std::string::size_type fromStart = f.find("_");
            std::string::size_type toStart = f.find("-");
            if (fromStart != string::npos && toStart != string::npos)
            {
                string nProcsName(f.substr(0, fromStart));
                string fromName(f.substr(fromStart+1, toStart-(fromStart+1)));
                string toName(f.substr(toStart+1));

                label groupEnd = -1;
                if
                (
                    Foam::read(fromName.c_str(), groupStart)
                 && Foam::read(toName.c_str(), groupEnd)
                 && Foam::read(nProcsName.c_str(), nProcs)
                )
                {
                    groupSize = groupEnd-groupStart+1;
                    return -1;
                }
            }
        }

        // Detect "processorsN"
        label n;
        if (Foam::read(f.c_str(), n))
        {
            nProcs = n;
        }
        return -1;
    }
    else
    {
        // Detect "processorN"
        label proci;
        if (Foam::read(f.c_str(), proci))
        {
            return proci;
        }
        else
        {
            return -1;
        }
    }
}


Foam::label Foam::fileOperation::detectProcessorPath(const fileName& fName)
{
    fileName path, pDir, local;
    label start, size, nProcs;
    return splitProcessorPath(fName, path, pDir, local, start, size, nProcs);
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

const Foam::fileOperation& Foam::fileHandler()
{
    if (!fileOperation::fileHandlerPtr_.valid())
    {
        word handler(getEnv("FOAM_FILEHANDLER"));

        if (!handler.size())
        {
            handler = fileOperation::defaultFileHandler;
        }

        fileOperation::fileHandlerPtr_ = fileOperation::New(handler, true);
    }
    return fileOperation::fileHandlerPtr_();
}


void Foam::fileHandler(autoPtr<fileOperation>& newHandlerPtr)
{
    if (fileOperation::fileHandlerPtr_.valid())
    {
        if
        (
            newHandlerPtr.valid()
         && newHandlerPtr->type() == fileOperation::fileHandlerPtr_->type()
        )
        {
            return;
        }
    }
    fileOperation::fileHandlerPtr_.clear();

    if (newHandlerPtr.valid())
    {
        fileOperation::fileHandlerPtr_ = newHandlerPtr;
    }
}


Foam::fileName Foam::search(const word& file, const fileName& directory)
{
    // Search the current directory for the file
    fileNameList files(fileHandler().readDir(directory));
    forAll(files, i)
    {
        if (files[i] == file)
        {
            return directory/file;
        }
    }

    // If not found search each of the sub-directories
    fileNameList dirs(fileHandler().readDir(directory, fileName::DIRECTORY));
    forAll(dirs, i)
    {
        fileName path = search(file, directory/dirs[i]);
        if (path != fileName::null)
        {
            return path;
        }
    }

    return fileName::null;
}


void Foam::cpFiles(const fileName& srcDir, const fileName& targetDir)
{
    mkDir(targetDir);

    const fileNameList srcFiles(readDir(srcDir, fileName::FILE, true));

    forAll(srcFiles, filei)
    {
        cp(srcDir/srcFiles[filei], targetDir);
    }
}


// ************************************************************************* //
