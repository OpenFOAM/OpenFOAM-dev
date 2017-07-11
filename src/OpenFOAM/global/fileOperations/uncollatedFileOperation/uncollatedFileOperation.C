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

#include "uncollatedFileOperation.H"
#include "Time.H"
#include "IFstream.H"
#include "OFstream.H"
#include "addToRunTimeSelectionTable.H"
#include "masterUncollatedFileOperation.H"
#include "decomposedBlockData.H"
#include "dummyISstream.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(uncollatedFileOperation, 0);
    addToRunTimeSelectionTable(fileOperation, uncollatedFileOperation, word);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::fileName Foam::fileOperations::uncollatedFileOperation::filePathInfo
(
    const bool checkGlobal,
    const bool isFile,
    const IOobject& io
) const
{
    if (io.instance().isAbsolute())
    {
        fileName objectPath = io.instance()/io.name();

        if (isFileOrDir(isFile, objectPath))
        {
            return objectPath;
        }
        else
        {
            return fileName::null;
        }
    }
    else
    {
        fileName path = io.path();
        fileName objectPath = path/io.name();

        if (isFileOrDir(isFile, objectPath))
        {
            return objectPath;
        }
        else
        {
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
                // Constant & system can come from global case

                fileName parentObjectPath =
                    io.rootPath()/io.time().globalCaseName()
                   /io.instance()/io.db().dbDir()/io.local()/io.name();

                if (isFileOrDir(isFile, parentObjectPath))
                {
                    return parentObjectPath;
                }
            }

            // Check if parallel "procesors" directory
            if (io.time().processorCase())
            {
                fileName path = fileOperations::masterUncollatedFileOperation::
                processorsPath
                (
                    io,
                    io.instance()
                );
                fileName objectPath = path/io.name();

                if (isFileOrDir(isFile, objectPath))
                {
                    return objectPath;
                }
            }


            // Check for approximately same time. E.g. if time = 1e-2 and
            // directory is 0.01 (due to different time formats)
            if (!Foam::isDir(path))
            {
                word newInstancePath = io.time().findInstancePath
                (
                    instant(io.instance())
                );

                if (newInstancePath.size())
                {
                    fileName fName
                    (
                        io.rootPath()/io.caseName()
                       /newInstancePath/io.db().dbDir()/io.local()/io.name()
                    );

                    if (isFileOrDir(isFile, fName))
                    {
                        return fName;
                    }
                }
            }
        }

        return fileName::null;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::uncollatedFileOperation::uncollatedFileOperation
(
    const bool verbose
)
{
    if (verbose)
    {
        Info<< "I/O    : " << typeName << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::uncollatedFileOperation::~uncollatedFileOperation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileOperations::uncollatedFileOperation::mkDir
(
    const fileName& dir,
    mode_t mode
) const
{
    return Foam::mkDir(dir, mode);
}


bool Foam::fileOperations::uncollatedFileOperation::chMod
(
    const fileName& fName,
    mode_t mode
) const
{
    return Foam::chMod(fName, mode);
}


mode_t Foam::fileOperations::uncollatedFileOperation::mode
(
    const fileName& fName,
    const bool followLink
) const
{
    return Foam::mode(fName, followLink);
}


Foam::fileName::Type Foam::fileOperations::uncollatedFileOperation::type
(
    const fileName& fName,
    const bool followLink
) const
{
    return Foam::type(fName, followLink);
}


bool Foam::fileOperations::uncollatedFileOperation::exists
(
    const fileName& fName,
    const bool checkGzip,
    const bool followLink
) const
{
    return Foam::exists(fName, checkGzip, followLink);
}


bool Foam::fileOperations::uncollatedFileOperation::isDir
(
    const fileName& fName,
    const bool followLink
) const
{
    return Foam::isDir(fName, followLink);
}


bool Foam::fileOperations::uncollatedFileOperation::isFile
(
    const fileName& fName,
    const bool checkGzip,
    const bool followLink
) const
{
    return Foam::isFile(fName, checkGzip, followLink);
}


off_t Foam::fileOperations::uncollatedFileOperation::fileSize
(
    const fileName& fName,
    const bool followLink
) const
{
    return Foam::fileSize(fName, followLink);
}


time_t Foam::fileOperations::uncollatedFileOperation::lastModified
(
    const fileName& fName,
    const bool followLink
) const
{
    return Foam::lastModified(fName, followLink);
}


double Foam::fileOperations::uncollatedFileOperation::highResLastModified
(
    const fileName& fName,
    const bool followLink
) const
{
    return Foam::highResLastModified(fName, followLink);
}


bool Foam::fileOperations::uncollatedFileOperation::mvBak
(
    const fileName& fName,
    const std::string& ext
) const
{
    return Foam::mvBak(fName, ext);
}


bool Foam::fileOperations::uncollatedFileOperation::rm
(
    const fileName& fName
) const
{
    return Foam::rm(fName);
}


bool Foam::fileOperations::uncollatedFileOperation::rmDir
(
    const fileName& dir
) const
{
    return Foam::rmDir(dir);
}


Foam::fileNameList Foam::fileOperations::uncollatedFileOperation::readDir
(
    const fileName& dir,
    const fileName::Type type,
    const bool filtergz,
    const bool followLink
) const
{
    return Foam::readDir(dir, type, filtergz, followLink);
}


bool Foam::fileOperations::uncollatedFileOperation::cp
(
    const fileName& src,
    const fileName& dst,
    const bool followLink
) const
{
    return Foam::cp(src, dst, followLink);
}


bool Foam::fileOperations::uncollatedFileOperation::ln
(
    const fileName& src,
    const fileName& dst
) const
{
    return Foam::ln(src, dst);
}


bool Foam::fileOperations::uncollatedFileOperation::mv
(
    const fileName& src,
    const fileName& dst,
    const bool followLink
) const
{
    return Foam::mv(src, dst, followLink);
}


Foam::fileName Foam::fileOperations::uncollatedFileOperation::filePath
(
    const bool checkGlobal,
    const IOobject& io,
    const word& typeName
) const
{
    if (debug)
    {
        Pout<< "uncollatedFileOperation::filePath :"
            << " objectPath:" << io.objectPath()
            << " checkGlobal:" << checkGlobal << endl;
    }

    fileName objPath(filePathInfo(checkGlobal, true, io));

    if (debug)
    {
        Pout<< "uncollatedFileOperation::filePath :"
            << " Returning from file searching:" << endl
            << "    objectPath:" << io.objectPath() << endl
            << "    filePath  :" << objPath << endl << endl;
    }
    return objPath;
}


Foam::fileName Foam::fileOperations::uncollatedFileOperation::dirPath
(
    const bool checkGlobal,
    const IOobject& io
) const
{
    if (debug)
    {
        Pout<< "uncollatedFileOperation::dirPath :"
            << " objectPath:" << io.objectPath()
            << " checkGlobal:" << checkGlobal << endl;
    }

    fileName objPath(filePathInfo(checkGlobal, false, io));

    if (debug)
    {
        Pout<< "uncollatedFileOperation::dirPath :"
            << " Returning from directory searching:" << endl
            << "    objectPath:" << io.objectPath() << endl
            << "    dirPath   :" << objPath << endl << endl;
    }
    return objPath;
}


Foam::fileNameList Foam::fileOperations::uncollatedFileOperation::readObjects
(
    const objectRegistry& db,
    const fileName& instance,
    const fileName& local,
    word& newInstance
) const
{
    if (debug)
    {
        Pout<< "uncollatedFileOperation::readObjects :"
            << " db:" << db.objectPath()
            << " instance:" << instance << endl;
    }

    //- Use non-time searching version
    fileNameList objectNames
    (
        fileOperation::readObjects(db, instance, local, newInstance)
    );

    if (newInstance.empty())
    {
        // Find similar time
        fileName newInst = db.time().findInstancePath(instant(instance));
        if (!newInst.empty() && newInst != instance)
        {
            // Try with new time
            objectNames = fileOperation::readObjects
            (
                db,
                newInst,
                local,
                newInstance
            );
        }
    }

    if (debug)
    {
        Pout<< "uncollatedFileOperation::readObjects :"
            << " newInstance:" << newInstance
            << " objectNames:" << objectNames << endl;
    }

    return objectNames;
}


bool Foam::fileOperations::uncollatedFileOperation::readHeader
(
    IOobject& io,
    const fileName& fName,
    const word& typeName
) const
{
    if (fName.empty())
    {
        if (IOobject::debug)
        {
            InfoInFunction
                << "file " << io.objectPath() << " could not be opened"
                << endl;
        }

        return false;
    }

    autoPtr<ISstream> isPtr(NewIFstream(fName));

    if (!isPtr.valid() || !isPtr->good())
    {
        return false;
    }

    bool ok = io.readHeader(isPtr());

    if (io.headerClassName() == decomposedBlockData::typeName)
    {
        // Read the header inside the container (master data)
        ok = decomposedBlockData::readMasterHeader(io, isPtr());
    }

    return ok;
}


Foam::autoPtr<Foam::ISstream>
Foam::fileOperations::uncollatedFileOperation::readStream
(
    regIOobject& io,
    const fileName& fName,
    const word& typeName,
    const bool valid
) const
{
    autoPtr<ISstream> isPtr;

    if (!valid)
    {
        isPtr = autoPtr<ISstream>(new dummyISstream());
        return isPtr;
    }

    if (fName.empty())
    {
        FatalErrorInFunction
            << "cannot find file " << io.objectPath()
            << exit(FatalError);
    }

    isPtr = NewIFstream(fName);

    if (!isPtr.valid() || !isPtr->good())
    {
        FatalIOError
        (
            "uncollatedFileOperation::readStream()",
            __FILE__,
            __LINE__,
            fName,
            0
        )   << "cannot open file"
            << exit(FatalIOError);
    }
    else if (!io.readHeader(isPtr()))
    {
        FatalIOErrorInFunction(isPtr())
            << "problem while reading header for object " << io.name()
            << exit(FatalIOError);
    }

    if (io.headerClassName() != decomposedBlockData::typeName)
    {
        return isPtr;
    }
    else
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
                << "could not detect processor number"
                << " from objectPath:" << io.objectPath()
                << exit(FatalIOError);
        }

        // Read data and return as stream
        return decomposedBlockData::readBlock(proci, isPtr(), io);
    }
}


bool Foam::fileOperations::uncollatedFileOperation::read
(
    regIOobject& io,
    const bool masterOnly,
    const IOstream::streamFormat format,
    const word& typeName
) const
{
    bool ok = true;
    if (Pstream::master() || !masterOnly)
    {
        if (debug)
        {
            Pout<< "uncollatedFileOperation::read() : "
                << "reading object " << io.objectPath()
                << " from file " << endl;
        }

        // Set flag for e.g. codeStream
        const bool oldGlobal = io.globalObject();
        io.globalObject() = masterOnly;
        // If codeStream originates from dictionary which is
        // not IOdictionary we have a problem so use global
        const bool oldFlag = regIOobject::masterOnlyReading;
        regIOobject::masterOnlyReading = masterOnly;

        // Read file
        ok = io.readData(io.readStream(typeName));
        io.close();

        // Restore flags
        io.globalObject() = oldGlobal;
        regIOobject::masterOnlyReading = oldFlag;
    }

    if (masterOnly && Pstream::parRun())
    {
        // Master reads headerclassname from file. Make sure this gets
        // transfered as well as contents.
        Pstream::scatter(io.headerClassName());
        Pstream::scatter(io.note());

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
    return ok;
}


Foam::autoPtr<Foam::ISstream>
Foam::fileOperations::uncollatedFileOperation::NewIFstream
(
    const fileName& filePath
) const
{
    return autoPtr<ISstream>(new IFstream(filePath));
}


Foam::autoPtr<Foam::Ostream>
Foam::fileOperations::uncollatedFileOperation::NewOFstream
(
    const fileName& pathName,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool valid
) const
{
    return autoPtr<Ostream>(new OFstream(pathName, fmt, ver, cmp));
}


// ************************************************************************* //
