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

Class
    Foam::fileOperation

\*---------------------------------------------------------------------------*/

#ifndef fileOperation_H
#define fileOperation_H

#include "ISstream.H"
#include "Ostream.H"
#include "autoPtr.H"
#include "fileNameList.H"
#include "instantList.H"
#include "fileMonitor.H"
#include "labelList.H"
#include "Switch.H"
#include "tmpNrc.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class IOobject;
class regIOobject;
class objectRegistry;
class Time;

/*---------------------------------------------------------------------------*\
                        Class fileOperation Declaration
\*---------------------------------------------------------------------------*/

class fileOperation
{
public:

        //- Enumeration for the location of an IOobject
        enum pathType
        {
            NOTFOUND,               // not found
            ABSOLUTE,               // instance is absolute directory
            OBJECT,                 // io.objectPath() exists
            WRITEOBJECT,            // write path exists
            PROCUNCOLLATED,         // objectPath exists in processor0
            PROCBASEOBJECT,         // objectPath exists in specified, constant
                                    // processorsDir (usually 'processorsDDD')
            PROCOBJECT,             // objectPath exists in locally differing
                                    // processorsDir (e.g. 'processorsDDD_0-1')
            PARENTOBJECT,           // parent of object path
            FINDINSTANCE,           // file found in time directory
            PROCUNCOLLATEDINSTANCE, // as PROCUNCOLLATED but with instance
            PROCBASEINSTANCE,       // as PROCBASEOBJECT but with instance
            PROCINSTANCE            // as PROCOBJECT but with instance
        };

        static const NamedEnum<pathType, 12> pathTypeNames_;

        //- Description of processor directory naming:
        //  - processor directory naming
        //  - whether directory contains a range
        //    (so differs on different processors)
        //  - index in range
        typedef Tuple2<fileName, Tuple2<pathType, label>> dirIndex;

        typedef List<dirIndex> dirIndexList;


protected:

    // Protected data

        //- Communicator to use
        const label comm_;

        //- Detected processors directories
        mutable HashTable<dirIndexList> procsDirs_;

        //- file-change monitor for all registered files
        mutable autoPtr<fileMonitor> monitorPtr_;


   // Protected Member Functions

        fileMonitor& monitor() const;

        //- Sort directory entries according to time value
        static instantList sortTimes(const fileNameList&, const word&);

        //- Merge two times
        static void mergeTimes
        (
            const instantList& extraTimes,
            const word& constantName,
            instantList& times
        );

        //- Helper: check for file (isFile) or directory (!isFile)
        static bool isFileOrDir(const bool isFile, const fileName&);

        //- Detect presence of processorsDDD
        void cacheProcessorsPath(const fileName& fName) const;

        //- Lookup name of processorsDDD using cache. Return empty fileName
        //  if not found
        tmpNrc<dirIndexList> lookupProcessorsPath(const fileName&) const;

        //- Does ioobject exist. Is either a directory (empty name()) or
        //  a file
        bool exists(IOobject& io) const;


public:

    // Static data

        //- Return the processors directory name (usually "processors")
        static word processorsBaseDir;

        //- Default fileHandler
        static word defaultFileHandler;


    // Public data types

        //- Runtime type information
        TypeName("fileOperation");


    //- Static fileOperation
    static autoPtr<fileOperation> fileHandlerPtr_;


    // Constructors

        //- Construct null
        fileOperation(const label comm);


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fileOperation,
            word,
            (
                const bool verbose
            ),
            (verbose)
        );


    // Selectors

        //- Select type
        static autoPtr<fileOperation> New(const word& type, const bool verbose);


    //- Destructor
    virtual ~fileOperation();


    // Member Functions

        // OSSpecific equivalents

            //- Make directory
            virtual bool mkDir(const fileName&, mode_t=0777) const = 0;

            //- Set the file mode
            virtual bool chMod(const fileName&, const mode_t) const = 0;

            //- Return the file mode
            virtual mode_t mode
            (
                const fileName&,
                const bool checkVariants = true,
                const bool followLink = true
            ) const = 0;

            //- Return the file type: directory, file or link
            virtual fileType type
            (
                const fileName&,
                const bool checkVariants = true,
                const bool followLink = true
            ) const = 0;

            //- Does the name exist (as directory or file) in the file system?
            //  Optionally enable/disable check for gzip file.
            virtual bool exists
            (
                const fileName&,
                const bool checkVariants = true,
                const bool followLink = true
            ) const = 0;

            //- Does the name exist as a directory in the file system?
            virtual bool isDir
            (
                const fileName&,
                const bool followLink = true
            ) const = 0;

            //- Does the name exist as a file in the file system?
            //  Optionally enable/disable check for gzip file.
            virtual bool isFile
            (
                const fileName&,
                const bool checkVariants = true,
                const bool followLink = true
            ) const = 0;

            //- Return size of file
            virtual off_t fileSize
            (
                const fileName&,
                const bool checkVariants = true,
                const bool followLink = true
            ) const = 0;

            //- Return time of last file modification
            virtual time_t lastModified
            (
                const fileName&,
                const bool checkVariants = true,
                const bool followLink = true
            ) const = 0;

            //- Return time of last file modification
            virtual double highResLastModified
            (
                const fileName&,
                const bool checkVariants = true,
                const bool followLink = true
            ) const = 0;

            //- Read a directory and return the entries as a string list
            virtual fileNameList readDir
            (
                const fileName&,
                const fileType = fileType::file,
                const bool filterVariants = true,
                const bool followLink = true
            ) const = 0;

            //- Copy, recursively if necessary, the source to the destination
            virtual bool cp
            (
                const fileName& src,
                const fileName& dst,
                const bool followLink = true
            ) const = 0;

            //- Create a softlink. dst should not exist. Returns true if
            //  successful.
            virtual bool ln(const fileName& src, const fileName& dst) const = 0;

            //- Rename src to dst
            virtual bool mv
            (
                const fileName& src,
                const fileName& dst,
                const bool followLink = false
            ) const = 0;

            //- Rename to a corresponding backup file
            //  If the backup file already exists, attempt with
            //  "01" .. "99" suffix
            virtual bool mvBak
            (
                const fileName&,
                const std::string& ext = "bak"
            ) const = 0;

            //- Remove a file, returning true if successful otherwise false
            virtual bool rm(const fileName&) const = 0;

            //- Remove a directory and its contents
            virtual bool rmDir(const fileName&) const = 0;

            // //- Open a shared library. Return handle to library. Print error
            // //  message if library cannot be loaded (check = true)
            // virtual void* dlOpen
            // (
            //     const fileName& lib,
            //     const bool check = true
            // ) const = 0;


        // (reg)IOobject functionality

            //- Generate disk file name for object. Opposite of filePath.
            //  Optional wanted typeName.
            virtual fileName objectPath
            (
                const IOobject& io,
                const word& typeName
            ) const;

            //- Search for an object. checkGlobal : also check undecomposed case
            //  Optional wanted typeName.
            virtual fileName filePath
            (
                const bool checkGlobal,
                const IOobject&,
                const word& typeName
            ) const = 0;

            //- Search for a directory. checkGlobal : also check undecomposed
            //  case
            virtual fileName dirPath
            (
                const bool checkGlobal,
                const IOobject&
            ) const = 0;

            //- Search directory for objects. Used in IOobjectList.
            virtual fileNameList readObjects
            (
                const objectRegistry& db,
                const fileName& instance,
                const fileName& local,
                word& newInstance
            ) const;

            //- Read object header from supplied file
            virtual bool readHeader
            (
                IOobject&,
                const fileName&,
                const word& typeName
            ) const = 0;

            //- Reads header for regIOobject and returns an ISstream
            //  to read the contents.
            virtual autoPtr<ISstream> readStream
            (
                regIOobject&,
                const fileName&,
                const word& typeName,
                const bool read = true
            ) const = 0;

            //- Top-level read
            //  Format set to ASCII if defaultFormat and file format are ASCII
            virtual bool read
            (
                regIOobject&,
                const bool masterOnly,
                const IOstream::streamFormat defaultFormat,
                const word& typeName
            ) const = 0;

            //- Writes a regIOobject (so header, contents and divider).
            //  Returns success state. Default action is to write to
            //  the objectPath using writeData. If !write the
            //  file does not need to be written (this is used e.g. to
            //  suppress empty local lagrangian data)
            virtual bool writeObject
            (
                const regIOobject&,
                IOstream::streamFormat format=IOstream::ASCII,
                IOstream::versionNumber version=IOstream::currentVersion,
                IOstream::compressionType compression=IOstream::UNCOMPRESSED,
                const bool write = true
            ) const;


        // Filename (not IOobject) operations

            //- Search for a file or directory. Use IOobject version in
            //  preference
            virtual fileName filePath(const fileName&) const;

            //- Generate an ISstream that reads a file
            virtual autoPtr<ISstream> NewIFstream
            (
                const fileName& filePath,
                IOstream::streamFormat format=IOstream::ASCII,
                IOstream::versionNumber version=IOstream::currentVersion
            ) const = 0;

            //- Generate an Ostream that writes a file
            virtual autoPtr<Ostream> NewOFstream
            (
                const fileName& filePath,
                IOstream::streamFormat format=IOstream::ASCII,
                IOstream::versionNumber version=IOstream::currentVersion,
                IOstream::compressionType compression=IOstream::UNCOMPRESSED,
                const bool write = true
            ) const = 0;


        // File modification checking

            //- Add watching of a file. Returns handle
            virtual label addWatch(const fileName&) const;

            //- Remove watch on a file (using handle)
            virtual bool removeWatch(const label) const;

            //- Find index (or -1) of file in list of handles
            virtual label findWatch
            (
                const labelList& watchIndices,
                const fileName&
            ) const;

            //- Helper: add watches for list of regIOobjects
            virtual void addWatches(regIOobject&, const fileNameList&) const;

            //- Get name of file being watched (using handle)
            virtual fileName getFile(const label) const;

            //- Update state of all files
            virtual void updateStates
            (
                const bool masterOnly,
                const bool syncPar
            ) const;

            //- Get current state of file (using handle)
            virtual fileMonitor::fileState getState(const label) const;

            //- Set current state of file (using handle) to unmodified
            virtual void setUnmodified(const label) const;


        // Other

            //- Actual name of processors dir (for use in mode PROCOBJECT,
            //  PROCINSTANCE)
            virtual word processorsDir(const IOobject& io) const
            {
                return processorsBaseDir;
            }

            //- Actual name of processors dir (for use in mode PROCOBJECT,
            //  PROCINSTANCE)
            virtual word processorsDir(const fileName&) const
            {
                return processorsBaseDir;
            }

            //- Set number of processor directories/results. Only used in
            //  decomposePar
            virtual void setNProcs(const label nProcs);

            //- Get number of processor directories/results. Used for e.g.
            //  reconstructPar, argList checking
            virtual label nProcs
            (
                const fileName& dir,
                const fileName& local = ""
            ) const;

            //- Get sorted list of times
            virtual instantList findTimes(const fileName&, const word&) const;

            //- Find instance where IOobject is. Fails if cannot be found
            //  and readOpt() is MUST_READ/MUST_READ_IF_MODIFIED. Otherwise
            //  returns stopInstance.
            virtual IOobject findInstance
            (
                const IOobject& io,
                const scalar startValue,
                const word& stopInstance
            ) const;

            //- Callback for time change
            virtual void setTime(const Time&) const
            {}

            //- Forcibly wait until all output done. Flush any cached data
            virtual void flush() const;

            //- Generate path (like io.path) from root+casename with any
            //  'processorXXX' replaced by procDir (usually 'processors')
            fileName processorsCasePath
            (
                const IOobject&,
                const word& procDir
            ) const;

            //- Generate path (like io.path) with provided instance and any
            //  'processorXXX' replaced by procDir (usually 'processors')
            fileName processorsPath
            (
                const IOobject&,
                const word& instance,
                const word& procDir
            ) const;

            //- Operating on fileName: replace processorXXX with procDir
            fileName processorsPath(const fileName&, const word& procDir) const;

            //- Split fileName into part before 'processor' and part after.
            //  Returns -1 or processor number and optionally number
            //  of processors. Use with care.
            //  - path/"processor"+Foam::name(proci)/local reconstructs input
            //  - path/"processors"+Foam::name(nProcs)/local reconstructs
            //    collated processors equivalence
            static label splitProcessorPath
            (
                const fileName&,
                fileName& path,
                fileName& procDir,
                fileName& local,
                label& groupStart,
                label& groupSize,
                label& nProcs
            );

            //- Detect processor number from '/aa/bb/processorDDD/cc'
            static label detectProcessorPath(const fileName&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Global declarations: defined in fileOperation.C

//- Get current file handler
const fileOperation& fileHandler();

//- Reset file handler
void fileHandler(autoPtr<fileOperation>&);

//- Recursively search the given directory for the file
//  returning the path relative to the directory or
//  fileName::null if not found
fileName search(const word& file, const fileName& directory);

//- Copy all the files from the source to the target directory
void cpFiles(const fileName& srcDir, const fileName& targetDir);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
