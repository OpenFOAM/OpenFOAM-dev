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
    Foam::fileOperations::masterUncollatedFileOperationInitialise

Description
    fileOperations that performs all file operations on the master processor.
    Requires the calls to be parallel synchronised!

    Limitations: - no /processor in filename
                 - no /uniform/ in the filename

    The main logic is in ::filePath which returns a
    - same path on all processors. This can either be a global file
      (system/controlDict, processorXXX/0/uniform/) or a collated file
      (processors/0/p)
    - same path on all processors of the local communicator
      (processors4_0-1/0/p)
    - different path on all processors (processor0/0/p)

    system/controlDict:
    filePath worldmaster: <globalRoot>/system/controlDict
             localmaster:           ,,
             slave      :           ,,

    processor0/uniform/time
    filePath worldmaster: <globalRoot>/processorXXX/uniform/time
             localmaster:           ,,
             slave      :           ,,

    processors0/0/p
    processors10/0/p
    processors10_2-4/0/p

\*---------------------------------------------------------------------------*/

#ifndef fileOperations_masterUncollatedFileOperation_H
#define fileOperations_masterUncollatedFileOperation_H

#include "fileOperation.H"
#include "HashPtrTable.H"
#include "unthreadedInitialise.H"
#include "boolList.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class PstreamBuffers;

namespace fileOperations
{

/*---------------------------------------------------------------------------*\
                Class masterUncollatedFileOperation Declaration
\*---------------------------------------------------------------------------*/

class masterUncollatedFileOperation
:
    public fileOperation
{
protected:

    // Protected data

        //- Any communicator allocated by me
        const label myComm_;

        //- Cached times for a given directory
        mutable HashPtrTable<instantList> times_;


    // Protected classes

        class mkDirOp
        {
            const mode_t mode_;
        public:
            mkDirOp(const mode_t mode)
            :
                mode_(mode)
            {}

            bool operator()(const fileName& fName) const
            {
                return Foam::mkDir(fName, mode_);
            }
        };

        class chModOp
        {
            const mode_t mode_;
        public:
            chModOp(const mode_t mode)
            :
                mode_(mode)
            {}

            bool operator()(const fileName& fName) const
            {
                return Foam::chMod(fName, mode_);
            }
        };

        class modeOp
        {
            const bool checkVariants_;
            const bool followLink_;
        public:
            modeOp(const bool checkVariants, const bool followLink)
            :
                checkVariants_(checkVariants),
                followLink_(followLink)
            {}

            mode_t operator()(const fileName& fName) const
            {
                return Foam::mode(fName, checkVariants_, followLink_);
            }
        };

        class typeOp
        {
            const bool checkVariants_;
            const bool followLink_;
        public:
            typeOp(const bool checkVariants, const bool followLink)
            :
                checkVariants_(checkVariants),
                followLink_(followLink)
            {}

            label operator()(const fileName& fName) const
            {
                return label(Foam::type(fName, checkVariants_, followLink_));
            }
        };

        class existsOp
        {
            const bool checkVariants_;
            const bool followLink_;
        public:
            existsOp(const bool checkVariants, const bool followLink)
            :
                checkVariants_(checkVariants),
                followLink_(followLink)
            {}

            bool operator()(const fileName& fName) const
            {
                return Foam::exists(fName, checkVariants_, followLink_);
            }
        };

        class isDirOp
        {
            const bool followLink_;
        public:
            isDirOp(const bool followLink)
            :
                followLink_(followLink)
            {}

        public:
            bool operator()(const fileName& fName) const
            {
                return Foam::isDir(fName, followLink_);
            }
        };

        class isFileOp
        {
            const bool checkVariants_;
            const bool followLink_;
        public:
            isFileOp(const bool checkVariants, const bool followLink)
            :
                checkVariants_(checkVariants),
                followLink_(followLink)
            {}
        public:
            bool operator()(const fileName& fName) const
            {
                return Foam::isFile(fName, checkVariants_, followLink_);
            }
        };

        class fileSizeOp
        {
            const bool checkVariants_;
            const bool followLink_;
        public:
            fileSizeOp(const bool checkVariants, const bool followLink)
            :
                checkVariants_(checkVariants),
                followLink_(followLink)
            {}

        public:
            off_t operator()(const fileName& fName) const
            {
                return Foam::fileSize(fName, checkVariants_, followLink_);
            }
        };

        class lastModifiedOp
        {
            const bool checkVariants_;
            const bool followLink_;
        public:
            lastModifiedOp(const bool checkVariants, const bool followLink)
            :
                checkVariants_(checkVariants),
                followLink_(followLink)
            {}

        public:
            time_t operator()(const fileName& fName) const
            {
                return Foam::lastModified(fName, checkVariants_, followLink_);
            }
        };

        class lastModifiedHROp
        {
            const bool checkVariants_;
            const bool followLink_;
        public:
            lastModifiedHROp
            (
                const bool checkVariants,
                const bool followLink
            )
            :
                checkVariants_(checkVariants),
                followLink_(followLink)
            {}

        public:
            double operator()(const fileName& fName) const
            {
                return
                    Foam::highResLastModified
                    (
                        fName,
                        checkVariants_,
                        followLink_
                    );
            }
        };

        class mvBakOp
        {
            std::string ext_;
        public:
            mvBakOp(const std::string& ext)
            :
                ext_(ext)
            {}

            bool operator()(const fileName& fName) const
            {
                return Foam::mvBak(fName, ext_);
            }
        };

        class rmOp
        {
        public:
            bool operator()(const fileName& fName) const
            {
                return Foam::rm(fName);
            }
        };

        class rmDirOp
        {
        public:
            bool operator()(const fileName& fName) const
            {
                return Foam::rmDir(fName);
            }
        };

        class cpOp
        {
            const bool followLink_;
        public:
            cpOp(const bool followLink)
            :
                followLink_(followLink)
            {}

        public:
            bool operator()(const fileName& src, const fileName& dest) const
            {
                return Foam::cp(src, dest, followLink_);
            }
        };

        class lnOp
        {
        public:
            bool operator()(const fileName& src, const fileName& dest) const
            {
                return Foam::ln(src, dest);
            }
        };

        class mvOp
        {
            const bool followLink_;
        public:
            mvOp(const bool followLink)
            :
                followLink_(followLink)
            {}

        public:
            bool operator()(const fileName& src, const fileName& dest) const
            {
                return Foam::mv(src, dest, followLink_);
            }
        };

        class fileOrNullOp
        {
            const bool isFile_;
        public:
            fileOrNullOp(const bool isFile)
            :
                isFile_(isFile)
            {}

            fileName operator()(const fileName& fName) const
            {
                return
                (
                    (isFile_ && Foam::isFile(fName))
                 || (!isFile_ && Foam::isDir(fName))
                  ? fName
                  : fileName::null
                );
            }
        };

        class readDirOp
        {
            const fileType type_;
            const bool filterVariants_;
            const bool followLink_;

        public:

            readDirOp
            (
                const fileType type,
                const bool filterVariants,
                const bool followLink
            )
            :
                type_(type),
                filterVariants_(filterVariants),
                followLink_(followLink)
            {}

            fileNameList operator()(const fileName& fName) const
            {
                return Foam::readDir
                (
                    fName,
                    type_,
                    filterVariants_,
                    followLink_
                );
            }
        };


    // Private Member Functions

        //- Get the list of processors that are part of this communicator
        static labelList subRanks(const label n);

        template<class Type>
        Type scatterList(const UList<Type>&, const int, const label comm) const;

        template<class Type, class fileOp>
        Type masterOp
        (
            const fileName&,
            const fileOp& fop,
            const int tag,
            const label comm
        ) const;

        template<class Type, class fileOp>
        Type masterOp
        (
            const fileName&,
            const fileName&,
            const fileOp& fop,
            const int tag,
            const label comm
        ) const;

        //- Equivalent of Time::findInstance
        static word findInstancePath
        (
            const instantList& timeDirs,
            const instant& t
        );

        //- Search (locally!) for object; return info on how it was found.
        //  Does not do any parallel communication.
        //    checkGlobal : also check undecomposed case
        //    isFile      : true:check for file  false:check for directory
        // searchType    : how was found
        // processorsDir : name of processor directory
        // instance      : instance
        virtual fileName filePathInfo
        (
            const bool checkGlobal,
            const bool isFile,
            const IOobject&,
            pathType& searchType,
            word& processorsDir,
            word& instance
        ) const;

        //- Construct filePath
        fileName localObjectPath
        (
            const IOobject&,
            const pathType& searchType,
            const word& processorsDir,
            const word& instancePath
        ) const;

        //- Detect file (possibly compressed), read file contents and send
        //  to processors
        static void readAndSend
        (
            const fileName& fName,
            const labelUList& procs,
            PstreamBuffers& pBufs
        );

        //- Read files on comms master
        static autoPtr<ISstream> read
        (
            IOobject& io,
            const label comm,
            const bool uniform,             // on comms master only
            const fileNameList& filePaths,  // on comms master only
            const boolList& read            // on comms master only
        );

        //- Helper: check IO for local existence. Like filePathInfo but
        //  without parent searchign and instance searching
        bool exists(const dirIndexList&, IOobject& io) const;


public:

        //- Runtime type information
        TypeName("masterUncollated");


    // Static data

        //- Max size of parallel communications. Switches from non-blocking
        //  to scheduled when reading/writing files. Read as float to enable
        //  easy specification of large sizes.
        static float maxMasterFileBufferSize;


    // Constructors

        //- Construct null
        masterUncollatedFileOperation(const bool verbose);

        //- Construct from communicator
        masterUncollatedFileOperation(const label comm, const bool verbose);


    //- Destructor
    virtual ~masterUncollatedFileOperation();


    // Member Functions

        // OSSpecific equivalents

            //- Make directory
            virtual bool mkDir(const fileName&, mode_t=0777) const;

            //- Set the file mode
            virtual bool chMod(const fileName&, const mode_t) const;

            //- Return the file mode
            virtual mode_t mode
            (
                const fileName&,
                const bool checkVariants = true,
                const bool followLink = true
            ) const;

            //- Return the file type: directory, file or link
            virtual fileType type
            (
                const fileName&,
                const bool checkVariants = true,
                const bool followLink = true
            ) const;

            //- Does the name exist (as directory or file) in the file system?
            //  Optionally enable/disable check for gzip file.
            virtual bool exists
            (
                const fileName&,
                const bool checkVariants = true,
                const bool followLink = true
            ) const;

            //- Does the name exist as a directory in the file system?
            virtual bool isDir
            (
                const fileName&,
                const bool followLink = true
            ) const;

            //- Does the name exist as a file in the file system?
            //  Optionally enable/disable check for gzip file.
            virtual bool isFile
            (
                const fileName&,
                const bool checkVariants = true,
                const bool followLink = true
            ) const;

            //- Return size of file
            virtual off_t fileSize
            (
                const fileName&,
                const bool checkVariants = true,
                const bool followLink = true
            ) const;

            //- Return time of last file modification
            virtual time_t lastModified
            (
                const fileName&,
                const bool checkVariants = true,
                const bool followLink = true
            ) const;

            //- Return time of last file modification
            virtual double highResLastModified
            (
                const fileName&,
                const bool checkVariants = true,
                const bool followLink = true
            ) const;

            //- Read a directory and return the entries as a string list
            virtual fileNameList readDir
            (
                const fileName&,
                const fileType = fileType::file,
                const bool filterVariants = true,
                const bool followLink = true
            ) const;

            //- Copy, recursively if necessary, the source to the destination
            virtual bool cp
            (
                const fileName& src,
                const fileName& dst,
                const bool followLink = true
            ) const;

            //- Create a softlink. dst should not exist. Returns true if
            //  successful.
            virtual bool ln(const fileName& src, const fileName& dst) const;

            //- Rename src to dst
            virtual bool mv
            (
                const fileName& src,
                const fileName& dst,
                const bool followLink = false
            ) const;

            //- Rename to a corresponding backup file
            //  If the backup file already exists, attempt with
            //  "01" .. "99" suffix
            virtual bool mvBak
            (
                const fileName&,
                const std::string& ext = "bak"
            ) const;

            //- Remove a file, returning true if successful otherwise false
            virtual bool rm(const fileName&) const;

            //- Remove a directory and its contents
            virtual bool rmDir(const fileName&) const;

//            //- Open a shared library. Return handle to library. Print error
//            //  message if library cannot be loaded (check = true)
//            virtual void* dlOpen
//            (
//                const fileName& lib,
//                const bool check = true
//            ) const;


        // (reg)IOobject functionality

            //- Search for an object. checkGlobal : also check undecomposed case
            virtual fileName filePath
            (
                const bool checkGlobal,
                const IOobject&,
                const word& typeName
            ) const;

            //- Search for a directory. checkGlobal : also check undecomposed
            //  case
            virtual fileName dirPath
            (
                const bool checkGlobal,
                const IOobject&
            ) const;

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
            ) const;

            //- Reads header for regIOobject and returns an ISstream
            //  to read the contents.
            virtual autoPtr<ISstream> readStream
            (
                regIOobject&,
                const fileName&,
                const word& typeName,
                const bool read = true
            ) const;

            //- Top-level read
            //  Format set to ASCII if defaultFormat and file format are ASCII
            virtual bool read
            (
                regIOobject&,
                const bool masterOnly,
                const IOstream::streamFormat defaultFormat,
                const word& typeName
            ) const;

            //- Writes a regIOobject (so header, contents and divider).
            //  Returns success state.
            virtual bool writeObject
            (
                const regIOobject&,
                IOstream::streamFormat format=IOstream::ASCII,
                IOstream::versionNumber version=IOstream::currentVersion,
                IOstream::compressionType compression=IOstream::UNCOMPRESSED,
                const bool write = true
            ) const;

            //- Generate an ISstream that reads a file
            virtual autoPtr<ISstream> NewIFstream
            (
                const fileName& filePath,
                IOstream::streamFormat format=IOstream::ASCII,
                IOstream::versionNumber version=IOstream::currentVersion
            ) const;

            //- Generate an Ostream that writes a file
            virtual autoPtr<Ostream> NewOFstream
            (
                const fileName& filePath,
                IOstream::streamFormat format=IOstream::ASCII,
                IOstream::versionNumber version=IOstream::currentVersion,
                IOstream::compressionType compression=IOstream::UNCOMPRESSED,
                const bool write = true
            ) const;


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

            //- Same file?
            static bool uniformFile(const fileNameList&);

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
            virtual void setTime(const Time&) const;

            //- Forcibly wait until all output done. Flush any cached data
            virtual void flush() const;

            //- Return cached times
            const HashPtrTable<instantList>& times() const
            {
                return times_;
            }
};


/*---------------------------------------------------------------------------*\
           Class masterUncollatedFileOperationInitialise Declaration
\*---------------------------------------------------------------------------*/

class masterUncollatedFileOperationInitialise
:
    public unthreadedInitialise
{
public:

    // Constructors

        //- Construct from components
        masterUncollatedFileOperationInitialise(int& argc, char**& argv);


    //- Destructor
    virtual ~masterUncollatedFileOperationInitialise()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileOperations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "masterUncollatedFileOperationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
