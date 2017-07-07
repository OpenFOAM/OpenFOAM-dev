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

Class
    Foam::fileOperations::masterUncollatedFileOperation

Description
    fileOperations that performs all file operations on the master processor.
    Requires the calls to be parallel synchronised!

\*---------------------------------------------------------------------------*/

#ifndef fileOperations_masterUncollatedFileOperation_H
#define fileOperations_masterUncollatedFileOperation_H

#include "fileOperation.H"
#include "OSspecific.H"
#include "HashPtrTable.H"

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
            const bool followLink_;
        public:
            modeOp(const bool followLink)
            :
                followLink_(followLink)
            {}

            mode_t operator()(const fileName& fName) const
            {
                return Foam::mode(fName, followLink_);
            }
        };

        class typeOp
        {
            const bool followLink_;
        public:
            typeOp(const bool followLink)
            :
                followLink_(followLink)
            {}

            label operator()(const fileName& fName) const
            {
                return Foam::type(fName, followLink_);
            }
        };

        class existsOp
        {
            const bool checkGzip_;
            const bool followLink_;
        public:
            existsOp(const bool checkGzip, const bool followLink)
            :
                checkGzip_(checkGzip),
                followLink_(followLink)
            {}

            bool operator()(const fileName& fName) const
            {
                return Foam::exists(fName, checkGzip_, followLink_);
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
            const bool checkGzip_;
            const bool followLink_;
        public:
            isFileOp(const bool checkGzip, const bool followLink)
            :
                checkGzip_(checkGzip),
                followLink_(followLink)
            {}
        public:
            bool operator()(const fileName& fName) const
            {
                return Foam::isFile(fName, checkGzip_, followLink_);
            }
        };

        class fileSizeOp
        {
            const bool followLink_;
        public:
            fileSizeOp(const bool followLink)
            :
                followLink_(followLink)
            {}

        public:
            off_t operator()(const fileName& fName) const
            {
                return Foam::fileSize(fName, followLink_);
            }
        };

        class lastModifiedOp
        {
            const bool followLink_;
        public:
            lastModifiedOp(const bool followLink)
            :
                followLink_(followLink)
            {}

        public:
            time_t operator()(const fileName& fName) const
            {
                return Foam::lastModified(fName, followLink_);
            }
        };

        class lastModifiedHROp
        {
            const bool followLink_;
        public:
            lastModifiedHROp(const bool followLink)
            :
                followLink_(followLink)
            {}

        public:
            double operator()(const fileName& fName) const
            {
                return Foam::highResLastModified(fName, followLink_);
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
            const fileName::Type type_;
            const bool filtergz_;
            const bool followLink_;
        public:
            readDirOp
            (
                const fileName::Type type,
                const bool filtergz,
                const bool followLink
            )
            :
                type_(type),
                filtergz_(filtergz),
                followLink_(followLink)
            {}

            fileNameList operator()(const fileName& fName) const
            {
                return Foam::readDir(fName, type_, filtergz_, followLink_);
            }
        };


    // Private Member Functions

        template<class Type>
        Type scatterList(const UList<Type>&) const;

        template<class Type, class fileOp>
        Type masterOp(const fileName&, const fileOp& fop) const;

        template<class Type, class fileOp>
        Type masterOp
        (
            const fileName&,
            const fileName&,
            const fileOp& fop
        ) const;

        //- Equivalent of Time::findInstance
        static word findInstancePath
        (
            const instantList& timeDirs,
            const instant& t
        );

        //- Search for object; return info on how it was found
        //    checkGlobal : also check undecomposed case
        //    isFile      : true:check for file  false:check for directory
        fileName filePathInfo
        (
            const bool checkGlobal,
            const bool isFile,
            const IOobject&,
            pathType&,
            word&
        ) const;

        //- Construct filePath
        static fileName objectPath
        (
            const IOobject&,
            const pathType&,
            const word&
        );

        //- Read file contents and send to processors
        static void readAndSend
        (
            const fileName& filePath,
            const IOstream::compressionType cmp,
            const labelUList& procs,
            PstreamBuffers& pBufs
        );


public:

        //- Runtime type information
        TypeName("masterUncollated");


    // Static data

        //- Max size of parallel communications. Switches from non-blocking
        //  to scheduled when reading/writing files. Read as float to enable
        //  easy specificiation of large sizes.
        static float maxMasterFileBufferSize;


    // Constructors

        //- Construct null
        masterUncollatedFileOperation(const bool verbose);


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
                const bool followLink = true
            ) const;

            //- Return the file type: DIRECTORY, FILE or LINK
            virtual fileName::Type type
            (
                const fileName&,
                const bool followLink = true
            ) const;

            //- Does the name exist (as DIRECTORY or FILE) in the file system?
            //  Optionally enable/disable check for gzip file.
            virtual bool exists
            (
                const fileName&,
                const bool checkGzip=true,
                const bool followLink = true
            ) const;

            //- Does the name exist as a DIRECTORY in the file system?
            virtual bool isDir
            (
                const fileName&,
                const bool followLink = true
            ) const;

            //- Does the name exist as a FILE in the file system?
            //  Optionally enable/disable check for gzip file.
            virtual bool isFile
            (
                const fileName&,
                const bool checkGzip=true,
                const bool followLink = true
            ) const;

            //- Return size of file
            virtual off_t fileSize
            (
                const fileName&,
                const bool followLink = true
            ) const;

            //- Return time of last file modification
            virtual time_t lastModified
            (
                const fileName&,
                const bool followLink = true
            ) const;

            //- Return time of last file modification
            virtual double highResLastModified
            (
                const fileName&,
                const bool followLink = true
            ) const;

            //- Read a directory and return the entries as a string list
            virtual fileNameList readDir
            (
                const fileName&,
                const fileName::Type=fileName::FILE,
                const bool filtergz=true,
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

            //- Remove a dirctory and its contents
            virtual bool rmDir(const fileName&) const;

//            //- Open a shared library. Return handle to library. Print error
//            //  message if library cannot be loaded (check = true)
//            virtual void* dlOpen
//            (
//                const fileName& lib,
//                const bool check = true
//            ) const;


        // (reg)IOobject functinality

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
                const bool valid = true
            ) const;

            //- Top-level read
            virtual bool read
            (
                regIOobject&,
                const bool masterOnly,
                const IOstream::streamFormat format,
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
                const bool valid = true
            ) const;

            //- Generate an ISstream that reads a file
            virtual autoPtr<ISstream> NewIFstream(const fileName&) const;

            //- Generate an Ostream that writes a file
            virtual autoPtr<Ostream> NewOFstream
            (
                const fileName& pathname,
                IOstream::streamFormat format=IOstream::ASCII,
                IOstream::versionNumber version=IOstream::currentVersion,
                IOstream::compressionType compression=IOstream::UNCOMPRESSED,
                const bool valid = true
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

            //- Callback for time change
            virtual void setTime(const Time&) const;

            //- root+casename with any 'processorXXX' replaced by 'processsors'
            static fileName processorsCasePath(const IOobject&);

            //- Like io.path with provided instance and any 'processorXXX'
            //  replaced by 'processsors'
            static fileName processorsPath(const IOobject&, const word&);

            //- Operating on fileName: replace processorXXX with processors
            static fileName processorsPath(const fileName&);

            //- Split fileName into part before processor and part after.
            //  Returns -1 or processor number. Use with care.
            //  - path/"processor" + Foam::name(proci)/local reconstructs input
            //  - path/"processors"/local reconstructs processors equivalence
            static label splitProcessorPath
            (
                const fileName&,
                fileName& path,
                fileName& local
            );

            //- Return cached times
            const HashPtrTable<instantList>& times() const
            {
                return times_;
            }
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
