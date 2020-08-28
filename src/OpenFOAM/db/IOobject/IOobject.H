/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::IOobject

Description
    IOobject defines the attributes of an object for which implicit
    objectRegistry management is supported, and provides the infrastructure
    for performing stream I/O.

    An IOobject is constructed with an object name, a class name, an instance
    path, a reference to a objectRegistry, and parameters determining its
    storage status.

    \par Read options

    Define what is done on object construction and explicit reads:
      - \par MUST_READ
        Object must be read from Istream on construction. \n
        Error if Istream does not exist or can't be read.
        Does not check timestamp or re-read.
      - \par MUST_READ_IF_MODIFIED
        Object must be read from Istream on construction. \n
        Error if Istream does not exist or can't be read. If object is
        registered its timestamp will be checked every timestep and possibly
        re-read.
      - \par READ_IF_PRESENT
        Read object from Istream if Istream exists, otherwise don't. \n
        Error only if Istream exists but can't be read.
        Does not check timestamp or re-read.
      - \par NO_READ
        Don't read

    \par Write options

    Define what is done on object destruction and explicit writes:
      - \par AUTO_WRITE
        Object is written automatically when requested to by the
        objectRegistry.
      - \par NO_WRITE
        No automatic write on destruction but can be written explicitly

SourceFiles
    IOobject.C
    IOobjectReadHeader.C
    IOobjectWriteHeader.C
    IOobjectPrint.C

\*---------------------------------------------------------------------------*/

#ifndef IOobject_H
#define IOobject_H

#include "fileName.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "InfoProxy.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Time;
class objectRegistry;

/*---------------------------------------------------------------------------*\
                          Class IOobject Declaration
\*---------------------------------------------------------------------------*/

class IOobject
{

public:

    //- Keyword for the FoamFile header sub-dictionary
    static constexpr const char* foamFile = "FoamFile";

    // Public data types

        //- Enumeration defining the valid states of an IOobject
        enum objectState
        {
            GOOD,
            BAD
        };

        //- Enumeration defining the read options
        enum readOption
        {
            MUST_READ,
            MUST_READ_IF_MODIFIED,
            READ_IF_PRESENT,
            NO_READ
        };

        //- Enumeration defining the write options
        enum writeOption
        {
            AUTO_WRITE = 0,
            NO_WRITE = 1
        };

        //- Enumeration defining the file checking options
        enum fileCheckTypes
        {
            timeStamp,
            timeStampMaster,
            inotify,
            inotifyMaster
        };

        static const NamedEnum<fileCheckTypes, 4> fileCheckTypesNames;

private:

    // Private Data

        //- Name
        word name_;

        //- Class name read from header
        word headerClassName_;

        //- Optional note
        string note_;

        //- Instance path component
        fileName instance_;

        //- Local path component
        fileName local_;

        //- objectRegistry reference
        const objectRegistry& db_;

        //- Read option
        readOption rOpt_;

        //- Write option
        writeOption wOpt_;

        //- Register object created from this IOobject with registry if true
        bool registerObject_;

        //- IOobject state
        objectState objState_;


protected:

    // Protected Member Functions

        //- Set the object state to bad
        void setBad(const string&);


public:

    //- Runtime type information
    TypeName("IOobject");


    // Static Member Functions

        //- Split path into instance, local, name components
        //  input               IOobject(instance, local, name)
        //  -----               ------
        //  "foo"               ("", "", "foo")
        //  "foo/bar"           ("foo", "", "bar")
        //  "/XXX/bar"          ("/XXX", "", "bar")
        //  "foo/bar/"          ERROR - no name
        //  "foo/xxx/bar"       ("foo", "xxx", "bar")
        //  "foo/xxx/yyy/bar"   ("foo", "xxx/yyy", "bar")
        static bool fileNameComponents
        (
            const fileName& path,
            fileName& instance,
            fileName& local,
            word& name
        );

        template<class Name>
        static inline word groupName(Name name, const word& group);

        //- Return group (extension part of name)
        static word group(const word& name);

        //- Return member (name without the extension)
        static word member(const word& name);

        //- Return the name of the object within the given model
        //  as <model>:<name>
        template<class Name>
        static inline word modelName(Name name, const word& model);

        //- Type of file modification checking
        static fileCheckTypes fileModificationChecking;


    // Constructors

        //- Construct from name, instance, registry, io options
        IOobject
        (
            const word& name,
            const fileName& instance,
            const objectRegistry& registry,
            readOption r=NO_READ,
            writeOption w=NO_WRITE,
            bool registerObject=true
        );

        //- Construct from name, instance, local, registry, io options
        IOobject
        (
            const word& name,
            const fileName& instance,
            const fileName& local,
            const objectRegistry& registry,
            readOption r=NO_READ,
            writeOption w=NO_WRITE,
            bool registerObject=true
        );

        //- Construct from path, registry, io options
        //  Uses fileNameComponents() to split path into components.
        IOobject
        (
            const fileName& path,
            const objectRegistry& registry,
            readOption r=NO_READ,
            writeOption w=NO_WRITE,
            bool registerObject=true
        );

        //- Construct from copy resetting registry
        IOobject
        (
            const IOobject& io,
            const objectRegistry& registry
        );

        //- Construct from copy resetting name
        IOobject
        (
            const IOobject& io,
            const word& name
        );

        //- Copy constructor
        IOobject(const IOobject& io) = default;

        //- Clone
        autoPtr<IOobject> clone() const
        {
            return autoPtr<IOobject>(new IOobject(*this));
        }

        //- Clone resetting registry
        autoPtr<IOobject> clone(const objectRegistry& registry) const
        {
            return autoPtr<IOobject>(new IOobject(*this, registry));
        }


    //- Destructor
    virtual ~IOobject();


    // Member Functions

        // General access

            //- Return time
            const Time& time() const;

            //- Return the local objectRegistry
            const objectRegistry& db() const;

            //- Return name
            const word& name() const
            {
                return name_;
            }

            //- Return name of the class name read from header
            const word& headerClassName() const
            {
                return headerClassName_;
            }

            //- Return name of the class name read from header
            word& headerClassName()
            {
                return headerClassName_;
            }

            //- Return non-constant access to the optional note
            string& note()
            {
                return note_;
            }

            //- Return the optional note
            const string& note() const
            {
                return note_;
            }

            //- Rename
            virtual void rename(const word& newName)
            {
                name_ = newName;
            }

            //- Register object created from this IOobject with registry if true
            bool& registerObject()
            {
                return registerObject_;
            }

            //- Register object created from this IOobject with registry if true
            bool registerObject() const
            {
                return registerObject_;
            }


        // Read/write options

            readOption readOpt() const
            {
                return rOpt_;
            }

            readOption& readOpt()
            {
                return rOpt_;
            }

            writeOption writeOpt() const
            {
                return wOpt_;
            }

            writeOption& writeOpt()
            {
                return wOpt_;
            }


        // Path components

            //- Return group (extension part of name)
            word group() const;

            //- Return member (name without the extension)
            word member() const;

            //- Return the name of the object within this model
            //  as <model>:<name>
            inline word modelName(const char* name) const;

            const fileName& rootPath() const;

            const fileName& caseName() const;

            const fileName& instance() const
            {
                return instance_;
            }

            fileName& instance()
            {
                return instance_;
            }

            const fileName& local() const
            {
                return local_;
            }

            //- Return complete path
            fileName path() const;

            //- Return complete path with alternative instance and local
            fileName path
            (
                const word& instance,
                const fileName& local = ""
            ) const;

            //- Return the path relative to the case
            fileName localPath() const;

            //- Return complete path + object name
            fileName objectPath() const
            {
                return path()/name();
            }

            //- Return complete localPath + object name
            fileName localObjectPath() const
            {
                return localPath()/name();
            }

            //- Helper for filePath that searches locally.
            fileName localFilePath(const word& typeName) const;

            //- Helper for filePath that searches up if in parallel
            fileName globalFilePath(const word& typeName) const;


        // Reading

            //- Read header
            bool readHeader(Istream&);

            //- Read header (uses typeFilePath to find file) and check header
            //  info. Optionally checks headerClassName against type
            template<class Type>
            bool typeHeaderOk(const bool checkType = true);

            //- Helper: warn that type does not support re-reading
            template<class Type>
            void warnNoRereading() const;

        // Writing

            //- Write the standard OpenFOAM file/dictionary banner
            //  Optionally without -*- C++ -*- editor hint (eg, for logs)
            template<class Stream>
            static inline Stream& writeBanner(Stream& os, bool noHint=false);

            //- Write the standard file section divider
            template<class Stream>
            static inline Stream& writeDivider(Stream& os);

            //- Write the standard end file divider
            template<class Stream>
            static inline Stream& writeEndDivider(Stream& os);

            //- Write header
            bool writeHeader(Ostream&) const;

            //- Write header. Allow override of type
            bool writeHeader(Ostream&, const word& objectType) const;


        // Error Handling

            bool good() const
            {
                return objState_ == GOOD;
            }

            bool bad() const
            {
                return objState_ == BAD;
            }


        // Info

            //- Return info proxy.
            //  Used to print token information to a stream
            InfoProxy<IOobject> info() const
            {
                return *this;
            }


    // Member Operators

        void operator=(const IOobject&);
};


template<>
Ostream& operator<<(Ostream& os, const InfoProxy<IOobject>& ip);


//- Template function for obtaining global status
template<class T>
inline bool typeGlobal()
{
    return false;
}

//- Template function for obtaining filePath
template<class T>
inline fileName typeFilePath(const IOobject& io)
{
    return
    (
        typeGlobal<T>()
      ? io.globalFilePath(T::typeName)
      : io.localFilePath(T::typeName)
    );
}

inline IOobject unregister(const IOobject& io)
{
    IOobject uio(io);
    uio.registerObject() = false;
    return uio;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IOobjectI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "IOobjectTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
