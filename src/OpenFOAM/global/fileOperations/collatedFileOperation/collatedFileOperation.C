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
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

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

    fileName prefix;
    fileName postfix;
    label proci = splitProcessorPath(io.objectPath(), prefix, postfix);

    if (debug)
    {
        Pout<< "writeObject:" << " : For local object : "
            << io.name()
            << " appending processor " << proci
            << " data to " << pathName << endl;
    }

    if (proci == -1)
    {
        FatalErrorInFunction
            << "Not a valid processor path " << pathName
            << exit(FatalError);
    }


    // Create string from all data to write
    string buf;
    {
        OStringStream os(fmt, ver);
        if (proci == 0)
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

        if (proci == 0)
        {
            IOobject::writeEndDivider(os);
        }

        buf = os.str();
    }


    bool append = (proci > 0);

    // Note: cannot do append + compression. This is a limitation
    // of ogzstream (or rather most compressed formats)

    OFstream os
    (
        pathName,
        IOstream::BINARY,
        ver,
        IOstream::UNCOMPRESSED, // no compression
        append
    );

    if (!os.good())
    {
        FatalIOErrorInFunction(os)
            << "Cannot open for appending"
            << exit(FatalIOError);
    }

    if (proci == 0)
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
    os << nl << "// Processor" << proci << nl << slice << nl;

    return os.good();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::collatedFileOperation::collatedFileOperation
(
    const bool verbose
)
:
    masterUncollatedFileOperation(false),
    writer_(maxThreadFileBufferSize)
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
                << "         Requires thread support enabled in MPI, "
                   "otherwise the simulation" << nl
                << "         may \"hang\".  If thread support cannot be "
                   "enabled, deactivate threading" << nl
                << "         by setting maxThreadFileBufferSize to 0 in "
                   "$FOAM_ETC/controlDict"
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
{}


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
        return masterUncollatedFileOperation::objectPath
        (
            io,
            fileOperation::PROCESSORSOBJECT,
            io.instance()
        );
    }
    else
    {
        return masterUncollatedFileOperation::objectPath
        (
            io,
            fileOperation::OBJECT,
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
            Pout<< "writeObject:"
                << " : For object : " << io.name()
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
        fileName path(processorsPath(io, inst));

        mkDir(path);
        fileName pathName(path/io.name());

        if (io.global())
        {
            if (debug)
            {
                Pout<< "writeObject:" << " : For global object : " << io.name()
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
            // processors/ file
            if (debug)
            {
                Pout<< "writeObject:"
                    << " : For object : " << io.name()
                    << " appending to " << pathName << endl;
            }

            return appendObject(io, pathName, fmt, ver, cmp);
        }
        else
        {
            if (debug)
            {
                Pout<< "writeObject:"
                    << " : For object : " << io.name()
                    << " starting collating output to " << pathName << endl;
            }

            threadedCollatedOFstream os(writer_, pathName, fmt, ver, cmp);

            // If any of these fail, return (leave error handling to Ostream
            // class)
            if (!os.good())
            {
                return false;
            }
            if (Pstream::master() && !io.writeHeader(os))
            {
                return false;
            }
            // Write the data to the Ostream
            if (!io.writeData(os))
            {
                return false;
            }
            if (Pstream::master())
            {
                IOobject::writeEndDivider(os);
            }

            return true;
        }
    }
}


// ************************************************************************* //
