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

\*---------------------------------------------------------------------------*/

#include "OFstream.H"
#include "OSspecific.H"
#include "gzstream.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(OFstream, 0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::OFstreamAllocator::OFstreamAllocator
(
    const fileName& filePath,
    IOstream::compressionType compression,
    const bool append
)
:
    ofPtr_(nullptr)
{
    if (filePath.empty())
    {
        if (OFstream::debug)
        {
            InfoInFunction << "Cannot open null file " << endl;
        }
    }
    ofstream::openmode mode(ofstream::out);
    if (append)
    {
        mode |= ofstream::app;
    }

    if (compression == IOstream::COMPRESSED)
    {
        // Get identically named uncompressed version out of the way
        fileType pathType = Foam::type(filePath, false, false);
        if (pathType == fileType::file || pathType == fileType::link)
        {
            rm(filePath);
        }
        fileName gzfilePath(filePath + ".gz");

        if (!append && Foam::type(gzfilePath) == fileType::link)
        {
            // Disallow writing into softlink to avoid any problems with
            // e.g. softlinked initial fields
            rm(gzfilePath);
        }

        ofPtr_ = new ogzstream(gzfilePath.c_str(), mode);
    }
    else
    {
        // get identically named compressed version out of the way
        fileName gzfilePath(filePath + ".gz");
        fileType gzType = Foam::type(gzfilePath, false, false);
        if (gzType == fileType::file || gzType == fileType::link)
        {
            rm(gzfilePath);
        }
        if
        (
            !append
         && Foam::type(filePath, false, false) == fileType::link
        )
        {
            // Disallow writing into softlink to avoid any problems with
            // e.g. softlinked initial fields
            rm(filePath);
        }

        ofPtr_ = new ofstream(filePath.c_str(), mode);
    }
}


Foam::OFstreamAllocator::~OFstreamAllocator()
{
    delete ofPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::OFstream::OFstream
(
    const fileName& filePath,
    streamFormat format,
    versionNumber version,
    compressionType compression,
    const bool append
)
:
    OFstreamAllocator(filePath, compression, append),
    OSstream(*ofPtr_, "OFstream.sinkFile_", format, version, compression),
    filePath_(filePath)
{
    setClosed();
    setState(ofPtr_->rdstate());

    if (!good())
    {
        if (debug)
        {
            InfoInFunction
                << "Could not open file " << filePath
                << "for input\n"
                   "in stream " << info() << Foam::endl;
        }

        setBad();
    }
    else
    {
        setOpened();
    }

    lineNumber_ = 1;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::OFstream::~OFstream()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

std::ostream& Foam::OFstream::stdStream()
{
    if (!ofPtr_)
    {
        FatalErrorInFunction
            << "No stream allocated." << abort(FatalError);
    }
    return *ofPtr_;
}


const std::ostream& Foam::OFstream::stdStream() const
{
    if (!ofPtr_)
    {
        FatalErrorInFunction
            << "No stream allocated." << abort(FatalError);
    }
    return *ofPtr_;
}


void Foam::OFstream::print(Ostream& os) const
{
    os  << "    OFstream: ";
    OSstream::print(os);
}


// ************************************************************************* //
