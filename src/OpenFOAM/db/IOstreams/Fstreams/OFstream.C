/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    const fileName& pathname,
    IOstream::compressionType compression
)
:
    ofPtr_(NULL)
{
    if (pathname.empty())
    {
        if (OFstream::debug)
        {
            Info<< "OFstreamAllocator::OFstreamAllocator(const fileName&) : "
                   "cannot open null file " << endl;
        }
    }

    if (compression == IOstream::COMPRESSED)
    {
        // get identically named uncompressed version out of the way
        if (isFile(pathname, false))
        {
            rm(pathname);
        }

        ofPtr_ = new ogzstream((pathname + ".gz").c_str());
    }
    else
    {
        // get identically named compressed version out of the way
        if (isFile(pathname + ".gz", false))
        {
            rm(pathname + ".gz");
        }

        ofPtr_ = new ofstream(pathname.c_str());
    }
}


Foam::OFstreamAllocator::~OFstreamAllocator()
{
    delete ofPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::OFstream::OFstream
(
    const fileName& pathname,
    streamFormat format,
    versionNumber version,
    compressionType compression
)
:
    OFstreamAllocator(pathname, compression),
    OSstream(*ofPtr_, "OFstream.sinkFile_", format, version, compression),
    pathname_(pathname)
{
    setClosed();
    setState(ofPtr_->rdstate());

    if (!good())
    {
        if (debug)
        {
            Info<< "OFstream::OFstream(const fileName&,"
                   "streamFormat, versionNumber, compressionType) : "
                   "could not open file " << pathname
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
        FatalErrorIn("OFstream::stdStream()")
            << "No stream allocated." << abort(FatalError);
    }
    return *ofPtr_;
}


const std::ostream& Foam::OFstream::stdStream() const
{
    if (!ofPtr_)
    {
        FatalErrorIn("OFstreamAllocator::stdStream() const")
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
