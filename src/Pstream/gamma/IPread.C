/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Description
    Read token and binary block from IPstream

\*---------------------------------------------------------------------------*/

#include "IPstream.H"
#include "long.H"
#include "PstreamGlobals.H"

extern "C"
{
#   include <linux/gamma/libgamma.h>
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

IPstream::IPstream
(
    const commsTypes commsType,
    const int fromProcNo,
    const label bufSize,
    streamFormat format,
    versionNumber version
)
:
    Pstream(commsType, bufSize),
    Istream(format, version),
    fromProcNo_(fromProcNo),
    messageSize_(0)
{
    // Blocking read.

    setOpened();
    setGood();

    if (Pstream::debug)
    {
        Pout<< "IPstream::IPstream : Starting receive from " << fromProcNo_
            << " recvIndex:" << PstreamGlobals::recvIndex[fromProcNo_]
            << Foam::endl;
    }

    PstreamGlobals::gammaWait(fromProcNo_);

    label ready = PstreamGlobals::consumeIndex[fromProcNo_];
    messageSize_ = PstreamGlobals::recvBufLen[ready][fromProcNo_];

    if (!bufSize)
    {
        if (Pstream::debug)
        {
            Pout<< "IPstream::IPstream : sizing buffer to " << messageSize_
                << endl;
        }

        buf_.setSize(messageSize_);
    }

    PstreamGlobals::copyReceive(fromProcNo_, buf_.begin(), buf_.size());

    if (Pstream::debug)
    {
        Pout<< "IPstream::IPstream : Received " << messageSize_
            << " from " << fromProcNo_
            << Foam::endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label IPstream::read
(
    const commsTypes commsType,
    const int fromProcNo,
    char* buf,
    const std::streamsize bufSize
)
{
    // Blocking read.
    label messageSize;

    if (Pstream::debug)
    {
        Pout<< "IPstream::read : Starting receive from " << fromProcNo
            << " recvIndex:" << PstreamGlobals::recvIndex[fromProcNo]
            << Foam::endl;
    }

    PstreamGlobals::gammaWait(fromProcNo);
    messageSize = PstreamGlobals::copyReceive(fromProcNo, buf, bufSize);

    if (Pstream::debug)
    {
        Pout<< "IPstream::read : Received " << messageSize
            << " from " << fromProcNo
            << Foam::endl;
    }

    return messageSize;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
