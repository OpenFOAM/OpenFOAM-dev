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
    Write primitive and binary block from OPstream gamma-mpi

\*---------------------------------------------------------------------------*/

#include "OPstream.H"
#include "long.H"
#include "PstreamGlobals.H"

extern "C" {

#include <linux/gamma/libgamma.h>

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Largest message sent so far. This tracks the size of the receive
// buffer on the receiving end. Done so we only send out resize messages
// if necessary
//! \cond fileScope
labelList maxSendSize;
//! \endcond


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

OPstream::~OPstream()
{
    if (Pstream::debug)
    {
        Pout<< "OPstream::~OPstream() to processor " << toProcNo_
            << Foam::endl;
    }

    if
    (
       !write
        (
            commsType_,
            toProcNo_,
            buf_.begin(),
            bufPosition_
        )
    )
    {
        FatalErrorIn("OPstream::~OPstream()")
            << "GAMMA cannot send outgoing message"
            << Foam::abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool OPstream::write
(
    const commsTypes commsType,
    const int toProcNo,
    const char* buf,
    const std::streamsize bufSize
)
{
    if (PstreamGlobals::getSizeFromHeader(buf, bufSize) != -1)
    {
        FatalErrorIn("OPstream::write")
            << "Problem: Trying to send message of size " << bufSize
            << " that corresponds to the special resizeMessage."
            << Foam::abort(FatalError);
    }

    if (maxSendSize.empty())
    {
        // Intialize maxSendSize to the initial size of the receive buffers.
        maxSendSize.setSize(Pstream::nProcs());
        maxSendSize = PstreamGlobals::initialBufferLen;
        maxSendSize[Pstream::myProcNo()] = 0;

        if (Pstream::debug)
        {
            forAll(maxSendSize, procNo)
            {
                Pout<< "OPstream::write() : for toProcNo:" << procNo
                    << " set maxSendSize to " << maxSendSize[procNo]
                    << Foam::endl;
            }
        }
    }

    if (Pstream::debug)
    {
        Pout<< "OPstream::write() : proc:" << toProcNo
            << " maxSendSize:" << maxSendSize[toProcNo]
            << Foam::endl;
    }

    if (bufSize > maxSendSize[toProcNo])
    {
        // Send resize message.
        if (Pstream::debug)
        {
            Pout<< "OPstream::write() : Sending resize message to proc "
            << toProcNo
            << " for size:" << bufSize
            << Foam::endl;
        }

        PstreamGlobals::setResizeMessage(bufSize);
        gamma_send_flowctl
        (
            toProcNo,
            reinterpret_cast<char*>(PstreamGlobals::resizeMessage),
            PstreamGlobals::resizeMessageLen*sizeof(uint64_t)
        );

        maxSendSize[toProcNo] = bufSize;
    }


    // Do normal send
    // ~~~~~~~~~~~~~~

    // Note: could be put into allocation of buf.
    //gamma_mlock(const_cast<char*>(buf), bufSize);

    if (Pstream::debug)
    {
        Pout<< "OPstream::write() : Sending to proc " << toProcNo
            << " bytes:" << bufSize << Foam::endl;
    }

    gamma_send_flowctl
    (
        toProcNo,
        const_cast<char*>(buf),
        bufSize
    );

    //gamma_munlock(const_cast<char*>(buf), bufSize);

    if (Pstream::debug)
    {
        Pout<< "OPstream::write() : Sent " << bufSize
            << " to proc " << toProcNo
            << Foam::endl;
    }


    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
