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

\*---------------------------------------------------------------------------*/

#include "PstreamGlobals.H"
#include "IOstreams.H"
#include "Pstream.H"

extern "C" {

#include <linux/gamma/libgamma.h>

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Receive buffers
FixedList<List<List<char> >, 4> PstreamGlobals::recvBuf;

// Length of receive buffers
FixedList<labelList, 4> PstreamGlobals::recvBufLen;

labelList PstreamGlobals::recvIndex;
labelList PstreamGlobals::consumeIndex;

// These are all signalling nans and probably different from the ones that
// the fpu might ever generate.
uint64_t PstreamGlobals::resizeMessage[PstreamGlobals::resizeMessageLen] =
{
    0x7ff7ffffffffffABllu,
    0x7ff7ffffffffffCDllu,
    0x7ff7ffffffffff12llu,
    0x7ff7ffffffffff30llu,
    0x7ff7ffffffffff19llu,
    0x0000000000000000llu       // this word gets overwritten with the length.
};


// Wrapper around gamma_wait
void PstreamGlobals::gammaWait(const label procNo)
{
    // Last request. Block.
    gamma_wait(procNo, 1);

    // Currently unconsumed received message
    label ready = PstreamGlobals::consumeIndex[procNo];

    // Check received length
    if (PstreamGlobals::recvBufLen[ready][procNo] == -2)
    {
        // Was resize message. Consume and rewait (is always followed by
        // real message)

        if (Pstream::debug)
        {
            Pout<< "PstreamGlobals::gammaWait : "
                << "Resize event. consumeIndex:" << ready
                << " Restarting receive from " << procNo << endl;
        }
        // Consume resize message
        PstreamGlobals::recvBufLen[ready][procNo] = -1;
        PstreamGlobals::consumeIndex[procNo] =
            PstreamGlobals::recvBuf.fcIndex(ready);
        // And rewait
        gamma_wait(procNo, 1);
    }
}


// Copies data from global receive buffer into buf.
label PstreamGlobals::copyReceive
(
    const label procNo,
    char* buf,
    const label bufSize
)
{
    // Get the ready buffer
    label ready = consumeIndex[procNo];

    // Actually received
    label receivedLen = recvBufLen[ready][procNo];

    if (Pstream::debug)
    {
        Pout<< "copyReceive : for proc " << procNo
            << " copying " << receivedLen << " bytes out of buffer " << ready
            << endl;
    }

    if (receivedLen < 0)
    {
        FatalErrorIn
        (
            "Pstream::copyReceive(const label, char*, const label)"
        )   << "Illegal message length "
            << receivedLen
            << " received from proc " << procNo << " into buffer " << ready
            << endl
            << "This is probably caused by receiving more than is actually"
            << " sent (e.g. gather without scatter)." << endl
            << abort(FatalError);
    }

    if (receivedLen > bufSize)
    {
        FatalErrorIn
        (
            "Pstream::copyReceive(const label, char*, const label)"
        )   << "buffer ("
            << bufSize
            << ") not large enough for incomming message ("
            << receivedLen << ')'
            << " received from proc " << procNo << " into buffer " << ready
            << abort(FatalError);
    }

    // Copy out of receive buffer
    memcpy
    (
        buf,
        recvBuf[ready][procNo].begin(),
        receivedLen
    );
    // Release receive buffer
    recvBufLen[ready][procNo] = -1;
    // Go to next buffer to consume
    consumeIndex[procNo] = recvBuf.fcIndex(ready);

    return receivedLen;
}


// Checks whether an incoming message is a resize message. If not returns -1,
// otherwise returns size read from header.
label PstreamGlobals::getSizeFromHeader(const char* buf, const label len)
{
    if (len != resizeMessageLen*sizeof(uint64_t))
    {
        return -1;
    }

    const uint64_t* dPtr = reinterpret_cast<const uint64_t*>(buf);

    // Check all but the last word
    for (label i = 0; i < resizeMessageLen-1; i++)
    {
        if (*dPtr++ != resizeMessage[i])
        {
            return -1;
        }
    }

    return *reinterpret_cast<const label*>(dPtr);
}


void PstreamGlobals::setResizeMessage(const label len)
{
    reinterpret_cast<label&>(resizeMessage[resizeMessageLen-1]) = len;
}


label PstreamGlobals::getMaxBufSize(const int procNo)
{
    label maxSz = 0;

    forAll(recvBuf, i)
    {
        maxSz = max(maxSz, recvBuf[i][procNo].size());
    }
    return maxSz;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
