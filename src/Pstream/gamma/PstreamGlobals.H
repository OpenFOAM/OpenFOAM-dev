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

Namespace
    Foam::PstreamGlobals

Description
    Global functions and variables for working with parallel streams,
    but principally for gamma/mpi

SourceFiles
    PstreamGlobals.C

\*---------------------------------------------------------------------------*/

#ifndef PstreamGlobals_H
#define PstreamGlobals_H

#include "FixedList.H"
#include "labelList.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class PstreamGlobals Declaration
\*---------------------------------------------------------------------------*/

namespace PstreamGlobals
{

//- Block wait for message on port procNo
void gammaWait(const label procNo);

//- Helper routine to copy out newly received data
label copyReceive
(
    const label procNo,
    char* buf,
    const label bufSize
);


//- Receive buffers
extern FixedList<List<List<char> >, 4> recvBuf;

//- Length of receive buffers
extern FixedList<labelList, 4> recvBufLen;

//- Currently active buffer in receiving
extern labelList recvIndex;
//- Receive buffer that has to be consumed
extern labelList consumeIndex;


//- Special message to signal resizing
const int resizeMessageLen = 6;
extern uint64_t resizeMessage[];
//- Initial buffer length. Should be able to contain the message comfortably.
const int initialBufferLen = 2*resizeMessageLen*sizeof(uint64_t);

//- Helper routine to check if a message is a resize message.
//  Returns -1 if not or the new size.
label getSizeFromHeader(const char* buf, const label len);
//- Change the resize message to contain the new length
void setResizeMessage(const label len);

//- Get max size of all receive buffers to procNo
label getMaxBufSize(const int procNo);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
