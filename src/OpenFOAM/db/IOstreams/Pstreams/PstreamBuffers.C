/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "PstreamBuffers.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{

    DynamicList<char> PstreamBuffers::nullBuf(0);
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::PstreamBuffers::PstreamBuffers
(
    const UPstream::commsTypes commsType,
    const int tag,
    const label comm,
    IOstream::streamFormat format,
    IOstream::versionNumber version
)
:
    commsType_(commsType),
    tag_(tag),
    comm_(comm),
    format_(format),
    version_(version),
    sendBuf_(UPstream::nProcs(comm)),
    recvBuf_(UPstream::nProcs(comm)),
    recvBufPos_(UPstream::nProcs(comm),  0),
    finishedSendsCalled_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PstreamBuffers::~PstreamBuffers()
{
    // Check that all data has been consumed.
    forAll(recvBufPos_, proci)
    {
        if (recvBufPos_[proci] < recvBuf_[proci].size())
        {
            FatalErrorInFunction
                << "Message from processor " << proci
                << " not fully consumed. messageSize:" << recvBuf_[proci].size()
                << " bytes of which only " << recvBufPos_[proci]
                << " consumed."
                << Foam::abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::PstreamBuffers::finishedSends(const bool block)
{
    finishedSendsCalled_ = true;

    if (commsType_ == UPstream::commsTypes::nonBlocking)
    {
        Pstream::exchange<DynamicList<char>, char>
        (
            sendBuf_,
            recvBuf_,
            tag_,
            comm_,
            block
        );
    }
}


void Foam::PstreamBuffers::finishedSends(labelList& recvSizes, const bool block)
{
    finishedSendsCalled_ = true;

    if (commsType_ == UPstream::commsTypes::nonBlocking)
    {
        Pstream::exchangeSizes(sendBuf_, recvSizes, comm_);

        Pstream::exchange<DynamicList<char>, char>
        (
            sendBuf_,
            recvSizes,
            recvBuf_,
            tag_,
            comm_,
            block
        );
    }
    else
    {
        FatalErrorInFunction
            << "Obtaining sizes not supported in "
            << UPstream::commsTypeNames[commsType_] << endl
            << " since transfers already in progress. Use non-blocking instead."
            << exit(FatalError);

        // Note: maybe possible only if using different tag from write started
        // by ~UOPstream. Needs some work.
    }
}


void Foam::PstreamBuffers::clear()
{
    forAll(sendBuf_, i)
    {
        sendBuf_[i].clear();
    }
    forAll(recvBuf_, i)
    {
        recvBuf_[i].clear();
    }
    recvBufPos_ = 0;
    finishedSendsCalled_ = false;
}


// ************************************************************************* //
