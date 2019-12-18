/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Exchange data.

\*---------------------------------------------------------------------------*/

#include "Pstream.H"
#include "contiguous.H"
#include "PstreamCombineReduceOps.H"
#include "UPstream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Container, class T>
void Foam::Pstream::exchange
(
    const UList<Container>& sendBufs,
    const labelUList& recvSizes,
    List<Container>& recvBufs,
    const int tag,
    const label comm,
    const bool block
)
{
    if (!contiguous<T>())
    {
        FatalErrorInFunction
            << "Continuous data only." << sizeof(T) << Foam::abort(FatalError);
    }

    if (sendBufs.size() != UPstream::nProcs(comm))
    {
        FatalErrorInFunction
            << "Size of list " << sendBufs.size()
            << " does not equal the number of processors "
            << UPstream::nProcs(comm)
            << Foam::abort(FatalError);
    }

    recvBufs.setSize(sendBufs.size());

    if (UPstream::parRun() && UPstream::nProcs(comm) > 1)
    {
        label startOfRequests = Pstream::nRequests();

        // Set up receives
        // ~~~~~~~~~~~~~~~

        forAll(recvSizes, proci)
        {
            label nRecv = recvSizes[proci];

            if (proci != Pstream::myProcNo(comm) && nRecv > 0)
            {
                recvBufs[proci].setSize(nRecv);
                UIPstream::read
                (
                    UPstream::commsTypes::nonBlocking,
                    proci,
                    reinterpret_cast<char*>(recvBufs[proci].begin()),
                    nRecv*sizeof(T),
                    tag,
                    comm
                );
            }
        }


        // Set up sends
        // ~~~~~~~~~~~~

        forAll(sendBufs, proci)
        {
            if (proci != Pstream::myProcNo(comm) && sendBufs[proci].size() > 0)
            {
                if
                (
                   !UOPstream::write
                    (
                        UPstream::commsTypes::nonBlocking,
                        proci,
                        reinterpret_cast<const char*>(sendBufs[proci].begin()),
                        sendBufs[proci].size()*sizeof(T),
                        tag,
                        comm
                    )
                )
                {
                    FatalErrorInFunction
                        << "Cannot send outgoing message. "
                        << "to:" << proci << " nBytes:"
                        << label(sendBufs[proci].size()*sizeof(T))
                        << Foam::abort(FatalError);
                }
            }
        }


        // Wait for all to finish
        // ~~~~~~~~~~~~~~~~~~~~~~

        if (block)
        {
            Pstream::waitRequests(startOfRequests);
        }
    }

    // Do myself
    recvBufs[Pstream::myProcNo(comm)] = sendBufs[Pstream::myProcNo(comm)];
}


template<class Container>
void Foam::Pstream::exchangeSizes
(
    const Container& sendBufs,
    labelList& recvSizes,
    const label comm
)
{
    if (sendBufs.size() != UPstream::nProcs(comm))
    {
        FatalErrorInFunction
            << "Size of container " << sendBufs.size()
            << " does not equal the number of processors "
            << UPstream::nProcs(comm)
            << Foam::abort(FatalError);
    }

    labelList sendSizes(sendBufs.size());
    forAll(sendBufs, proci)
    {
        sendSizes[proci] = sendBufs[proci].size();
    }
    recvSizes.setSize(sendSizes.size());
    allToAll(sendSizes, recvSizes, comm);
}


template<class Container, class T>
void Foam::Pstream::exchange
(
    const UList<Container>& sendBufs,
    List<Container>& recvBufs,
    const int tag,
    const label comm,
    const bool block
)
{
    labelList recvSizes;
    exchangeSizes(sendBufs, recvSizes, comm);

    exchange<Container, T>(sendBufs, recvSizes, recvBufs, tag, comm, block);
}


// ************************************************************************* //
