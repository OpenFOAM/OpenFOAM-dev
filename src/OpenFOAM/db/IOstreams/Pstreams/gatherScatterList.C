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

Description
    Gather data from all processors onto single processor according to some
    communication schedule (usually linear-to-master or tree-to-master).
    The gathered data will be a list with element procID the data from processor
    procID. Before calling every processor should insert its value into
    Values[UPstream::myProcNo(comm)].
    Note: after gather every processor only knows its own data and that of the
    processors below it. Only the 'master' of the communication schedule holds
    a fully filled List. Use scatter to distribute the data.

\*---------------------------------------------------------------------------*/

#include "IPstream.H"
#include "OPstream.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Pstream::gatherList
(
    const List<UPstream::commsStruct>& comms,
    List<T>& Values,
    const int tag,
    const label comm
)
{
    if (UPstream::nProcs(comm) > 1)
    {
        if (Values.size() != UPstream::nProcs(comm))
        {
            FatalErrorIn
            (
                "UPstream::gatherList(const List<UPstream::commsStruct>&"
                ", List<T>)"
            )   << "Size of list:" << Values.size()
                << " does not equal the number of processors:"
                << UPstream::nProcs(comm)
                << Foam::abort(FatalError);
        }

        // Get my communication order
        const commsStruct& myComm = comms[UPstream::myProcNo(comm)];

        // Receive from my downstairs neighbours
        forAll(myComm.below(), belowI)
        {
            label belowID = myComm.below()[belowI];
            const labelList& belowLeaves = comms[belowID].allBelow();

            if (contiguous<T>())
            {
                List<T> receivedValues(belowLeaves.size() + 1);

                UIPstream::read
                (
                    UPstream::scheduled,
                    belowID,
                    reinterpret_cast<char*>(receivedValues.begin()),
                    receivedValues.byteSize(),
                    tag,
                    comm
                );

                Values[belowID] = receivedValues[0];

                forAll(belowLeaves, leafI)
                {
                    Values[belowLeaves[leafI]] = receivedValues[leafI + 1];
                }
            }
            else
            {
                IPstream fromBelow(UPstream::scheduled, belowID, 0, tag, comm);
                fromBelow >> Values[belowID];

                if (debug & 2)
                {
                    Pout<< " received through "
                        << belowID << " data from:" << belowID
                        << " data:" << Values[belowID] << endl;
                }

                // Receive from all other processors below belowID
                forAll(belowLeaves, leafI)
                {
                    label leafID = belowLeaves[leafI];
                    fromBelow >> Values[leafID];

                    if (debug & 2)
                    {
                        Pout<< " received through "
                            << belowID << " data from:" << leafID
                            << " data:" << Values[leafID] << endl;
                    }
                }
            }
        }

        // Send up from Values:
        // - my own value first
        // - all belowLeaves next
        if (myComm.above() != -1)
        {
            const labelList& belowLeaves = myComm.allBelow();

            if (debug & 2)
            {
                Pout<< " sending to " << myComm.above()
                    << " data from me:" << UPstream::myProcNo(comm)
                    << " data:" << Values[UPstream::myProcNo(comm)] << endl;
            }

            if (contiguous<T>())
            {
                List<T> sendingValues(belowLeaves.size() + 1);
                sendingValues[0] = Values[UPstream::myProcNo(comm)];

                forAll(belowLeaves, leafI)
                {
                    sendingValues[leafI + 1] = Values[belowLeaves[leafI]];
                }

                OPstream::write
                (
                    UPstream::scheduled,
                    myComm.above(),
                    reinterpret_cast<const char*>(sendingValues.begin()),
                    sendingValues.byteSize(),
                    tag,
                    comm
                );
            }
            else
            {
                OPstream toAbove
                (
                    UPstream::scheduled,
                    myComm.above(),
                    0,
                    tag,
                    comm
                );
                toAbove << Values[UPstream::myProcNo(comm)];

                forAll(belowLeaves, leafI)
                {
                    label leafID = belowLeaves[leafI];

                    if (debug & 2)
                    {
                        Pout<< " sending to "
                            << myComm.above() << " data from:" << leafID
                            << " data:" << Values[leafID] << endl;
                    }
                    toAbove << Values[leafID];
                }
            }
        }
    }
}


template <class T>
void Pstream::gatherList(List<T>& Values, const int tag, const label comm)
{
    if (UPstream::nProcs(comm) < UPstream::nProcsSimpleSum)
    {
        gatherList(UPstream::linearCommunication(comm), Values, tag, comm);
    }
    else
    {
        gatherList(UPstream::treeCommunication(comm), Values, tag, comm);
    }
}


template<class T>
void Pstream::scatterList
(
    const List<UPstream::commsStruct>& comms,
    List<T>& Values,
    const int tag,
    const label comm
)
{
    if (UPstream::nProcs(comm) > 1)
    {
        if (Values.size() != UPstream::nProcs(comm))
        {
            FatalErrorIn
            (
                "UPstream::scatterList(const List<UPstream::commsStruct>&"
                ", List<T>)"
            )   << "Size of list:" << Values.size()
                << " does not equal the number of processors:"
                << UPstream::nProcs(comm)
                << Foam::abort(FatalError);
        }

        // Get my communication order
        const commsStruct& myComm = comms[UPstream::myProcNo(comm)];

        // Reveive from up
        if (myComm.above() != -1)
        {
            const labelList& notBelowLeaves = myComm.allNotBelow();

            if (contiguous<T>())
            {
                List<T> receivedValues(notBelowLeaves.size());

                UIPstream::read
                (
                    UPstream::scheduled,
                    myComm.above(),
                    reinterpret_cast<char*>(receivedValues.begin()),
                    receivedValues.byteSize(),
                    tag,
                    comm
                );

                forAll(notBelowLeaves, leafI)
                {
                    Values[notBelowLeaves[leafI]] = receivedValues[leafI];
                }
            }
            else
            {
                IPstream fromAbove
                (
                    UPstream::scheduled,
                    myComm.above(),
                    0,
                    tag,
                    comm
                );

                forAll(notBelowLeaves, leafI)
                {
                    label leafID = notBelowLeaves[leafI];
                    fromAbove >> Values[leafID];

                    if (debug)
                    {
                        Pout<< " received through "
                            << myComm.above() << " data for:" << leafID
                            << " data:" << Values[leafID] << endl;
                    }
                }
            }
        }

        // Send to my downstairs neighbours
        forAll(myComm.below(), belowI)
        {
            label belowID = myComm.below()[belowI];
            const labelList& notBelowLeaves = comms[belowID].allNotBelow();

            if (contiguous<T>())
            {
                List<T> sendingValues(notBelowLeaves.size());

                forAll(notBelowLeaves, leafI)
                {
                    sendingValues[leafI] = Values[notBelowLeaves[leafI]];
                }

                OPstream::write
                (
                    UPstream::scheduled,
                    belowID,
                    reinterpret_cast<const char*>(sendingValues.begin()),
                    sendingValues.byteSize(),
                    tag,
                    comm
                );
            }
            else
            {
                OPstream toBelow(UPstream::scheduled, belowID, 0, tag, comm);

                // Send data destined for all other processors below belowID
                forAll(notBelowLeaves, leafI)
                {
                    label leafID = notBelowLeaves[leafI];
                    toBelow << Values[leafID];

                    if (debug)
                    {
                        Pout<< " sent through "
                            << belowID << " data for:" << leafID
                            << " data:" << Values[leafID] << endl;
                    }
                }
            }
        }
    }
}


template <class T>
void Pstream::scatterList(List<T>& Values, const int tag, const label comm)
{
    if (UPstream::nProcs(comm) < UPstream::nProcsSimpleSum)
    {
        scatterList(UPstream::linearCommunication(comm), Values, tag, comm);
    }
    else
    {
        scatterList(UPstream::treeCommunication(comm), Values, tag, comm);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
