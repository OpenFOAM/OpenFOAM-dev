/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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

#include "mapDistributeBase.H"
#include "commSchedule.H"
#include "HashSet.H"
#include "globalIndex.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mapDistributeBase, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::List<Foam::labelPair> Foam::mapDistributeBase::schedule
(
    const labelListList& subMap,
    const labelListList& constructMap,
    const int tag
)
{
    // Communications: send and receive processor
    List<labelPair> allComms;

    {
        HashSet<labelPair, labelPair::Hash<>> commsSet(Pstream::nProcs());

        // Find what communication is required
        forAll(subMap, proci)
        {
            if (proci != Pstream::myProcNo())
            {
                if (subMap[proci].size())
                {
                    // I need to send to proci
                    commsSet.insert(labelPair(Pstream::myProcNo(), proci));
                }
                if (constructMap[proci].size())
                {
                    // I need to receive from proci
                    commsSet.insert(labelPair(proci, Pstream::myProcNo()));
                }
            }
        }
        allComms = commsSet.toc();
    }


    // Reduce
    if (Pstream::master())
    {
        // Receive and merge
        for
        (
            int slave=Pstream::firstSlave();
            slave<=Pstream::lastSlave();
            slave++
        )
        {
            IPstream fromSlave(Pstream::commsTypes::scheduled, slave, 0, tag);
            List<labelPair> nbrData(fromSlave);

            forAll(nbrData, i)
            {
                if (findIndex(allComms, nbrData[i]) == -1)
                {
                    label sz = allComms.size();
                    allComms.setSize(sz+1);
                    allComms[sz] = nbrData[i];
                }
            }
        }
        // Send back
        for
        (
            int slave=Pstream::firstSlave();
            slave<=Pstream::lastSlave();
            slave++
        )
        {
            OPstream toSlave(Pstream::commsTypes::scheduled, slave, 0, tag);
            toSlave << allComms;
        }
    }
    else
    {
        {
            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo(),
                0,
                tag
            );
            toMaster << allComms;
        }
        {
            IPstream fromMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo(),
                0,
                tag
            );
            fromMaster >> allComms;
        }
    }


    // Determine my schedule.
    labelList mySchedule
    (
        commSchedule
        (
            Pstream::nProcs(),
            allComms
        ).procSchedule()[Pstream::myProcNo()]
    );

    // Processors involved in my schedule
    return List<labelPair>(UIndirectList<labelPair>(allComms, mySchedule));


    // if (debug)
    //{
    //    Pout<< "I need to:" << endl;
    //    const List<labelPair>& comms = schedule();
    //    forAll(comms, i)
    //    {
    //        const labelPair& twoProcs = comms[i];
    //        label sendProc = twoProcs[0];
    //        label recvProc = twoProcs[1];
    //
    //        if (recvProc == Pstream::myProcNo())
    //        {
    //            Pout<< "    receive from " << sendProc << endl;
    //        }
    //        else
    //        {
    //            Pout<< "    send to " << recvProc << endl;
    //        }
    //    }
    //}
}


const Foam::List<Foam::labelPair>& Foam::mapDistributeBase::schedule() const
{
    if (schedulePtr_.empty())
    {
        schedulePtr_.reset
        (
            new List<labelPair>
            (
                schedule(subMap_, constructMap_, Pstream::msgType())
            )
        );
    }
    return schedulePtr_();
}


void Foam::mapDistributeBase::checkReceivedSize
(
    const label proci,
    const label expectedSize,
    const label receivedSize
)
{
    if (receivedSize != expectedSize)
    {
        FatalErrorInFunction
            << "Expected from processor " << proci
            << " " << expectedSize << " but received "
            << receivedSize << " elements."
            << abort(FatalError);
    }
}


void Foam::mapDistributeBase::printLayout(Ostream& os) const
{
    // Determine offsets of remote data.
    labelList minIndex(Pstream::nProcs(), labelMax);
    labelList maxIndex(Pstream::nProcs(), labelMin);
    forAll(constructMap_, proci)
    {
        const labelList& construct = constructMap_[proci];
        if (constructHasFlip_)
        {
            forAll(construct, i)
            {
                label index = mag(construct[i])-1;
                minIndex[proci] = min(minIndex[proci], index);
                maxIndex[proci] = max(maxIndex[proci], index);
            }
        }
        else
        {
            forAll(construct, i)
            {
                label index = construct[i];
                minIndex[proci] = min(minIndex[proci], index);
                maxIndex[proci] = max(maxIndex[proci], index);
            }
        }
    }

    label localSize;
    if (maxIndex[Pstream::myProcNo()] == labelMin)
    {
        localSize = 0;
    }
    else
    {
        localSize = maxIndex[Pstream::myProcNo()]+1;
    }

    os  << "Layout: (constructSize:" << constructSize_
        << " subHasFlip:" << subHasFlip_
        << " constructHasFlip:" << constructHasFlip_
        << ")" << endl
        << "local (processor " << Pstream::myProcNo() << "):" << endl
        << "    start : 0" << endl
        << "    size  : " << localSize << endl;

    label offset = localSize;
    forAll(minIndex, proci)
    {
        if (proci != Pstream::myProcNo())
        {
            if (constructMap_[proci].size() > 0)
            {
                if (minIndex[proci] != offset)
                {
                    FatalErrorInFunction
                        << "offset:" << offset
                        << " proci:" << proci
                        << " minIndex:" << minIndex[proci]
                        << abort(FatalError);
                }

                label size = maxIndex[proci]-minIndex[proci]+1;
                os  << "processor " << proci << ':' << endl
                    << "    start : " << offset << endl
                    << "    size  : " << size << endl;

                offset += size;
            }
        }
    }
}


void Foam::mapDistributeBase::calcCompactAddressing
(
    const globalIndex& globalNumbering,
    const labelList& elements,
    List<Map<label>>& compactMap
) const
{
    compactMap.setSize(Pstream::nProcs());

    // Count all (non-local) elements needed. Just for presizing map.
    labelList nNonLocal(Pstream::nProcs(), 0);

    forAll(elements, i)
    {
        label globalIndex = elements[i];

        if (globalIndex != -1 && !globalNumbering.isLocal(globalIndex))
        {
            label proci = globalNumbering.whichProcID(globalIndex);
            nNonLocal[proci]++;
        }
    }

    forAll(compactMap, proci)
    {
        compactMap[proci].clear();
        if (proci != Pstream::myProcNo())
        {
            compactMap[proci].resize(2*nNonLocal[proci]);
        }
    }


    // Collect all (non-local) elements needed.
    forAll(elements, i)
    {
        label globalIndex = elements[i];

        if (globalIndex != -1 && !globalNumbering.isLocal(globalIndex))
        {
            label proci = globalNumbering.whichProcID(globalIndex);
            label index = globalNumbering.toLocal(proci, globalIndex);
            label nCompact = compactMap[proci].size();
            compactMap[proci].insert(index, nCompact);
        }
    }
}


void Foam::mapDistributeBase::calcCompactAddressing
(
    const globalIndex& globalNumbering,
    const labelListList& cellCells,
    List<Map<label>>& compactMap
) const
{
    compactMap.setSize(Pstream::nProcs());

    // Count all (non-local) elements needed. Just for presizing map.
    labelList nNonLocal(Pstream::nProcs(), 0);

    forAll(cellCells, cellI)
    {
        const labelList& cCells = cellCells[cellI];

        forAll(cCells, i)
        {
            label globalIndex = cCells[i];

            if (globalIndex != -1 && !globalNumbering.isLocal(globalIndex))
            {
                label proci = globalNumbering.whichProcID(globalIndex);
                nNonLocal[proci]++;
            }
        }
    }

    forAll(compactMap, proci)
    {
        compactMap[proci].clear();
        if (proci != Pstream::myProcNo())
        {
            compactMap[proci].resize(2*nNonLocal[proci]);
        }
    }


    // Collect all (non-local) elements needed.
    forAll(cellCells, cellI)
    {
        const labelList& cCells = cellCells[cellI];

        forAll(cCells, i)
        {
            label globalIndex = cCells[i];

            if (globalIndex != -1 && !globalNumbering.isLocal(globalIndex))
            {
                label proci = globalNumbering.whichProcID(globalIndex);
                label index = globalNumbering.toLocal(proci, globalIndex);
                label nCompact = compactMap[proci].size();
                compactMap[proci].insert(index, nCompact);
            }
        }
    }
}


void Foam::mapDistributeBase::exchangeAddressing
(
    const int tag,
    const globalIndex& globalNumbering,
    labelList& elements,
    List<Map<label>>& compactMap,
    labelList& compactStart
)
{
    // The overall compact addressing is
    // - myProcNo data first (uncompacted)
    // - all other processors consecutively

    compactStart.setSize(Pstream::nProcs());
    compactStart[Pstream::myProcNo()] = 0;
    constructSize_ = globalNumbering.localSize();
    forAll(compactStart, proci)
    {
        if (proci != Pstream::myProcNo())
        {
            compactStart[proci] = constructSize_;
            constructSize_ += compactMap[proci].size();
        }
    }



    // Find out what to receive/send in compact addressing.

    // What I want to receive is what others have to send
    labelListList wantedRemoteElements(Pstream::nProcs());
    // Compact addressing for received data
    constructMap_.setSize(Pstream::nProcs());
    forAll(compactMap, proci)
    {
        if (proci == Pstream::myProcNo())
        {
            // All my own elements are used
            label nLocal = globalNumbering.localSize();
            wantedRemoteElements[proci] = identity(nLocal);
            constructMap_[proci] = identity(nLocal);
        }
        else
        {
            // Remote elements wanted from processor proci
            labelList& remoteElem = wantedRemoteElements[proci];
            labelList& localElem = constructMap_[proci];
            remoteElem.setSize(compactMap[proci].size());
            localElem.setSize(compactMap[proci].size());
            label i = 0;
            forAllIter(Map<label>, compactMap[proci], iter)
            {
                const label compactI = compactStart[proci] + iter();
                remoteElem[i] = iter.key();
                localElem[i]  = compactI;
                iter() = compactI;
                i++;
            }
        }
    }

    subMap_.setSize(Pstream::nProcs());
    Pstream::exchange<labelList, label>
    (
        wantedRemoteElements,
        subMap_,
        tag,
        Pstream::worldComm  // TBD
    );

    // Renumber elements
    forAll(elements, i)
    {
        elements[i] = renumber(globalNumbering, compactMap, elements[i]);
    }
}


void Foam::mapDistributeBase::exchangeAddressing
(
    const int tag,
    const globalIndex& globalNumbering,
    labelListList& cellCells,
    List<Map<label>>& compactMap,
    labelList& compactStart
)
{
    // The overall compact addressing is
    // - myProcNo data first (uncompacted)
    // - all other processors consecutively

    compactStart.setSize(Pstream::nProcs());
    compactStart[Pstream::myProcNo()] = 0;
    constructSize_ = globalNumbering.localSize();
    forAll(compactStart, proci)
    {
        if (proci != Pstream::myProcNo())
        {
            compactStart[proci] = constructSize_;
            constructSize_ += compactMap[proci].size();
        }
    }



    // Find out what to receive/send in compact addressing.

    // What I want to receive is what others have to send
    labelListList wantedRemoteElements(Pstream::nProcs());
    // Compact addressing for received data
    constructMap_.setSize(Pstream::nProcs());
    forAll(compactMap, proci)
    {
        if (proci == Pstream::myProcNo())
        {
            // All my own elements are used
            label nLocal = globalNumbering.localSize();
            wantedRemoteElements[proci] = identity(nLocal);
            constructMap_[proci] = identity(nLocal);
        }
        else
        {
            // Remote elements wanted from processor proci
            labelList& remoteElem = wantedRemoteElements[proci];
            labelList& localElem = constructMap_[proci];
            remoteElem.setSize(compactMap[proci].size());
            localElem.setSize(compactMap[proci].size());
            label i = 0;
            forAllIter(Map<label>, compactMap[proci], iter)
            {
                const label compactI = compactStart[proci] + iter();
                remoteElem[i] = iter.key();
                localElem[i]  = compactI;
                iter() = compactI;
                i++;
            }
        }
    }

    subMap_.setSize(Pstream::nProcs());
    Pstream::exchange<labelList, label>
    (
        wantedRemoteElements,
        subMap_,
        tag,
        Pstream::worldComm      // TBD
    );

    // Renumber elements
    forAll(cellCells, cellI)
    {
        labelList& cCells = cellCells[cellI];

        forAll(cCells, i)
        {
            cCells[i] = renumber(globalNumbering, compactMap, cCells[i]);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mapDistributeBase::mapDistributeBase()
:
    constructSize_(0),
    subHasFlip_(false),
    constructHasFlip_(false),
    schedulePtr_()
{}


Foam::mapDistributeBase::mapDistributeBase
(
    const label constructSize,
    const labelListList&& subMap,
    const labelListList&& constructMap,
    const bool subHasFlip,
    const bool constructHasFlip
)
:
    constructSize_(constructSize),
    subMap_(move(subMap)),
    constructMap_(move(constructMap)),
    subHasFlip_(subHasFlip),
    constructHasFlip_(constructHasFlip),
    schedulePtr_()
{}


Foam::mapDistributeBase::mapDistributeBase
(
    const labelList& sendProcs,
    const labelList& recvProcs
)
:
    constructSize_(0),
    subHasFlip_(false),
    constructHasFlip_(false),
    schedulePtr_()
{
    if (sendProcs.size() != recvProcs.size())
    {
        FatalErrorInFunction
            << "The send and receive data is not the same length. sendProcs:"
            << sendProcs.size() << " recvProcs:" << recvProcs.size()
            << abort(FatalError);
    }

    // Per processor the number of samples we have to send/receive.
    labelList nSend(Pstream::nProcs(), 0);
    labelList nRecv(Pstream::nProcs(), 0);

    forAll(sendProcs, sampleI)
    {
        label sendProc = sendProcs[sampleI];
        label recvProc = recvProcs[sampleI];

        // Note that also need to include local communication (both
        // RecvProc and sendProc on local processor)

        if (Pstream::myProcNo() == sendProc)
        {
            // I am the sender. Count destination processor.
            nSend[recvProc]++;
        }
        if (Pstream::myProcNo() == recvProc)
        {
            // I am the receiver.
            nRecv[sendProc]++;
        }
    }

    subMap_.setSize(Pstream::nProcs());
    constructMap_.setSize(Pstream::nProcs());
    forAll(nSend, proci)
    {
        subMap_[proci].setSize(nSend[proci]);
        constructMap_[proci].setSize(nRecv[proci]);
    }
    nSend = 0;
    nRecv = 0;

    forAll(sendProcs, sampleI)
    {
        label sendProc = sendProcs[sampleI];
        label recvProc = recvProcs[sampleI];

        if (Pstream::myProcNo() == sendProc)
        {
            // I am the sender. Store index I need to send.
            subMap_[recvProc][nSend[recvProc]++] = sampleI;
        }
        if (Pstream::myProcNo() == recvProc)
        {
            // I am the receiver.
            constructMap_[sendProc][nRecv[sendProc]++] = sampleI;
            // Largest entry inside constructMap
            constructSize_ = sampleI+1;
        }
    }
}


Foam::mapDistributeBase::mapDistributeBase
(
    const globalIndex& globalNumbering,
    labelList& elements,
    List<Map<label>>& compactMap,
    const int tag
)
:
    constructSize_(0),
    subHasFlip_(false),
    constructHasFlip_(false),
    schedulePtr_()
{
    // Construct per processor compact addressing of the global elements
    // needed. The ones from the local processor are not included since
    // these are always all needed.
    calcCompactAddressing
    (
        globalNumbering,
        elements,
        compactMap
    );

    //// Sort remote elements needed (not really necessary)
    // forAll(compactMap, proci)
    //{
    //    if (proci != Pstream::myProcNo())
    //    {
    //        Map<label>& globalMap = compactMap[proci];
    //
    //        SortableList<label> sorted(move(globalMap.toc()));
    //
    //        forAll(sorted, i)
    //        {
    //            Map<label>::iterator iter = globalMap.find(sorted[i]);
    //            iter() = i;
    //        }
    //    }
    //}


    // Exchange what I need with processor that supplies it. Renumber elements
    // into compact numbering
    labelList compactStart;
    exchangeAddressing
    (
        tag,
        globalNumbering,
        elements,
        compactMap,
        compactStart
    );

    if (debug)
    {
        printLayout(Pout);
    }
}


Foam::mapDistributeBase::mapDistributeBase
(
    const globalIndex& globalNumbering,
    labelListList& cellCells,
    List<Map<label>>& compactMap,
    const int tag
)
:
    constructSize_(0),
    subHasFlip_(false),
    constructHasFlip_(false),
    schedulePtr_()
{
    // Construct per processor compact addressing of the global elements
    // needed. The ones from the local processor are not included since
    // these are always all needed.
    calcCompactAddressing
    (
        globalNumbering,
        cellCells,
        compactMap
    );

    //// Sort remote elements needed (not really necessary)
    // forAll(compactMap, proci)
    //{
    //    if (proci != Pstream::myProcNo())
    //    {
    //        Map<label>& globalMap = compactMap[proci];
    //
    //        SortableList<label> sorted(move(globalMap.toc()));
    //
    //        forAll(sorted, i)
    //        {
    //            Map<label>::iterator iter = globalMap.find(sorted[i]);
    //            iter() = i;
    //        }
    //    }
    //}


    // Exchange what I need with processor that supplies it. Renumber elements
    // into compact numbering
    labelList compactStart;
    exchangeAddressing
    (
        tag,
        globalNumbering,
        cellCells,
        compactMap,
        compactStart
    );

    if (debug)
    {
        printLayout(Pout);
    }
}


Foam::mapDistributeBase::mapDistributeBase(const mapDistributeBase& map)
:
    constructSize_(map.constructSize_),
    subMap_(map.subMap_),
    constructMap_(map.constructMap_),
    subHasFlip_(map.subHasFlip_),
    constructHasFlip_(map.constructHasFlip_),
    schedulePtr_()
{}


Foam::mapDistributeBase::mapDistributeBase(mapDistributeBase&& map)
:
    constructSize_(map.constructSize_),
    subMap_(move(map.subMap_)),
    constructMap_(move(map.constructMap_)),
    subHasFlip_(map.subHasFlip_),
    constructHasFlip_(map.constructHasFlip_),
    schedulePtr_()
{}


Foam::mapDistributeBase::mapDistributeBase(Istream& is)
{
    is >> *this;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::mapDistributeBase::transfer(mapDistributeBase& rhs)
{
    constructSize_ = rhs.constructSize_;
    subMap_.transfer(rhs.subMap_);
    constructMap_.transfer(rhs.constructMap_);
    subHasFlip_ = rhs.subHasFlip_;
    constructHasFlip_ = rhs.constructHasFlip_;
    schedulePtr_.clear();
}


Foam::label Foam::mapDistributeBase::renumber
(
    const globalIndex& globalNumbering,
    const List<Map<label>>& compactMap,
    const label globalI
)
{
    if (globalI == -1)
    {
        return globalI;
    }
    if (globalNumbering.isLocal(globalI))
    {
        return globalNumbering.toLocal(globalI);
    }
    else
    {
        label proci = globalNumbering.whichProcID(globalI);
        label index = globalNumbering.toLocal(proci, globalI);
        return compactMap[proci][index];
    }
}


void Foam::mapDistributeBase::compact(const boolList& elemIsUsed, const int tag)
{
    // 1. send back to sender. Have sender delete the corresponding element
    //    from the submap and do the same to the constructMap locally
    //    (and in same order).

    // Send elemIsUsed field to neighbour. Use nonblocking code from
    // mapDistributeBase but in reverse order.
    if (Pstream::parRun())
    {
        label startOfRequests = Pstream::nRequests();

        // Set up receives from neighbours

        List<boolList> recvFields(Pstream::nProcs());

        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            const labelList& map = subMap_[domain];

            if (domain != Pstream::myProcNo() && map.size())
            {
                recvFields[domain].setSize(map.size());
                IPstream::read
                (
                    Pstream::commsTypes::nonBlocking,
                    domain,
                    reinterpret_cast<char*>(recvFields[domain].begin()),
                    recvFields[domain].size()*sizeof(bool),
                    tag
                );
            }
        }


        List<boolList> sendFields(Pstream::nProcs());

        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            const labelList& map = constructMap_[domain];

            if (domain != Pstream::myProcNo() && map.size())
            {
                boolList& subField = sendFields[domain];
                subField.setSize(map.size());
                forAll(map, i)
                {
                    subField[i] = accessAndFlip
                    (
                        elemIsUsed,
                        map[i],
                        constructHasFlip_,
                        noOp()          // do not flip elemIsUsed value
                    );
                }

                OPstream::write
                (
                    Pstream::commsTypes::nonBlocking,
                    domain,
                    reinterpret_cast<const char*>(subField.begin()),
                    subField.size()*sizeof(bool),
                    tag
                );
            }
        }



        // Set up 'send' to myself - write directly into recvFields

        {
            const labelList& map = constructMap_[Pstream::myProcNo()];

            recvFields[Pstream::myProcNo()].setSize(map.size());
            forAll(map, i)
            {
                recvFields[Pstream::myProcNo()][i] = accessAndFlip
                (
                    elemIsUsed,
                    map[i],
                    constructHasFlip_,
                    noOp()              // do not flip elemIsUsed value
                );
            }
        }


        // Wait for all to finish

        Pstream::waitRequests(startOfRequests);


        // Compact out all submap entries that are referring to unused elements
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            const labelList& map = subMap_[domain];

            labelList newMap(map.size());
            label newI = 0;

            forAll(map, i)
            {
                if (recvFields[domain][i])
                {
                    // So element is used on destination side
                    newMap[newI++] = map[i];
                }
            }
            if (newI < map.size())
            {
                newMap.setSize(newI);
                subMap_[domain].transfer(newMap);
            }
        }
    }


    // 2. remove from construct map - since end-result (element in elemIsUsed)
    //    not used.

    label maxConstructIndex = -1;

    for (label domain = 0; domain < Pstream::nProcs(); domain++)
    {
        const labelList& map = constructMap_[domain];

        labelList newMap(map.size());
        label newI = 0;

        forAll(map, i)
        {
            label destinationI = map[i];
            if (constructHasFlip_)
            {
                destinationI = mag(destinationI)-1;
            }

            // Is element is used on destination side
            if (elemIsUsed[destinationI])
            {
                maxConstructIndex = max(maxConstructIndex, destinationI);

                newMap[newI++] = map[i];
            }
        }
        if (newI < map.size())
        {
            newMap.setSize(newI);
            constructMap_[domain].transfer(newMap);
        }
    }

    constructSize_ = maxConstructIndex+1;

    // Clear the schedule (note:not necessary if nothing changed)
    schedulePtr_.clear();
}


void Foam::mapDistributeBase::compact
(
    const boolList& elemIsUsed,
    const label localSize,            // max index for subMap
    labelList& oldToNewSub,
    labelList& oldToNewConstruct,
    const int tag
)
{
    // 1. send back to sender. Have sender delete the corresponding element
    //    from the submap and do the same to the constructMap locally
    //    (and in same order).

    // Send elemIsUsed field to neighbour. Use nonblocking code from
    // mapDistributeBase but in reverse order.
    if (Pstream::parRun())
    {
        label startOfRequests = Pstream::nRequests();

        // Set up receives from neighbours

        List<boolList> recvFields(Pstream::nProcs());

        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            const labelList& map = subMap_[domain];

            if (domain != Pstream::myProcNo() && map.size())
            {
                recvFields[domain].setSize(map.size());
                IPstream::read
                (
                    Pstream::commsTypes::nonBlocking,
                    domain,
                    reinterpret_cast<char*>(recvFields[domain].begin()),
                    recvFields[domain].size()*sizeof(bool),
                    tag
                );
            }
        }


        List<boolList> sendFields(Pstream::nProcs());

        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            const labelList& map = constructMap_[domain];

            if (domain != Pstream::myProcNo() && map.size())
            {
                boolList& subField = sendFields[domain];
                subField.setSize(map.size());
                forAll(map, i)
                {
                    label index = map[i];
                    if (constructHasFlip_)
                    {
                        index = mag(index)-1;
                    }
                    subField[i] = elemIsUsed[index];
                }

                OPstream::write
                (
                    Pstream::commsTypes::nonBlocking,
                    domain,
                    reinterpret_cast<const char*>(subField.begin()),
                    subField.size()*sizeof(bool),
                    tag
                );
            }
        }



        // Set up 'send' to myself - write directly into recvFields

        {
            const labelList& map = constructMap_[Pstream::myProcNo()];

            recvFields[Pstream::myProcNo()].setSize(map.size());
            forAll(map, i)
            {
                label index = map[i];
                if (constructHasFlip_)
                {
                    index = mag(index)-1;
                }
                recvFields[Pstream::myProcNo()][i] = elemIsUsed[index];
            }
        }


        // Wait for all to finish

        Pstream::waitRequests(startOfRequests);




        // Work out which elements on the sending side are needed
        {
            oldToNewSub.setSize(localSize, -1);

            boolList sendElemIsUsed(localSize, false);

            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = subMap_[domain];
                forAll(map, i)
                {
                    if (recvFields[domain][i])
                    {
                        label index = map[i];
                        if (subHasFlip_)
                        {
                            index = mag(index)-1;
                        }
                        sendElemIsUsed[index] = true;
                    }
                }
            }

            label newI = 0;
            forAll(sendElemIsUsed, i)
            {
                if (sendElemIsUsed[i])
                {
                    oldToNewSub[i] = newI++;
                }
            }
        }


        // Compact out all submap entries that are referring to unused elements
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            const labelList& map = subMap_[domain];

            labelList newMap(map.size());
            label newI = 0;

            forAll(map, i)
            {
                if (recvFields[domain][i])
                {
                    // So element is used on destination side
                    label index = map[i];
                    label sign = 1;
                    if (subHasFlip_)
                    {
                        if (index < 0)
                        {
                            sign = -1;
                        }
                        index = mag(index)-1;
                    }
                    label newIndex = oldToNewSub[index];
                    if (subHasFlip_)
                    {
                        newIndex = sign*(newIndex+1);
                    }
                    newMap[newI++] = newIndex;
                }
            }
            newMap.setSize(newI);
            subMap_[domain].transfer(newMap);
        }
    }


    // 2. remove from construct map - since end-result (element in elemIsUsed)
    //    not used.


    oldToNewConstruct.setSize(elemIsUsed.size(), -1);
    constructSize_ = 0;
    forAll(elemIsUsed, i)
    {
        if (elemIsUsed[i])
        {
            oldToNewConstruct[i] = constructSize_++;
        }
    }

    for (label domain = 0; domain < Pstream::nProcs(); domain++)
    {
        const labelList& map = constructMap_[domain];

        labelList newMap(map.size());
        label newI = 0;

        forAll(map, i)
        {
            label destinationI = map[i];
            label sign = 1;
            if (constructHasFlip_)
            {
                if (destinationI < 0)
                {
                    sign = -1;
                }
                destinationI = mag(destinationI)-1;
            }

            // Is element is used on destination side
            if (elemIsUsed[destinationI])
            {
                label newIndex = oldToNewConstruct[destinationI];
                if (constructHasFlip_)
                {
                    newIndex = sign*(newIndex+1);
                }
                newMap[newI++] = newIndex;
            }
        }
        newMap.setSize(newI);
        constructMap_[domain].transfer(newMap);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::mapDistributeBase::operator=(const mapDistributeBase& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorInFunction
            << "Attempted assignment to self"
            << abort(FatalError);
    }
    constructSize_ = rhs.constructSize_;
    subMap_ = rhs.subMap_;
    constructMap_ = rhs.constructMap_;
    subHasFlip_ = rhs.subHasFlip_;
    constructHasFlip_ = rhs.constructHasFlip_;
    schedulePtr_.clear();
}


// * * * * * * * * * * * * * * Istream Operator  * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, mapDistributeBase& map)
{
    is.fatalCheck("operator>>(Istream&, mapDistributeBase&)");

    is  >> map.constructSize_ >> map.subMap_ >> map.constructMap_
        >> map.subHasFlip_ >> map.constructHasFlip_;

    return is;
}


// * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const mapDistributeBase& map)
{
    os  << map.constructSize_ << token::NL
        << map.subMap_ << token::NL
        << map.constructMap_ << token::NL
        << map.subHasFlip_ << token::SPACE << map.constructHasFlip_
        << token::NL;

    return os;
}


// ************************************************************************* //
