/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "UPstream.H"
#include "debug.H"
#include "registerSwitch.H"
#include "dictionary.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(UPstream, 0);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::UPstream::commsTypes,
        3
    >::names[] =
    {
        "blocking",
        "scheduled",
        "nonBlocking"
    };
}


const Foam::NamedEnum<Foam::UPstream::commsTypes, 3>
    Foam::UPstream::commsTypeNames;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::UPstream::setParRun(const label nProcs, const bool haveThreads)
{
    if (nProcs == 0)
    {
        parRun_ = false;
        haveThreads_ = haveThreads;

        freeCommunicator(UPstream::worldComm);
        label comm = allocateCommunicator(-1, labelList(1, label(0)), false);
        if (comm != UPstream::worldComm)
        {
            FatalErrorIn("UPstream::setParRun(const label)")
                << "problem : comm:" << comm
                << "  UPstream::worldComm:" << UPstream::worldComm
                << Foam::exit(FatalError);
        }

        Pout.prefix() = "";
        Perr.prefix() = "";
    }
    else
    {
        parRun_ = true;
        haveThreads_ = haveThreads;

        // Redo worldComm communicator (this has been created at static
        // initialisation time)
        freeCommunicator(UPstream::worldComm);
        label comm = allocateCommunicator(-1, identity(nProcs), true);
        if (comm != UPstream::worldComm)
        {
            FatalErrorInFunction
                << "problem : comm:" << comm
                << "  UPstream::worldComm:" << UPstream::worldComm
                << Foam::exit(FatalError);
        }

        Pout.prefix() = '[' +  name(myProcNo(Pstream::worldComm)) + "] ";
        Perr.prefix() = '[' +  name(myProcNo(Pstream::worldComm)) + "] ";
    }
}


Foam::List<Foam::UPstream::commsStruct> Foam::UPstream::calcLinearComm
(
    const label nProcs
)
{
    List<commsStruct> linearCommunication(nProcs);

    // Master
    labelList belowIDs(nProcs - 1);
    forAll(belowIDs, i)
    {
        belowIDs[i] = i + 1;
    }

    linearCommunication[0] = commsStruct
    (
        nProcs,
        0,
        -1,
        belowIDs,
        labelList(0)
    );

    // Slaves. Have no below processors, only communicate up to master
    for (label procID = 1; procID < nProcs; procID++)
    {
        linearCommunication[procID] = commsStruct
        (
            nProcs,
            procID,
            0,
            labelList(0),
            labelList(0)
        );
    }
    return linearCommunication;
}


void Foam::UPstream::collectReceives
(
    const label procID,
    const List<DynamicList<label>>& receives,
    DynamicList<label>& allReceives
)
{
    // Append my children (and my children children etc.) to allReceives.

    const DynamicList<label>& myChildren = receives[procID];

    forAll(myChildren, childI)
    {
        allReceives.append(myChildren[childI]);
        collectReceives(myChildren[childI], receives, allReceives);
    }
}


Foam::List<Foam::UPstream::commsStruct> Foam::UPstream::calcTreeComm
(
    label nProcs
)
{
    // Tree like schedule. For 8 procs:
    // (level 0)
    //      0 receives from 1
    //      2 receives from 3
    //      4 receives from 5
    //      6 receives from 7
    // (level 1)
    //      0 receives from 2
    //      4 receives from 6
    // (level 2)
    //      0 receives from 4
    //
    // The sends/receives for all levels are collected per processor
    //  (one send per processor; multiple receives possible) creating a table:
    //
    // So per processor:
    // proc     receives from   sends to
    // ----     -------------   --------
    //  0       1,2,4           -
    //  1       -               0
    //  2       3               0
    //  3       -               2
    //  4       5               0
    //  5       -               4
    //  6       7               4
    //  7       -               6

    label nLevels = 1;
    while ((1 << nLevels) < nProcs)
    {
        nLevels++;
    }

    List<DynamicList<label>> receives(nProcs);
    labelList sends(nProcs, -1);

    // Info<< "Using " << nLevels << " communication levels" << endl;

    label offset = 2;
    label childOffset = offset/2;

    for (label level = 0; level < nLevels; level++)
    {
        label receiveID = 0;
        while (receiveID < nProcs)
        {
            // Determine processor that sends and we receive from
            label sendID = receiveID + childOffset;

            if (sendID < nProcs)
            {
                receives[receiveID].append(sendID);
                sends[sendID] = receiveID;
            }

            receiveID += offset;
        }

        offset <<= 1;
        childOffset <<= 1;
    }

    // For all processors find the processors it receives data from
    // (and the processors they receive data from etc.)
    List<DynamicList<label>> allReceives(nProcs);
    for (label procID = 0; procID < nProcs; procID++)
    {
        collectReceives(procID, receives, allReceives[procID]);
    }


    List<commsStruct> treeCommunication(nProcs);

    for (label procID = 0; procID < nProcs; procID++)
    {
        treeCommunication[procID] = commsStruct
        (
            nProcs,
            procID,
            sends[procID],
            receives[procID].shrink(),
            allReceives[procID].shrink()
        );
    }
    return treeCommunication;
}


Foam::label Foam::UPstream::allocateCommunicator
(
    const label parentIndex,
    const labelList& subRanks,
    const bool doPstream
)
{
    label index;
    if (!freeComms_.empty())
    {
        index = freeComms_.pop();
    }
    else
    {
        // Extend storage
        index = parentCommunicator_.size();

        myProcNo_.append(-1);
        procIDs_.append(List<int>(0));
        parentCommunicator_.append(-1);
        linearCommunication_.append(List<commsStruct>(0));
        treeCommunication_.append(List<commsStruct>(0));
    }

    if (debug)
    {
        Pout<< "Communicators : Allocating communicator " << index << endl
            << "    parent : " << parentIndex << endl
            << "    procs  : " << subRanks << endl
            << endl;
    }

    // Initialise; overwritten by allocatePstreamCommunicator
    myProcNo_[index] = 0;

    // Convert from label to int
    procIDs_[index].setSize(subRanks.size());
    forAll(procIDs_[index], i)
    {
        procIDs_[index][i] = subRanks[i];

        // Enforce incremental order (so index is rank in next communicator)
        if (i >= 1 && subRanks[i] <= subRanks[i-1])
        {
            FatalErrorInFunction
                << "subranks not sorted : " << subRanks
                << " when allocating subcommunicator from parent "
                << parentIndex
                << Foam::abort(FatalError);
        }
    }
    parentCommunicator_[index] = parentIndex;

    linearCommunication_[index] = calcLinearComm(procIDs_[index].size());
    treeCommunication_[index] = calcTreeComm(procIDs_[index].size());


    if (doPstream && parRun())
    {
        allocatePstreamCommunicator(parentIndex, index);
    }

    return index;
}


void Foam::UPstream::freeCommunicator
(
    const label communicator,
    const bool doPstream
)
{
    if (debug)
    {
        Pout<< "Communicators : Freeing communicator " << communicator << endl
            << "    parent   : " << parentCommunicator_[communicator] << endl
            << "    myProcNo : " << myProcNo_[communicator] << endl
            << endl;
    }

    if (doPstream && parRun())
    {
        freePstreamCommunicator(communicator);
    }
    myProcNo_[communicator] = -1;
    // procIDs_[communicator].clear();
    parentCommunicator_[communicator] = -1;
    linearCommunication_[communicator].clear();
    treeCommunication_[communicator].clear();

    freeComms_.push(communicator);
}


void Foam::UPstream::freeCommunicators(const bool doPstream)
{
    forAll(myProcNo_, communicator)
    {
        if (myProcNo_[communicator] != -1)
        {
            freeCommunicator(communicator, doPstream);
        }
    }
}


int Foam::UPstream::baseProcNo(const label myComm, const int myProcID)
{
    int procID = myProcID;
    label comm = myComm;

    while (parent(comm) != -1)
    {
        const List<int>& parentRanks = UPstream::procID(comm);
        procID = parentRanks[procID];
        comm = UPstream::parent(comm);
    }

    return procID;
}


Foam::label Foam::UPstream::procNo(const label myComm, const int baseProcID)
{
    const List<int>& parentRanks = procID(myComm);
    label parentComm = parent(myComm);

    if (parentComm == -1)
    {
        return findIndex(parentRanks, baseProcID);
    }
    else
    {
        label parentRank = procNo(parentComm, baseProcID);
        return findIndex(parentRanks, parentRank);
    }
}


Foam::label Foam::UPstream::procNo
(
    const label myComm,
    const label currentComm,
    const int currentProcID
)
{
    label physProcID = UPstream::baseProcNo(currentComm, currentProcID);
    return procNo(myComm, physProcID);
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::UPstream::parRun_(false);

bool Foam::UPstream::haveThreads_(false);

Foam::LIFOStack<Foam::label> Foam::UPstream::freeComms_;

Foam::DynamicList<int> Foam::UPstream::myProcNo_(10);

Foam::DynamicList<Foam::List<int>> Foam::UPstream::procIDs_(10);

Foam::DynamicList<Foam::label> Foam::UPstream::parentCommunicator_(10);

int Foam::UPstream::msgType_(1);


Foam::DynamicList<Foam::List<Foam::UPstream::commsStruct>>
Foam::UPstream::linearCommunication_(10);

Foam::DynamicList<Foam::List<Foam::UPstream::commsStruct>>
Foam::UPstream::treeCommunication_(10);


// Allocate a serial communicator. This gets overwritten in parallel mode
// (by UPstream::setParRun())
Foam::UPstream::communicator serialComm
(
    -1,
    Foam::labelList(1, Foam::label(0)),
    false
);


bool Foam::UPstream::floatTransfer
(
    Foam::debug::optimisationSwitch("floatTransfer", 0)
);
registerOptSwitch
(
    "floatTransfer",
    bool,
    Foam::UPstream::floatTransfer
);

int Foam::UPstream::nProcsSimpleSum
(
    Foam::debug::optimisationSwitch("nProcsSimpleSum", 16)
);
registerOptSwitch
(
    "nProcsSimpleSum",
    int,
    Foam::UPstream::nProcsSimpleSum
);

Foam::UPstream::commsTypes Foam::UPstream::defaultCommsType
(
    commsTypeNames.read(Foam::debug::optimisationSwitches().lookup("commsType"))
);

namespace Foam
{
    // Register re-reader
    class addcommsTypeToOpt
    :
        public ::Foam::simpleRegIOobject
    {
    public:

        addcommsTypeToOpt(const char* name)
        :
            ::Foam::simpleRegIOobject(Foam::debug::addOptimisationObject, name)
        {}

        virtual ~addcommsTypeToOpt()
        {}

        virtual void readData(Foam::Istream& is)
        {
            UPstream::defaultCommsType = UPstream::commsTypeNames.read
            (
                is
            );
        }

        virtual void writeData(Foam::Ostream& os) const
        {
            os << UPstream::commsTypeNames[UPstream::defaultCommsType];
        }
    };

    addcommsTypeToOpt addcommsTypeToOpt_("commsType");
}

Foam::label Foam::UPstream::worldComm(0);

Foam::label Foam::UPstream::warnComm(-1);

int Foam::UPstream::nPollProcInterfaces
(
    Foam::debug::optimisationSwitch("nPollProcInterfaces", 0)
);
registerOptSwitch
(
    "nPollProcInterfaces",
    int,
    Foam::UPstream::nPollProcInterfaces
);


// ************************************************************************* //
