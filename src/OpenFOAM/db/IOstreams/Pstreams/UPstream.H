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

Class
    Foam::UPstream

Description
    Inter-processor communications stream

SourceFiles
    UPstream.C
    UPstreamCommsStruct.C
    gatherScatter.C
    combineGatherScatter.C
    gatherScatterList.C

\*---------------------------------------------------------------------------*/

#ifndef UPstream_H
#define UPstream_H

#include "labelList.H"
#include "DynamicList.H"
#include "HashTable.H"
#include "string.H"
#include "NamedEnum.H"
#include "ListOps.H"
#include "LIFOStack.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class UPstream Declaration
\*---------------------------------------------------------------------------*/

class UPstream
{

public:

    //- Types of communications
    enum class commsTypes
    {
        blocking,
        scheduled,
        nonBlocking
    };

    static const NamedEnum<commsTypes, 3> commsTypeNames;

    // Public classes

        //- Structure for communicating between processors
        class commsStruct
        {
            // Private data

                //- procID of above processor
                label above_;

                //- procIDs of processors directly below me
                labelList below_;

                //- procIDs of all processors below (so not just directly below)
                labelList allBelow_;

                //- procIDs of all processors not below. (inverse set of
                //  allBelow_ and minus myProcNo)
                labelList allNotBelow_;


        public:

            // Constructors

                //- Construct null
                commsStruct();

                //- Construct from components
                commsStruct
                (
                    const label,
                    const labelList&,
                    const labelList&,
                    const labelList&
                );

                //- Construct from components; construct allNotBelow_
                commsStruct
                (
                    const label nProcs,
                    const label myProcID,
                    const label,
                    const labelList&,
                    const labelList&
                );


            // Member Functions

                // Access

                    label above() const
                    {
                        return above_;
                    }

                    const labelList& below() const
                    {
                        return below_;
                    }

                    const labelList& allBelow() const
                    {
                        return allBelow_;
                    }

                    const labelList& allNotBelow() const
                    {
                        return allNotBelow_;
                    }


            // Member operators

                bool operator==(const commsStruct&) const;

                bool operator!=(const commsStruct&) const;


             // Ostream Operator

                friend Ostream& operator<<(Ostream&, const commsStruct&);
        };


        //- combineReduce operator for lists. Used for counting.
        class listEq
        {

        public:

            template<class T>
            void operator()(T& x, const T& y) const
            {
                forAll(y, i)
                {
                    if (y[i].size())
                    {
                        x[i] = y[i];
                    }
                }
            }
        };


private:

    // Private data

        //- By default this is not a parallel run
        static bool parRun_;

        //- Have support for threads?
        static bool haveThreads_;

        //- Standard transfer message type
        static int msgType_;

        // Communicator specific data

        //- Free communicators
        static LIFOStack<label> freeComms_;

        //- My processor number
        static DynamicList<int> myProcNo_;

        //- List of process IDs
        static DynamicList<List<int>> procIDs_;

        //- Parent communicator
        static DynamicList<label> parentCommunicator_;

        //- Linear communication schedule
        static DynamicList<List<commsStruct>> linearCommunication_;

        //- Multi level communication schedule
        static DynamicList<List<commsStruct>> treeCommunication_;


    // Private Member Functions

        //- Set data for parallel running
        static void setParRun(const label nProcs, const bool haveThreads);

        //- Calculate linear communication schedule
        static List<commsStruct> calcLinearComm(const label nProcs);

        //- Calculate tree communication schedule
        static List<commsStruct> calcTreeComm(const label nProcs);

        //- Helper function for tree communication schedule determination
        //  Collects all processorIDs below a processor
        static void collectReceives
        (
            const label procID,
            const List<DynamicList<label>>& receives,
            DynamicList<label>& allReceives
        );

        //- Allocate a communicator with index
        static void allocatePstreamCommunicator
        (
            const label parentIndex,
            const label index
        );

        //- Free a communicator
        static void freePstreamCommunicator
        (
            const label index
        );


protected:

    // Protected data

        //- Communications type of this stream
        commsTypes commsType_;

public:

    // Declare name of the class and its debug switch
    ClassName("UPstream");


    // Static data

        //- Should compact transfer be used in which floats replace doubles
        //  reducing the bandwidth requirement at the expense of some loss
        //  in accuracy
        static bool floatTransfer;

        //- Number of processors at which the sum algorithm changes from linear
        //  to tree
        static int nProcsSimpleSum;

        //- Default commsType
        static commsTypes defaultCommsType;

        //- Number of polling cycles in processor updates
        static int nPollProcInterfaces;

        //- Default communicator (all processors)
        static label worldComm;

        //- Debugging: warn for use of any communicator differing from warnComm
        static label warnComm;


    // Constructors

        //- Construct given optional buffer size
        UPstream(const commsTypes commsType)
        :
            commsType_(commsType)
        {}


    // Member functions

        //- Allocate a new communicator
        static label allocateCommunicator
        (
            const label parent,
            const labelList& subRanks,
            const bool doPstream = true
        );

        //- Free a previously allocated communicator
        static void freeCommunicator
        (
            const label communicator,
            const bool doPstream = true
        );

        //- Free all communicators
        static void freeCommunicators(const bool doPstream);

        //- Helper class for allocating/freeing communicators
        class communicator
        {
            label comm_;

            //- Disallow copy and assignment
            communicator(const communicator&);
            void operator=(const communicator&);

        public:

            communicator
            (
                const label parent,
                const labelList& subRanks,
                const bool doPstream
            )
            :
                comm_(allocateCommunicator(parent, subRanks, doPstream))
            {}

            ~communicator()
            {
                freeCommunicator(comm_);
            }

            operator label() const
            {
                return comm_;
            }
        };

        //- Return physical processor number (i.e. processor number in
        //  worldComm) given communicator and procssor
        static int baseProcNo(const label myComm, const int procID);

        //- Return processor number in communicator (given physical processor
        //  number) (= reverse of baseProcNo)
        static label procNo(const label comm, const int baseProcID);

        //- Return processor number in communicator (given processor number
        //  and communicator)
        static label procNo
        (
            const label myComm,
            const label currentComm,
            const int currentProcID
        );

        //- Add the valid option this type of communications library
        //  adds/requires on the command line
        static void addValidParOptions(HashTable<string>& validParOptions);

        //- Initialisation function called from main
        //  Spawns slave processes and initialises inter-communication
        static bool init(int& argc, char**& argv, const bool needsThread);

        // Non-blocking comms

            //- Get number of outstanding requests
            static label nRequests();

            //- Truncate number of outstanding requests
            static void resetRequests(const label sz);

            //- Wait until all requests (from start onwards) have finished.
            static void waitRequests(const label start = 0);

            //- Wait until request i has finished.
            static void waitRequest(const label i);

            //- Non-blocking comms: has request i finished?
            static bool finishedRequest(const label i);

            static int allocateTag(const char*);

            static int allocateTag(const word&);

            static void freeTag(const char*, const int tag);

            static void freeTag(const word&, const int tag);


        //- Is this a parallel run?
        static bool& parRun()
        {
            return parRun_;
        }

        //- Have support for threads
        static bool haveThreads()
        {
            return haveThreads_;
        }

        //- Number of processes in parallel run
        static label nProcs(const label communicator = 0)
        {
            return procIDs_[communicator].size();
        }

        //- Process index of the master
        static int masterNo()
        {
            return 0;
        }

        //- Am I the master process
        static bool master(const label communicator = 0)
        {
            return myProcNo_[communicator] == masterNo();
        }

        //- Number of this process (starting from masterNo() = 0)
        static int myProcNo(const label communicator = 0)
        {
            return myProcNo_[communicator];
        }

        static label parent(const label communicator)
        {
            return parentCommunicator_(communicator);
        }

        //- Process ID of given process index
        static List<int>& procID(label communicator)
        {
            return procIDs_[communicator];
        }

        //- Process index of first slave
        static int firstSlave()
        {
            return 1;
        }

        //- Process index of last slave
        static int lastSlave(const label communicator = 0)
        {
            return nProcs(communicator) - 1;
        }

        //- Communication schedule for linear all-to-master (proc 0)
        static const List<commsStruct>& linearCommunication
        (
            const label communicator = 0
        )
        {
            return linearCommunication_[communicator];
        }

        //- Communication schedule for tree all-to-master (proc 0)
        static const List<commsStruct>& treeCommunication
        (
            const label communicator = 0
        )
        {
            return treeCommunication_[communicator];
        }

        //- Message tag of standard messages
        static int& msgType()
        {
            return msgType_;
        }


            //- Get the communications type of the stream
            commsTypes commsType() const
            {
                return commsType_;
            }

            //- Set the communications type of the stream
            commsTypes commsType(const commsTypes ct)
            {
                commsTypes oldCommsType = commsType_;
                commsType_ = ct;
                return oldCommsType;
            }


        //- Exit program
        static void exit(int errnum = 1);

        //- Abort program
        static void abort();

        //- Exchange label with all processors (in the communicator).
        //  sendData[proci] is the label to send to proci.
        //  After return recvData contains the data from the other processors.
        static void allToAll
        (
            const labelUList& sendData,
            labelUList& recvData,
            const label communicator = 0
        );

        //- Exchange data with all processors (in the communicator)
        //  sendSizes, sendOffsets give (per processor) the slice of
        //  sendData to send, similarly recvSizes, recvOffsets give the slice
        //  of recvData to receive
        static void allToAll
        (
            const char* sendData,
            const UList<int>& sendSizes,
            const UList<int>& sendOffsets,

            char* recvData,
            const UList<int>& recvSizes,
            const UList<int>& recvOffsets,

            const label communicator = 0
        );

        //- Receive data from all processors on the master
        static void gather
        (
            const char* sendData,
            int sendSize,

            char* recvData,
            const UList<int>& recvSizes,
            const UList<int>& recvOffsets,
            const label communicator = 0
        );

        //- Send data to all processors from the root of the communicator
        static void scatter
        (
            const char* sendData,
            const UList<int>& sendSizes,
            const UList<int>& sendOffsets,

            char* recvData,
            int recvSize,
            const label communicator = 0
        );
};


Ostream& operator<<(Ostream&, const UPstream::commsStruct&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
