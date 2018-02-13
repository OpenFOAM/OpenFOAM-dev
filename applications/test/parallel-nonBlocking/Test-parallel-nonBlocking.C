/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

Application
    Test-parallel-nonBlocking

Description
    Test for various non-blocking parallel routines.

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "mapDistribute.H"
#include "argList.H"
#include "Time.H"
#include "IPstream.H"
#include "OPstream.H"
#include "vector.H"
#include "IOstreams.H"
#include "Random.H"
#include "Tuple2.H"
#include "PstreamBuffers.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "setRootCase.H"
    #include "createTime.H"


    // Test PstreamBuffers
    // ~~~~~~~~~~~~~~~~~~~
    if (false)
    {
        Perr<< "\nStarting transfers\n" << endl;

        vector data
        (
            Pstream::myProcNo(),
            Pstream::myProcNo(),
            Pstream::myProcNo()
        );

        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        if (Pstream::myProcNo() != Pstream::masterNo())
        {
            Perr<< "slave sending to master "
                << Pstream::masterNo() << endl;
            UOPstream toMaster(Pstream::masterNo(), pBufs);
            toMaster << data;
        }

        // Start sending and receiving and block
        pBufs.finishedSends();

        // Consume
        DynamicList<vector> allData;
        if (Pstream::myProcNo() == Pstream::masterNo())
        {
            // Collect my own data
            allData.append(data);

            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                Perr << "master receiving from slave " << slave << endl;
                UIPstream fromSlave(slave, pBufs);
                allData.append(vector(fromSlave));
            }
        }


        // Send allData back
        PstreamBuffers pBufs2(Pstream::commsTypes::nonBlocking);
        if (Pstream::myProcNo() == Pstream::masterNo())
        {
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                Perr << "master sending to slave " << slave << endl;
                UOPstream toSlave(slave, pBufs2);
                toSlave << allData;
            }
        }

        // Start sending and receiving and block
        pBufs2.finishedSends();

        // Consume
        if (Pstream::myProcNo() != Pstream::masterNo())
        {
            Perr<< "slave receiving from master "
                << Pstream::masterNo() << endl;
            UIPstream fromMaster(Pstream::masterNo(), pBufs2);
            fromMaster >> allData;
            Perr<< allData << endl;
        }
    }


    // Test non-blocking reductions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    scalar data1 = 1.0;
    label request1 = -1;
    {
        Foam::reduce(data1, sumOp<scalar>(), Pstream::msgType(), request1);
    }

    scalar data2 = 0.1;
    label request2 = -1;
    {
        Foam::reduce(data2, sumOp<scalar>(), Pstream::msgType(), request2);
    }


    // Do a non-blocking send in between
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            UOPstream toProc(proci, pBufs);
            toProc << Pstream::myProcNo();
        }

        // Start sending and receiving and block
        pBufs.finishedSends();

        // Consume
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            UIPstream fromProc(proci, pBufs);
            label data;
            fromProc >> data;

            if (data != proci)
            {
                FatalErrorInFunction
                    << "From processor " << proci << " received " << data
                    << " but expected " << proci
                    << exit(FatalError);
            }
        }
    }


    if (request1 != -1)
    {
        Pout<< "Waiting for non-blocking reduce with request " << request1
            << endl;
        Pstream::waitRequest(request1);
    }
    Info<< "Reduced data1:" << data1 << endl;

    if (request2 != -1)
    {
        Pout<< "Waiting for non-blocking reduce with request " << request1
            << endl;
        Pstream::waitRequest(request2);
    }
    Info<< "Reduced data2:" << data2 << endl;


    // Clear any outstanding requests
    Pstream::resetRequests(0);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
