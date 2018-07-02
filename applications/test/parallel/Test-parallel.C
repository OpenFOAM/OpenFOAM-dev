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

Application
    parallelTest

Description
    Test for various parallel routines.

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

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "setRootCase.H"
    #include "createTime.H"


    // Test mapDistribute
    // ~~~~~~~~~~~~~~~~~~

    if (true)
    {
        Random rndGen(43544*Pstream::myProcNo());

        // Generate random data.
        List<Tuple2<label, List<scalar>>> complexData(100);
        forAll(complexData, i)
        {
            complexData[i].first() =
                rndGen.sampleAB<label>(0, Pstream::nProcs());
            complexData[i].second().setSize(3);
            complexData[i].second()[0] = 1;
            complexData[i].second()[1] = 2;
            complexData[i].second()[2] = 3;
        }

        // Send all ones to processor indicated by .first()


        // Count how many to send
        labelList nSend(Pstream::nProcs(), 0);
        forAll(complexData, i)
        {
            label proci = complexData[i].first();
            nSend[proci]++;
        }

        // Collect items to be sent
        labelListList sendMap(Pstream::nProcs());
        forAll(sendMap, proci)
        {
            sendMap[proci].setSize(nSend[proci]);
        }
        nSend = 0;
        forAll(complexData, i)
        {
            label proci = complexData[i].first();
            sendMap[proci][nSend[proci]++] = i;
        }

        // Sync how many to send
        labelList nRecv;
        Pstream::exchangeSizes(sendMap, nRecv);

        // Collect items to be received
        labelListList recvMap(Pstream::nProcs());
        forAll(recvMap, proci)
        {
            recvMap[proci].setSize(nRecv[proci]);
        }

        label constructSize = 0;
        // Construct with my own elements first
        forAll(recvMap[Pstream::myProcNo()], i)
        {
            recvMap[Pstream::myProcNo()][i] = constructSize++;
        }
        // Construct from other processors
        forAll(recvMap, proci)
        {
            if (proci != Pstream::myProcNo())
            {
                forAll(recvMap[proci], i)
                {
                    recvMap[proci][i] = constructSize++;
                }
            }
        }



        // Construct distribute map (destructively)
        mapDistribute map(constructSize, sendMap.xfer(), recvMap.xfer());

        // Distribute complexData
        map.distribute(complexData);

        Pout<< "complexData:" << complexData << endl;
    }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Perr<< "\nStarting transfers\n" << endl;

    vector data(0, 1, 2);

    if (Pstream::parRun())
    {
        if (Pstream::myProcNo() != Pstream::masterNo())
        {
            {
                Perr<< "slave sending to master "
                    << Pstream::masterNo() << endl;
                OPstream toMaster
                (
                    Pstream::commsTypes::blocking,
                    Pstream::masterNo()
                );
                toMaster << data;
            }

            Perr<< "slave receiving from master "
                << Pstream::masterNo() << endl;
            IPstream fromMaster
            (
                Pstream::commsTypes::blocking,
                Pstream::masterNo()
            );
            fromMaster >> data;

            Perr<< data << endl;
        }
        else
        {
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                Perr << "master receiving from slave " << slave << endl;
                IPstream fromSlave(Pstream::commsTypes::blocking, slave);
                fromSlave >> data;

                Perr<< data << endl;
            }

            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                Perr << "master sending to slave " << slave << endl;
                OPstream toSlave(Pstream::commsTypes::blocking, slave);
                toSlave << data;
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
