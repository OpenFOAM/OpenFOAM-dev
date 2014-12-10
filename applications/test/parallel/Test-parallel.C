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

#   include "setRootCase.H"
#   include "createTime.H"


    // Test mapDistribute
    // ~~~~~~~~~~~~~~~~~~

    if (false)
    {
        Random rndGen(43544*Pstream::myProcNo());

        // Generate random data.
        List<Tuple2<label, List<scalar> > > complexData(100);
        forAll(complexData, i)
        {
            complexData[i].first() = rndGen.integer(0, Pstream::nProcs()-1);
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
            label procI = complexData[i].first();
            nSend[procI]++;
        }

        // Sync how many to send
        labelListList allNTrans(Pstream::nProcs());
        allNTrans[Pstream::myProcNo()] = nSend;
        combineReduce(allNTrans, UPstream::listEq());

        // Collect items to be sent
        labelListList sendMap(Pstream::nProcs());
        forAll(sendMap, procI)
        {
            sendMap[procI].setSize(nSend[procI]);
        }
        nSend = 0;
        forAll(complexData, i)
        {
            label procI = complexData[i].first();
            sendMap[procI][nSend[procI]++] = i;
        }

        // Collect items to be received
        labelListList recvMap(Pstream::nProcs());
        forAll(recvMap, procI)
        {
            recvMap[procI].setSize(allNTrans[procI][Pstream::myProcNo()]);
        }

        label constructSize = 0;
        // Construct with my own elements first
        forAll(recvMap[Pstream::myProcNo()], i)
        {
            recvMap[Pstream::myProcNo()][i] = constructSize++;
        }
        // Construct from other processors
        forAll(recvMap, procI)
        {
            if (procI != Pstream::myProcNo())
            {
                forAll(recvMap[procI], i)
                {
                    recvMap[procI][i] = constructSize++;
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
                OPstream toMaster(Pstream::blocking, Pstream::masterNo());
                toMaster << data;
            }

            Perr<< "slave receiving from master "
                << Pstream::masterNo() << endl;
            IPstream fromMaster(Pstream::blocking, Pstream::masterNo());
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
                IPstream fromSlave(Pstream::blocking, slave);
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
                OPstream toSlave(Pstream::blocking, slave);
                toSlave << data;
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
