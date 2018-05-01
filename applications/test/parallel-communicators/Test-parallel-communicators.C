/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Test-parallel-communicators

Description
    Checks communication using user-defined communicators

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IPstream.H"
#include "OPstream.H"
#include "vector.H"
#include "IOstreams.H"
#include "PstreamReduceOps.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar sumReduce
(
    const label comm,
    const scalar localValue
)
{
    scalar sum = 0;
    if (Pstream::parRun())
    {
        if (UPstream::master(comm))
        {
            // Add master value and all slaves
            sum = localValue;

            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave(comm);
                slave++
            )
            {
                scalar slaveValue;
                UIPstream::read
                (
                    Pstream::commsTypes::blocking,
                    slave,
                    reinterpret_cast<char*>(&slaveValue),
                    sizeof(scalar),
                    UPstream::msgType(),    // tag
                    comm                    // communicator
                );

                sum += slaveValue;
            }

            // Send back to slaves

            for
            (
                int slave=UPstream::firstSlave();
                slave<=UPstream::lastSlave(comm);
                slave++
            )
            {
                UOPstream::write
                (
                    UPstream::commsTypes::blocking,
                    slave,
                    reinterpret_cast<const char*>(&sum),
                    sizeof(scalar),
                    UPstream::msgType(),    // tag
                    comm                    // communicator
                );
            }
        }
        else
        {
            {
                UOPstream::write
                (
                    UPstream::commsTypes::blocking,
                    UPstream::masterNo(),
                    reinterpret_cast<const char*>(&localValue),
                    sizeof(scalar),
                    UPstream::msgType(),    // tag
                    comm                    // communicator
                );
            }

            {
                UIPstream::read
                (
                    UPstream::commsTypes::blocking,
                    UPstream::masterNo(),
                    reinterpret_cast<char*>(&sum),
                    sizeof(scalar),
                    UPstream::msgType(),    // tag
                    comm                    // communicator
                );
            }
        }
    }
    return sum;
}


int main(int argc, char *argv[])
{

    #include "setRootCase.H"
    #include "createTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Allocate a communicator
    label n = Pstream::nProcs(UPstream::worldComm);

    DynamicList<label> bottom;
    DynamicList<label> top;

    for (label i = 0; i < n/2; i++)
    {
        bottom.append(i);
    }
    for (label i = n/2; i < n; i++)
    {
        top.append(i);
    }

    // Pout<< "bottom:" << bottom << endl;
    Pout<< "top             :" << top << endl;


    scalar localValue = 111*UPstream::myProcNo(UPstream::worldComm);
    Pout<< "localValue      :" << localValue << endl;


    label comm = Pstream::allocateCommunicator
    (
        UPstream::worldComm,
        top
    );

    Pout<< "allocated comm  :" << comm << endl;
    Pout<< "comm myproc     :" << Pstream::myProcNo(comm)
        << endl;


    if (Pstream::myProcNo(comm) != -1)
    {
        // scalar sum = sumReduce(comm, localValue);
        // scalar sum = localValue;
        // reduce
        //(
        //    UPstream::treeCommunication(comm),
        //    sum,
        //    sumOp<scalar>(),
        //    Pstream::msgType(),
        //    comm
        //);
        scalar sum = returnReduce
        (
            localValue,
            sumOp<scalar>(),
            Pstream::msgType(),
            comm
        );
        Pout<< "sum             :" << sum << endl;
    }

    Pstream::freeCommunicator(comm);


    Pout<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
