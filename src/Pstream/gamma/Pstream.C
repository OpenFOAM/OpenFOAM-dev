/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Pstream for GAMMA

    GAMMA has a (polling) receive handler which gets called every time a
    received message is complete. Ours stores the length of the currently
    received message and sets up the next buffer to store the next message
    in.
    Note that the pattern between two processors can be
    - send
    - receive
    - receive
    - send
    since the first swap might belong to a local exchange and the second to
    a reduce. Since gamma has to have the receive buffers already set up we
    have to allocate them big enough. To prevent excessive amounts needed we
    dynamically resize them (never shrink) by sending special 'resize' messages
    before sending a largish message.

    Because of this we actually need four receive buffers:
    - send
    - receive resize message
    - receive normal message
    - receive resize message
    - receive normal message
    - send

    The special resize message is a message with a special header which
    (hopefully) should never appear in normal exchanges (it actually checks
    for this in the OPstream::send)

\*---------------------------------------------------------------------------*/

#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "OSspecific.H"
#include "PstreamGlobals.H"

#include <cstring>
#include <cstdlib>
#include <csignal>

extern "C"
{
#   include <linux/gamma/libgamma.h>
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Receive handler to copy out received message length and switch buffers.
static void handler(void)
{
    label current = PstreamGlobals::recvIndex[gamma_active_port];

    List<char>& buf = PstreamGlobals::recvBuf[current][gamma_active_port];
    label bufLen = PstreamGlobals::recvBufLen[current][gamma_active_port];

    if (bufLen != -1)
    {
        FatalErrorIn("Pstream::handler(void)")
            << "Buffer length not reset : "
            << bufLen
            << " when receiving message of size " << gamma_msglen
            << " from processor " << gamma_active_port << endl
            << "This means that the existing data has not been consumed yet"
            << " (by IPstream::read) and means your communication pattern"
            << " is probably not balanced (a receive for every send)"
            << endl
            << "This can happen if you have e.g. gather without scatter."
            << endl
            << "A workaround is to increase the depth of the circular"
            << " receive buffers in PstreamGlobals.H"
            << abort(FatalError);
    }


    // Some checks
    if
    (
        gamma_msglen < 0
     || gamma_msglen > buf.size()
    )
    {
        FatalErrorIn("Pstream::handler(void)")
            << "Received message of size " << gamma_msglen
            << " from processor " << gamma_active_port
            << Foam::endl
            << "but global receive buffer is only of size "
            << buf.size()
            << abort(FatalError);
    }

    // Check for resize message
    label resizeLen = PstreamGlobals::getSizeFromHeader
    (
        buf.begin(),
        gamma_msglen
    );

    if (resizeLen != -1)
    {
        if (Pstream::debug)
        {
            Pout<< "Pstream::handler : Resize message:" << resizeLen
                << " from proc " << gamma_active_port
                << " current size:"
                << PstreamGlobals::getMaxBufSize(gamma_active_port)
                << Foam::endl;
        }

        // Saved current buffer.
        List<char> savedBuf;

        if (resizeLen > PstreamGlobals::getMaxBufSize(gamma_active_port))
        {
            if (Pstream::debug)
            {
                Pout<< "Pstream::handler :"
                    << " resizing receive buffer for processor "
                    << gamma_active_port
                    << " from "
                    << PstreamGlobals::getMaxBufSize(gamma_active_port)
                    << " to " << resizeLen << Foam::endl;
            }

            // Save the pointer (that gamma knows about) so we can safely
            // gamma_switch_to_buffer with a valid pointer.
            // Not sure if necessary but do anyway.
            savedBuf.transfer(buf);

            // Resize all the buffers
            forAll(PstreamGlobals::recvBuf, i)
            {
                List<char>& chars =
                    PstreamGlobals::recvBuf[i][gamma_active_port];

//                gamma_munlock(chars.begin(), chars.size());
                chars.setSize(resizeLen);
//                gamma_mlock(chars.begin(), chars.size());
            }
        }

        // Update length with special value to denote resize was done.
        PstreamGlobals::recvBufLen[current][gamma_active_port] = -2;
    }
    else
    {
        // Update length with actual message length
        PstreamGlobals::recvBufLen[current][gamma_active_port] = gamma_msglen;
    }

    // Go to next buffer.
    label next = PstreamGlobals::recvBuf.fcIndex(current);
    PstreamGlobals::recvIndex[gamma_active_port] = next;

//    gamma_switch_to_buffer
    gamma_post_recv
    (
        gamma_active_port,
        PstreamGlobals::recvBuf[next][gamma_active_port].begin(),
        PstreamGlobals::recvBuf[next][gamma_active_port].size()
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Pstream::addValidParOptions(HashTable<string>& validParOptions)
{
    validParOptions.insert("np", "");
    validParOptions.insert("p4pg", "PI file");
    validParOptions.insert("p4wd", "directory");
    validParOptions.insert("p4amslave", "");
    validParOptions.insert("p4yourname", "hostname");

    validParOptions.insert("machinefile", "machine file");
    validParOptions.insert("GAMMANP", "numProcs");
    validParOptions.insert("GAMMAHOME", "gamma cwd");
    validParOptions.insert("GAMMA", "1(enable) or 0(disable)");
}


bool Pstream::init(int& argc, char**& argv)
{
    int numprocs = 0;

    string npString("-GAMMANP");

    for (label i = 0; i < argc; i++)
    {
        if (argv[i] == npString)
        {
            if (i+1 < argc)
            {
                numprocs = atoi(argv[i+1]);
                break;
            }
        }
    }

    // Initialize GAMMA
    unsigned char smallNumprocs = numprocs;

    gamma_init(smallNumprocs, argc, argv);

    myProcNo_ = gamma_my_node();

    // Make sure printing with prefix.
    setParRun();

    procIDs_.setSize(numprocs);

    forAll(procIDs_, procNo)
    {
        procIDs_[procNo] = procNo;
    }


    // Allocate receive buffers.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    // Make sure each receive buffer is at least large enough to receive
    // the resize message.

    // Current active buffer
    PstreamGlobals::recvIndex.setSize(numprocs);
    PstreamGlobals::recvIndex = 0;
    PstreamGlobals::consumeIndex.setSize(numprocs);
    PstreamGlobals::consumeIndex = 0;

    forAll(PstreamGlobals::recvBuf, i)
    {
        PstreamGlobals::recvBufLen[i].setSize(numprocs);
        PstreamGlobals::recvBufLen[i] = -1;

        List<List<char> >& buffers = PstreamGlobals::recvBuf[i];

        buffers.setSize(numprocs);
        forAll(buffers, procNo)
        {
            if (procNo != myProcNo_)
            {
                buffers[procNo].setSize(PstreamGlobals::initialBufferLen);

                // Acc. to gamma sources all buffers need to be in memory.
                // Either locked or "write touched".
//                gamma_mlock
 //               (
  //                  buffers[procNo].begin(),
   //                 buffers[procNo].size()
    //            );
            }
        }
    }


    // Lock the special resize message
    //    gamma_mlock
    //    (
    //       reinterpret_cast<char*>(PstreamGlobals::resizeMessage),
    //      PstreamGlobals::resizeMessageLen*sizeof(uint64_t)
    // );


    // Attach current receive buffers
    forAll(procIDs_, procNo)
    {
        if (procNo != myProcNo_)
        {
            // Buffer index (always 0 at this point)
            label current = PstreamGlobals::recvIndex[procNo];

            // Current buffer for this processor.
            List<char>& buf = PstreamGlobals::recvBuf[current][procNo];

            gamma_set_active_port
            (
                procNo,             //unsigned short port,
                procNo,             //unsigned short dest_node,
                gamma_my_par_pid(), //unsigned char dest_par_pid,
                myProcNo_,          //unsigned short dest_port,
                handler,            //callback
                procNo,             //unsigned short semaphore,
                GO_BACK,            //unsigned char buffer_kind,
                buf.begin(),
                buf.size()
            );
        }
    }


    // Make sure all have allocated the ports (so set the receive buffers)
    gamma_sync();

    Info<< "GAMMA Pstream initialized with:" << nl
        << "    floatTransfer         : " << floatTransfer << nl
        << "    nProcsSimpleSum       : " << nProcsSimpleSum << nl
        << "    scheduledTransfer     : " << Pstream::scheduledTransfer << nl
        << Foam::endl;

    // Now that nprocs is known construct communication tables.
    initCommunicationSchedule();

    return true;
}


void Pstream::exit(int errnum)
{
    //    gamma_munlockall();
    gamma_exit();
    //gamma_abort();
}


void Pstream::abort()
{
    Pout<< "**Pstream::abort()**" << endl;
    // gamma_munlockall();
    gamma_abort();
}


void reduce(scalar& Value, const sumOp<scalar>& bop)
{
    if (!Pstream::parRun())
    {
        return;
    }

    if (Pstream::debug)
    {
        Pout<< "**entering Pstream::reduce for " << Value << Foam::endl;
    }


    if (Pstream::master())
    {
        for
        (
            int slave=Pstream::firstSlave();
            slave<=Pstream::lastSlave();
            slave++
        )
        {
            scalar value;

            if
            (
               !IPstream::read
                (
                    slave,
                    reinterpret_cast<char*>(&value),    // buf
                    sizeof(Value)                       // bufSize
                )
            )
            {
                FatalErrorIn
                (
                    "reduce(scalar& Value, const sumOp<scalar>& sumOp)"
                )   << "IPstream::read failed"
                    << Foam::abort(FatalError);
            }

            Value = bop(Value, value);
        }
    }
    else
    {
        if
        (
           !OPstream::write
            (
                Pstream::masterNo(),
                reinterpret_cast<const char*>(&Value),  // buf
                sizeof(Value),                          // bufSize
                false                                   // non-buffered
            )
        )
        {
            FatalErrorIn
            (
                "reduce(scalar& Value, const sumOp<scalar>& sumOp)"
            )   << "OPstream::write failed"
                << Foam::abort(FatalError);
        }
    }

    if (Pstream::master())
    {
        for
        (
            int slave=Pstream::firstSlave();
            slave<=Pstream::lastSlave();
            slave++
        )
        {
            if
            (
               !OPstream::write
                (
                    slave,
                    reinterpret_cast<const char*>(&Value),  // buf
                    sizeof(Value),                          // bufSize,
                    false                                   // non-buffered
                )
            )
            {
                FatalErrorIn
                (
                    "reduce(scalar& Value, const sumOp<scalar>& sumOp)"
                )   << "OPstream::write failed"
                    << Foam::abort(FatalError);
            }
        }
    }
    else
    {
        if
        (
           !IPstream::read
            (
                Pstream::masterNo(),
                reinterpret_cast<char*>(&Value),    // buf
                sizeof(Value)                       // bufSize
            )
        )
        {
            FatalErrorIn
            (
                "reduce(scalar& Value, const sumOp<scalar>& sumOp)"
            )   << "IPstream::read failed"
                << Foam::abort(FatalError);
        }
    }

    if (Pstream::debug)
    {
        Pout<< "**exiting Pstream::reduce with " << Value << Foam::endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
