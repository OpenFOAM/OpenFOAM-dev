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

\*---------------------------------------------------------------------------*/

#include "commSchedule.H"
#include "SortableList.H"
#include "boolList.H"
#include "IOstreams.H"
#include "IOmanip.H"
#include "OStringStream.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(commSchedule, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::commSchedule::outstandingComms
(
    const labelList& commToSchedule,
    DynamicList<label>& procComms
) const
{
    label nOutstanding = 0;

    forAll(procComms, i)
    {
        if (commToSchedule[procComms[i]] == -1)
        {
            nOutstanding++;
        }
    }
    return nOutstanding;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from separate addressing
Foam::commSchedule::commSchedule
(
    const label nProcs,
    const List<labelPair>& comms
)
:
    schedule_(comms.size()),
    procSchedule_(nProcs)
{
    // Determine comms per processor.
    List<DynamicList<label> > procToComms(nProcs);

    forAll(comms, commI)
    {
        label proc0 = comms[commI][0];
        label proc1 = comms[commI][1];

        if (proc0 < 0 || proc0 >= nProcs || proc1 < 0 || proc1 >= nProcs)
        {
            FatalErrorIn
            (
                "commSchedule::commSchedule"
                "(const label, const List<labelPair>&)"
            )   << "Illegal processor " << comms[commI] << abort(FatalError);
        }

        procToComms[proc0].append(commI);
        procToComms[proc1].append(commI);
    }
    // Note: no need to shrink procToComms. Are small.

    if (debug && Pstream::master())
    {
        Pout<< "commSchedule::commSchedule : Wanted communication:" << endl;

        forAll(comms, i)
        {
            const labelPair& twoProcs = comms[i];

            Pout<< i << ": "
                << twoProcs[0] << " with " << twoProcs[1] << endl;
        }
        Pout<< endl;


        Pout<< "commSchedule::commSchedule : Schedule:" << endl;

        // Print header. Use buffered output to prevent parallel output messing
        // up.
        {
            OStringStream os;
            os  << "iter|";
            for (int i = 0; i < nProcs; i++)
            {
                os  << setw(3) << i;
            }
            Pout<< os.str().c_str() << endl;
        }
        {
            OStringStream os;
            os  << "----+";
            for (int i = 0; i < nProcs; i++)
            {
                os  << "---";
            }
            Pout<< os.str().c_str() << endl;
        }
    }

    // Schedule all. Note: crap scheduler. Assumes all communication takes
    // equally long.

    label nScheduled = 0;

    label iter = 0;

    // Per index into comms the time when it was scheduled
    labelList commToSchedule(comms.size(), -1);

    while (nScheduled < comms.size())
    {
        label oldNScheduled = nScheduled;

        // Find unscheduled comms. This is the comms where the two processors
        // still have the most unscheduled comms.

        boolList busy(nProcs, false);

        while (true)
        {
            label maxCommI = -1;
            label maxNeed = labelMin;

            forAll(comms, commI)
            {
                label proc0 = comms[commI][0];
                label proc1 = comms[commI][1];

                if
                (
                    commToSchedule[commI] == -1             // unscheduled
                && !busy[proc0]
                && !busy[proc1]
                )
                {
                    label need =
                        outstandingComms(commToSchedule, procToComms[proc0])
                      + outstandingComms(commToSchedule, procToComms[proc1]);

                    if (need > maxNeed)
                    {
                        maxNeed = need;
                        maxCommI = commI;
                    }
                }
            }


            if (maxCommI == -1)
            {
                // Found no unscheduled procs.
                break;
            }

            // Schedule commI in this iteration
            commToSchedule[maxCommI] = nScheduled++;
            busy[comms[maxCommI][0]] = true;
            busy[comms[maxCommI][1]] = true;
        }

        if (debug && Pstream::master())
        {
            label nIterComms = nScheduled-oldNScheduled;

            if (nIterComms > 0)
            {
                labelList procToComm(nProcs, -1);

                forAll(commToSchedule, commI)
                {
                    label sched = commToSchedule[commI];

                    if (sched >= oldNScheduled && sched < nScheduled)
                    {
                        label proc0 = comms[commI][0];
                        procToComm[proc0] = commI;
                        label proc1 = comms[commI][1];
                        procToComm[proc1] = commI;
                    }
                }

                // Print it
                OStringStream os;
                os  << setw(3) << iter << " |";
                forAll(procToComm, procI)
                {
                    if (procToComm[procI] == -1)
                    {
                        os  << "   ";
                    }
                    else
                    {
                        os  << setw(3) << procToComm[procI];
                    }
                }
                Pout<< os.str().c_str() << endl;
            }
        }

        iter++;
    }

    if (debug && Pstream::master())
    {
        Pout<< endl;
    }


    // Sort commToSchedule and obtain order in comms
    schedule_ = SortableList<label>(commToSchedule).indices();

    // Sort schedule_ by processor

    labelList nProcScheduled(nProcs, 0);

    // Count
    forAll(schedule_, i)
    {
        label commI = schedule_[i];
        const labelPair& twoProcs = comms[commI];

        nProcScheduled[twoProcs[0]]++;
        nProcScheduled[twoProcs[1]]++;
    }
    // Allocate
    forAll(procSchedule_, procI)
    {
        procSchedule_[procI].setSize(nProcScheduled[procI]);
    }
    nProcScheduled = 0;
    // Fill
    forAll(schedule_, i)
    {
        label commI = schedule_[i];
        const labelPair& twoProcs = comms[commI];

        label proc0 = twoProcs[0];
        procSchedule_[proc0][nProcScheduled[proc0]++] = commI;

        label proc1 = twoProcs[1];
        procSchedule_[proc1][nProcScheduled[proc1]++] = commI;
    }

    if (debug && Pstream::master())
    {
        Pout<< "commSchedule::commSchedule : Per processor:" << endl;

        forAll(procSchedule_, procI)
        {
            const labelList& procComms = procSchedule_[procI];

            Pout<< "Processor " << procI << " talks to processors:" << endl;

            forAll(procComms, i)
            {
                const labelPair& twoProcs = comms[procComms[i]];

                label nbr = (twoProcs[1] == procI ? twoProcs[0] : twoProcs[1]);

                Pout<< "    " << nbr << endl;
            }
        }
        Pout<< endl;
    }
}


// ************************************************************************* //
