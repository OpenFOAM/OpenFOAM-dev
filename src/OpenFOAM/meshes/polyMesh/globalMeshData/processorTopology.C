/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "processorTopology.H"
#include "polyBoundaryMesh.H"
#include "processorPolyPatch.H"
#include "commSchedule.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::processorTopology::procNeighbours
(
    const label nProcs,
    const polyBoundaryMesh& patches
)
{
    // Determine number of processor neighbours and max neighbour id.

    label nNeighbours = 0;

    label maxNb = 0;

    boolList isNeighbourProc(nProcs, false);

    forAll(patches, patchi)
    {
        const polyPatch& patch = patches[patchi];

        if (isA<processorPolyPatch>(patch))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patch);

            label pNeighbProcNo = procPatch.neighbProcNo();

            if (!isNeighbourProc[pNeighbProcNo])
            {
                nNeighbours++;

                maxNb = max(maxNb, procPatch.neighbProcNo());

                isNeighbourProc[pNeighbProcNo] = true;
            }
        }
    }

    labelList neighbours(nNeighbours, -1);

    nNeighbours = 0;

    forAll(isNeighbourProc, proci)
    {
        if (isNeighbourProc[proci])
        {
            neighbours[nNeighbours++] = proci;
        }
    }

    procPatchMap_.setSize(maxNb + 1);
    procPatchMap_ = -1;

    forAll(patches, patchi)
    {
        const polyPatch& patch = patches[patchi];

        if (isA<processorPolyPatch>(patch))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patch);

            // Construct reverse map
            procPatchMap_[procPatch.neighbProcNo()] = patchi;
        }
    }

    return neighbours;
}


Foam::lduSchedule Foam::processorTopology::nonBlockingSchedule
(
    const polyBoundaryMesh& patches
)
{
    lduSchedule patchSchedule(2*patches.size());

    label patchEvali = 0;

    // 1. All non-processor patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Have evaluate directly after initEvaluate. Could have them separated
    // as long as they're not intermingled with processor patches since
    // then e.g. any reduce parallel traffic would interfere with the
    // processor swaps.

    forAll(patches, patchi)
    {
        if (!isA<processorPolyPatch>(patches[patchi]))
        {
            patchSchedule[patchEvali].patch = patchi;
            patchSchedule[patchEvali++].init = true;
            patchSchedule[patchEvali].patch = patchi;
            patchSchedule[patchEvali++].init = false;
        }
    }

    // 2. All processor patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    // 2a. initEvaluate
    forAll(patches, patchi)
    {
        if (isA<processorPolyPatch>(patches[patchi]))
        {
            patchSchedule[patchEvali].patch = patchi;
            patchSchedule[patchEvali++].init = true;
        }
    }

    // 2b. evaluate
    forAll(patches, patchi)
    {
        if (isA<processorPolyPatch>(patches[patchi]))
        {
            patchSchedule[patchEvali].patch = patchi;
            patchSchedule[patchEvali++].init = false;
        }
    }

    return patchSchedule;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorTopology::processorTopology
(
    const polyBoundaryMesh& patches,
    const label comm
)
:
    procNbrProcs_(Pstream::nProcs(comm)),
    procPatchMap_(),
    patchSchedule_(2*patches.size())
{
    if (Pstream::parRun())
    {
        // Fill my 'slot' with my neighbours
        procNbrProcs_[Pstream::myProcNo(comm)] =
            procNeighbours(procNbrProcs_.size(), patches);

        // Distribute to all processors
        Pstream::gatherList(procNbrProcs_, Pstream::msgType(), comm);
        Pstream::scatterList(procNbrProcs_, Pstream::msgType(), comm);
    }

    if
    (
        Pstream::parRun()
     && Pstream::defaultCommsType == Pstream::commsTypes::scheduled
    )
    {
        label patchEvali = 0;

        // 1. All non-processor patches
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        forAll(patches, patchi)
        {
            if (!isA<processorPolyPatch>(patches[patchi]))
            {
                patchSchedule_[patchEvali].patch = patchi;
                patchSchedule_[patchEvali++].init = true;
                patchSchedule_[patchEvali].patch = patchi;
                patchSchedule_[patchEvali++].init = false;
            }
        }

        // 2. All processor patches
        // ~~~~~~~~~~~~~~~~~~~~~~~~

        // Determine the schedule for all. Insert processor pair once
        // to determine the schedule. Each processor pair stands for both
        // send and receive.
        label nComms = 0;
        forAll(procNbrProcs_, proci)
        {
            nComms += procNbrProcs_[proci].size();
        }
        DynamicList<labelPair> comms(nComms);

        forAll(procNbrProcs_, proci)
        {
            const labelList& nbrs = procNbrProcs_[proci];

            forAll(nbrs, i)
            {
                if (proci < nbrs[i])
                {
                    comms.append(labelPair(proci, nbrs[i]));
                }
            }
        }
        comms.shrink();

        // Determine a schedule.
        labelList mySchedule
        (
            commSchedule
            (
                Pstream::nProcs(comm),
                comms
            ).procSchedule()[Pstream::myProcNo(comm)]
        );

        forAll(mySchedule, iter)
        {
            label commI = mySchedule[iter];

            // Get the other processor
            label nb = comms[commI][0];
            if (nb == Pstream::myProcNo(comm))
            {
                nb = comms[commI][1];
            }
            label patchi = procPatchMap_[nb];

            if (Pstream::myProcNo(comm) > nb)
            {
                patchSchedule_[patchEvali].patch = patchi;
                patchSchedule_[patchEvali++].init = true;
                patchSchedule_[patchEvali].patch = patchi;
                patchSchedule_[patchEvali++].init = false;
            }
            else
            {
                patchSchedule_[patchEvali].patch = patchi;
                patchSchedule_[patchEvali++].init = false;
                patchSchedule_[patchEvali].patch = patchi;
                patchSchedule_[patchEvali++].init = true;
            }
        }
    }
    else
    {
        patchSchedule_ = nonBlockingSchedule(patches);
    }
}


// ************************************************************************* //
