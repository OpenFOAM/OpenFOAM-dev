/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "regionSplitBase.H"
#include "Map.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionSplitBase, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::globalIndex>
Foam::regionSplitBase::localGlobalIndex(const label localSize)
{
    labelList offsets(Pstream::nProcs() + 1, 0);

    for (label i = Pstream::myProcNo() + 1; i <= Pstream::nProcs(); i ++)
    {
        offsets[i] = localSize;
    }

    return autoPtr<globalIndex>(new globalIndex(move(offsets)));
}


Foam::label Foam::regionSplitBase::compactLocalRegionSplit
(
    labelList& cellRegion
)
{
    Map<label> nonCompactToCompactRegion(cellRegion.size()/8);

    forAll(cellRegion, celli)
    {
        const label nonCompactRegioni = cellRegion[celli];

        Map<label>::const_iterator iter =
            nonCompactToCompactRegion.find(nonCompactRegioni);

        if (iter == nonCompactToCompactRegion.end())
        {
            // First time encountering this non-compact region. Assign it a new
            // compact region index.
            const label compactRegioni = nonCompactToCompactRegion.size();
            cellRegion[celli] = compactRegioni;
            nonCompactToCompactRegion.insert(nonCompactRegioni, compactRegioni);
        }
        else
        {
            cellRegion[celli] = iter();
        }
    }

    return nonCompactToCompactRegion.size();
}


Foam::label Foam::regionSplitBase::compactGlobalRegionSplit
(
    const globalIndex& globalNonCompactRegions,
    labelList& cellRegion
)
{
    // Count per originating processor the number of regions
    labelList nOriginating(Pstream::nProcs(), 0);
    {
        labelHashSet haveRegion(cellRegion.size()/8);

        forAll(cellRegion, celli)
        {
            const label region = cellRegion[celli];

            // Count originating processor. Use isLocal as efficiency since
            // most cells are locally originating.
            if (globalNonCompactRegions.isLocal(region))
            {
                if (haveRegion.insert(region))
                {
                    nOriginating[Pstream::myProcNo()] ++;
                }
            }
            else
            {
                const label proci = globalNonCompactRegions.whichProcID(region);
                if (haveRegion.insert(region))
                {
                    nOriginating[proci] ++;
                }
            }
        }
    }

    // Global numbering for compacted local regions
    const globalIndex globalCompact(nOriginating[Pstream::myProcNo()]);

    // Renumber into compact indices. Note that since we've already made
    // all regions global we now need a Map to store the compacting information
    // instead of a labelList - otherwise we could have used a straight
    // labelList.

    // Local compaction map
    Map<label> globalToCompact(2*nOriginating[Pstream::myProcNo()]);

    // Remote regions we want the compact number for
    List<labelHashSet> nonLocal(Pstream::nProcs());
    forAll(nonLocal, proci)
    {
        if (proci != Pstream::myProcNo())
        {
            nonLocal[proci].resize(2*nOriginating[proci]);
        }
    }

    forAll(cellRegion, celli)
    {
        label region = cellRegion[celli];
        if (globalNonCompactRegions.isLocal(region))
        {
            // Insert new compact region (if not yet present)
            globalToCompact.insert
            (
                region,
                globalCompact.toGlobal(globalToCompact.size())
            );
        }
        else
        {
            const label proci = globalNonCompactRegions.whichProcID(region);
            nonLocal[proci].insert(region);
        }
    }

    // Now we have all the local regions compacted. Now we need to get the
    // non-local ones from the processors to whom they are local.
    // Convert the nonLocal (labelHashSets) to labelLists.
    labelListList sendNonLocal(Pstream::nProcs());
    forAll(sendNonLocal, proci)
    {
        sendNonLocal[proci] = nonLocal[proci].toc();
    }

    // Get the wanted region labels into recvNonLocal
    labelListList recvNonLocal;
    Pstream::exchange<labelList, label>(sendNonLocal, recvNonLocal);

    // Now we have the wanted compact region labels that proci wants in
    // recvNonLocal[proci]. Construct corresponding list of compact
    // region labels to send back.
    labelListList sendWantedLocal(Pstream::nProcs());
    forAll(recvNonLocal, proci)
    {
        const labelList& nonLocal = recvNonLocal[proci];
        sendWantedLocal[proci].setSize(nonLocal.size());

        forAll(nonLocal, i)
        {
            sendWantedLocal[proci][i] = globalToCompact[nonLocal[i]];
        }
    }

    // Send back (into recvNonLocal)
    recvNonLocal.clear();
    Pstream::exchange<labelList, label>(sendWantedLocal, recvNonLocal);
    sendWantedLocal.clear();

    // Now recvNonLocal contains for every element in setNonLocal the
    // corresponding compact number. Insert these into the local compaction
    // map.
    forAll(recvNonLocal, proci)
    {
        const labelList& wantedRegions = sendNonLocal[proci];
        const labelList& compactRegions = recvNonLocal[proci];

        forAll(wantedRegions, i)
        {
            globalToCompact.insert(wantedRegions[i], compactRegions[i]);
        }
    }

    // Finally renumber the regions
    forAll(cellRegion, celli)
    {
        cellRegion[celli] = globalToCompact[cellRegion[celli]];
    }

    return globalCompact.size();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSplitBase::regionSplitBase(const label size)
:
    labelList(size, -1),
    nRegions_(-1)
{}


// ************************************************************************* //
