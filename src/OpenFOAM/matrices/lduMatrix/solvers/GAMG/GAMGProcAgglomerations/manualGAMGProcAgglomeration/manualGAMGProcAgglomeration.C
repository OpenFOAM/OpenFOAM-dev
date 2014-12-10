/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "manualGAMGProcAgglomeration.H"
#include "addToRunTimeSelectionTable.H"
#include "GAMGAgglomeration.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(manualGAMGProcAgglomeration, 0);

    addToRunTimeSelectionTable
    (
        GAMGProcAgglomeration,
        manualGAMGProcAgglomeration,
        GAMGAgglomeration
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::manualGAMGProcAgglomeration::manualGAMGProcAgglomeration
(
    GAMGAgglomeration& agglom,
    const dictionary& controlDict
)
:
    GAMGProcAgglomeration(agglom, controlDict),
    procAgglomMaps_(controlDict.lookup("processorAgglomeration"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::manualGAMGProcAgglomeration::
~manualGAMGProcAgglomeration()
{
    forAllReverse(comms_, i)
    {
        if (comms_[i] != -1)
        {
            UPstream::freeCommunicator(comms_[i]);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::manualGAMGProcAgglomeration::agglomerate()
{
    if (debug)
    {
        Pout<< nl << "Starting mesh overview" << endl;
        printStats(Pout, agglom_);
    }

    if (agglom_.size() >= 1)
    {
        forAll(procAgglomMaps_, i)
        {
            const label fineLevelIndex = procAgglomMaps_[i].first();

            if (fineLevelIndex >= agglom_.size())
            {
                WarningIn("manualGAMGProcAgglomeration::agglomerate()")
                    << "Ignoring specification for level " << fineLevelIndex
                    << " since outside agglomeration." << endl;

                continue;
            }

            if (agglom_.hasMeshLevel(fineLevelIndex))
            {
                // Get the fine mesh
                const lduMesh& levelMesh = agglom_.meshLevel(fineLevelIndex);
                label nProcs = UPstream::nProcs(levelMesh.comm());

                if (nProcs > 1)
                {
                    // My processor id
                    const label myProcID = Pstream::myProcNo(levelMesh.comm());

                    const List<labelList>& clusters =
                        procAgglomMaps_[i].second();

                    // Coarse to fine master processor
                    labelList coarseToMaster(clusters.size());

                    // Fine to coarse map
                    labelList procAgglomMap(nProcs, -1);

                    // Cluster for my processor (with master index first)
                    labelList agglomProcIDs;



                    forAll(clusters, coarseI)
                    {
                        const labelList& cluster = clusters[coarseI];
                        coarseToMaster[coarseI] = cluster[0];

                        forAll(cluster, i)
                        {
                            procAgglomMap[cluster[i]] = coarseI;
                        }

                        label masterIndex = findIndex
                        (
                            cluster,
                            coarseToMaster[coarseI]
                        );

                        if (masterIndex == -1)
                        {
                            FatalErrorIn
                            (
                                "manualGAMGProcAgglomeration::agglomerate()"
                            )   << "At level " << fineLevelIndex
                                << " the master processor "
                                << coarseToMaster[coarseI]
                                << " is not in the cluster "
                                << cluster
                                << exit(FatalError);
                        }

                        if (findIndex(cluster, myProcID) != -1)
                        {
                            // This is my cluster. Make sure master index is
                            // first
                            agglomProcIDs = cluster;
                            Swap(agglomProcIDs[0], agglomProcIDs[masterIndex]);
                        }
                    }


                    // Check that we've done all processors
                    if (findIndex(procAgglomMap, -1) != -1)
                    {
                        FatalErrorIn
                        (
                            "manualGAMGProcAgglomeration::agglomerate()"
                        )   << "At level " << fineLevelIndex
                            << " processor "
                            << findIndex(procAgglomMap, -1)
                            << " is not in any cluster"
                            << exit(FatalError);
                    }


                    // Allocate a communicator for the processor-agglomerated
                    // matrix
                    comms_.append
                    (
                        UPstream::allocateCommunicator
                        (
                            levelMesh.comm(),
                            coarseToMaster
                        )
                    );

                    // Use procesor agglomeration maps to do the actual
                    // collecting
                    if (Pstream::myProcNo(levelMesh.comm()) != -1)
                    {
                        GAMGProcAgglomeration::agglomerate
                        (
                            fineLevelIndex,
                            procAgglomMap,
                            coarseToMaster,
                            agglomProcIDs,
                            comms_.last()
                        );
                    }
                }
            }
        }

        // Print a bit
        if (debug)
        {
            Pout<< nl << "Agglomerated mesh overview" << endl;
            printStats(Pout, agglom_);
        }
    }

    return true;
}


// ************************************************************************* //
