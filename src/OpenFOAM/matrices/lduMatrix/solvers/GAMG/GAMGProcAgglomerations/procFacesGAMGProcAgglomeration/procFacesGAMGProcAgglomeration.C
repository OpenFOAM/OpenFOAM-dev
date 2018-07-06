/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

\*---------------------------------------------------------------------------*/

#include "procFacesGAMGProcAgglomeration.H"
#include "addToRunTimeSelectionTable.H"
#include "GAMGAgglomeration.H"
#include "Random.H"
#include "lduMesh.H"
#include "processorLduInterface.H"
#include "processorGAMGInterface.H"
#include "pairGAMGAgglomeration.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(procFacesGAMGProcAgglomeration, 0);

    addToRunTimeSelectionTable
    (
        GAMGProcAgglomeration,
        procFacesGAMGProcAgglomeration,
        GAMGAgglomeration
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::lduPrimitiveMesh>
Foam::procFacesGAMGProcAgglomeration::singleCellMesh
(
    const label singleCellMeshComm,
    const lduMesh& mesh,
    scalarField& faceWeights
) const
{
    // Count number of faces per processor
    List<Map<label>> procFaces(UPstream::nProcs(mesh.comm()));
    Map<label>& myNeighbours = procFaces[UPstream::myProcNo(mesh.comm())];

    {
        const lduInterfacePtrsList interfaces(mesh.interfaces());
        forAll(interfaces, intI)
        {
            if (interfaces.set(intI))
            {
                const processorLduInterface& pp =
                    refCast<const processorLduInterface>
                    (
                        interfaces[intI]
                    );
                label size = interfaces[intI].faceCells().size();
                myNeighbours.insert(pp.neighbProcNo(), size);
            }
        }
    }

    Pstream::gatherList(procFaces, Pstream::msgType(), mesh.comm());
    Pstream::scatterList(procFaces, Pstream::msgType(), mesh.comm());

    autoPtr<lduPrimitiveMesh> singleCellMeshPtr;

    if (Pstream::master(mesh.comm()))
    {
        // I am master
        label nCells = Pstream::nProcs(mesh.comm());

        DynamicList<label> l(3*nCells);
        DynamicList<label> u(3*nCells);
        DynamicList<scalar> weight(3*nCells);

        DynamicList<label> nbrs;
        DynamicList<scalar> weights;

        forAll(procFaces, proci)
        {
            const Map<label>& neighbours = procFaces[proci];

            // Add all the higher processors
            nbrs.clear();
            weights.clear();
            forAllConstIter(Map<label>, neighbours, iter)
            {
                if (iter.key() > proci)
                {
                    nbrs.append(iter.key());
                    weights.append(iter());
                }
                sort(nbrs);
                forAll(nbrs, i)
                {
                    l.append(proci);
                    u.append(nbrs[i]);
                    weight.append(weights[i]);
                }
            }
        }

        faceWeights.transfer(weight);

        PtrList<const lduInterface> primitiveInterfaces(0);
        const lduSchedule ps(0);

        singleCellMeshPtr.reset
        (
            new lduPrimitiveMesh
            (
                nCells,
                l,
                u,
                primitiveInterfaces,
                ps,
                singleCellMeshComm
            )
        );
    }
    return singleCellMeshPtr;
}


Foam::tmp<Foam::labelField>
Foam::procFacesGAMGProcAgglomeration::processorAgglomeration
(
    const lduMesh& mesh
) const
{
    label singleCellMeshComm = UPstream::allocateCommunicator
    (
        mesh.comm(),
        labelList(1, label(0))            // only processor 0
    );

    scalarField faceWeights;
    autoPtr<lduPrimitiveMesh> singleCellMeshPtr
    (
        singleCellMesh
        (
            singleCellMeshComm,
            mesh,
            faceWeights
        )
    );

    tmp<labelField> tfineToCoarse(new labelField(0));
    labelField& fineToCoarse = tfineToCoarse.ref();

    if (singleCellMeshPtr.valid())
    {
        // On master call the agglomerator
        const lduPrimitiveMesh& singleCellMesh = singleCellMeshPtr();

        label nCoarseProcs;
        fineToCoarse = pairGAMGAgglomeration::agglomerate
        (
            nCoarseProcs,
            singleCellMesh,
            faceWeights
        );

        labelList coarseToMaster(nCoarseProcs, labelMax);
        forAll(fineToCoarse, celli)
        {
            label coarseI = fineToCoarse[celli];
            coarseToMaster[coarseI] = min(coarseToMaster[coarseI], celli);
        }

        // Sort according to master and redo restriction
        labelList newToOld;
        sortedOrder(coarseToMaster, newToOld);
        labelList oldToNew(invert(newToOld.size(), newToOld));

        fineToCoarse = UIndirectList<label>(oldToNew, fineToCoarse)();
    }

    Pstream::scatter(fineToCoarse, Pstream::msgType(), mesh.comm());
    UPstream::freeCommunicator(singleCellMeshComm);

    return tfineToCoarse;
}


bool Foam::procFacesGAMGProcAgglomeration::doProcessorAgglomeration
(
    const lduMesh& mesh
) const
{
    // Check the need for further agglomeration on all processors
    bool doAgg = mesh.lduAddr().size() < nAgglomeratingCells_;
    mesh.reduce(doAgg, orOp<bool>());
    return doAgg;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::procFacesGAMGProcAgglomeration::procFacesGAMGProcAgglomeration
(
    GAMGAgglomeration& agglom,
    const dictionary& controlDict
)
:
    GAMGProcAgglomeration(agglom, controlDict),
    nAgglomeratingCells_(readLabel(controlDict.lookup("nAgglomeratingCells")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::procFacesGAMGProcAgglomeration::~procFacesGAMGProcAgglomeration()
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

bool Foam::procFacesGAMGProcAgglomeration::agglomerate()
{
    if (debug)
    {
        Pout<< nl << "Starting mesh overview" << endl;
        printStats(Pout, agglom_);
    }

    if (agglom_.size() >= 1)
    {
        Random rndGen(0);

        for
        (
            label fineLevelIndex = 2;
            fineLevelIndex < agglom_.size();
            fineLevelIndex++
        )
        {
            if (agglom_.hasMeshLevel(fineLevelIndex))
            {
                // Get the fine mesh
                const lduMesh& levelMesh = agglom_.meshLevel(fineLevelIndex);

                label levelComm = levelMesh.comm();
                label nProcs = UPstream::nProcs(levelComm);

                if (nProcs > 1 && doProcessorAgglomeration(levelMesh))
                {
                    tmp<labelField> tprocAgglomMap
                    (
                        processorAgglomeration(levelMesh)
                    );
                    const labelField& procAgglomMap = tprocAgglomMap();

                    // Master processor
                    labelList masterProcs;
                    // Local processors that agglomerate. agglomProcIDs[0] is in
                    // masterProc.
                    List<label> agglomProcIDs;
                    GAMGAgglomeration::calculateRegionMaster
                    (
                        levelComm,
                        procAgglomMap,
                        masterProcs,
                        agglomProcIDs
                    );

                    // Allocate a communicator for the processor-agglomerated
                    // matrix
                    comms_.append
                    (
                        UPstream::allocateCommunicator
                        (
                            levelComm,
                            masterProcs
                        )
                    );


                    // Use procesor agglomeration maps to do the actual
                    // collecting.
                    GAMGProcAgglomeration::agglomerate
                    (
                        fineLevelIndex,
                        procAgglomMap,
                        masterProcs,
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

    return true;
}


// ************************************************************************* //
