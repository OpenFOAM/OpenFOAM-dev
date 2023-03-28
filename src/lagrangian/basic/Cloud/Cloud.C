/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "Cloud.H"
#include "processorPolyPatch.H"
#include "globalMeshData.H"
#include "meshToMesh.H"
#include "PstreamCombineReduceOps.H"
#include "polyTopoChangeMap.H"
#include "polyMeshMap.H"
#include "polyDistributionMap.H"
#include "Time.H"
#include "OFstream.H"
#include "wallPolyPatch.H"
#include "nonConformalCyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
struct IDLListAppendEqOp
{
    void operator()(IDLList<Type>& x, const IDLList<Type>& y) const
    {
        if (y.size())
        {
            forAllConstIter(typename IDLList<Type>, y, iter)
            {
                x.append(new Type(iter()));
            }
        }
    }
};

}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class ParticleType>
Foam::labelList Foam::Cloud<ParticleType>::patchNbrProc
(
    const polyMesh& pMesh
)
{
    const polyBoundaryMesh& pbm = pMesh.boundaryMesh();

    labelList result(pbm.size(), -1);

    if (Pstream::parRun())
    {
        forAll(pbm, patchi)
        {
            if (isA<processorPolyPatch>(pbm[patchi]))
            {
                const processorPolyPatch& ppp =
                    refCast<const processorPolyPatch>(pbm[patchi]);

                result[patchi] = ppp.neighbProcNo();
            }
        }
    }

    return result;
}


template<class ParticleType>
Foam::labelList Foam::Cloud<ParticleType>::patchNbrProcPatch
(
    const polyMesh& pMesh
)
{
    const polyBoundaryMesh& pbm = pMesh.boundaryMesh();

    labelList result(pbm.size(), -1);

    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        forAll(pbm, patchi)
        {
            if (isA<processorPolyPatch>(pbm[patchi]))
            {
                const processorPolyPatch& ppp =
                    refCast<const processorPolyPatch>(pbm[patchi]);

                UOPstream(ppp.neighbProcNo(), pBufs)()
                  << ppp.index();
            }
        }

        pBufs.finishedSends();

        forAll(pbm, patchi)
        {
            if (isA<processorPolyPatch>(pbm[patchi]))
            {
                const processorPolyPatch& ppp =
                    refCast<const processorPolyPatch>(pbm[patchi]);

                UIPstream(ppp.neighbProcNo(), pBufs)()
                  >> result[patchi];
            }
        }
    }

    return result;
}


template<class ParticleType>
Foam::labelListList Foam::Cloud<ParticleType>::patchNonConformalCyclicPatches
(
    const polyMesh& pMesh
)
{
    const polyBoundaryMesh& pbm = pMesh.boundaryMesh();

    labelListList result(pbm.size(), labelList());

    forAll(pbm, patchi)
    {
        if (isA<nonConformalCyclicPolyPatch>(pbm[patchi]))
        {
            const nonConformalCyclicPolyPatch& nccPp =
                refCast<const nonConformalCyclicPolyPatch>(pbm[patchi]);

            result[nccPp.origPatchID()].append(patchi);
        }
    }

    return result;
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::storeRays() const
{
    const polyBoundaryMesh& pbm = pMesh_.boundaryMesh();

    forAll(patchNonConformalCyclicPatches_, patchi)
    {
        forAll(patchNonConformalCyclicPatches_[patchi], i)
        {
            const label nccPatchi =
                patchNonConformalCyclicPatches_[patchi][i];

            const nonConformalCyclicPolyPatch& nccPp =
                refCast<const nonConformalCyclicPolyPatch>(pbm[nccPatchi]);

            if (nccPp.owner()) nccPp.rays();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::Cloud<ParticleType>::Cloud
(
    const polyMesh& pMesh,
    const word& cloudName,
    const IDLList<ParticleType>& particles
)
:
    cloud(pMesh, cloudName),
    IDLList<ParticleType>(),
    pMesh_(pMesh),
    patchNbrProc_(patchNbrProc(pMesh)),
    patchNbrProcPatch_(patchNbrProcPatch(pMesh)),
    patchNonConformalCyclicPatches_(patchNonConformalCyclicPatches(pMesh)),
    globalPositionsPtr_(),
    timeIndex_(-1)
{
    // Ask for the tetBasePtIs and oldCellCentres to trigger all processors to
    // build them, otherwise, if some processors have no particles then there
    // is a comms mismatch.
    pMesh_.tetBasePtIs();
    pMesh_.oldCellCentres();

    if (particles.size())
    {
        IDLList<ParticleType>::operator=(particles);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
void Foam::Cloud<ParticleType>::addParticle(ParticleType* pPtr)
{
    this->append(pPtr);
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::deleteParticle(ParticleType& p)
{
    delete(this->remove(&p));
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::deleteLostParticles()
{
    forAllIter(typename Cloud<ParticleType>, *this, pIter)
    {
        ParticleType& p = pIter();

        if (p.cell() == -1)
        {
            WarningInFunction
                << "deleting lost particle at position "
                << p.position(pMesh_) << endl;

            deleteParticle(p);
        }
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::cloudReset(const Cloud<ParticleType>& c)
{
    // Reset particle count and particles only
    // - not changing the cloud object registry or reference to the polyMesh
    ParticleType::particleCount_ = 0;
    IDLList<ParticleType>::operator=(c);
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::changeTimeStep()
{
    forAllIter(typename Cloud<ParticleType>, *this, pIter)
    {
        pIter().reset(0);
    }

    timeIndex_ = pMesh_.time().timeIndex();
}


template<class ParticleType>
template<class TrackCloudType>
void Foam::Cloud<ParticleType>::move
(
    TrackCloudType& cloud,
    typename ParticleType::trackingData& td
)
{
    // If the time has changed, modify the particles accordingly
    if (timeIndex_ != pMesh_.time().timeIndex())
    {
        changeTimeStep();
    }

    // Clear the global positions as these are about to change
    globalPositionsPtr_.clear();

    // Ensure rays are available for non conformal transfers
    storeRays();

    // Create transfer buffers
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // Create lists of particles and patch indices to transfer
    List<IDLList<ParticleType>> sendParticles(Pstream::nProcs());
    List<DynamicList<label>> sendPatchIndices(Pstream::nProcs());

    // While there are particles to transfer
    while (true)
    {
        // Clear the transfer lists
        forAll(sendParticles, proci)
        {
            sendParticles[proci].clear();
            sendPatchIndices[proci].clear();
        }

        // Loop over all particles
        forAllIter(typename Cloud<ParticleType>, *this, pIter)
        {
            ParticleType& p = pIter();

            // Move the particle
            const bool keepParticle = p.move(cloud, td);

            // If the particle is to be kept
            if (keepParticle)
            {
                if (td.sendToProc != -1)
                {
                    #ifdef FULLDEBUG
                    if (!Pstream::parRun() || !p.onBoundaryFace(pMesh_))
                    {
                        FatalErrorInFunction
                            << "Switch processor flag is true when no parallel "
                            << "transfer is possible. This is a bug."
                            << exit(FatalError);
                    }
                    #endif

                    p.prepareForParallelTransfer(cloud, td);

                    sendParticles[td.sendToProc].append(this->remove(&p));

                    sendPatchIndices[td.sendToProc].append(td.sendToPatch);
                }
            }
            else
            {
                deleteParticle(p);
            }
        }

        // If running in serial then everything has been moved, so finish
        if (!Pstream::parRun())
        {
            break;
        }

        // Clear transfer buffers
        pBufs.clear();

        // Stream into send buffers
        forAll(sendParticles, proci)
        {
            if (sendParticles[proci].size())
            {
                UOPstream particleStream(proci, pBufs);

                particleStream
                    << sendPatchIndices[proci]
                    << sendParticles[proci];
            }
        }

        // Start sending. Sets number of bytes transferred.
        labelList receiveSizes(Pstream::nProcs());
        pBufs.finishedSends(receiveSizes);

        // Determine if any particles were transferred. If not, then finish.
        bool transferred = false;
        forAll(receiveSizes, proci)
        {
            if (receiveSizes[proci])
            {
                transferred = true;
                break;
            }
        }
        reduce(transferred, orOp<bool>());
        if (!transferred)
        {
            break;
        }

        // Retrieve from receive buffers and add into the cloud
        forAll(receiveSizes, proci)
        {
            if (receiveSizes[proci])
            {
                UIPstream particleStream(proci, pBufs);

                const labelList receivePatchIndices(particleStream);

                IDLList<ParticleType> newParticles(particleStream);

                label i = 0;

                forAllIter(typename Cloud<ParticleType>, newParticles, iter)
                {
                    const label patchi = receivePatchIndices[i ++];

                    ParticleType& p = iter();

                    td.sendToPatch = patchi;

                    p.correctAfterParallelTransfer(cloud, td);

                    addParticle(newParticles.remove(&p));
                }
            }
        }
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::topoChange(const polyTopoChangeMap& map)
{
    if (map.reverseCellMap().empty()) return;

    // Ask for the tetBasePtIs to trigger all processors to build
    // them, otherwise, if some processors have no particles then
    // there is a comms mismatch.
    pMesh_.tetBasePtIs();
    pMesh_.oldCellCentres();

    if (!globalPositionsPtr_.valid())
    {
        FatalErrorInFunction
            << "Global positions are not available. "
            << "Cloud::storeGlobalPositions has not been called."
            << exit(FatalError);
    }

    const vectorField& positions = globalPositionsPtr_();

    label particlei = 0;
    forAllIter(typename Cloud<ParticleType>, *this, iter)
    {
        const label celli = map.reverseCellMap()[iter().cell()];
        iter().map(pMesh_, positions[particlei++], celli);
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::mapMesh(const polyMeshMap& map)
{
    // Ask for the tetBasePtIs to trigger all processors to build
    // them, otherwise, if some processors have no particles then
    // there is a comms mismatch.
    pMesh_.tetBasePtIs();
    pMesh_.oldCellCentres();

    // Update cached mesh indexing
    patchNbrProc_ = patchNbrProc(pMesh_);
    patchNbrProcPatch_ = patchNbrProcPatch(pMesh_);
    patchNonConformalCyclicPatches_ = patchNonConformalCyclicPatches(pMesh_);

    if (!globalPositionsPtr_.valid())
    {
        FatalErrorInFunction
            << "Global positions are not available. "
            << "Cloud::storeGlobalPositions has not been called."
            << exit(FatalError);
    }

    const vectorField& positions = globalPositionsPtr_();

    // Loop the particles. Map those that remain on this processor, and
    // transfer others into send arrays.
    List<DynamicList<label>> sendCellIndices(Pstream::nProcs());
    List<DynamicList<point>> sendPositions(Pstream::nProcs());
    List<IDLList<ParticleType>> sendParticles(Pstream::nProcs());
    {
        label particlei = 0;
        forAllIter(typename Cloud<ParticleType>, *this, iter)
        {
            const point& pos = positions[particlei ++];

            const remote tgtProcCell =
                map.mapper().srcToTgtPoint(iter().cell(), pos);

            if (tgtProcCell == remote())
            {
                WarningInFunction
                    << "Particle at " << pos << " mapped to a location outside "
                    << "of the new mesh. This particle will be removed." << nl;
            }
            else if (tgtProcCell.proci == Pstream::myProcNo())
            {
                iter().map(pMesh_, pos, tgtProcCell.elementi);
            }
            else
            {
                sendCellIndices[tgtProcCell.proci].append(tgtProcCell.elementi);
                sendPositions[tgtProcCell.proci].append(pos);
                sendParticles[tgtProcCell.proci].append(this->remove(iter));
            }
        }
    }

    // If serial then there is nothing more to do
    if (!Pstream::parRun())
    {
        return;
    }

    // Create transfer buffers
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // Stream into send buffers
    forAll(sendParticles, proci)
    {
        if (sendParticles[proci].size())
        {
            UOPstream particleStream(proci, pBufs);

            particleStream
                << sendCellIndices[proci]
                << sendPositions[proci]
                << sendParticles[proci];
        }
    }

    // Finish sending
    labelList receiveSizes(Pstream::nProcs());
    pBufs.finishedSends(receiveSizes);

    // Retrieve from receive buffers and map into the new mesh
    forAll(sendParticles, proci)
    {
        if (receiveSizes[proci])
        {
            UIPstream particleStream(proci, pBufs);

            const labelList receiveCellIndices(particleStream);
            const List<point> receivePositions(particleStream);
            IDLList<ParticleType> receiveParticles(particleStream);

            label particlei = 0;
            forAllIter(typename Cloud<ParticleType>, receiveParticles, iter)
            {
                const label celli = receiveCellIndices[particlei];
                const vector& pos = receivePositions[particlei ++];

                iter().map(pMesh_, pos, celli);
                this->append(receiveParticles.remove(iter));
            }
        }
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::distribute(const polyDistributionMap& map)
{
    // Ask for the tetBasePtIs to trigger all processors to build
    // them, otherwise, if some processors have no particles then
    // there is a comms mismatch.
    pMesh_.tetBasePtIs();
    pMesh_.oldCellCentres();

    // Update cached mesh indexing
    patchNbrProc_ = patchNbrProc(pMesh_);
    patchNbrProcPatch_ = patchNbrProcPatch(pMesh_);
    patchNonConformalCyclicPatches_ = patchNonConformalCyclicPatches(pMesh_);

    if (!globalPositionsPtr_.valid())
    {
        FatalErrorInFunction
            << "Global positions are not available. "
            << "Cloud::storeGlobalPositions has not been called."
            << exit(FatalError);
    }

    const vectorField& positions = globalPositionsPtr_();

    // Distribute the global positions
    List<List<point>> cellParticlePositions(map.nOldCells());
    {
        labelList cellParticleis(map.nOldCells(), 0);
        forAllIter(typename Cloud<ParticleType>, *this, iter)
        {
            cellParticleis[iter().cell()] ++;
        }
        forAll(cellParticlePositions, celli)
        {
            cellParticlePositions[celli].resize(cellParticleis[celli]);
        }

        label particlei = 0;
        cellParticleis = 0;
        forAllIter(typename Cloud<ParticleType>, *this, iter)
        {
            const label celli = iter().cell();
            label& cellParticlei = cellParticleis[celli];

            cellParticlePositions[celli][cellParticlei ++] =
                positions[particlei ++];
        }
    }
    distributionMapBase::distribute
    (
        Pstream::commsTypes::nonBlocking,
        List<labelPair>(),
        pMesh_.nCells(),
        map.cellMap().subMap(),
        false,
        map.cellMap().constructMap(),
        false,
        cellParticlePositions,
        ListAppendEqOp<point>(),
        flipOp(),
        List<point>()
    );

    // Distribute the particles
    List<IDLList<ParticleType>> cellParticles(map.nOldCells());
    forAllIter(typename Cloud<ParticleType>, *this, iter)
    {
        cellParticles[iter().cell()].append(this->remove(iter));
    }
    distributionMapBase::distribute
    (
        Pstream::commsTypes::nonBlocking,
        List<labelPair>(),
        pMesh_.nCells(),
        map.cellMap().subMap(),
        false,
        map.cellMap().constructMap(),
        false,
        cellParticles,
        IDLListAppendEqOp<ParticleType>(),
        flipOp(),
        IDLList<ParticleType>()
    );

    // Locate the particles within the new mesh
    forAll(cellParticles, celli)
    {
        label cellParticlei = 0;
        forAllIter(typename IDLList<ParticleType>, cellParticles[celli], iter)
        {
            const point& pos = cellParticlePositions[celli][cellParticlei++];

            iter().map(pMesh_, pos, celli);

            this->append(cellParticles[celli].remove(iter));
        }
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::writePositions() const
{
    OFstream pObj
    (
        this->db().time().path()/this->name() + "_positions.obj"
    );

    forAllConstIter(typename Cloud<ParticleType>, *this, pIter)
    {
        const point pos = pIter().position(pMesh_);

        pObj<< "v " << pos.x() << " " << pos.y() << " " << pos.z() << nl;
    }

    pObj.flush();
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::storeGlobalPositions() const
{
    // Store the global positions for later use by mapping functions. It would
    // be preferable not to need this. If the objects passed to had a copy of
    // the old mesh then the global positions could be recovered within the
    // mapping functions, and this pre-processing would not be necessary.

    globalPositionsPtr_.reset(new vectorField(this->size()));

    vectorField& positions = globalPositionsPtr_();

    label particlei = 0;
    forAllConstIter(typename Cloud<ParticleType>, *this, iter)
    {
        positions[particlei++] = iter().position(pMesh_);
    }
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "CloudIO.C"

// ************************************************************************* //
