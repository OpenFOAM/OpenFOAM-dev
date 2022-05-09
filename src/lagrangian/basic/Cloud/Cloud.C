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

#include "Cloud.H"
#include "processorPolyPatch.H"
#include "globalMeshData.H"
#include "PstreamCombineReduceOps.H"
#include "polyTopoChangeMap.H"
#include "Time.H"
#include "OFstream.H"
#include "wallPolyPatch.H"
#include "cyclicAMIPolyPatch.H"
#include "nonConformalCyclicPolyPatch.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class ParticleType>
void Foam::Cloud<ParticleType>::checkPatches() const
{
    const polyBoundaryMesh& pbm = polyMesh_.boundaryMesh();
    bool ok = true;
    forAll(pbm, patchi)
    {
        if (isA<cyclicAMIPolyPatch>(pbm[patchi]))
        {
            const cyclicAMIPolyPatch& cami =
                refCast<const cyclicAMIPolyPatch>(pbm[patchi]);

            ok = ok && cami.singlePatchProc() != -1;
        }
    }

    if (!ok)
    {
        FatalErrorInFunction
            << "Particle tracking across AMI patches is only currently "
            << "supported for cases where the AMI patches reside on a "
            << "single processor" << abort(FatalError);
    }
}


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
    const polyBoundaryMesh& pbm = polyMesh_.boundaryMesh();

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
    polyMesh_(pMesh),
    patchNbrProc_(patchNbrProc(pMesh)),
    patchNbrProcPatch_(patchNbrProcPatch(pMesh)),
    patchNonConformalCyclicPatches_(patchNonConformalCyclicPatches(pMesh)),
    globalPositionsPtr_()
{
    checkPatches();

    // Ask for the tetBasePtIs and oldCellCentres to trigger all processors to
    // build them, otherwise, if some processors have no particles then there
    // is a comms mismatch.
    polyMesh_.tetBasePtIs();
    polyMesh_.oldCellCentres();

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
                << "deleting lost particle at position " << p.position()
                << endl;

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
template<class TrackCloudType>
void Foam::Cloud<ParticleType>::move
(
    TrackCloudType& cloud,
    typename ParticleType::trackingData& td,
    const scalar trackTime
)
{
    // Clear the global positions as these are about to change
    globalPositionsPtr_.clear();

    // Ensure rays are available for non conformal transfers
    storeRays();

    // Initialise the stepFraction moved for the particles
    forAllIter(typename Cloud<ParticleType>, *this, pIter)
    {
        pIter().reset(0);
    }

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
            bool keepParticle = p.move(cloud, td, trackTime);

            // If the particle is to be kept
            if (keepParticle)
            {
                if (td.sendToProc != -1)
                {
                    #ifdef FULLDEBUG
                    if (!Pstream::parRun() || !p.onBoundaryFace())
                    {
                        FatalErrorInFunction
                            << "Switch processor flag is true when no parallel "
                            << "transfer is possible. This is a bug."
                            << exit(FatalError);
                    }
                    #endif

                    p.prepareForParallelTransfer(td);

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

                IDLList<ParticleType> newParticles
                (
                    particleStream,
                    typename ParticleType::iNew(polyMesh_)
                );

                label i = 0;

                forAllIter(typename Cloud<ParticleType>, newParticles, iter)
                {
                    const label patchi = receivePatchIndices[i ++];

                    ParticleType& p = iter();

                    td.sendToPatch = patchi;

                    p.correctAfterParallelTransfer(td);

                    addParticle(newParticles.remove(&p));
                }
            }
        }
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::autoMap(const polyTopoChangeMap& mapper)
{
    if (!globalPositionsPtr_.valid())
    {
        FatalErrorInFunction
            << "Global positions are not available. "
            << "Cloud::storeGlobalPositions has not been called."
            << exit(FatalError);
    }

    // Ask for the tetBasePtIs to trigger all processors to build
    // them, otherwise, if some processors have no particles then
    // there is a comms mismatch.
    polyMesh_.tetBasePtIs();
    polyMesh_.oldCellCentres();

    const vectorField& positions = globalPositionsPtr_();

    label i = 0;
    forAllIter(typename Cloud<ParticleType>, *this, iter)
    {
        iter().autoMap(positions[i], mapper);
        ++ i;
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
        const ParticleType& p = pIter();
        pObj<< "v " << p.position().x() << " " << p.position().y() << " "
            << p.position().z() << nl;
    }

    pObj.flush();
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::storeGlobalPositions() const
{
    // Store the global positions for later use by autoMap. It would be
    // preferable not to need this. If the polyTopoChangeMap object passed to
    // autoMap had a copy of the old mesh then the global positions could be
    // recovered within autoMap, and this pre-processing would not be necessary.

    globalPositionsPtr_.reset(new vectorField(this->size()));

    vectorField& positions = globalPositionsPtr_();

    label i = 0;
    forAllConstIter(typename Cloud<ParticleType>, *this, iter)
    {
        positions[i] = iter().position();
        ++ i;
    }
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "CloudIO.C"

// ************************************************************************* //
