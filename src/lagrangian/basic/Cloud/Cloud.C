/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
#include "mapPolyMesh.H"
#include "Time.H"
#include "OFstream.H"
#include "wallPolyPatch.H"
#include "cyclicAMIPolyPatch.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class ParticleType>
void Foam::Cloud<ParticleType>::checkPatches() const
{
    const polyBoundaryMesh& pbm = polyMesh_.boundaryMesh();
    bool ok = true;
    forAll(pbm, patchI)
    {
        if (isA<cyclicAMIPolyPatch>(pbm[patchI]))
        {
            const cyclicAMIPolyPatch& cami =
                refCast<const cyclicAMIPolyPatch>(pbm[patchI]);

            if (cami.owner())
            {
                ok = ok && (cami.AMI().singlePatchProc() != -1);
            }
        }
    }

    if (!ok)
    {
        FatalErrorIn("void Foam::Cloud<ParticleType>::initCloud(const bool)")
            << "Particle tracking across AMI patches is only currently "
            << "supported for cases where the AMI patches reside on a "
            << "single processor" << abort(FatalError);
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::calcCellWallFaces() const
{
    cellWallFacesPtr_.reset(new PackedBoolList(pMesh().nCells(), false));

    PackedBoolList& cellWallFaces = cellWallFacesPtr_();

    const polyBoundaryMesh& patches = polyMesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        if (isA<wallPolyPatch>(patches[patchI]))
        {
            const polyPatch& patch = patches[patchI];

            const labelList& pFaceCells = patch.faceCells();

            forAll(pFaceCells, pFCI)
            {
                cellWallFaces[pFaceCells[pFCI]] = true;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::Cloud<ParticleType>::Cloud
(
    const polyMesh& pMesh,
    const IDLList<ParticleType>& particles
)
:
    cloud(pMesh),
    IDLList<ParticleType>(),
    polyMesh_(pMesh),
    labels_(),
    nTrackingRescues_(),
    cellWallFacesPtr_()
{
    checkPatches();

    // Ask for the tetBasePtIs to trigger all processors to build
    // them, otherwise, if some processors have no particles then
    // there is a comms mismatch.
    polyMesh_.tetBasePtIs();

    IDLList<ParticleType>::operator=(particles);
}


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
    labels_(),
    nTrackingRescues_(),
    cellWallFacesPtr_()
{
    checkPatches();

    // Ask for the tetBasePtIs to trigger all processors to build
    // them, otherwise, if some processors have no particles then
    // there is a comms mismatch.
    polyMesh_.tetBasePtIs();

    IDLList<ParticleType>::operator=(particles);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
const Foam::PackedBoolList& Foam::Cloud<ParticleType>::cellHasWallFaces()
const
{
    if (!cellWallFacesPtr_.valid())
    {
        calcCellWallFaces();
    }

    return cellWallFacesPtr_();
}


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
void Foam::Cloud<ParticleType>::cloudReset(const Cloud<ParticleType>& c)
{
    // Reset particle cound and particles only
    // - not changing the cloud object registry or reference to the polyMesh
    ParticleType::particleCount_ = 0;
    IDLList<ParticleType>::operator=(c);
}


template<class ParticleType>
template<class TrackData>
void Foam::Cloud<ParticleType>::move(TrackData& td, const scalar trackTime)
{
    const polyBoundaryMesh& pbm = pMesh().boundaryMesh();
    const globalMeshData& pData = polyMesh_.globalData();

    // Which patches are processor patches
    const labelList& procPatches = pData.processorPatches();

    // Indexing of patches into the procPatches list
    const labelList& procPatchIndices = pData.processorPatchIndices();

    // Indexing of equivalent patch on neighbour processor into the
    // procPatches list on the neighbour
    const labelList& procPatchNeighbours = pData.processorPatchNeighbours();

    // Which processors this processor is connected to
    const labelList& neighbourProcs = pData[Pstream::myProcNo()];

    // Indexing from the processor number into the neighbourProcs list
    labelList neighbourProcIndices(Pstream::nProcs(), -1);

    forAll(neighbourProcs, i)
    {
        neighbourProcIndices[neighbourProcs[i]] = i;
    }

    // Initialise the stepFraction moved for the particles
    forAllIter(typename Cloud<ParticleType>, *this, pIter)
    {
        pIter().stepFraction() = 0;
    }

    // Reset nTrackingRescues
    nTrackingRescues_ = 0;


    // List of lists of particles to be transfered for all of the
    // neighbour processors
    List<IDLList<ParticleType> > particleTransferLists
    (
        neighbourProcs.size()
    );

    // List of destination processorPatches indices for all of the
    // neighbour processors
    List<DynamicList<label> > patchIndexTransferLists
    (
        neighbourProcs.size()
    );

    // Allocate transfer buffers
    PstreamBuffers pBufs(Pstream::nonBlocking);


    // While there are particles to transfer
    while (true)
    {
        particleTransferLists = IDLList<ParticleType>();
        forAll(patchIndexTransferLists, i)
        {
            patchIndexTransferLists[i].clear();
        }

        // Loop over all particles
        forAllIter(typename Cloud<ParticleType>, *this, pIter)
        {
            ParticleType& p = pIter();

            // Move the particle
            bool keepParticle = p.move(td, trackTime);

            // If the particle is to be kept
            // (i.e. it hasn't passed through an inlet or outlet)
            if (keepParticle)
            {
                // If we are running in parallel and the particle is on a
                // boundary face
                if (Pstream::parRun() && p.face() >= pMesh().nInternalFaces())
                {
                    label patchI = pbm.whichPatch(p.face());

                    // ... and the face is on a processor patch
                    // prepare it for transfer
                    if (procPatchIndices[patchI] != -1)
                    {
                        label n = neighbourProcIndices
                        [
                            refCast<const processorPolyPatch>
                            (
                                pbm[patchI]
                            ).neighbProcNo()
                        ];

                        p.prepareForParallelTransfer(patchI, td);

                        particleTransferLists[n].append(this->remove(&p));

                        patchIndexTransferLists[n].append
                        (
                            procPatchNeighbours[patchI]
                        );
                    }
                }
            }
            else
            {
                deleteParticle(p);
            }
        }

        if (!Pstream::parRun())
        {
            break;
        }


        // Clear transfer buffers
        pBufs.clear();

        // Stream into send buffers
        forAll(particleTransferLists, i)
        {
            if (particleTransferLists[i].size())
            {
                UOPstream particleStream
                (
                    neighbourProcs[i],
                    pBufs
                );

                particleStream
                    << patchIndexTransferLists[i]
                    << particleTransferLists[i];
            }
        }


        // Start sending. Sets number of bytes transferred
        labelListList allNTrans(Pstream::nProcs());
        pBufs.finishedSends(allNTrans);


        bool transfered = false;

        forAll(allNTrans, i)
        {
            forAll(allNTrans[i], j)
            {
                if (allNTrans[i][j])
                {
                    transfered = true;
                    break;
                }
            }
        }

        if (!transfered)
        {
            break;
        }

        // Retrieve from receive buffers
        forAll(neighbourProcs, i)
        {
            label neighbProci = neighbourProcs[i];

            label nRec = allNTrans[neighbProci][Pstream::myProcNo()];

            if (nRec)
            {
                UIPstream particleStream(neighbProci, pBufs);

                labelList receivePatchIndex(particleStream);

                IDLList<ParticleType> newParticles
                (
                    particleStream,
                    typename ParticleType::iNew(polyMesh_)
                );

                label pI = 0;

                forAllIter(typename Cloud<ParticleType>, newParticles, newpIter)
                {
                    ParticleType& newp = newpIter();

                    label patchI = procPatches[receivePatchIndex[pI++]];

                    newp.correctAfterParallelTransfer(patchI, td);

                    addParticle(newParticles.remove(&newp));
                }
            }
        }
    }

    if (cloud::debug)
    {
        reduce(nTrackingRescues_, sumOp<label>());

        if (nTrackingRescues_ > 0)
        {
            Info<< nTrackingRescues_ << " tracking rescue corrections" << endl;
        }
    }
}


template<class ParticleType>
template<class TrackData>
void Foam::Cloud<ParticleType>::autoMap
(
    TrackData& td,
    const mapPolyMesh& mapper
)
{
    if (cloud::debug)
    {
        Info<< "Cloud<ParticleType>::autoMap(TrackData&, const mapPolyMesh&) "
            << "for lagrangian cloud " << cloud::name() << endl;
    }

    const labelList& reverseCellMap = mapper.reverseCellMap();
    const labelList& reverseFaceMap = mapper.reverseFaceMap();

    // Reset stored data that relies on the mesh
//    polyMesh_.clearCellTree();
    cellWallFacesPtr_.clear();

    // Ask for the tetBasePtIs to trigger all processors to build
    // them, otherwise, if some processors have no particles then
    // there is a comms mismatch.
    polyMesh_.tetBasePtIs();


    forAllIter(typename Cloud<ParticleType>, *this, pIter)
    {
        ParticleType& p = pIter();

        if (reverseCellMap[p.cell()] >= 0)
        {
            p.cell() = reverseCellMap[p.cell()];

            if (p.face() >= 0 && reverseFaceMap[p.face()] >= 0)
            {
                p.face() = reverseFaceMap[p.face()];
            }
            else
            {
                p.face() = -1;
            }

            p.initCellFacePt();
        }
        else
        {
            label trackStartCell = mapper.mergedCell(p.cell());

            if (trackStartCell < 0)
            {
                trackStartCell = 0;
                p.cell() = 0;
            }
            else
            {
                p.cell() = trackStartCell;
            }

            vector pos = p.position();

            const_cast<vector&>(p.position()) =
                polyMesh_.cellCentres()[trackStartCell];

            p.stepFraction() = 0;

            p.initCellFacePt();

            p.track(pos, td);
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
        const ParticleType& p = pIter();
        pObj<< "v " << p.position().x() << " " << p.position().y() << " "
            << p.position().z() << nl;
    }

    pObj.flush();
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "CloudIO.C"

// ************************************************************************* //
