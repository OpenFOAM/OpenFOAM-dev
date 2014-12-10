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

#include "FaceCellWave.H"
#include "polyMesh.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "cyclicAMIPolyPatch.H"
#include "OPstream.H"
#include "IPstream.H"
#include "PstreamReduceOps.H"
#include "debug.H"
#include "typeInfo.H"
#include "SubField.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type, class TrackingData>
const Foam::scalar Foam::FaceCellWave<Type, TrackingData>::geomTol_ = 1e-6;

template<class Type, class TrackingData>
Foam::scalar Foam::FaceCellWave<Type, TrackingData>::propagationTol_ = 0.01;

template<class Type, class TrackingData>
int Foam::FaceCellWave<Type, TrackingData>::dummyTrackData_ = 12345;

namespace Foam
{
    //- Combine operator for AMIInterpolation
    template<class Type, class TrackingData>
    class combine
    {
        FaceCellWave<Type, TrackingData>& solver_;

        const cyclicAMIPolyPatch& patch_;

        public:

            combine
            (
                FaceCellWave<Type, TrackingData>& solver,
                const cyclicAMIPolyPatch& patch
            )
            :
                solver_(solver),
                patch_(patch)
            {}


            void operator()
            (
                Type& x,
                const label faceI,
                const Type& y,
                const scalar weight
            ) const
            {
                if (y.valid(solver_.data()))
                {
                    label meshFaceI = -1;
                    if (patch_.owner())
                    {
                        meshFaceI = patch_.start() + faceI;
                    }
                    else
                    {
                        meshFaceI = patch_.neighbPatch().start() + faceI;
                    }
                    x.updateFace
                    (
                        solver_.mesh(),
                        meshFaceI,
                        y,
                        solver_.propagationTol(),
                        solver_.data()
                    );
                }
            }
    };
}



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Update info for cellI, at position pt, with information from
// neighbouring face/cell.
// Updates:
//      - changedCell_, changedCells_, nChangedCells_,
//      - statistics: nEvals_, nUnvisitedCells_
template<class Type, class TrackingData>
bool Foam::FaceCellWave<Type, TrackingData>::updateCell
(
    const label cellI,
    const label neighbourFaceI,
    const Type& neighbourInfo,
    const scalar tol,
    Type& cellInfo
)
{
    nEvals_++;

    bool wasValid = cellInfo.valid(td_);

    bool propagate =
        cellInfo.updateCell
        (
            mesh_,
            cellI,
            neighbourFaceI,
            neighbourInfo,
            tol,
            td_
        );

    if (propagate)
    {
        if (!changedCell_[cellI])
        {
            changedCell_[cellI] = true;
            changedCells_[nChangedCells_++] = cellI;
        }
    }

    if (!wasValid && cellInfo.valid(td_))
    {
        --nUnvisitedCells_;
    }

    return propagate;
}


// Update info for faceI, at position pt, with information from
// neighbouring face/cell.
// Updates:
//      - changedFace_, changedFaces_, nChangedFaces_,
//      - statistics: nEvals_, nUnvisitedFaces_
template<class Type, class TrackingData>
bool Foam::FaceCellWave<Type, TrackingData>::updateFace
(
    const label faceI,
    const label neighbourCellI,
    const Type& neighbourInfo,
    const scalar tol,
    Type& faceInfo
)
{
    nEvals_++;

    bool wasValid = faceInfo.valid(td_);

    bool propagate =
        faceInfo.updateFace
        (
            mesh_,
            faceI,
            neighbourCellI,
            neighbourInfo,
            tol,
            td_
        );

    if (propagate)
    {
        if (!changedFace_[faceI])
        {
            changedFace_[faceI] = true;
            changedFaces_[nChangedFaces_++] = faceI;
        }
    }

    if (!wasValid && faceInfo.valid(td_))
    {
        --nUnvisitedFaces_;
    }

    return propagate;
}


// Update info for faceI, at position pt, with information from
// same face.
// Updates:
//      - changedFace_, changedFaces_, nChangedFaces_,
//      - statistics: nEvals_, nUnvisitedFaces_
template<class Type, class TrackingData>
bool Foam::FaceCellWave<Type, TrackingData>::updateFace
(
    const label faceI,
    const Type& neighbourInfo,
    const scalar tol,
    Type& faceInfo
)
{
    nEvals_++;

    bool wasValid = faceInfo.valid(td_);

    bool propagate =
        faceInfo.updateFace
        (
            mesh_,
            faceI,
            neighbourInfo,
            tol,
            td_
        );

    if (propagate)
    {
        if (!changedFace_[faceI])
        {
            changedFace_[faceI] = true;
            changedFaces_[nChangedFaces_++] = faceI;
        }
    }

    if (!wasValid && faceInfo.valid(td_))
    {
        --nUnvisitedFaces_;
    }

    return propagate;
}


// For debugging: check status on both sides of cyclic
template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::checkCyclic
(
    const polyPatch& patch
) const
{
    const cyclicPolyPatch& nbrPatch =
        refCast<const cyclicPolyPatch>(patch).neighbPatch();

    forAll(patch, patchFaceI)
    {
        label i1 = patch.start() + patchFaceI;
        label i2 = nbrPatch.start() + patchFaceI;

        if
        (
           !allFaceInfo_[i1].sameGeometry
            (
                mesh_,
                allFaceInfo_[i2],
                geomTol_,
                td_
            )
        )
        {
            FatalErrorIn
            (
                "FaceCellWave<Type, TrackingData>"
                "::checkCyclic(const polyPatch&)"
            )   << "problem: i:" << i1 << "  otheri:" << i2
                << "   faceInfo:" << allFaceInfo_[i1]
                << "   otherfaceInfo:" << allFaceInfo_[i2]
                << abort(FatalError);
        }

        if (changedFace_[i1] != changedFace_[i2])
        {
            FatalErrorIn
            (
                "FaceCellWave<Type, TrackingData>"
                "::checkCyclic(const polyPatch&)"
            )   << " problem: i:" << i1 << "  otheri:" << i2
                << "   faceInfo:" << allFaceInfo_[i1]
                << "   otherfaceInfo:" << allFaceInfo_[i2]
                << "   changedFace:" << changedFace_[i1]
                << "   otherchangedFace:" << changedFace_[i2]
                << abort(FatalError);
        }
    }
}


// Check if has cyclic patches
template<class Type, class TrackingData>
template<class PatchType>
bool Foam::FaceCellWave<Type, TrackingData>::hasPatch() const
{
    forAll(mesh_.boundaryMesh(), patchI)
    {
        if (isA<PatchType>(mesh_.boundaryMesh()[patchI]))
        {
            return true;
        }
    }
    return false;
}


// Copy face information into member data
template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::setFaceInfo
(
    const labelList& changedFaces,
    const List<Type>& changedFacesInfo
)
{
    forAll(changedFaces, changedFaceI)
    {
        label faceI = changedFaces[changedFaceI];

        bool wasValid = allFaceInfo_[faceI].valid(td_);

        // Copy info for faceI
        allFaceInfo_[faceI] = changedFacesInfo[changedFaceI];

        // Maintain count of unset faces
        if (!wasValid && allFaceInfo_[faceI].valid(td_))
        {
            --nUnvisitedFaces_;
        }

        // Mark faceI as changed, both on list and on face itself.

        changedFace_[faceI] = true;
        changedFaces_[nChangedFaces_++] = faceI;
    }
}


// Merge face information into member data
template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::mergeFaceInfo
(
    const polyPatch& patch,
    const label nFaces,
    const labelList& changedFaces,
    const List<Type>& changedFacesInfo
)
{
    for (label changedFaceI = 0; changedFaceI < nFaces; changedFaceI++)
    {
        const Type& neighbourWallInfo = changedFacesInfo[changedFaceI];
        label patchFaceI = changedFaces[changedFaceI];

        label meshFaceI = patch.start() + patchFaceI;

        Type& currentWallInfo = allFaceInfo_[meshFaceI];

        if (!currentWallInfo.equal(neighbourWallInfo, td_))
        {
            updateFace
            (
                meshFaceI,
                neighbourWallInfo,
                propagationTol_,
                currentWallInfo
            );
        }
    }
}


// Construct compact patchFace change arrays for a (slice of a) single patch.
// changedPatchFaces in local patch numbering.
// Return length of arrays.
template<class Type, class TrackingData>
Foam::label Foam::FaceCellWave<Type, TrackingData>::getChangedPatchFaces
(
    const polyPatch& patch,
    const label startFaceI,
    const label nFaces,
    labelList& changedPatchFaces,
    List<Type>& changedPatchFacesInfo
) const
{
    label nChangedPatchFaces = 0;

    for (label i = 0; i < nFaces; i++)
    {
        label patchFaceI = i + startFaceI;

        label meshFaceI = patch.start() + patchFaceI;

        if (changedFace_[meshFaceI])
        {
            changedPatchFaces[nChangedPatchFaces] = patchFaceI;
            changedPatchFacesInfo[nChangedPatchFaces] = allFaceInfo_[meshFaceI];
            nChangedPatchFaces++;
        }
    }
    return nChangedPatchFaces;
}


// Handle leaving domain. Implementation referred to Type
template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::leaveDomain
(
    const polyPatch& patch,
    const label nFaces,
    const labelList& faceLabels,
    List<Type>& faceInfo
) const
{
    const vectorField& fc = mesh_.faceCentres();

    for (label i = 0; i < nFaces; i++)
    {
        label patchFaceI = faceLabels[i];

        label meshFaceI = patch.start() + patchFaceI;
        faceInfo[i].leaveDomain(mesh_, patch, patchFaceI, fc[meshFaceI], td_);
    }
}


// Handle entering domain. Implementation referred to Type
template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::enterDomain
(
    const polyPatch& patch,
    const label nFaces,
    const labelList& faceLabels,
    List<Type>& faceInfo
) const
{
    const vectorField& fc = mesh_.faceCentres();

    for (label i = 0; i < nFaces; i++)
    {
        label patchFaceI = faceLabels[i];

        label meshFaceI = patch.start() + patchFaceI;
        faceInfo[i].enterDomain(mesh_, patch, patchFaceI, fc[meshFaceI], td_);
    }
}


// Transform. Implementation referred to Type
template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::transform
(
    const tensorField& rotTensor,
    const label nFaces,
    List<Type>& faceInfo
)
{
    if (rotTensor.size() == 1)
    {
        const tensor& T = rotTensor[0];

        for (label faceI = 0; faceI < nFaces; faceI++)
        {
            faceInfo[faceI].transform(mesh_, T, td_);
        }
    }
    else
    {
        for (label faceI = 0; faceI < nFaces; faceI++)
        {
            faceInfo[faceI].transform(mesh_, rotTensor[faceI], td_);
        }
    }
}


// Offset mesh face. Used for transferring from one cyclic half to the other.
template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::offset
(
    const polyPatch&,
    const label cycOffset,
    const label nFaces,
    labelList& faces
)
{
    for (label faceI = 0; faceI < nFaces; faceI++)
    {
        faces[faceI] += cycOffset;
    }
}


// Tranfer all the information to/from neighbouring processors
template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::handleProcPatches()
{
    const globalMeshData& pData = mesh_.globalData();

    // Which patches are processor patches
    const labelList& procPatches = pData.processorPatches();

    // Send all

    PstreamBuffers pBufs(Pstream::nonBlocking);

    forAll(procPatches, i)
    {
        label patchI = procPatches[i];

        const processorPolyPatch& procPatch =
            refCast<const processorPolyPatch>(mesh_.boundaryMesh()[patchI]);

        // Allocate buffers
        label nSendFaces;
        labelList sendFaces(procPatch.size());
        List<Type> sendFacesInfo(procPatch.size());

        // Determine which faces changed on current patch
        nSendFaces = getChangedPatchFaces
        (
            procPatch,
            0,
            procPatch.size(),
            sendFaces,
            sendFacesInfo
        );

        // Adapt wallInfo for leaving domain
        leaveDomain
        (
            procPatch,
            nSendFaces,
            sendFaces,
            sendFacesInfo
        );

        if (debug & 2)
        {
            Pout<< " Processor patch " << patchI << ' ' << procPatch.name()
                << " communicating with " << procPatch.neighbProcNo()
                << "  Sending:" << nSendFaces
                << endl;
        }

        UOPstream toNeighbour(procPatch.neighbProcNo(), pBufs);
        //writeFaces(nSendFaces, sendFaces, sendFacesInfo, toNeighbour);
        toNeighbour
            << SubList<label>(sendFaces, nSendFaces)
            << SubList<Type>(sendFacesInfo, nSendFaces);
    }

    pBufs.finishedSends();

    // Receive all

    forAll(procPatches, i)
    {
        label patchI = procPatches[i];

        const processorPolyPatch& procPatch =
            refCast<const processorPolyPatch>(mesh_.boundaryMesh()[patchI]);

        // Allocate buffers
        labelList receiveFaces;
        List<Type> receiveFacesInfo;

        {
            UIPstream fromNeighbour(procPatch.neighbProcNo(), pBufs);
            fromNeighbour >> receiveFaces >> receiveFacesInfo;
        }

        if (debug & 2)
        {
            Pout<< " Processor patch " << patchI << ' ' << procPatch.name()
                << " communicating with " << procPatch.neighbProcNo()
                << "  Receiving:" << receiveFaces.size()
                << endl;
        }

        // Apply transform to received data for non-parallel planes
        if (!procPatch.parallel())
        {
            transform
            (
                procPatch.forwardT(),
                receiveFaces.size(),
                receiveFacesInfo
            );
        }

        // Adapt wallInfo for entering domain
        enterDomain
        (
            procPatch,
            receiveFaces.size(),
            receiveFaces,
            receiveFacesInfo
        );

        // Merge received info
        mergeFaceInfo
        (
            procPatch,
            receiveFaces.size(),
            receiveFaces,
            receiveFacesInfo
        );
    }
}


// Transfer information across cyclic halves.
template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::handleCyclicPatches()
{
    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchI];

        if (isA<cyclicPolyPatch>(patch))
        {
            const cyclicPolyPatch& nbrPatch =
                refCast<const cyclicPolyPatch>(patch).neighbPatch();

            // Allocate buffers
            label nReceiveFaces;
            labelList receiveFaces(patch.size());
            List<Type> receiveFacesInfo(patch.size());

            // Determine which faces changed
            nReceiveFaces = getChangedPatchFaces
            (
                nbrPatch,
                0,
                nbrPatch.size(),
                receiveFaces,
                receiveFacesInfo
            );

            // Adapt wallInfo for leaving domain
            leaveDomain
            (
                nbrPatch,
                nReceiveFaces,
                receiveFaces,
                receiveFacesInfo
            );

            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patch);

            if (!cycPatch.parallel())
            {
                // received data from other half
                transform
                (
                    cycPatch.forwardT(),
                    nReceiveFaces,
                    receiveFacesInfo
                );
            }

            if (debug & 2)
            {
                Pout<< " Cyclic patch " << patchI << ' ' << cycPatch.name()
                    << "  Changed : " << nReceiveFaces
                    << endl;
            }

            // Half2: Adapt wallInfo for entering domain
            enterDomain
            (
                cycPatch,
                nReceiveFaces,
                receiveFaces,
                receiveFacesInfo
            );

            // Merge into global storage
            mergeFaceInfo
            (
                cycPatch,
                nReceiveFaces,
                receiveFaces,
                receiveFacesInfo
            );

            if (debug)
            {
                checkCyclic(cycPatch);
            }
        }
    }
}


// Transfer information across cyclic halves.
template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::handleAMICyclicPatches()
{
    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchI];

        if (isA<cyclicAMIPolyPatch>(patch))
        {
            const cyclicAMIPolyPatch& cycPatch =
                refCast<const cyclicAMIPolyPatch>(patch);

            List<Type> receiveInfo;

            {
                const cyclicAMIPolyPatch& nbrPatch =
                    refCast<const cyclicAMIPolyPatch>(patch).neighbPatch();

                // Get nbrPatch data (so not just changed faces)
                typename List<Type>::subList sendInfo
                (
                    nbrPatch.patchSlice
                    (
                        allFaceInfo_
                    )
                );

                if (!nbrPatch.parallel() || nbrPatch.separated())
                {
                    // Adapt sendInfo for leaving domain
                    const vectorField::subField fc = nbrPatch.faceCentres();
                    forAll(sendInfo, i)
                    {
                        sendInfo[i].leaveDomain(mesh_, nbrPatch, i, fc[i], td_);
                    }
                }

                // Transfer sendInfo to cycPatch
                combine<Type, TrackingData> cmb(*this, cycPatch);

                if (cycPatch.applyLowWeightCorrection())
                {
                    List<Type> defVals
                    (
                        cycPatch.patchInternalList(allCellInfo_)
                    );

                    cycPatch.interpolate(sendInfo, cmb, receiveInfo, defVals);
                }
                else
                {
                    cycPatch.interpolate(sendInfo, cmb, receiveInfo);
                }
            }

            // Apply transform to received data for non-parallel planes
            if (!cycPatch.parallel())
            {
                transform
                (
                    cycPatch.forwardT(),
                    receiveInfo.size(),
                    receiveInfo
                );
            }

            if (!cycPatch.parallel() || cycPatch.separated())
            {
                // Adapt receiveInfo for entering domain
                const vectorField::subField fc = cycPatch.faceCentres();
                forAll(receiveInfo, i)
                {
                    receiveInfo[i].enterDomain(mesh_, cycPatch, i, fc[i], td_);
                }
            }

            // Merge into global storage
            forAll(receiveInfo, i)
            {
                label meshFaceI = cycPatch.start()+i;

                Type& currentWallInfo = allFaceInfo_[meshFaceI];

                if
                (
                    receiveInfo[i].valid(td_)
                && !currentWallInfo.equal(receiveInfo[i], td_)
                )
                {
                    updateFace
                    (
                        meshFaceI,
                        receiveInfo[i],
                        propagationTol_,
                        currentWallInfo
                    );
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Set up only. Use setFaceInfo and iterate() to do actual calculation.
template<class Type, class TrackingData>
Foam::FaceCellWave<Type, TrackingData>::FaceCellWave
(
    const polyMesh& mesh,
    UList<Type>& allFaceInfo,
    UList<Type>& allCellInfo,
    TrackingData& td
)
:
    mesh_(mesh),
    allFaceInfo_(allFaceInfo),
    allCellInfo_(allCellInfo),
    td_(td),
    changedFace_(mesh_.nFaces(), false),
    changedFaces_(mesh_.nFaces()),
    nChangedFaces_(0),
    changedCell_(mesh_.nCells(), false),
    changedCells_(mesh_.nCells()),
    nChangedCells_(0),
    hasCyclicPatches_(hasPatch<cyclicPolyPatch>()),
    hasCyclicAMIPatches_
    (
        returnReduce(hasPatch<cyclicAMIPolyPatch>(), orOp<bool>())
    ),
    nEvals_(0),
    nUnvisitedCells_(mesh_.nCells()),
    nUnvisitedFaces_(mesh_.nFaces())
{
    if
    (
        allFaceInfo.size() != mesh_.nFaces()
     || allCellInfo.size() != mesh_.nCells()
    )
    {
        FatalErrorIn
        (
            "FaceCellWave<Type, TrackingData>::FaceCellWave"
            "(const polyMesh&, const labelList&, const List<Type>,"
            " UList<Type>&, UList<Type>&, const label maxIter)"
        )   << "face and cell storage not the size of mesh faces, cells:"
            << endl
            << "    allFaceInfo   :" << allFaceInfo.size() << endl
            << "    mesh_.nFaces():" << mesh_.nFaces() << endl
            << "    allCellInfo   :" << allCellInfo.size() << endl
            << "    mesh_.nCells():" << mesh_.nCells()
            << exit(FatalError);
    }
}


// Iterate, propagating changedFacesInfo across mesh, until no change (or
// maxIter reached). Initial cell values specified.
template<class Type, class TrackingData>
Foam::FaceCellWave<Type, TrackingData>::FaceCellWave
(
    const polyMesh& mesh,
    const labelList& changedFaces,
    const List<Type>& changedFacesInfo,
    UList<Type>& allFaceInfo,
    UList<Type>& allCellInfo,
    const label maxIter,
    TrackingData& td
)
:
    mesh_(mesh),
    allFaceInfo_(allFaceInfo),
    allCellInfo_(allCellInfo),
    td_(td),
    changedFace_(mesh_.nFaces(), false),
    changedFaces_(mesh_.nFaces()),
    nChangedFaces_(0),
    changedCell_(mesh_.nCells(), false),
    changedCells_(mesh_.nCells()),
    nChangedCells_(0),
    hasCyclicPatches_(hasPatch<cyclicPolyPatch>()),
    hasCyclicAMIPatches_
    (
        returnReduce(hasPatch<cyclicAMIPolyPatch>(), orOp<bool>())
    ),
    nEvals_(0),
    nUnvisitedCells_(mesh_.nCells()),
    nUnvisitedFaces_(mesh_.nFaces())
{
    if
    (
        allFaceInfo.size() != mesh_.nFaces()
     || allCellInfo.size() != mesh_.nCells()
    )
    {
        FatalErrorIn
        (
            "FaceCellWave<Type, TrackingData>::FaceCellWave"
            "(const polyMesh&, const labelList&, const List<Type>,"
            " UList<Type>&, UList<Type>&, const label maxIter)"
        )   << "face and cell storage not the size of mesh faces, cells:"
            << endl
            << "    allFaceInfo   :" << allFaceInfo.size() << endl
            << "    mesh_.nFaces():" << mesh_.nFaces() << endl
            << "    allCellInfo   :" << allCellInfo.size() << endl
            << "    mesh_.nCells():" << mesh_.nCells()
            << exit(FatalError);
    }

    // Copy initial changed faces data
    setFaceInfo(changedFaces, changedFacesInfo);

    // Iterate until nothing changes
    label iter = iterate(maxIter);

    if ((maxIter > 0) && (iter >= maxIter))
    {
        FatalErrorIn
        (
            "FaceCellWave<Type, TrackingData>::FaceCellWave"
            "(const polyMesh&, const labelList&, const List<Type>,"
            " UList<Type>&, UList<Type>&, const label maxIter)"
        )
            << "Maximum number of iterations reached. Increase maxIter." << endl
            << "    maxIter:" << maxIter << endl
            << "    nChangedCells:" << nChangedCells_ << endl
            << "    nChangedFaces:" << nChangedFaces_ << endl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type, class TrackingData>
Foam::label Foam::FaceCellWave<Type, TrackingData>::getUnsetCells() const
{
    return nUnvisitedCells_;
}


template<class Type, class TrackingData>
Foam::label Foam::FaceCellWave<Type, TrackingData>::getUnsetFaces() const
{
    return nUnvisitedFaces_;
}



// Propagate cell to face
template<class Type, class TrackingData>
Foam::label Foam::FaceCellWave<Type, TrackingData>::faceToCell()
{
    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();
    label nInternalFaces = mesh_.nInternalFaces();

    for
    (
        label changedFaceI = 0;
        changedFaceI < nChangedFaces_;
        changedFaceI++
    )
    {
        label faceI = changedFaces_[changedFaceI];
        if (!changedFace_[faceI])
        {
            FatalErrorIn("FaceCellWave<Type, TrackingData>::faceToCell()")
                << "Face " << faceI
                << " not marked as having been changed"
                << abort(FatalError);
        }


        const Type& neighbourWallInfo = allFaceInfo_[faceI];

        // Evaluate all connected cells

        // Owner
        label cellI = owner[faceI];
        Type& currentWallInfo = allCellInfo_[cellI];

        if (!currentWallInfo.equal(neighbourWallInfo, td_))
        {
            updateCell
            (
                cellI,
                faceI,
                neighbourWallInfo,
                propagationTol_,
                currentWallInfo
            );
        }

        // Neighbour.
        if (faceI < nInternalFaces)
        {
            cellI = neighbour[faceI];
            Type& currentWallInfo2 = allCellInfo_[cellI];

            if (!currentWallInfo2.equal(neighbourWallInfo, td_))
            {
                updateCell
                (
                    cellI,
                    faceI,
                    neighbourWallInfo,
                    propagationTol_,
                    currentWallInfo2
                );
            }
        }

        // Reset status of face
        changedFace_[faceI] = false;
    }

    // Handled all changed faces by now
    nChangedFaces_ = 0;

    if (debug & 2)
    {
        Pout<< " Changed cells            : " << nChangedCells_ << endl;
    }

    // Sum nChangedCells over all procs
    label totNChanged = nChangedCells_;

    reduce(totNChanged, sumOp<label>());

    return totNChanged;
}


// Propagate cell to face
template<class Type, class TrackingData>
Foam::label Foam::FaceCellWave<Type, TrackingData>::cellToFace()
{
    const cellList& cells = mesh_.cells();

    for
    (
        label changedCellI = 0;
        changedCellI < nChangedCells_;
        changedCellI++
    )
    {
        label cellI = changedCells_[changedCellI];
        if (!changedCell_[cellI])
        {
            FatalErrorIn("FaceCellWave<Type, TrackingData>::cellToFace()")
                << "Cell " << cellI << " not marked as having been changed"
                << abort(FatalError);
        }

        const Type& neighbourWallInfo = allCellInfo_[cellI];

        // Evaluate all connected faces

        const labelList& faceLabels = cells[cellI];
        forAll(faceLabels, faceLabelI)
        {
            label faceI = faceLabels[faceLabelI];
            Type& currentWallInfo = allFaceInfo_[faceI];

            if (!currentWallInfo.equal(neighbourWallInfo, td_))
            {
                updateFace
                (
                    faceI,
                    cellI,
                    neighbourWallInfo,
                    propagationTol_,
                    currentWallInfo
                );
            }
        }

        // Reset status of cell
        changedCell_[cellI] = false;
    }

    // Handled all changed cells by now
    nChangedCells_ = 0;

    if (hasCyclicPatches_)
    {
        // Transfer changed faces across cyclic halves
        handleCyclicPatches();
    }

    if (hasCyclicAMIPatches_)
    {
        handleAMICyclicPatches();
    }

    if (Pstream::parRun())
    {
        // Transfer changed faces from neighbouring processors.
        handleProcPatches();
    }

    if (debug & 2)
    {
        Pout<< " Changed faces            : " << nChangedFaces_ << endl;
    }

    // Sum nChangedFaces over all procs
    label totNChanged = nChangedFaces_;

    reduce(totNChanged, sumOp<label>());

    return totNChanged;
}


// Iterate
template<class Type, class TrackingData>
Foam::label Foam::FaceCellWave<Type, TrackingData>::iterate(const label maxIter)
{
    if (hasCyclicPatches_)
    {
        // Transfer changed faces across cyclic halves
        handleCyclicPatches();
    }

    if (hasCyclicAMIPatches_)
    {
        handleAMICyclicPatches();
    }

    if (Pstream::parRun())
    {
        // Transfer changed faces from neighbouring processors.
        handleProcPatches();
    }

    label iter = 0;

    while (iter < maxIter)
    {
        if (debug)
        {
            Info<< " Iteration " << iter << endl;
        }

        nEvals_ = 0;

        label nCells = faceToCell();

        if (debug)
        {
            Info<< " Total changed cells      : " << nCells << endl;
        }

        if (nCells == 0)
        {
            break;
        }

        label nFaces = cellToFace();

        if (debug)
        {
            Info<< " Total changed faces      : " << nFaces << nl
                << " Total evaluations        : " << nEvals_ << nl
                << " Remaining unvisited cells: " << nUnvisitedCells_ << nl
                << " Remaining unvisited faces: " << nUnvisitedFaces_ << endl;
        }

        if (nFaces == 0)
        {
            break;
        }

        ++iter;
    }

    return iter;
}


// ************************************************************************* //
