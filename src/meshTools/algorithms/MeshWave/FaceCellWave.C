/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    template<class Type, class TrackingData>
    class combine
    {
        //- Combine operator for AMIInterpolation

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
                const label facei,
                const Type& y,
                const scalar weight
            ) const
            {
                if (y.valid(solver_.data()))
                {
                    label meshFacei = -1;
                    if (patch_.owner())
                    {
                        meshFacei = patch_.start() + facei;
                    }
                    else
                    {
                        meshFacei = patch_.neighbPatch().start() + facei;
                    }
                    x.updateFace
                    (
                        solver_.mesh(),
                        meshFacei,
                        y,
                        solver_.propagationTol(),
                        solver_.data()
                    );
                }
            }
    };
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type, class TrackingData>
bool Foam::FaceCellWave<Type, TrackingData>::updateCell
(
    const label celli,
    const label neighbourFacei,
    const Type& neighbourInfo,
    const scalar tol,
    Type& cellInfo
)
{
    // Update info for celli, at position pt, with information from
    // neighbouring face/cell.
    // Updates:
    //      - changedCell_, changedCells_
    //      - statistics: nEvals_, nUnvisitedCells_

    nEvals_++;

    bool wasValid = cellInfo.valid(td_);

    bool propagate =
        cellInfo.updateCell
        (
            mesh_,
            celli,
            neighbourFacei,
            neighbourInfo,
            tol,
            td_
        );

    if (propagate)
    {
        if (!changedCell_[celli])
        {
            changedCell_[celli] = true;
            changedCells_.append(celli);
        }
    }

    if (!wasValid && cellInfo.valid(td_))
    {
        --nUnvisitedCells_;
    }

    return propagate;
}


template<class Type, class TrackingData>
bool Foam::FaceCellWave<Type, TrackingData>::updateFace
(
    const label facei,
    const label neighbourCelli,
    const Type& neighbourInfo,
    const scalar tol,
    Type& faceInfo
)
{
    // Update info for facei, at position pt, with information from
    // neighbouring face/cell.
    // Updates:
    //      - changedFace_, changedFaces_,
    //      - statistics: nEvals_, nUnvisitedFaces_

    nEvals_++;

    bool wasValid = faceInfo.valid(td_);

    bool propagate =
        faceInfo.updateFace
        (
            mesh_,
            facei,
            neighbourCelli,
            neighbourInfo,
            tol,
            td_
        );

    if (propagate)
    {
        if (!changedFace_[facei])
        {
            changedFace_[facei] = true;
            changedFaces_.append(facei);
        }
    }

    if (!wasValid && faceInfo.valid(td_))
    {
        --nUnvisitedFaces_;
    }

    return propagate;
}


template<class Type, class TrackingData>
bool Foam::FaceCellWave<Type, TrackingData>::updateFace
(
    const label facei,
    const Type& neighbourInfo,
    const scalar tol,
    Type& faceInfo
)
{
    // Update info for facei, at position pt, with information from
    // same face.
    // Updates:
    //      - changedFace_, changedFaces_,
    //      - statistics: nEvals_, nUnvisitedFaces_

    nEvals_++;

    bool wasValid = faceInfo.valid(td_);

    bool propagate =
        faceInfo.updateFace
        (
            mesh_,
            facei,
            neighbourInfo,
            tol,
            td_
        );

    if (propagate)
    {
        if (!changedFace_[facei])
        {
            changedFace_[facei] = true;
            changedFaces_.append(facei);
        }
    }

    if (!wasValid && faceInfo.valid(td_))
    {
        --nUnvisitedFaces_;
    }

    return propagate;
}


template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::checkCyclic
(
    const polyPatch& patch
) const
{
    // For debugging: check status on both sides of cyclic

    const cyclicPolyPatch& nbrPatch =
        refCast<const cyclicPolyPatch>(patch).neighbPatch();

    forAll(patch, patchFacei)
    {
        label i1 = patch.start() + patchFacei;
        label i2 = nbrPatch.start() + patchFacei;

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
            FatalErrorInFunction
                << "   faceInfo:" << allFaceInfo_[i1]
                << "   otherfaceInfo:" << allFaceInfo_[i2]
                << abort(FatalError);
        }

        if (changedFace_[i1] != changedFace_[i2])
        {
            FatalErrorInFunction
                << "   faceInfo:" << allFaceInfo_[i1]
                << "   otherfaceInfo:" << allFaceInfo_[i2]
                << "   changedFace:" << changedFace_[i1]
                << "   otherchangedFace:" << changedFace_[i2]
                << abort(FatalError);
        }
    }
}


template<class Type, class TrackingData>
template<class PatchType>
bool Foam::FaceCellWave<Type, TrackingData>::hasPatch() const
{
    forAll(mesh_.boundaryMesh(), patchi)
    {
        if (isA<PatchType>(mesh_.boundaryMesh()[patchi]))
        {
            return true;
        }
    }
    return false;
}


template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::setFaceInfo
(
    const labelList& changedFaces,
    const List<Type>& changedFacesInfo
)
{
    forAll(changedFaces, changedFacei)
    {
        label facei = changedFaces[changedFacei];

        bool wasValid = allFaceInfo_[facei].valid(td_);

        // Copy info for facei
        allFaceInfo_[facei] = changedFacesInfo[changedFacei];

        // Maintain count of unset faces
        if (!wasValid && allFaceInfo_[facei].valid(td_))
        {
            --nUnvisitedFaces_;
        }

        // Mark facei as changed, both on list and on face itself.

        changedFace_[facei] = true;
        changedFaces_.append(facei);
    }
}


template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::mergeFaceInfo
(
    const polyPatch& patch,
    const label nFaces,
    const labelList& changedFaces,
    const List<Type>& changedFacesInfo
)
{
    // Merge face information into member data

    for (label changedFacei = 0; changedFacei < nFaces; changedFacei++)
    {
        const Type& neighbourWallInfo = changedFacesInfo[changedFacei];
        label patchFacei = changedFaces[changedFacei];

        label meshFacei = patch.start() + patchFacei;

        Type& currentWallInfo = allFaceInfo_[meshFacei];

        if (!currentWallInfo.equal(neighbourWallInfo, td_))
        {
            updateFace
            (
                meshFacei,
                neighbourWallInfo,
                propagationTol_,
                currentWallInfo
            );
        }
    }
}


template<class Type, class TrackingData>
Foam::label Foam::FaceCellWave<Type, TrackingData>::getChangedPatchFaces
(
    const polyPatch& patch,
    const label startFacei,
    const label nFaces,
    labelList& changedPatchFaces,
    List<Type>& changedPatchFacesInfo
) const
{
    // Construct compact patchFace change arrays for a (slice of a) single
    // patch. changedPatchFaces in local patch numbering.
    // Return length of arrays.
    label nChangedPatchFaces = 0;

    for (label i = 0; i < nFaces; i++)
    {
        label patchFacei = i + startFacei;

        label meshFacei = patch.start() + patchFacei;

        if (changedFace_[meshFacei])
        {
            changedPatchFaces[nChangedPatchFaces] = patchFacei;
            changedPatchFacesInfo[nChangedPatchFaces] = allFaceInfo_[meshFacei];
            nChangedPatchFaces++;
        }
    }
    return nChangedPatchFaces;
}


template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::leaveDomain
(
    const polyPatch& patch,
    const label nFaces,
    const labelList& faceLabels,
    List<Type>& faceInfo
) const
{
    // Handle leaving domain. Implementation referred to Type

    const vectorField& fc = mesh_.faceCentres();

    for (label i = 0; i < nFaces; i++)
    {
        label patchFacei = faceLabels[i];

        label meshFacei = patch.start() + patchFacei;
        faceInfo[i].leaveDomain(mesh_, patch, patchFacei, fc[meshFacei], td_);
    }
}


template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::enterDomain
(
    const polyPatch& patch,
    const label nFaces,
    const labelList& faceLabels,
    List<Type>& faceInfo
) const
{
    // Handle entering domain. Implementation referred to Type

    const vectorField& fc = mesh_.faceCentres();

    for (label i = 0; i < nFaces; i++)
    {
        label patchFacei = faceLabels[i];

        label meshFacei = patch.start() + patchFacei;
        faceInfo[i].enterDomain(mesh_, patch, patchFacei, fc[meshFacei], td_);
    }
}


template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::transform
(
    const tensorField& rotTensor,
    const label nFaces,
    List<Type>& faceInfo
)
{
    // Transform. Implementation referred to Type

    if (rotTensor.size() == 1)
    {
        const tensor& T = rotTensor[0];

        for (label facei = 0; facei < nFaces; facei++)
        {
            faceInfo[facei].transform(mesh_, T, td_);
        }
    }
    else
    {
        for (label facei = 0; facei < nFaces; facei++)
        {
            faceInfo[facei].transform(mesh_, rotTensor[facei], td_);
        }
    }
}


template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::transform
(
    const vectorTensorTransform& trans,
    const label nFaces,
    List<Type>& faceInfo
)
{
    // Transform. Implementation referred to Type

    if (trans.hasR())
    {
        for (label facei = 0; facei < nFaces; facei++)
        {
            faceInfo[facei].transform(mesh_, trans.R(), td_);
        }
    }
}


template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::offset
(
    const polyPatch&,
    const label cycOffset,
    const label nFaces,
    labelList& faces
)
{
    // Offset mesh face. Used for transferring from one cyclic half to the
    // other.

    for (label facei = 0; facei < nFaces; facei++)
    {
        faces[facei] += cycOffset;
    }
}


template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::handleProcPatches()
{
    // Transfer all the information to/from neighbouring processors

    const globalMeshData& pData = mesh_.globalData();

    // Which patches are processor patches
    const labelList& procPatches = pData.processorPatches();

    // Send all

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    forAll(procPatches, i)
    {
        label patchi = procPatches[i];

        const processorPolyPatch& procPatch =
            refCast<const processorPolyPatch>(mesh_.boundaryMesh()[patchi]);

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
            Pout<< " Processor patch " << patchi << ' ' << procPatch.name()
                << " communicating with " << procPatch.neighbProcNo()
                << "  Sending:" << nSendFaces
                << endl;
        }

        UOPstream toNeighbour(procPatch.neighbProcNo(), pBufs);
        // writeFaces(nSendFaces, sendFaces, sendFacesInfo, toNeighbour);
        toNeighbour
            << SubList<label>(sendFaces, nSendFaces)
            << SubList<Type>(sendFacesInfo, nSendFaces);
    }

    pBufs.finishedSends();

    // Receive all

    forAll(procPatches, i)
    {
        label patchi = procPatches[i];

        const processorPolyPatch& procPatch =
            refCast<const processorPolyPatch>(mesh_.boundaryMesh()[patchi]);

        // Allocate buffers
        labelList receiveFaces;
        List<Type> receiveFacesInfo;

        {
            UIPstream fromNeighbour(procPatch.neighbProcNo(), pBufs);
            fromNeighbour >> receiveFaces >> receiveFacesInfo;
        }

        if (debug & 2)
        {
            Pout<< " Processor patch " << patchi << ' ' << procPatch.name()
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


template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::handleCyclicPatches()
{
    // Transfer information across cyclic halves.

    forAll(mesh_.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchi];

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
                Pout<< " Cyclic patch " << patchi << ' ' << cycPatch.name()
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


template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::handleAMICyclicPatches()
{
    // Transfer information across cyclicAMI halves.

    forAll(mesh_.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchi];

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

                List<Type> defVals;
                if (cycPatch.applyLowWeightCorrection())
                {
                    defVals = cycPatch.patchInternalList(allCellInfo_);
                }

                if (cycPatch.owner())
                {
                    forAll(cycPatch.AMIs(), i)
                    {
                        List<Type> sendInfoT(sendInfo);
                        transform
                        (
                            cycPatch.AMITransforms()[i],
                            sendInfoT.size(),
                            sendInfoT
                        );
                        cycPatch.AMIs()[i].interpolateToSource
                        (
                            sendInfoT,
                            cmb,
                            receiveInfo,
                            defVals
                        );
                    }
                }
                else
                {
                    forAll(cycPatch.neighbPatch().AMIs(), i)
                    {
                        List<Type> sendInfoT(sendInfo);
                        transform
                        (
                            cycPatch.neighbPatch().AMITransforms()[i],
                            sendInfoT.size(),
                            sendInfoT
                        );
                        cycPatch.neighbPatch().AMIs()[i].interpolateToTarget
                        (
                            sendInfoT,
                            cmb,
                            receiveInfo,
                            defVals
                        );
                    }
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
                label meshFacei = cycPatch.start()+i;

                Type& currentWallInfo = allFaceInfo_[meshFacei];

                if
                (
                    receiveInfo[i].valid(td_)
                && !currentWallInfo.equal(receiveInfo[i], td_)
                )
                {
                    updateFace
                    (
                        meshFacei,
                        receiveInfo[i],
                        propagationTol_,
                        currentWallInfo
                    );
                }
            }
        }
    }
}


template<class Type, class TrackingData>
void Foam::FaceCellWave<Type, TrackingData>::handleExplicitConnections()
{
    // Collect changed information

    DynamicList<label> f0Baffle(explicitConnections_.size());
    DynamicList<Type> f0Info(explicitConnections_.size());

    DynamicList<label> f1Baffle(explicitConnections_.size());
    DynamicList<Type> f1Info(explicitConnections_.size());

    forAll(explicitConnections_, connI)
    {
        const labelPair& baffle = explicitConnections_[connI];

        label f0 = baffle[0];
        if (changedFace_[f0])
        {
            f0Baffle.append(connI);
            f0Info.append(allFaceInfo_[f0]);
        }

        label f1 = baffle[1];
        if (changedFace_[f1])
        {
            f1Baffle.append(connI);
            f1Info.append(allFaceInfo_[f1]);
        }
    }


    // Update other side with changed information

    forAll(f1Info, index)
    {
        const labelPair& baffle = explicitConnections_[f1Baffle[index]];

        label f0 = baffle[0];
        Type& currentWallInfo = allFaceInfo_[f0];

        if (!currentWallInfo.equal(f1Info[index], td_))
        {
            updateFace
            (
                f0,
                f1Info[index],
                propagationTol_,
                currentWallInfo
            );
        }
    }

    forAll(f0Info, index)
    {
        const labelPair& baffle = explicitConnections_[f0Baffle[index]];

        label f1 = baffle[1];
        Type& currentWallInfo = allFaceInfo_[f1];

        if (!currentWallInfo.equal(f0Info[index], td_))
        {
            updateFace
            (
                f1,
                f0Info[index],
                propagationTol_,
                currentWallInfo
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

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
    explicitConnections_(0),
    allFaceInfo_(allFaceInfo),
    allCellInfo_(allCellInfo),
    td_(td),
    changedFace_(mesh_.nFaces(), false),
    changedFaces_(mesh_.nFaces()),
    changedCell_(mesh_.nCells(), false),
    changedCells_(mesh_.nCells()),
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
        FatalErrorInFunction
            << "face and cell storage not the size of mesh faces, cells:"
            << endl
            << "    allFaceInfo   :" << allFaceInfo.size() << endl
            << "    mesh_.nFaces():" << mesh_.nFaces() << endl
            << "    allCellInfo   :" << allCellInfo.size() << endl
            << "    mesh_.nCells():" << mesh_.nCells()
            << exit(FatalError);
    }
}


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
    explicitConnections_(0),
    allFaceInfo_(allFaceInfo),
    allCellInfo_(allCellInfo),
    td_(td),
    changedFace_(mesh_.nFaces(), false),
    changedFaces_(mesh_.nFaces()),
    changedCell_(mesh_.nCells(), false),
    changedCells_(mesh_.nCells()),
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
        FatalErrorInFunction
            << "face and cell storage not the size of mesh faces, cells:"
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
        FatalErrorInFunction
            << "Maximum number of iterations reached. Increase maxIter." << endl
            << "    maxIter:" << maxIter << endl
            << "    nChangedCells:" << changedCells_.size() << endl
            << "    nChangedFaces:" << changedFaces_.size() << endl
            << exit(FatalError);
    }
}


template<class Type, class TrackingData>
Foam::FaceCellWave<Type, TrackingData>::FaceCellWave
(
    const polyMesh& mesh,
    const labelPairList& explicitConnections,
    const bool handleCyclicAMI,
    const labelList& changedFaces,
    const List<Type>& changedFacesInfo,
    UList<Type>& allFaceInfo,
    UList<Type>& allCellInfo,
    const label maxIter,
    TrackingData& td
)
:
    mesh_(mesh),
    explicitConnections_(explicitConnections),
    allFaceInfo_(allFaceInfo),
    allCellInfo_(allCellInfo),
    td_(td),
    changedFace_(mesh_.nFaces(), false),
    changedFaces_(mesh_.nFaces()),
    changedCell_(mesh_.nCells(), false),
    changedCells_(mesh_.nCells()),
    hasCyclicPatches_(hasPatch<cyclicPolyPatch>()),
    hasCyclicAMIPatches_
    (
        handleCyclicAMI
     && returnReduce(hasPatch<cyclicAMIPolyPatch>(), orOp<bool>())
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
        FatalErrorInFunction
            << "face and cell storage not the size of mesh faces, cells:"
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
        FatalErrorInFunction
            << "Maximum number of iterations reached. Increase maxIter." << endl
            << "    maxIter:" << maxIter << endl
            << "    nChangedCells:" << changedCells_.size() << endl
            << "    nChangedFaces:" << changedFaces_.size() << endl
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


template<class Type, class TrackingData>
Foam::label Foam::FaceCellWave<Type, TrackingData>::faceToCell()
{
    // Propagate face to cell

    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();
    label nInternalFaces = mesh_.nInternalFaces();

    forAll(changedFaces_, changedFacei)
    {
        label facei = changedFaces_[changedFacei];
        if (!changedFace_[facei])
        {
            FatalErrorInFunction
                << "Face " << facei
                << " not marked as having been changed"
                << abort(FatalError);
        }


        const Type& neighbourWallInfo = allFaceInfo_[facei];

        // Evaluate all connected cells

        // Owner
        label celli = owner[facei];
        Type& currentWallInfo = allCellInfo_[celli];

        if (!currentWallInfo.equal(neighbourWallInfo, td_))
        {
            updateCell
            (
                celli,
                facei,
                neighbourWallInfo,
                propagationTol_,
                currentWallInfo
            );
        }

        // Neighbour.
        if (facei < nInternalFaces)
        {
            celli = neighbour[facei];
            Type& currentWallInfo2 = allCellInfo_[celli];

            if (!currentWallInfo2.equal(neighbourWallInfo, td_))
            {
                updateCell
                (
                    celli,
                    facei,
                    neighbourWallInfo,
                    propagationTol_,
                    currentWallInfo2
                );
            }
        }

        // Reset status of face
        changedFace_[facei] = false;
    }

    // Handled all changed faces by now
    changedFaces_.clear();

    if (debug & 2)
    {
        Pout<< " Changed cells            : " << changedCells_.size() << endl;
    }

    // Sum changedCells over all procs
    label totNChanged = changedCells_.size();

    reduce(totNChanged, sumOp<label>());

    return totNChanged;
}


template<class Type, class TrackingData>
Foam::label Foam::FaceCellWave<Type, TrackingData>::cellToFace()
{
    // Propagate cell to face

    const cellList& cells = mesh_.cells();

    forAll(changedCells_, changedCelli)
    {
        label celli = changedCells_[changedCelli];
        if (!changedCell_[celli])
        {
            FatalErrorInFunction
                << "Cell " << celli << " not marked as having been changed"
                << abort(FatalError);
        }

        const Type& neighbourWallInfo = allCellInfo_[celli];

        // Evaluate all connected faces

        const labelList& faceLabels = cells[celli];
        forAll(faceLabels, faceLabelI)
        {
            label facei = faceLabels[faceLabelI];
            Type& currentWallInfo = allFaceInfo_[facei];

            if (!currentWallInfo.equal(neighbourWallInfo, td_))
            {
                updateFace
                (
                    facei,
                    celli,
                    neighbourWallInfo,
                    propagationTol_,
                    currentWallInfo
                );
            }
        }

        // Reset status of cell
        changedCell_[celli] = false;
    }

    // Handled all changed cells by now
    changedCells_.clear();


    // Transfer across any explicitly provided internal connections
    handleExplicitConnections();

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
        Pout<< " Changed faces            : " << changedFaces_.size() << endl;
    }

    // Sum nChangedFaces over all procs
    label totNChanged = changedFaces_.size();

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
