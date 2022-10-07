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

#include "FvFaceCellWave.H"
#include "processorFvPatch.H"
#include "cyclicFvPatch.H"
#include "CompactListList.H"
#include "OPstream.H"
#include "IPstream.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type, class TrackingData>
const Foam::scalar Foam::FvFaceCellWave<Type, TrackingData>::geomTol_ = 1e-6;

template<class Type, class TrackingData>
Foam::scalar Foam::FvFaceCellWave<Type, TrackingData>::propagationTol_ = 1e-2;

template<class Type, class TrackingData>
int Foam::FvFaceCellWave<Type, TrackingData>::defaultTrackingData_ = -1;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type, class TrackingData>
template<class ListList>
Foam::labelList
Foam::FvFaceCellWave<Type, TrackingData>::listListSizes
(
    const ListList& ll
)
{
    labelList s(ll.size());
    forAll(ll, i)
    {
        s[i] = ll[i].size();
    }
    return s;
}


template<class Type, class TrackingData>
template<class ListList, class Value>
ListList Foam::FvFaceCellWave<Type, TrackingData>::sizesListList
(
    const labelList& s,
    const Value& value
)
{
    ListList ll(s.size());
    forAll(s, i)
    {
        ll[i] = typename ListList::value_type(s[i], value);
    }
    return ll;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type, class TrackingData>
const Type& Foam::FvFaceCellWave<Type, TrackingData>::faceInfo
(
    const labelPair& patchAndFacei
) const
{
    return
        patchAndFacei.first() == -1
      ? internalFaceInfo_[patchAndFacei.second()]
      : patchFaceInfo_[patchAndFacei.first()][patchAndFacei.second()];
}


template<class Type, class TrackingData>
Type& Foam::FvFaceCellWave<Type, TrackingData>::faceInfo
(
    const labelPair& patchAndFacei
)
{
    return
        patchAndFacei.first() == -1
      ? internalFaceInfo_[patchAndFacei.second()]
      : patchFaceInfo_[patchAndFacei.first()][patchAndFacei.second()];
}


template<class Type, class TrackingData>
bool Foam::FvFaceCellWave<Type, TrackingData>::faceChanged
(
    const labelPair& patchAndFacei
)
const
{
    return
        patchAndFacei.first() == -1
      ? internalFaceChanged_[patchAndFacei.second()]
      : patchFaceChanged_[patchAndFacei.first()][patchAndFacei.second()];
}


template<class Type, class TrackingData>
Foam::PackedBoolList::iteratorBase
Foam::FvFaceCellWave<Type, TrackingData>::faceChanged
(
    const labelPair& patchAndFacei
)
{
    return
        patchAndFacei.first() == -1
      ? internalFaceChanged_[patchAndFacei.second()]
      : patchFaceChanged_[patchAndFacei.first()][patchAndFacei.second()];
}


template<class Type, class TrackingData>
bool Foam::FvFaceCellWave<Type, TrackingData>::updateCell
(
    const label celli,
    const labelPair& neighbourPatchAndFacei,
    const Type& neighbourInfo,
    const scalar tol,
    Type& cellInfo
)
{
    bool propagate =
        cellInfo.updateCell
        (
            mesh_,
            celli,
            neighbourPatchAndFacei,
            neighbourInfo,
            tol,
            td_
        );

    if (propagate)
    {
        if (!cellChanged_[celli])
        {
            cellChanged_[celli] = true;
            changedCells_.append(celli);
        }
    }

    return propagate;
}


template<class Type, class TrackingData>
bool Foam::FvFaceCellWave<Type, TrackingData>::updateFace
(
    const labelPair& patchAndFacei,
    const label neighbourCelli,
    const Type& neighbourInfo,
    const scalar tol,
    Type& faceInfo
)
{
    bool propagate =
        faceInfo.updateFace
        (
            mesh_,
            patchAndFacei,
            neighbourCelli,
            neighbourInfo,
            tol,
            td_
        );

    if (propagate)
    {
        PackedBoolList::iteratorBase changed = faceChanged(patchAndFacei);

        if (!changed)
        {
            changed = true;
            changedPatchAndFaces_.append(patchAndFacei);
        }
    }

    return propagate;
}


template<class Type, class TrackingData>
bool Foam::FvFaceCellWave<Type, TrackingData>::updateFace
(
    const labelPair& patchAndFacei,
    const Type& neighbourInfo,
    const scalar tol,
    Type& faceInfo
)
{
    bool propagate =
        faceInfo.updateFace
        (
            mesh_,
            patchAndFacei,
            neighbourInfo,
            tol,
            td_
        );

    if (propagate)
    {
        PackedBoolList::iteratorBase changed = faceChanged(patchAndFacei);

        if (!changed)
        {
            changed = true;
            changedPatchAndFaces_.append(patchAndFacei);
        }
    }

    return propagate;
}


template<class Type, class TrackingData>
void Foam::FvFaceCellWave<Type, TrackingData>::checkCyclic
(
    const fvPatch& patch
) const
{
    const cyclicFvPatch& nbrPatch =
        refCast<const cyclicFvPatch>(patch).nbrPatch();

    forAll(patch, patchFacei)
    {
        const Type& info = faceInfo({patch.index(), patchFacei});
        const Type& nbrInfo = faceInfo({nbrPatch.index(), patchFacei});

        const bool changed = faceChanged({patch.index(), patchFacei});
        const bool nbrChanged = faceChanged({nbrPatch.index(), patchFacei});

        if (!info.sameGeometry(mesh_, nbrInfo, geomTol_, td_))
        {
            FatalErrorInFunction
                << "   faceInfo:" << info
                << "   otherfaceInfo:" << nbrInfo
                << abort(FatalError);
        }

        if (changed != nbrChanged)
        {
            FatalErrorInFunction
                << "   faceInfo:" << info
                << "   otherfaceInfo:" << nbrInfo
                << "   changedFace:" << changed
                << "   otherchangedFace:" << nbrChanged
                << abort(FatalError);
        }
    }
}


template<class Type, class TrackingData>
template<class PatchType>
bool Foam::FvFaceCellWave<Type, TrackingData>::hasPatch() const
{
    forAll(mesh_.boundary(), patchi)
    {
        if (isA<PatchType>(mesh_.boundary()[patchi]))
        {
            return true;
        }
    }

    return false;
}


template<class Type, class TrackingData>
void Foam::FvFaceCellWave<Type, TrackingData>::mergeFaceInfo
(
    const fvPatch& patch,
    const label nFaces,
    const labelList& changedPatchFaces,
    const List<Type>& changedPatchFacesInfo
)
{
    for (label changedFacei = 0; changedFacei < nFaces; changedFacei ++)
    {
        const label patchFacei = changedPatchFaces[changedFacei];

        Type& info = faceInfo({patch.index(), patchFacei});
        const Type& neighbourInfo = changedPatchFacesInfo[changedFacei];

        if (!info.equal(neighbourInfo, td_))
        {
            updateFace
            (
                {patch.index(), patchFacei},
                neighbourInfo,
                propagationTol_,
                info
            );
        }
    }
}


template<class Type, class TrackingData>
Foam::label
Foam::FvFaceCellWave<Type, TrackingData>::getChangedPatchFaces
(
    const fvPatch& patch,
    labelList& changedPatchFaces,
    List<Type>& changedPatchFacesInfo
) const
{
    // Construct compact patchFace change arrays for a single patch.
    // changedPatchFaces in local patch numbering. Return length of arrays.

    label i = 0;

    forAll(patch, patchFacei)
    {
        if (patchFaceChanged_[patch.index()][patchFacei])
        {
            changedPatchFaces[i] = patchFacei;
            changedPatchFacesInfo[i] = faceInfo({patch.index(), patchFacei});
            i ++;
        }
    }

    return i;
}


template<class Type, class TrackingData>
void Foam::FvFaceCellWave<Type, TrackingData>::transform
(
    const fvPatch& patch,
    const label nFaces,
    const labelList& patchFaces,
    const transformer& transform,
    List<Type>& faceInfo
)
{
    for (label i = 0; i < nFaces; i++)
    {
        faceInfo[i].transform(patch, patchFaces[i], transform, td_);
    }
}


template<class Type, class TrackingData>
void Foam::FvFaceCellWave<Type, TrackingData>::handleProcPatches()
{
    // Transfer information to/from neighbouring processors

    // Which patches are processor patches
    const labelList& procPatches = mesh_.globalData().processorPatches();

    // Send all
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    forAll(procPatches, i)
    {
        const label patchi = procPatches[i];

        const processorFvPatch& procPatch =
            refCast<const processorFvPatch>(mesh_.boundary()[patchi]);

        // Allocate buffers
        label nSendFaces;
        labelList sendFaces(procPatch.size());
        List<Type> sendFacesInfo(procPatch.size());

        // Determine which faces changed on current patch
        nSendFaces =
            getChangedPatchFaces
            (
                procPatch,
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

        // Send
        UOPstream toNeighbour(procPatch.neighbProcNo(), pBufs);
        toNeighbour
            << SubList<label>(sendFaces, nSendFaces)
            << SubList<Type>(sendFacesInfo, nSendFaces);
    }

    pBufs.finishedSends();

    // Receive all
    forAll(procPatches, i)
    {
        const label patchi = procPatches[i];

        const processorFvPatch& procPatch =
            refCast<const processorFvPatch>(mesh_.boundary()[patchi]);

        // Allocate buffers
        labelList receiveFaces;
        List<Type> receiveFacesInfo;

        // Receive
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

        // Transform info across the interface
        transform
        (
            procPatch,
            receiveFaces.size(),
            receiveFaces,
            procPatch.transform(),
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
void
Foam::FvFaceCellWave<Type, TrackingData>::handleCyclicPatches()
{
    // Transfer information across cyclics

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];

        if (!isA<cyclicFvPatch>(patch)) continue;

        const cyclicFvPatch& cycPatch = refCast<const cyclicFvPatch>(patch);
        const cyclicFvPatch& nbrCycPatch = cycPatch.nbrPatch();

        // Allocate buffers
        label nReceiveFaces;
        labelList receiveFaces(patch.size());
        List<Type> receiveFacesInfo(patch.size());

        // Determine which faces changed
        nReceiveFaces =
            getChangedPatchFaces
            (
                nbrCycPatch,
                receiveFaces,
                receiveFacesInfo
            );

        if (debug & 2)
        {
            Pout<< " Cyclic patch " << patchi << ' ' << cycPatch.name()
                << "  Changed : " << nReceiveFaces
                << endl;
        }

        // Transform info across the interface
        transform
        (
            cycPatch,
            nReceiveFaces,
            receiveFaces,
            cycPatch.transform(),
            receiveFacesInfo
        );

        // Merge into existing faces
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class TrackingData>
Foam::FvFaceCellWave<Type, TrackingData>::FvFaceCellWave
(
    const fvMesh& mesh,
    List<Type>& internalFaceInfo,
    List<List<Type>>& patchFaceInfo,
    List<Type>& cellInfo,
    TrackingData& td
)
:
    mesh_(mesh),
    internalFaceInfo_(internalFaceInfo),
    patchFaceInfo_(patchFaceInfo),
    cellInfo_(cellInfo),
    td_(td),
    internalFaceChanged_
    (
        mesh_.nInternalFaces(),
        false
    ),
    patchFaceChanged_
    (
        sizesListList<List<PackedBoolList>, bool>
        (
            listListSizes(mesh_.boundary()),
            false
        )
    ),
    cellChanged_
    (
        mesh_.nCells(),
        false
    ),
    changedPatchAndFaces_(mesh_.nInternalFaces()),
    changedCells_(mesh_.nCells()),
    hasCyclicPatches_(hasPatch<cyclicFvPatch>())
{
    if
    (
        internalFaceInfo.size() != mesh_.nInternalFaces()
     || listListSizes(patchFaceInfo) != listListSizes(mesh_.boundary())
     || cellInfo.size() != mesh_.nCells()
    )
    {
        FatalErrorInFunction
            << "Storage does not match the size of the mesh:" << endl
            << "          internalFaceInfo:" << internalFaceInfo.size() << endl
            << "    mesh_.nInternalFaces():" << mesh_.nInternalFaces() << endl
            << "          patchFaceInfo[i]:" << listListSizes(patchFaceInfo)
            << endl
            << "mesh_.boundary()[i].size():" << listListSizes(mesh_.boundary())
            << endl
            << "                  cellInfo:" << cellInfo.size() << endl
            << "            mesh_.nCells():" << mesh_.nCells()
            << exit(FatalError);
    }
}


template<class Type, class TrackingData>
Foam::FvFaceCellWave<Type, TrackingData>::FvFaceCellWave
(
    const fvMesh& mesh,
    const List<labelPair>& initialChangedPatchAndFaces,
    const List<Type>& initialChangedFacesInfo,
    List<Type>& internalFaceInfo,
    List<List<Type>>& patchFaceInfo,
    List<Type>& cellInfo,
    const label maxIter,
    TrackingData& td
)
:
    FvFaceCellWave
    (
        mesh,
        internalFaceInfo,
        patchFaceInfo,
        cellInfo,
        td
    )
{
    // Copy initial changed faces data
    setFaceInfo(initialChangedPatchAndFaces, initialChangedFacesInfo);

    // Iterate until nothing changes
    const label iter = iterate(maxIter);

    // Error if incomplete
    if ((maxIter > 0) && (iter >= maxIter))
    {
        FatalErrorInFunction
            << "Maximum number of iterations reached. Increase maxIter." << endl
            << "                     maxIter:" << maxIter << endl
            << "        changedCells_.size():" << changedCells_.size() << endl
            << "changedPatchAndFaces_.size():" << changedPatchAndFaces_.size()
            << endl << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class TrackingData>
void Foam::FvFaceCellWave<Type, TrackingData>::setFaceInfo
(
    const List<labelPair>& changedPatchAndFaces,
    const List<Type>& changedFacesInfo
)
{
    forAll(changedPatchAndFaces, i)
    {
        const labelPair& patchAndFacei = changedPatchAndFaces[i];

        Type& info = faceInfo(patchAndFacei);
        PackedBoolList::iteratorBase changed = faceChanged(patchAndFacei);

        // Copy info for this face
        info = changedFacesInfo[i];

        // Mark this face as changed
        changed = true;
        changedPatchAndFaces_.append(patchAndFacei);
    }
}


template<class Type, class TrackingData>
Foam::label Foam::FvFaceCellWave<Type, TrackingData>::faceToCell()
{
    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();

    forAll(changedPatchAndFaces_, changedFacei)
    {
        const labelPair& patchAndFacei = changedPatchAndFaces_[changedFacei];
        const label patchi = patchAndFacei.first();
        const label facei = patchAndFacei.second();

        if (!faceChanged(patchAndFacei))
        {
            FatalErrorInFunction
                << "Patch and face " << patchAndFacei
                << " not marked as having been changed"
                << abort(FatalError);
        }

        const Type& info = faceInfo(patchAndFacei);

        // Propagate to the owner
        {
            const label ownerCelli =
                patchi == -1
              ? owner[facei]
              : mesh_.boundary()[patchi].faceCells()[facei];

            Type& ownerInfo = cellInfo_[ownerCelli];

            if (!ownerInfo.equal(info, td_))
            {
                updateCell
                (
                    ownerCelli,
                    {patchi, facei},
                    info,
                    propagationTol_,
                    ownerInfo
                );
            }
        }

        // Propagate to the neighbour
        if (patchi == -1)
        {
            const label neighbourCelli = neighbour[facei];

            Type& neighbourInfo = cellInfo_[neighbourCelli];

            if (!neighbourInfo.equal(info, td_))
            {
                updateCell
                (
                    neighbourCelli,
                    {patchi, facei},
                    info,
                    propagationTol_,
                    neighbourInfo
                );
            }
        }

        // Reset status of face
        faceChanged(patchAndFacei) = false;
    }

    // Handled all changed faces by now
    changedPatchAndFaces_.clear();

    if (debug & 2)
    {
        Pout<< " Changed cells            : " << changedCells_.size() << endl;
    }

    return returnReduce(changedCells_.size(), sumOp<label>());
}


template<class Type, class TrackingData>
Foam::label Foam::FvFaceCellWave<Type, TrackingData>::cellToFace()
{
    const cellList& cells = mesh_.cells();

    forAll(changedCells_, changedCelli)
    {
        const label celli = changedCells_[changedCelli];

        if (!cellChanged_[celli])
        {
            FatalErrorInFunction
                << "Cell " << celli << " not marked as having been changed"
                << abort(FatalError);
        }

        const Type& info = cellInfo_[celli];

        // Propagate to connected faces
        forAll(cells[celli], cellFacei)
        {
            // Get the patch (if any) and face index
            label polyFacei = cells[celli][cellFacei];

            // Get the FV patches and faces associated with this poly face
            labelUList patches, faces;
            if (polyFacei < mesh_.nInternalFaces())
            {
                static label noPatchi = -1;
                patches.shallowCopy(labelUList(&noPatchi, 1));
                faces.shallowCopy(labelUList(&polyFacei, 1));
            }
            else
            {
                const label polyBFacei = polyFacei - mesh_.nInternalFaces();
                patches.shallowCopy(mesh_.polyBFacePatches()[polyBFacei]);
                faces.shallowCopy(mesh_.polyBFacePatchFaces()[polyBFacei]);
            }

            // Propagate into the connected FV faces
            forAll(patches, i)
            {
                Type& connectedInfo = faceInfo({patches[i], faces[i]});

                if (!connectedInfo.equal(info, td_))
                {
                    updateFace
                    (
                        {patches[i], faces[i]},
                        celli,
                        info,
                        propagationTol_,
                        connectedInfo
                    );
                }
            }
        }

        // Reset status of cell
        cellChanged_[celli] = false;
    }

    // Handled all changed cells by now
    changedCells_.clear();

    if (hasCyclicPatches_)
    {
        // Transfer changed faces across cyclics
        handleCyclicPatches();
    }

    if (Pstream::parRun())
    {
        // Transfer changed faces from neighbouring processors.
        handleProcPatches();
    }

    if (debug & 2)
    {
        Pout<< " Changed faces            : "
            << changedPatchAndFaces_.size() << endl;
    }

    return returnReduce(changedPatchAndFaces_.size(), sumOp<label>());
}


template<class Type, class TrackingData>
Foam::label Foam::FvFaceCellWave<Type, TrackingData>::iterate
(
    const label maxIter
)
{
    if (hasCyclicPatches_)
    {
        // Transfer changed faces across cyclics
        handleCyclicPatches();
    }

    if (Pstream::parRun())
    {
        // Transfer changed faces from neighbouring processors.
        handleProcPatches();
    }

    label iter = 0;

    while (iter < maxIter)
    {
        if (debug) Info<< " Iteration " << iter << endl;

        label nCells = faceToCell();

        if (debug) Info<< " Total changed cells      : " << nCells << endl;

        if (nCells == 0)
        {
            break;
        }

        label nFaces = cellToFace();

        if (debug) Info<< " Total changed faces      : " << nFaces << nl;

        if (nFaces == 0)
        {
            break;
        }

        ++iter;
    }

    return iter;
}


// ************************************************************************* //
