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

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "slicedVolFields.H"
#include "slicedSurfaceFields.H"
#include "SubField.H"
#include "demandDrivenData.H"
#include "fvMeshLduAddressing.H"
#include "fvMeshTopoChanger.H"
#include "fvMeshDistributor.H"
#include "fvMeshMover.H"
#include "fvMeshStitcher.H"
#include "nonConformalFvPatch.H"
#include "polyFacesFvsPatchLabelField.H"
#include "nonConformalPolyFacesFvsPatchLabelField.H"
#include "polyTopoChangeMap.H"
#include "MapFvFields.H"
#include "fvMeshMapper.H"
#include "pointMesh.H"
#include "pointMeshMapper.H"
#include "MapPointField.H"
#include "meshObjects.H"
#include "HashPtrTable.H"
#include "CompactListList.H"

#include "fvcSurfaceIntegrate.H"
#include "fvcReconstruct.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMesh, 0);
}

const Foam::HashSet<Foam::word> Foam::fvMesh::geometryFields
{
    "Vc",
    "Vc0",
    "Vc00",
    "Sf",
    "magSf",
    "Cc",
    "Cf",
    "meshPhi"
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMesh::clearGeomNotOldVol()
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << "clearGeomNotOldVol" << endl;
    }

    meshObjects::clearUpto
    <
        fvMesh,
        GeometricMeshObject,
        MoveableMeshObject
    >(*this);

    meshObjects::clearUpto
    <
        lduMesh,
        GeometricMeshObject,
        MoveableMeshObject
    >(*this);

    deleteDemandDrivenData(VPtr_);
    deleteDemandDrivenData(SfSlicePtr_);
    deleteDemandDrivenData(SfPtr_);
    deleteDemandDrivenData(magSfSlicePtr_);
    deleteDemandDrivenData(magSfPtr_);
    deleteDemandDrivenData(CSlicePtr_);
    deleteDemandDrivenData(CPtr_);
    deleteDemandDrivenData(CfSlicePtr_);
    deleteDemandDrivenData(CfPtr_);
}


void Foam::fvMesh::updateGeomNotOldVol()
{
    bool haveV = (VPtr_ != nullptr);
    bool haveSf = (SfSlicePtr_ != nullptr || SfPtr_ != nullptr);
    bool haveMagSf = (magSfSlicePtr_ != nullptr || magSfPtr_ != nullptr);
    bool haveCP = (CSlicePtr_ != nullptr || CPtr_ != nullptr);
    bool haveCf = (CfSlicePtr_ != nullptr || CfPtr_ != nullptr);

    clearGeomNotOldVol();

    // Now recreate the fields
    if (haveV)
    {
        (void)V();
    }
    if (haveSf)
    {
        (void)Sf();
    }
    if (haveMagSf)
    {
        (void)magSf();
    }
    if (haveCP)
    {
        (void)C();
    }
    if (haveCf)
    {
        (void)Cf();
    }
}


void Foam::fvMesh::clearGeom()
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << "Clearing geometric data" << endl;
    }

    clearGeomNotOldVol();

    deleteDemandDrivenData(phiPtr_);
    deleteDemandDrivenData(V0Ptr_);
    deleteDemandDrivenData(V00Ptr_);
}


void Foam::fvMesh::clearAddressing(const bool isMeshUpdate)
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << "isMeshUpdate: " << isMeshUpdate << endl;
    }

    if (isMeshUpdate)
    {
        // Part of a mesh update. Keep meshObjects that have an topoChange
        // callback
        meshObjects::clearUpto
        <
            fvMesh,
            TopologicalMeshObject,
            UpdateableMeshObject
        >
        (
            *this
        );
        meshObjects::clearUpto
        <
            lduMesh,
            TopologicalMeshObject,
            UpdateableMeshObject
        >
        (
            *this
        );
    }
    else
    {
        meshObjects::clear<fvMesh, TopologicalMeshObject>(*this);
        meshObjects::clear<lduMesh, TopologicalMeshObject>(*this);
    }

    deleteDemandDrivenData(lduPtr_);
    deleteDemandDrivenData(polyFacesBfPtr_);
    deleteDemandDrivenData(polyBFaceOffsetsPtr_);
    deleteDemandDrivenData(polyBFaceOffsetPatchesPtr_);
    deleteDemandDrivenData(polyBFaceOffsetPatchFacesPtr_);
    deleteDemandDrivenData(polyBFacePatchesPtr_);
    deleteDemandDrivenData(polyBFacePatchFacesPtr_);
}


void Foam::fvMesh::storeOldVol(const scalarField& V)
{
    if (curTimeIndex_ < time().timeIndex())
    {
        if (debug)
        {
            Pout<< FUNCTION_NAME
                << " Storing old time volumes since from time " << curTimeIndex_
                << " and time now " << time().timeIndex()
                << " V:" << V.size()
                << endl;
        }

        if (V00Ptr_ && V0Ptr_)
        {
            // Copy V0 into V00 storage
            *V00Ptr_ = *V0Ptr_;
        }

        if (V0Ptr_)
        {
            // Copy V into V0 storage
            V0Ptr_->scalarField::operator=(V);
        }
        else
        {
            // Allocate V0 storage, fill with V
            V0Ptr_ = new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    "Vc0",
                    time().name(),
                    *this,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    true
                ),
                *this,
                dimVolume
            );
            scalarField& V0 = *V0Ptr_;
            // Note: V0 now sized with current mesh, not with (potentially
            //       different size) V.
            V0.setSize(V.size());
            V0 = V;
        }

        if (debug)
        {
            Pout<< FUNCTION_NAME
                << " Stored old time volumes V0:" << V0Ptr_->size()
                << endl;
            if (V00Ptr_)
            {
                Pout<< FUNCTION_NAME
                    << " Stored oldold time volumes V00:" << V00Ptr_->size()
                    << endl;
            }
        }
    }
}


void Foam::fvMesh::storeOldTimeFields()
{
    storeOldTimeFields<PointField>();
    storeOldTimeFields<VolField>();
    storeOldTimeFields<SurfaceField>();
}


void Foam::fvMesh::nullOldestTimeFields()
{
    nullOldestTimeFields<PointField>();
    nullOldestTimeFields<VolField>();
    nullOldestTimeFields<SurfaceField>();
}


void Foam::fvMesh::clearOut()
{
    clearGeom();

    surfaceInterpolation::clearOut();

    clearAddressing();

    polyMesh::clearOut();
}


Foam::wordList Foam::fvMesh::polyFacesPatchTypes() const
{
    wordList wantedPatchTypes
    (
        boundary().size(),
        polyFacesFvsPatchLabelField::typeName
    );

    forAll(boundary(), patchi)
    {
        const fvPatch& fvp = boundary()[patchi];

        if (isA<nonConformalFvPatch>(fvp))
        {
            wantedPatchTypes[patchi] =
                nonConformalPolyFacesFvsPatchLabelField::typeName;
        }
    }

    return wantedPatchTypes;
}


Foam::surfaceLabelField::Boundary& Foam::fvMesh::polyFacesBfRef()
{
    if (!polyFacesBfPtr_)
    {
        polyFacesBf();
    }

    return *polyFacesBfPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMesh::fvMesh
(
    const IOobject& io,
    const bool changers,
    const stitchType stitch
)
:
    polyMesh(io),
    surfaceInterpolation(*this),
    data(static_cast<const objectRegistry&>(*this)),
    boundary_(*this, boundaryMesh()),
    stitcher_(fvMeshStitcher::New(*this, changers).ptr()),
    topoChanger_(nullptr),
    distributor_(nullptr),
    mover_(nullptr),
    lduPtr_(nullptr),
    polyFacesBfPtr_(nullptr),
    polyBFaceOffsetsPtr_(nullptr),
    polyBFaceOffsetPatchesPtr_(nullptr),
    polyBFaceOffsetPatchFacesPtr_(nullptr),
    polyBFacePatchesPtr_(nullptr),
    polyBFacePatchFacesPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    VPtr_(nullptr),
    V0Ptr_(nullptr),
    V00Ptr_(nullptr),
    SfSlicePtr_(nullptr),
    SfPtr_(nullptr),
    magSfSlicePtr_(nullptr),
    magSfPtr_(nullptr),
    CSlicePtr_(nullptr),
    CPtr_(nullptr),
    CfSlicePtr_(nullptr),
    CfPtr_(nullptr),
    phiPtr_(nullptr)
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << "Constructing fvMesh from IOobject" << endl;
    }

    // Stitch or Re-stitch if necessary
    if (stitch != stitchType::none)
    {
        stitcher_->connect(false, stitch == stitchType::geometric, true);
    }

    // Construct changers
    if (changers)
    {
        topoChanger_.set(fvMeshTopoChanger::New(*this).ptr());
        distributor_.set(fvMeshDistributor::New(*this).ptr());
        mover_.set(fvMeshMover::New(*this).ptr());

        // Check the existence of the cell volumes and read if present
        // and set the storage of V00
        if (fileHandler().isFile(time().timePath()/"Vc0"))
        {
            V0Ptr_ = new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    "Vc0",
                    time().name(),
                    *this,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    true
                ),
                *this
            );

            V00();
        }

        // Check the existence of the mesh fluxes and read if present
        if (fileHandler().isFile(time().timePath()/"meshPhi"))
        {
            phiPtr_ = new surfaceScalarField
            (
                IOobject
                (
                    "meshPhi",
                    time().name(),
                    *this,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    true
                ),
                *this
            );
        }
    }
}


Foam::fvMesh::fvMesh
(
    const IOobject& io,
    pointField&& points,
    const cellShapeList& shapes,
    const faceListList& boundaryFaces,
    const wordList& boundaryPatchNames,
    const PtrList<dictionary>& boundaryDicts,
    const word& defaultBoundaryPatchName,
    const word& defaultBoundaryPatchType,
    const bool syncPar
)
:
    polyMesh
    (
        io,
        std::move(points),
        shapes,
        boundaryFaces,
        boundaryPatchNames,
        boundaryDicts,
        defaultBoundaryPatchName,
        defaultBoundaryPatchType,
        syncPar
    ),
    surfaceInterpolation(*this),
    data(static_cast<const objectRegistry&>(*this)),
    boundary_(*this, boundaryMesh()),
    stitcher_(nullptr),
    topoChanger_(nullptr),
    distributor_(nullptr),
    mover_(nullptr),
    lduPtr_(nullptr),
    polyFacesBfPtr_(nullptr),
    polyBFaceOffsetsPtr_(nullptr),
    polyBFaceOffsetPatchesPtr_(nullptr),
    polyBFaceOffsetPatchFacesPtr_(nullptr),
    polyBFacePatchesPtr_(nullptr),
    polyBFacePatchFacesPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    VPtr_(nullptr),
    V0Ptr_(nullptr),
    V00Ptr_(nullptr),
    SfSlicePtr_(nullptr),
    SfPtr_(nullptr),
    magSfSlicePtr_(nullptr),
    magSfPtr_(nullptr),
    CSlicePtr_(nullptr),
    CPtr_(nullptr),
    CfSlicePtr_(nullptr),
    CfPtr_(nullptr),
    phiPtr_(nullptr)
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << "Constructing fvMesh from cellShapes" << endl;
    }
}


Foam::fvMesh::fvMesh
(
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    labelList&& allOwner,
    labelList&& allNeighbour,
    const bool syncPar
)
:
    polyMesh
    (
        io,
        std::move(points),
        std::move(faces),
        std::move(allOwner),
        std::move(allNeighbour),
        syncPar
    ),
    surfaceInterpolation(*this),
    data(static_cast<const objectRegistry&>(*this)),
    boundary_(*this, boundaryMesh()),
    stitcher_(nullptr),
    topoChanger_(nullptr),
    distributor_(nullptr),
    mover_(nullptr),
    lduPtr_(nullptr),
    polyFacesBfPtr_(nullptr),
    polyBFaceOffsetsPtr_(nullptr),
    polyBFaceOffsetPatchesPtr_(nullptr),
    polyBFaceOffsetPatchFacesPtr_(nullptr),
    polyBFacePatchesPtr_(nullptr),
    polyBFacePatchFacesPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    VPtr_(nullptr),
    V0Ptr_(nullptr),
    V00Ptr_(nullptr),
    SfSlicePtr_(nullptr),
    SfPtr_(nullptr),
    magSfSlicePtr_(nullptr),
    magSfPtr_(nullptr),
    CSlicePtr_(nullptr),
    CPtr_(nullptr),
    CfSlicePtr_(nullptr),
    CfPtr_(nullptr),
    phiPtr_(nullptr)
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << "Constructing fvMesh from components" << endl;
    }
}


Foam::fvMesh::fvMesh
(
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    cellList&& cells,
    const bool syncPar
)
:
    polyMesh
    (
        io,
        std::move(points),
        std::move(faces),
        std::move(cells),
        syncPar
    ),
    surfaceInterpolation(*this),
    data(static_cast<const objectRegistry&>(*this)),
    boundary_(*this),
    stitcher_(nullptr),
    topoChanger_(nullptr),
    distributor_(nullptr),
    mover_(nullptr),
    lduPtr_(nullptr),
    polyFacesBfPtr_(nullptr),
    polyBFaceOffsetsPtr_(nullptr),
    polyBFaceOffsetPatchesPtr_(nullptr),
    polyBFaceOffsetPatchFacesPtr_(nullptr),
    polyBFacePatchesPtr_(nullptr),
    polyBFacePatchFacesPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    VPtr_(nullptr),
    V0Ptr_(nullptr),
    V00Ptr_(nullptr),
    SfSlicePtr_(nullptr),
    SfPtr_(nullptr),
    magSfSlicePtr_(nullptr),
    magSfPtr_(nullptr),
    CSlicePtr_(nullptr),
    CPtr_(nullptr),
    CfSlicePtr_(nullptr),
    CfPtr_(nullptr),
    phiPtr_(nullptr)
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << "Constructing fvMesh from components" << endl;
    }
}


Foam::fvMesh::fvMesh(fvMesh&& mesh)
:
    polyMesh(Foam::move(mesh)),
    surfaceInterpolation(Foam::move(mesh)),
    data(static_cast<const objectRegistry&>(*this)),
    boundary_(Foam::move(mesh.boundary_)),
    stitcher_(Foam::move(mesh.stitcher_)),
    topoChanger_(Foam::move(mesh.topoChanger_)),
    distributor_(Foam::move(mesh.distributor_)),
    mover_(Foam::move(mesh.mover_)),
    lduPtr_(Foam::move(mesh.lduPtr_)),
    polyFacesBfPtr_(Foam::move(mesh.polyFacesBfPtr_)),
    polyBFaceOffsetsPtr_(Foam::move(mesh.polyBFaceOffsetsPtr_)),
    polyBFaceOffsetPatchesPtr_(Foam::move(mesh.polyBFaceOffsetPatchesPtr_)),
    polyBFaceOffsetPatchFacesPtr_
    (
        Foam::move(mesh.polyBFaceOffsetPatchFacesPtr_)
    ),
    polyBFacePatchesPtr_(Foam::move(mesh.polyBFacePatchesPtr_)),
    polyBFacePatchFacesPtr_(Foam::move(mesh.polyBFacePatchFacesPtr_)),
    curTimeIndex_(mesh.curTimeIndex_),
    VPtr_(Foam::move(mesh.VPtr_)),
    V0Ptr_(Foam::move(mesh.V0Ptr_)),
    V00Ptr_(Foam::move(mesh.V00Ptr_)),
    SfSlicePtr_(Foam::move(mesh.SfSlicePtr_)),
    SfPtr_(Foam::move(mesh.SfPtr_)),
    magSfSlicePtr_(Foam::move(mesh.magSfSlicePtr_)),
    magSfPtr_(Foam::move(mesh.magSfPtr_)),
    CSlicePtr_(Foam::move(mesh.CSlicePtr_)),
    CPtr_(Foam::move(mesh.CPtr_)),
    CfSlicePtr_(Foam::move(mesh.CfSlicePtr_)),
    CfPtr_(Foam::move(mesh.CfPtr_)),
    phiPtr_(Foam::move(mesh.phiPtr_))
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << "Moving fvMesh" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMesh::~fvMesh()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMesh::topoChanging() const
{
    return topoChanger_.valid() && topoChanger_->dynamic();
}


bool Foam::fvMesh::dynamic() const
{
    return
        (topoChanger_.valid() && topoChanger_->dynamic())
     || (mover_.valid() && mover_->dynamic());
}


bool Foam::fvMesh::update()
{
    if
    (
        stitcher_->stitches()
     || topoChanger_->dynamic()
     || distributor_->dynamic()
    )
    {
        nullOldestTimeFields();
    }

    if (!conformal()) stitcher_->disconnect(true, true);

    const bool hasV00 = V00Ptr_;
    deleteDemandDrivenData(V00Ptr_);

    if (!hasV00)
    {
        deleteDemandDrivenData(V0Ptr_);
    }

    // Set topoChanged_ false before any mesh change
    topoChanged_ = false;
    bool updated = topoChanger_->update();
    topoChanged_ = updated;

    updated = distributor_->update() || updated;

    if (hasV00)
    {
        // If V00 had been set reset to the mapped V0 prior to mesh-motion
        V00();
    }

    return updated;
}


bool Foam::fvMesh::move()
{
    if (!conformal()) stitcher_->disconnect(true, true);

    if (curTimeIndex_ < time().timeIndex() && stitcher_->stitches())
    {
        // Store all old-time fields. If we don't do this then we risk
        // triggering a store in the middle of mapping and potentially
        // overwriting a mapped old-time field with a not-yet-mapped new-time
        // field.
        storeOldTimeFields();
    }

    // Do not set moving false
    // Once the mesh starts moving it is considered to be moving
    // for the rest of the run
    const bool moved = mover_->update();

    curTimeIndex_ = time().timeIndex();

    stitcher_->connect(true, true, false);

    return moved;
}


void Foam::fvMesh::addFvPatches
(
    const List<polyPatch*> & p,
    const bool validBoundary
)
{
    if (boundary().size())
    {
        FatalErrorInFunction
            << " boundary already exists"
            << abort(FatalError);
    }

    // first add polyPatches
    addPatches(p, validBoundary);
    boundary_.addPatches(boundaryMesh());
}


void Foam::fvMesh::removeFvBoundary()
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << "Removing boundary patches." << endl;
    }

    // Remove fvBoundaryMesh data first.
    boundary_.clear();
    boundary_.setSize(0);
    polyMesh::removeBoundary();

    clearOut();
}


void Foam::fvMesh::swap(fvMesh& otherMesh)
{
    // Clear the sliced fields
    clearGeom();

    // Clear the current volume and other geometry factors
    surfaceInterpolation::clearOut();

    // Clear any non-updateable addressing
    clearAddressing(true);

    polyMesh::swap(otherMesh);

    auto updatePatches = []
    (
        const polyPatchList& patches,
        fvBoundaryMesh& boundaryMesh
    )
    {
        boundaryMesh.setSize(patches.size());

        forAll(patches, patchi)
        {
            // Construct new processor patches, as the decomposition may have
            // changed. Leave other patches as is.

            if (isA<processorPolyPatch>(patches[patchi]))
            {
                boundaryMesh.set
                (
                    patchi,
                    fvPatch::New
                    (
                        patches[patchi],
                        boundaryMesh
                    )
                );
            }
        }
    };

    updatePatches(boundaryMesh(), boundary_);
    updatePatches(otherMesh.boundaryMesh(), otherMesh.boundary_);
}


Foam::polyMesh::readUpdateState Foam::fvMesh::readUpdate
(
    const stitchType stitch
)
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << "Updating fvMesh.  ";
    }

    polyMesh::readUpdateState state = polyMesh::readUpdate();

    if (state == polyMesh::TOPO_PATCH_CHANGE)
    {
        boundary_.readUpdate(boundaryMesh());
    }

    if
    (
        stitcher_.valid()
     && stitch != stitchType::none
     && state != polyMesh::UNCHANGED
    )
    {
        stitcher_->disconnect(false, stitch == stitchType::geometric);
    }

    if (state == polyMesh::TOPO_PATCH_CHANGE)
    {
        if (debug)
        {
            Info<< "Boundary and topological update" << endl;
        }

        clearOut();
    }
    else if (state == polyMesh::TOPO_CHANGE)
    {
        if (debug)
        {
            Info<< "Topological update" << endl;
        }

        clearOut();
    }
    else if (state == polyMesh::POINTS_MOVED)
    {
        if (debug)
        {
            Info<< "Point motion update" << endl;
        }

        clearGeom();
    }
    else
    {
        if (debug)
        {
            Info<< "No update" << endl;
        }
    }

    if
    (
        stitcher_.valid()
     && stitch != stitchType::none
     && state != polyMesh::UNCHANGED
    )
    {
        stitcher_->connect(false, stitch == stitchType::geometric, true);
    }

    // If the mesh has been re-stitched with different geometry, then the
    // finite-volume topology has changed
    if
    (
        stitcher_.valid()
     && stitcher_->stitches()
     && state == polyMesh::POINTS_MOVED
    )
    {
        state = polyMesh::TOPO_CHANGE;
    }

    return state;
}


const Foam::fvBoundaryMesh& Foam::fvMesh::boundary() const
{
    return boundary_;
}


const Foam::lduAddressing& Foam::fvMesh::lduAddr() const
{
    if (!lduPtr_)
    {
        lduPtr_ = new fvMeshLduAddressing(*this);
    }

    return *lduPtr_;
}


bool Foam::fvMesh::conformal() const
{
    return !(polyFacesBfPtr_ && SfPtr_);
}


Foam::IOobject Foam::fvMesh::polyFacesBfIO(const IOobject::readOption r) const
{
    return
        IOobject
        (
            "polyFaces",
            pointsInstance(),
            typeName,
            *this,
            r,
            IOobject::NO_WRITE,
            false
        );
}


const Foam::surfaceLabelField::Boundary& Foam::fvMesh::polyFacesBf() const
{
    if (!polyFacesBfPtr_)
    {
        polyFacesBfPtr_ =
            new surfaceLabelField::Boundary
            (
                boundary(),
                surfaceLabelField::null(),
                polyFacesPatchTypes(),
                boundaryMesh().types()
            );
    }

    return *polyFacesBfPtr_;
}


const Foam::UCompactListList<Foam::label>&
Foam::fvMesh::polyBFacePatches() const
{
    if (!polyBFacePatchesPtr_)
    {
        const label nPolyBFaces = nFaces() - nInternalFaces();

        // Count face-poly-bFaces to get the offsets
        polyBFaceOffsetsPtr_ = new labelList(nPolyBFaces + 1, 0);
        labelList& offsets = *polyBFaceOffsetsPtr_;
        forAll(boundary(), patchi)
        {
            forAll(boundary()[patchi], patchFacei)
            {
                const label polyBFacei =
                    (
                        polyFacesBfPtr_
                      ? (*polyFacesBfPtr_)[patchi][patchFacei]
                      : boundary()[patchi].start() + patchFacei
                    )
                  - nInternalFaces();

                offsets[polyBFacei + 1] ++;
            }
        }
        for (label polyBFacei = 0; polyBFacei < nPolyBFaces; ++ polyBFacei)
        {
            offsets[polyBFacei + 1] += offsets[polyBFacei];
        }

        // Set the poly-bFace patches and patch-faces, using the offsets as
        // counters
        polyBFaceOffsetPatchesPtr_ = new labelList(offsets.last());
        polyBFaceOffsetPatchFacesPtr_ = new labelList(offsets.last());
        labelUList& patches = *polyBFaceOffsetPatchesPtr_;
        labelUList& patchFaces = *polyBFaceOffsetPatchFacesPtr_;
        forAll(boundary(), patchi)
        {
            forAll(boundary()[patchi], patchFacei)
            {
                const label polyBFacei =
                    (
                        polyFacesBfPtr_
                      ? (*polyFacesBfPtr_)[patchi][patchFacei]
                      : boundary()[patchi].start() + patchFacei
                    )
                  - nInternalFaces();
                patches[offsets[polyBFacei]] = patchi;
                patchFaces[offsets[polyBFacei]] = patchFacei;
                offsets[polyBFacei] ++;
            }
        }

        // Restore the offsets by removing the count
        for
        (
            label polyBFacei = nPolyBFaces - 1;
            polyBFacei >= 0;
            -- polyBFacei
        )
        {
            offsets[polyBFacei + 1] = offsets[polyBFacei];
        }
        offsets[0] = 0;

        // List-lists
        polyBFacePatchesPtr_ =
            new UCompactListList<label>(offsets, patches);
        polyBFacePatchFacesPtr_ =
            new UCompactListList<label>(offsets, patchFaces);
    }

    return *polyBFacePatchesPtr_;
}


const Foam::UCompactListList<Foam::label>&
Foam::fvMesh::polyBFacePatchFaces() const
{
    if (!polyBFacePatchFacesPtr_)
    {
        polyBFacePatches();
    }

    return *polyBFacePatchFacesPtr_;
}


const Foam::fvMeshStitcher& Foam::fvMesh::stitcher() const
{
    return stitcher_();
}


const Foam::fvMeshTopoChanger& Foam::fvMesh::topoChanger() const
{
    return topoChanger_();
}


const Foam::fvMeshDistributor& Foam::fvMesh::distributor() const
{
    return distributor_();
}


const Foam::fvMeshMover& Foam::fvMesh::mover() const
{
    return mover_();
}


void Foam::fvMesh::mapFields(const polyTopoChangeMap& map)
{
    if (debug)
    {
        Pout<< FUNCTION_NAME
            << " nOldCells:" << map.nOldCells()
            << " nCells:" << nCells()
            << " nOldFaces:" << map.nOldFaces()
            << " nFaces:" << nFaces()
            << endl;
    }


    // We require geometric properties valid for the old mesh
    if
    (
        map.cellMap().size() != nCells()
     || map.faceMap().size() != nFaces()
    )
    {
        FatalErrorInFunction
            << "polyTopoChangeMap does not correspond to the old mesh."
            << " nCells:" << nCells()
            << " cellMap:" << map.cellMap().size()
            << " nOldCells:" << map.nOldCells()
            << " nFaces:" << nFaces()
            << " faceMap:" << map.faceMap().size()
            << " nOldFaces:" << map.nOldFaces()
            << exit(FatalError);
    }

    // Create a fv mapper
    const fvMeshMapper fvMap(*this, map);

    // Map all the volFields in the objectRegistry
    #define mapVolFieldType(Type, nullArg)                                     \
        MapGeometricFields<Type, fvPatchField, fvMeshMapper, volMesh>(fvMap);
    FOR_ALL_FIELD_TYPES(mapVolFieldType);

    // Map all the surfaceFields in the objectRegistry
    #define mapSurfaceFieldType(Type, nullArg)                                 \
        MapGeometricFields<Type, fvsPatchField, fvMeshMapper, surfaceMesh>     \
        (fvMap);
    FOR_ALL_FIELD_TYPES(mapSurfaceFieldType);

    // Map all the dimensionedFields in the objectRegistry
    #define mapVolInternalFieldType(Type, nullArg)                             \
        MapDimensionedFields<Type, fvMeshMapper, volMesh>(fvMap);
    FOR_ALL_FIELD_TYPES(mapVolInternalFieldType);

    if (pointMesh::found(*this))
    {
        // Create the pointMesh mapper
        const pointMeshMapper mapper(pointMesh::New(*this), map);

        #define mapPointFieldType(Type, nullArg)                               \
            MapGeometricFields                                                 \
            <                                                                  \
                Type,                                                          \
                pointPatchField,                                               \
                pointMeshMapper,                                               \
                pointMesh                                                      \
            >                                                                  \
            (mapper);
        FOR_ALL_FIELD_TYPES(mapPointFieldType);
    }
}


void Foam::fvMesh::setPoints(const pointField& p)
{
    polyMesh::setPoints(p);

    clearGeom();

    // Update other local data
    boundary_.movePoints();
    surfaceInterpolation::movePoints();

    meshObjects::movePoints<fvMesh>(*this);
    meshObjects::movePoints<lduMesh>(*this);

    const_cast<Time&>(time()).functionObjects().movePoints(*this);
}


Foam::tmp<Foam::scalarField> Foam::fvMesh::movePoints(const pointField& p)
{
    // Set moving_ true
    // Note: once set it remains true for the rest of the run
    moving_ = true;

    // Store old time volumes if the time has been incremented
    // This will update V0, V00
    if (curTimeIndex_ < time().timeIndex())
    {
        storeOldVol(V());
    }

    if (!phiPtr_)
    {
        // Create mesh motion flux
        phiPtr_ = new surfaceScalarField
        (
            IOobject
            (
                "meshPhi",
                this->time().name(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                true
            ),
            *this,
            dimVolume/dimTime
        );
    }
    else
    {
        // Store old time mesh motion fluxes if the time has been incremented
        if (!topoChanging() && phiPtr_->timeIndex() != time().timeIndex())
        {
            phiPtr_->oldTime();
        }
    }

    surfaceScalarField& phi = *phiPtr_;

    // Move the polyMesh and set the mesh motion fluxes to the swept-volumes

    scalar rDeltaT = 1.0/time().deltaTValue();

    tmp<scalarField> tsweptVols = polyMesh::movePoints(p);

    scalarField& sweptVols = tsweptVols.ref();

    phi.primitiveFieldRef() =
        scalarField::subField(sweptVols, nInternalFaces());
    phi.primitiveFieldRef() *= rDeltaT;

    const fvPatchList& patches = boundary();

    surfaceScalarField::Boundary& phibf = phi.boundaryFieldRef();

    forAll(patches, patchi)
    {
        phibf[patchi] = patches[patchi].patchSlice(sweptVols);
        phibf[patchi] *= rDeltaT;
    }

    // Update or delete the local geometric properties as early as possible so
    // they can be used if necessary. These get recreated here instead of
    // demand driven since they might do parallel transfers which can conflict
    // with when they're actually being used.
    // Note that between above "polyMesh::movePoints(p)" and here nothing
    // should use the local geometric properties.
    updateGeomNotOldVol();


    // Update other local data
    boundary_.movePoints();
    surfaceInterpolation::movePoints();

    meshObjects::movePoints<fvMesh>(*this);
    meshObjects::movePoints<lduMesh>(*this);

    const_cast<Time&>(time()).functionObjects().movePoints(*this);

    return tsweptVols;
}


void Foam::fvMesh::topoChange(const polyTopoChangeMap& map)
{
    // Update polyMesh. This needs to keep volume existent!
    polyMesh::topoChange(map);

    if (VPtr_)
    {
        // Cache old time volumes if they exist and the time has been
        // incremented
        if (V0Ptr_ && !V0Ptr_->registered())
        {
            storeOldVol(map.oldCellVolumes());
        }

        // Few checks
        if (VPtr_ && (VPtr_->size() != map.nOldCells()))
        {
            FatalErrorInFunction
                << "Vc:" << V().size()
                << " not equal to the number of old cells "
                << map.nOldCells()
                << exit(FatalError);
        }

        if (V0Ptr_ && (V0Ptr_->size() != map.nOldCells()))
        {
            FatalErrorInFunction
                << "Vc0:" << V0Ptr_->size()
                << " not equal to the number of old cells "
                << map.nOldCells()
                << exit(FatalError);
        }
    }

    // Clear the sliced fields
    clearGeomNotOldVol();

    // Map the old volume. Just map to new cell labels.
    if (V0Ptr_ && !V0Ptr_->registered())
    {
        const labelList& cellMap = map.cellMap();

        scalarField& V0 = *V0Ptr_;

        scalarField savedV0(V0);
        V0.setSize(nCells());

        forAll(V0, i)
        {
            if (cellMap[i] > -1)
            {
                V0[i] = savedV0[cellMap[i]];
            }
            else
            {
                V0[i] = 0.0;
            }
        }

        // Inject volume of merged cells
        label nMerged = 0;
        forAll(map.reverseCellMap(), oldCelli)
        {
            label index = map.reverseCellMap()[oldCelli];

            if (index < -1)
            {
                label celli = -index-2;

                V0[celli] += savedV0[oldCelli];

                nMerged++;
            }
        }

        if (debug)
        {
            Info<< "Mapping old time volume V0. Merged "
                << nMerged << " out of " << nCells() << " cells" << endl;
        }
    }

    // Map all fields
    mapFields(map);

    // Clear the current volume and other geometry factors
    surfaceInterpolation::clearOut();

    // Clear any non-updateable addressing
    clearAddressing(true);

    meshObjects::topoChange<fvMesh>(*this, map);
    meshObjects::topoChange<lduMesh>(*this, map);

    const_cast<Time&>(time()).functionObjects().topoChange(map);

    if (topoChanger_.valid())
    {
        topoChanger_->topoChange(map);
    }

    if (distributor_.valid())
    {
        distributor_->topoChange(map);
    }

    if (mover_.valid())
    {
        mover_->topoChange(map);
    }
}


void Foam::fvMesh::mapMesh(const polyMeshMap& map)
{
    // Distribute polyMesh data
    polyMesh::mapMesh(map);

    // Clear the sliced fields
    clearGeomNotOldVol();

    // Clear the current volume and other geometry factors
    surfaceInterpolation::clearOut();

    // Clear any non-updateable addressing
    clearAddressing(true);

    meshObjects::mapMesh<fvMesh>(*this, map);
    meshObjects::mapMesh<lduMesh>(*this, map);

    const_cast<Time&>(time()).functionObjects().mapMesh(map);

    topoChanger_->mapMesh(map);
    distributor_->mapMesh(map);
    mover_->mapMesh(map);
}


void Foam::fvMesh::distribute(const polyDistributionMap& map)
{
    // Distribute polyMesh data
    polyMesh::distribute(map);

    // Clear the sliced fields
    clearGeomNotOldVol();

    // Clear the current volume and other geometry factors
    surfaceInterpolation::clearOut();

    // Clear any non-updateable addressing
    clearAddressing(true);

    meshObjects::distribute<fvMesh>(*this, map);
    meshObjects::distribute<lduMesh>(*this, map);

    const_cast<Time&>(time()).functionObjects().distribute(map);

    topoChanger_->distribute(map);
    distributor_->distribute(map);
    mover_->distribute(map);
}


void Foam::fvMesh::conform(const surfaceScalarField& phi)
{
    // Clear the geometry fields
    clearGeomNotOldVol();

    // Clear the current volume and other geometry factors
    surfaceInterpolation::clearOut();

    // Clear any non-updateable addressing
    clearAddressing(true);

    // Modify the mesh fluxes, if necessary
    if (notNull(phi) && phiPtr_)
    {
        for (label i = 0; i <= phi.nOldTimes(); ++ i)
        {
            phiRef().oldTime(i) = phi.oldTime(i);
        }
    }
}


void Foam::fvMesh::unconform
(
    const surfaceLabelField::Boundary& polyFacesBf,
    const surfaceVectorField& Sf,
    const surfaceVectorField& Cf,
    const surfaceScalarField& phi,
    const bool sync
)
{
    // Clear the geometry fields
    clearGeomNotOldVol();

    // Clear the current volume and other geometry factors
    surfaceInterpolation::clearOut();

    // Clear any non-updateable addressing
    clearAddressing(true);

    // Create non-sliced copies of geometry fields
    SfRef();
    magSfRef();
    CRef();
    CfRef();

    // Set the topology
    polyFacesBfRef() == polyFacesBf;

    // Set the face geometry
    SfRef() == Sf;
    magSfRef() == max(mag(Sf), dimensionedScalar(dimArea, rootVSmall));
    CRef().boundaryFieldRef() == Cf.boundaryField();
    CfRef() == Cf;

    // Communicate processor-coupled cell geometry. Cell-centre processor patch
    // fields must contain the (transformed) cell-centre locations on the other
    // side of the coupling. This is so that non-conformal patches can
    // construct weights and deltas without reference to the poly mesh
    // geometry.
    //
    // Note that the initEvaluate/evaluate communication does a transformation,
    // but it is wrong in this instance. A vector field gets transformed as if
    // it were a displacement, but the cell-centres need a positional
    // transformation. That's why there's the un-transform and re-transform bit
    // below just after the evaluate call.
    //
    // This transform handling is a bit of a hack. It would be nicer to have a
    // field attribute which identifies a field as needing a positional
    // transformation, and for it to apply automatically within the coupled
    // patch field. However, at the moment, the cell centres field is the only
    // vol-field containing an absolute position, so the hack is functionally
    // sufficient for now.
    if (sync && (Pstream::parRun() || !time().processorCase()))
    {
        volVectorField::Boundary& CBf = CRef().boundaryFieldRef();

        const label nReq = Pstream::nRequests();

        forAll(CBf, patchi)
        {
            if (isA<processorFvPatch>(CBf[patchi].patch()))
            {
                CBf[patchi].initEvaluate(Pstream::defaultCommsType);
            }
        }

        if
        (
            Pstream::parRun()
         && Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
        )
        {
            Pstream::waitRequests(nReq);
        }

        forAll(CBf, patchi)
        {
            if (isA<processorFvPatch>(CBf[patchi].patch()))
            {
                CBf[patchi].evaluate(Pstream::defaultCommsType);

                const transformer& t =
                    refCast<const processorFvPatch>(CBf[patchi].patch())
                   .transform();

                t.invTransform(CBf[patchi], CBf[patchi]);
                t.transformPosition(CBf[patchi], CBf[patchi]);
            }
        }
    }

    // Modify the mesh fluxes, if necessary
    if (notNull(phi) && phiPtr_)
    {
        for (label i = 0; i <= phi.nOldTimes(); ++ i)
        {
            phiRef().oldTime(i) = phi.oldTime(i);
        }
    }
}


void Foam::fvMesh::addPatch
(
    const label insertPatchi,
    const polyPatch& patch,
    const dictionary& patchFieldDict,
    const word& defaultPatchFieldType,
    const bool validBoundary
)
{
    // Remove my local data (see topoChange)
    // Clear mesh motion flux
    deleteDemandDrivenData(phiPtr_);

    // Clear the sliced fields
    clearGeomNotOldVol();

    // Clear the current volume and other geometry factors
    surfaceInterpolation::clearOut();

    // Clear any non-updateable addressing
    clearAddressing(true);


    const label sz = boundary_.size();

    polyMesh::addPatch
    (
        insertPatchi,
        patch,
        patchFieldDict,
        defaultPatchFieldType,
        validBoundary
    );

    boundary_.setSize(sz+1);
    boundary_.set
    (
        insertPatchi,
        fvPatch::New
        (
            boundaryMesh()[insertPatchi],
            boundary_
        )
    );

    objectRegistry& db = const_cast<objectRegistry&>(thisDb());
    AddPatchFields<volScalarField>
    (
        db,
        insertPatchi,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    AddPatchFields<volVectorField>
    (
        db,
        insertPatchi,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    AddPatchFields<volSphericalTensorField>
    (
        db,
        insertPatchi,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    AddPatchFields<volSymmTensorField>
    (
        db,
        insertPatchi,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    AddPatchFields<volTensorField>
    (
        db,
        insertPatchi,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );

    // Surface fields

    AddPatchFields<surfaceScalarField>
    (
        db,
        insertPatchi,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    AddPatchFields<surfaceVectorField>
    (
        db,
        insertPatchi,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    AddPatchFields<surfaceSphericalTensorField>
    (
        db,
        insertPatchi,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    AddPatchFields<surfaceSymmTensorField>
    (
        db,
        insertPatchi,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    AddPatchFields<surfaceTensorField>
    (
        db,
        insertPatchi,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
}


void Foam::fvMesh::reorderPatches
(
    const labelUList& newToOld,
    const bool validBoundary
)
{
    polyMesh::reorderPatches(newToOld, validBoundary);

    boundary_.shuffle(newToOld, validBoundary);

    objectRegistry& db = const_cast<objectRegistry&>(thisDb());
    ReorderPatchFields<volScalarField>(db, newToOld);
    ReorderPatchFields<volVectorField>(db, newToOld);
    ReorderPatchFields<volSphericalTensorField>(db, newToOld);
    ReorderPatchFields<volSymmTensorField>(db, newToOld);
    ReorderPatchFields<volTensorField>(db, newToOld);

    ReorderPatchFields<surfaceScalarField>(db, newToOld);
    ReorderPatchFields<surfaceVectorField>(db, newToOld);
    ReorderPatchFields<surfaceSphericalTensorField>(db, newToOld);
    ReorderPatchFields<surfaceSymmTensorField>(db, newToOld);
    ReorderPatchFields<surfaceTensorField>(db, newToOld);
}


bool Foam::fvMesh::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    bool ok = true;

    if (!conformal() && pointsWriteOpt() == IOobject::AUTO_WRITE)
    {
        // Create a full surface field with the polyFacesBf boundary field to
        // write to disk. Make the internal field uniform to save disk space.

        surfaceLabelField polyFaces
        (
            polyFacesBfIO(IOobject::NO_READ),
            *this,
            dimless,
            labelField(nInternalFaces(), -1),
            *polyFacesBfPtr_
        );

        ok = ok & polyFaces.write(write);
    }

    if (phiPtr_)
    {
        ok = ok && phiPtr_->write(write);
    }

    // Write V0 only if V00 exists
    if (V00Ptr_)
    {
        ok = ok && V0Ptr_->write(write);
    }

    if (topoChanger_.valid())
    {
        topoChanger_->write(write);
    }

    if (distributor_.valid())
    {
        distributor_->write(write);
    }

    if (mover_.valid())
    {
        mover_->write(write);
    }

    return ok && polyMesh::writeObject(fmt, ver, cmp, write);
}


bool Foam::fvMesh::write(const bool write) const
{
    return polyMesh::write(write);
}


template<>
typename Foam::pTraits<Foam::sphericalTensor>::labelType
Foam::fvMesh::validComponents<Foam::sphericalTensor>() const
{
    return Foam::pTraits<Foam::sphericalTensor>::labelType(1);
}


const Foam::fvSchemes& Foam::fvMesh::schemes() const
{
    if (!fvSchemes_.valid())
    {
        fvSchemes_ = new fvSchemes(*this);
    }

    return fvSchemes_;
}


const Foam::fvSolution& Foam::fvMesh::solution() const
{
    if (!fvSolution_.valid())
    {
        fvSolution_ = new fvSolution(*this);
    }

    return fvSolution_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::fvMesh::operator!=(const fvMesh& bm) const
{
    return &bm != this;
}


bool Foam::fvMesh::operator==(const fvMesh& bm) const
{
    return &bm == this;
}


// ************************************************************************* //
