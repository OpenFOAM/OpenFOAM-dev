/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
#include "slicedVolFields.H"
#include "slicedSurfaceFields.H"
#include "SubField.H"
#include "demandDrivenData.H"
#include "fvMeshLduAddressing.H"
#include "mapPolyMesh.H"
#include "MapFvFields.H"
#include "fvMeshMapper.H"
#include "mapClouds.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMesh, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMesh::clearGeomNotOldVol()
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << "clearGeomNotOldVol" << endl;
    }

    meshObject::clearUpto
    <
        fvMesh,
        GeometricMeshObject,
        MoveableMeshObject
    >(*this);

    meshObject::clearUpto
    <
        lduMesh,
        GeometricMeshObject,
        MoveableMeshObject
    >(*this);

    slicedVolScalarField::Internal* VPtr =
        static_cast<slicedVolScalarField::Internal*>(VPtr_);
    deleteDemandDrivenData(VPtr);
    VPtr_ = nullptr;

    deleteDemandDrivenData(SfPtr_);
    deleteDemandDrivenData(magSfPtr_);
    deleteDemandDrivenData(CPtr_);
    deleteDemandDrivenData(CfPtr_);
}


void Foam::fvMesh::updateGeomNotOldVol()
{
    bool haveV = (VPtr_ != nullptr);
    bool haveSf = (SfPtr_ != nullptr);
    bool haveMagSf = (magSfPtr_ != nullptr);
    bool haveCP = (CPtr_ != nullptr);
    bool haveCf = (CfPtr_ != nullptr);

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

    deleteDemandDrivenData(V0Ptr_);
    deleteDemandDrivenData(V00Ptr_);

    // Mesh motion flux cannot be deleted here because the old-time flux
    // needs to be saved.
}


void Foam::fvMesh::clearAddressing(const bool isMeshUpdate)
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << "isMeshUpdate: " << isMeshUpdate << endl;
    }

    if (isMeshUpdate)
    {
        // Part of a mesh update. Keep meshObjects that have an updateMesh
        // callback
        meshObject::clearUpto
        <
            fvMesh,
            TopologicalMeshObject,
            UpdateableMeshObject
        >
        (
            *this
        );
        meshObject::clearUpto
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
        meshObject::clear<fvMesh, TopologicalMeshObject>(*this);
        meshObject::clear<lduMesh, TopologicalMeshObject>(*this);
    }
    deleteDemandDrivenData(lduPtr_);
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
                    "V0",
                    time().timeName(),
                    *this,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
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

        curTimeIndex_ = time().timeIndex();

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


void Foam::fvMesh::clearOut()
{
    clearGeom();
    surfaceInterpolation::clearOut();

    clearAddressing();

    // Clear mesh motion flux
    deleteDemandDrivenData(phiPtr_);

    polyMesh::clearOut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMesh::fvMesh(const IOobject& io)
:
    polyMesh(io),
    surfaceInterpolation(*this),
    fvSchemes(static_cast<const objectRegistry&>(*this)),
    fvSolution(static_cast<const objectRegistry&>(*this)),
    data(static_cast<const objectRegistry&>(*this)),
    boundary_(*this, boundaryMesh()),
    lduPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    VPtr_(nullptr),
    V0Ptr_(nullptr),
    V00Ptr_(nullptr),
    SfPtr_(nullptr),
    magSfPtr_(nullptr),
    CPtr_(nullptr),
    CfPtr_(nullptr),
    phiPtr_(nullptr)
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << "Constructing fvMesh from IOobject" << endl;
    }

    // Check the existence of the cell volumes and read if present
    // and set the storage of V00
    if (fileHandler().isFile(time().timePath()/"V0"))
    {
        V0Ptr_ = new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "V0",
                time().timeName(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
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
                time().timeName(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                true
            ),
            *this
        );
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
        move(points),
        shapes,
        boundaryFaces,
        boundaryPatchNames,
        boundaryDicts,
        defaultBoundaryPatchName,
        defaultBoundaryPatchType,
        syncPar
    ),
    surfaceInterpolation(*this),
    fvSchemes(static_cast<const objectRegistry&>(*this)),
    fvSolution(static_cast<const objectRegistry&>(*this)),
    data(static_cast<const objectRegistry&>(*this)),
    boundary_(*this, boundaryMesh()),
    lduPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    VPtr_(nullptr),
    V0Ptr_(nullptr),
    V00Ptr_(nullptr),
    SfPtr_(nullptr),
    magSfPtr_(nullptr),
    CPtr_(nullptr),
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
        move(points),
        move(faces),
        move(allOwner),
        move(allNeighbour),
        syncPar
    ),
    surfaceInterpolation(*this),
    fvSchemes(static_cast<const objectRegistry&>(*this)),
    fvSolution(static_cast<const objectRegistry&>(*this)),
    data(static_cast<const objectRegistry&>(*this)),
    boundary_(*this, boundaryMesh()),
    lduPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    VPtr_(nullptr),
    V0Ptr_(nullptr),
    V00Ptr_(nullptr),
    SfPtr_(nullptr),
    magSfPtr_(nullptr),
    CPtr_(nullptr),
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
    polyMesh(io, move(points), move(faces), move(cells), syncPar),
    surfaceInterpolation(*this),
    fvSchemes(static_cast<const objectRegistry&>(*this)),
    fvSolution(static_cast<const objectRegistry&>(*this)),
    data(static_cast<const objectRegistry&>(*this)),
    boundary_(*this),
    lduPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    VPtr_(nullptr),
    V0Ptr_(nullptr),
    V00Ptr_(nullptr),
    SfPtr_(nullptr),
    magSfPtr_(nullptr),
    CPtr_(nullptr),
    CfPtr_(nullptr),
    phiPtr_(nullptr)
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << "Constructing fvMesh from components" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMesh::~fvMesh()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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


Foam::polyMesh::readUpdateState Foam::fvMesh::readUpdate()
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << "Updating fvMesh.  ";
    }

    polyMesh::readUpdateState state = polyMesh::readUpdate();

    if (state == polyMesh::TOPO_PATCH_CHANGE)
    {
        if (debug)
        {
            Info<< "Boundary and topological update" << endl;
        }

        boundary_.readUpdate(boundaryMesh());

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


void Foam::fvMesh::mapFields(const mapPolyMesh& meshMap)
{
    if (debug)
    {
        Pout<< FUNCTION_NAME
            << " nOldCells:" << meshMap.nOldCells()
            << " nCells:" << nCells()
            << " nOldFaces:" << meshMap.nOldFaces()
            << " nFaces:" << nFaces()
            << endl;
    }


    // We require geometric properties valid for the old mesh
    if
    (
        meshMap.cellMap().size() != nCells()
     || meshMap.faceMap().size() != nFaces()
    )
    {
        FatalErrorInFunction
            << "mapPolyMesh does not correspond to the old mesh."
            << " nCells:" << nCells()
            << " cellMap:" << meshMap.cellMap().size()
            << " nOldCells:" << meshMap.nOldCells()
            << " nFaces:" << nFaces()
            << " faceMap:" << meshMap.faceMap().size()
            << " nOldFaces:" << meshMap.nOldFaces()
            << exit(FatalError);
    }

    // Create a mapper
    const fvMeshMapper mapper(*this, meshMap);

    // Map all the volFields in the objectRegistry
    MapGeometricFields<scalar, fvPatchField, fvMeshMapper, volMesh>
    (mapper);
    MapGeometricFields<vector, fvPatchField, fvMeshMapper, volMesh>
    (mapper);
    MapGeometricFields<sphericalTensor, fvPatchField, fvMeshMapper, volMesh>
    (mapper);
    MapGeometricFields<symmTensor, fvPatchField, fvMeshMapper, volMesh>
    (mapper);
    MapGeometricFields<tensor, fvPatchField, fvMeshMapper, volMesh>
    (mapper);

    // Map all the surfaceFields in the objectRegistry
    MapGeometricFields<scalar, fvsPatchField, fvMeshMapper, surfaceMesh>
    (mapper);
    MapGeometricFields<vector, fvsPatchField, fvMeshMapper, surfaceMesh>
    (mapper);
    MapGeometricFields<symmTensor, fvsPatchField, fvMeshMapper, surfaceMesh>
    (mapper);
    MapGeometricFields<symmTensor, fvsPatchField, fvMeshMapper, surfaceMesh>
    (mapper);
    MapGeometricFields<tensor, fvsPatchField, fvMeshMapper, surfaceMesh>
    (mapper);

    // Map all the dimensionedFields in the objectRegistry
    MapDimensionedFields<scalar, fvMeshMapper, volMesh>(mapper);
    MapDimensionedFields<vector, fvMeshMapper, volMesh>(mapper);
    MapDimensionedFields<sphericalTensor, fvMeshMapper, volMesh>(mapper);
    MapDimensionedFields<symmTensor, fvMeshMapper, volMesh>(mapper);
    MapDimensionedFields<tensor, fvMeshMapper, volMesh>(mapper);

    // Map all the clouds in the objectRegistry
    mapClouds(*this, meshMap);


    const labelList& cellMap = meshMap.cellMap();

    // Map the old volume. Just map to new cell labels.
    if (V0Ptr_)
    {
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
        forAll(meshMap.reverseCellMap(), oldCelli)
        {
            label index = meshMap.reverseCellMap()[oldCelli];

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


    // Map the old-old volume. Just map to new cell labels.
    if (V00Ptr_)
    {
        scalarField& V00 = *V00Ptr_;

        scalarField savedV00(V00);
        V00.setSize(nCells());

        forAll(V00, i)
        {
            if (cellMap[i] > -1)
            {
                V00[i] = savedV00[cellMap[i]];
            }
            else
            {
                V00[i] = 0.0;
            }
        }

        // Inject volume of merged cells
        label nMerged = 0;
        forAll(meshMap.reverseCellMap(), oldCelli)
        {
            label index = meshMap.reverseCellMap()[oldCelli];

            if (index < -1)
            {
                label celli = -index-2;

                V00[celli] += savedV00[oldCelli];
                nMerged++;
            }
        }

        if (debug)
        {
            Info<< "Mapping old time volume V00. Merged "
                << nMerged << " out of " << nCells() << " cells" << endl;
        }
    }
}


Foam::tmp<Foam::scalarField> Foam::fvMesh::movePoints(const pointField& p)
{
    // Grab old time volumes if the time has been incremented
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
                this->time().timeName(),
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
        // Grab old time mesh motion fluxes if the time has been incremented
        if (phiPtr_->timeIndex() != time().timeIndex())
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

    meshObject::movePoints<fvMesh>(*this);
    meshObject::movePoints<lduMesh>(*this);

    return tsweptVols;
}


void Foam::fvMesh::updateMesh(const mapPolyMesh& mpm)
{
    // Update polyMesh. This needs to keep volume existent!
    polyMesh::updateMesh(mpm);

    if (VPtr_)
    {
        // Cache old time volumes if they exist and the time has been
        // incremented.  This will update V0, V00
        if (V0Ptr_)
        {
            storeOldVol(mpm.oldCellVolumes());
        }

        // Few checks
        if (VPtr_ && (V().size() != mpm.nOldCells()))
        {
            FatalErrorInFunction
                << "V:" << V().size()
                << " not equal to the number of old cells "
                << mpm.nOldCells()
                << exit(FatalError);
        }
        if (V0Ptr_ && (V0Ptr_->size() != mpm.nOldCells()))
        {
            FatalErrorInFunction
                << "V0:" << V0Ptr_->size()
                << " not equal to the number of old cells "
                << mpm.nOldCells()
                << exit(FatalError);
        }
        if (V00Ptr_ && (V00Ptr_->size() != mpm.nOldCells()))
        {
            FatalErrorInFunction
                << "V0:" << V00Ptr_->size()
                << " not equal to the number of old cells "
                << mpm.nOldCells()
                << exit(FatalError);
        }
    }

    // Clear the sliced fields
    clearGeomNotOldVol();

    // Map all fields
    mapFields(mpm);

    // Clear the current volume and other geometry factors
    surfaceInterpolation::clearOut();

    // Clear any non-updateable addressing
    clearAddressing(true);

    meshObject::updateMesh<fvMesh>(*this, mpm);
    meshObject::updateMesh<lduMesh>(*this, mpm);
}


void Foam::fvMesh::updateMesh(const mapDistributePolyMesh& mdpm)
{
    // Update polyMesh. This needs to keep volume existent!
    // polyMesh::updateMesh(mdpm);

    // if (VPtr_)
    // {
    //     // Cache old time volumes if they exist and the time has been
    //     // incremented.  This will update V0, V00
    //     if (V0Ptr_)
    //     {
    //         storeOldVol(mpm.oldCellVolumes());
    //     }
    //
    //     // Few checks
    //     if (VPtr_ && (V().size() != mdpm.nOldCells()))
    //     {
    //         FatalErrorInFunction
    //             << "V:" << V().size()
    //             << " not equal to the number of old cells "
    //             << mdpm.nOldCells()
    //             << exit(FatalError);
    //     }
    //     if (V0Ptr_ && (V0Ptr_->size() != mdpm.nOldCells()))
    //     {
    //         FatalErrorInFunction
    //             << "V0:" << V0Ptr_->size()
    //             << " not equal to the number of old cells "
    //             << mdpm.nOldCells()
    //             << exit(FatalError);
    //     }
    //     if (V00Ptr_ && (V00Ptr_->size() != mdpm.nOldCells()))
    //     {
    //         FatalErrorInFunction
    //             << "V0:" << V00Ptr_->size()
    //             << " not equal to the number of old cells "
    //             << mdpm.nOldCells()
    //             << exit(FatalError);
    //     }
    // }

    // Clear the sliced fields
    // clearGeomNotOldVol();
    clearGeom();

    // Map all fields
    // mapFields(mdpm);

    // Clear the current volume and other geometry factors
    surfaceInterpolation::clearOut();

    // Clear any non-updateable addressing
    clearAddressing(true);

    // meshObject::updateMesh<fvMesh>(*this, mdpm);
    // meshObject::updateMesh<lduMesh>(*this, mdpm);
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
    // Remove my local data (see updateMesh)
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
    if (phiPtr_)
    {
        ok = phiPtr_->write(write);
    }

    // Write V0 only if V00 exists
    if (V00Ptr_)
    {
        ok = ok && V0Ptr_->write(write);
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
