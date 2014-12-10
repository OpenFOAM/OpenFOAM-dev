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

#include "polyMesh.H"
#include "Time.H"
#include "cellIOList.H"
#include "wedgePolyPatch.H"
#include "emptyPolyPatch.H"
#include "globalMeshData.H"
#include "processorPolyPatch.H"
#include "polyMeshTetDecomposition.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "MeshObject.H"
#include "pointMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyMesh, 0);

    word polyMesh::defaultRegion = "region0";
    word polyMesh::meshSubDir = "polyMesh";
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::polyMesh::calcDirections() const
{
    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        solutionD_[cmpt] = 1;
    }

    // Knock out empty and wedge directions. Note:they will be present on all
    // domains.

    label nEmptyPatches = 0;
    label nWedgePatches = 0;

    vector emptyDirVec = vector::zero;
    vector wedgeDirVec = vector::zero;

    forAll(boundaryMesh(), patchi)
    {
        if (boundaryMesh()[patchi].size())
        {
            if (isA<emptyPolyPatch>(boundaryMesh()[patchi]))
            {
                nEmptyPatches++;
                emptyDirVec += sum(cmptMag(boundaryMesh()[patchi].faceAreas()));
            }
            else if (isA<wedgePolyPatch>(boundaryMesh()[patchi]))
            {
                const wedgePolyPatch& wpp = refCast<const wedgePolyPatch>
                (
                    boundaryMesh()[patchi]
                );

                nWedgePatches++;
                wedgeDirVec += cmptMag(wpp.centreNormal());
            }
        }
    }

    reduce(nEmptyPatches, maxOp<label>());
    reduce(nWedgePatches, maxOp<label>());

    if (nEmptyPatches)
    {
        reduce(emptyDirVec, sumOp<vector>());

        emptyDirVec /= mag(emptyDirVec);

        for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
        {
            if (emptyDirVec[cmpt] > 1e-6)
            {
                solutionD_[cmpt] = -1;
            }
            else
            {
                solutionD_[cmpt] = 1;
            }
        }
    }


    // Knock out wedge directions

    geometricD_ = solutionD_;

    if (nWedgePatches)
    {
        reduce(wedgeDirVec, sumOp<vector>());

        wedgeDirVec /= mag(wedgeDirVec);

        for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
        {
            if (wedgeDirVec[cmpt] > 1e-6)
            {
                geometricD_[cmpt] = -1;
            }
            else
            {
                geometricD_[cmpt] = 1;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyMesh::polyMesh(const IOobject& io)
:
    objectRegistry(io),
    primitiveMesh(),
    points_
    (
        IOobject
        (
            "points",
            time().findInstance(meshDir(), "points"),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    faces_
    (
        IOobject
        (
            "faces",
            time().findInstance(meshDir(), "faces"),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    owner_
    (
        IOobject
        (
            "owner",
            faces_.instance(),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    neighbour_
    (
        IOobject
        (
            "neighbour",
            faces_.instance(),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    clearedPrimitives_(false),
    boundary_
    (
        IOobject
        (
            "boundary",
            time().findInstance(meshDir(), "boundary"),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        *this
    ),
    bounds_(points_),
    comm_(UPstream::worldComm),
    geometricD_(Vector<label>::zero),
    solutionD_(Vector<label>::zero),
    tetBasePtIsPtr_(NULL),
    cellTreePtr_(NULL),
    pointZones_
    (
        IOobject
        (
            "pointZones",
            time().findInstance
            (
                meshDir(),
                "pointZones",
                IOobject::READ_IF_PRESENT
            ),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *this
    ),
    faceZones_
    (
        IOobject
        (
            "faceZones",
            time().findInstance
            (
                meshDir(),
                "faceZones",
                IOobject::READ_IF_PRESENT
            ),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *this
    ),
    cellZones_
    (
        IOobject
        (
            "cellZones",
            time().findInstance
            (
                meshDir(),
                "cellZones",
                IOobject::READ_IF_PRESENT
            ),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *this
    ),
    globalMeshDataPtr_(NULL),
    moving_(false),
    topoChanging_(false),
    curMotionTimeIndex_(time().timeIndex()),
    oldPointsPtr_(NULL)
{
    if (exists(owner_.objectPath()))
    {
        initMesh();
    }
    else
    {
        cellCompactIOList cLst
        (
            IOobject
            (
                "cells",
                time().findInstance(meshDir(), "cells"),
                meshSubDir,
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        // Set the primitive mesh
        initMesh(cLst);

        owner_.write();
        neighbour_.write();
    }

    // Calculate topology for the patches (processor-processor comms etc.)
    boundary_.updateMesh();

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();

    // Warn if global empty mesh
    if (returnReduce(nPoints(), sumOp<label>()) == 0)
    {
        WarningIn("polyMesh(const IOobject&)")
            << "no points in mesh" << endl;
    }
    if (returnReduce(nCells(), sumOp<label>()) == 0)
    {
        WarningIn("polyMesh(const IOobject&)")
            << "no cells in mesh" << endl;
    }

    // Initialise demand-driven data
    calcDirections();
}


Foam::polyMesh::polyMesh
(
    const IOobject& io,
    const Xfer<pointField>& points,
    const Xfer<faceList>& faces,
    const Xfer<labelList>& owner,
    const Xfer<labelList>& neighbour,
    const bool syncPar
)
:
    objectRegistry(io),
    primitiveMesh(),
    points_
    (
        IOobject
        (
            "points",
            instance(),
            meshSubDir,
            *this,
            io.readOpt(),
            IOobject::AUTO_WRITE
        ),
        points
    ),
    faces_
    (
        IOobject
        (
            "faces",
            instance(),
            meshSubDir,
            *this,
            io.readOpt(),
            IOobject::AUTO_WRITE
        ),
        faces
    ),
    owner_
    (
        IOobject
        (
            "owner",
            instance(),
            meshSubDir,
            *this,
            io.readOpt(),
            IOobject::AUTO_WRITE
        ),
        owner
    ),
    neighbour_
    (
        IOobject
        (
            "neighbour",
            instance(),
            meshSubDir,
            *this,
            io.readOpt(),
            IOobject::AUTO_WRITE
        ),
        neighbour
    ),
    clearedPrimitives_(false),
    boundary_
    (
        IOobject
        (
            "boundary",
            instance(),
            meshSubDir,
            *this,
            io.readOpt(),
            IOobject::AUTO_WRITE
        ),
        *this,
        polyPatchList()
    ),
    bounds_(points_, syncPar),
    comm_(UPstream::worldComm),
    geometricD_(Vector<label>::zero),
    solutionD_(Vector<label>::zero),
    tetBasePtIsPtr_(NULL),
    cellTreePtr_(NULL),
    pointZones_
    (
        IOobject
        (
            "pointZones",
            instance(),
            meshSubDir,
            *this,
            io.readOpt(),
            IOobject::NO_WRITE
        ),
        *this,
        PtrList<pointZone>()
    ),
    faceZones_
    (
        IOobject
        (
            "faceZones",
            instance(),
            meshSubDir,
            *this,
            io.readOpt(),
            IOobject::NO_WRITE
        ),
        *this,
        PtrList<faceZone>()
    ),
    cellZones_
    (
        IOobject
        (
            "cellZones",
            instance(),
            meshSubDir,
            *this,
            io.readOpt(),
            IOobject::NO_WRITE
        ),
        *this,
        PtrList<cellZone>()
    ),
    globalMeshDataPtr_(NULL),
    moving_(false),
    topoChanging_(false),
    curMotionTimeIndex_(time().timeIndex()),
    oldPointsPtr_(NULL)
{
    // Check if the faces and cells are valid
    forAll(faces_, faceI)
    {
        const face& curFace = faces_[faceI];

        if (min(curFace) < 0 || max(curFace) > points_.size())
        {
            FatalErrorIn
            (
                "polyMesh::polyMesh\n"
                "(\n"
                "    const IOobject& io,\n"
                "    const pointField& points,\n"
                "    const faceList& faces,\n"
                "    const cellList& cells\n"
                ")\n"
            )   << "Face " << faceI << "contains vertex labels out of range: "
                << curFace << " Max point index = " << points_.size()
                << abort(FatalError);
        }
    }

    // Set the primitive mesh
    initMesh();
}


Foam::polyMesh::polyMesh
(
    const IOobject& io,
    const Xfer<pointField>& points,
    const Xfer<faceList>& faces,
    const Xfer<cellList>& cells,
    const bool syncPar
)
:
    objectRegistry(io),
    primitiveMesh(),
    points_
    (
        IOobject
        (
            "points",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        points
    ),
    faces_
    (
        IOobject
        (
            "faces",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        faces
    ),
    owner_
    (
        IOobject
        (
            "owner",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0
    ),
    neighbour_
    (
        IOobject
        (
            "neighbour",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0
    ),
    clearedPrimitives_(false),
    boundary_
    (
        IOobject
        (
            "boundary",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *this,
        0
    ),
    bounds_(points_, syncPar),
    comm_(UPstream::worldComm),
    geometricD_(Vector<label>::zero),
    solutionD_(Vector<label>::zero),
    tetBasePtIsPtr_(NULL),
    cellTreePtr_(NULL),
    pointZones_
    (
        IOobject
        (
            "pointZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    faceZones_
    (
        IOobject
        (
            "faceZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    cellZones_
    (
        IOobject
        (
            "cellZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    globalMeshDataPtr_(NULL),
    moving_(false),
    topoChanging_(false),
    curMotionTimeIndex_(time().timeIndex()),
    oldPointsPtr_(NULL)
{
    // Check if faces are valid
    forAll(faces_, faceI)
    {
        const face& curFace = faces_[faceI];

        if (min(curFace) < 0 || max(curFace) > points_.size())
        {
            FatalErrorIn
            (
                "polyMesh::polyMesh\n"
                "(\n"
                "    const IOobject&,\n"
                "    const Xfer<pointField>&,\n"
                "    const Xfer<faceList>&,\n"
                "    const Xfer<cellList>&\n"
                ")\n"
            )   << "Face " << faceI << "contains vertex labels out of range: "
                << curFace << " Max point index = " << points_.size()
                << abort(FatalError);
        }
    }

    // transfer in cell list
    cellList cLst(cells);

    // Check if cells are valid
    forAll(cLst, cellI)
    {
        const cell& curCell = cLst[cellI];

        if (min(curCell) < 0 || max(curCell) > faces_.size())
        {
            FatalErrorIn
            (
                "polyMesh::polyMesh\n"
                "(\n"
                "    const IOobject&,\n"
                "    const Xfer<pointField>&,\n"
                "    const Xfer<faceList>&,\n"
                "    const Xfer<cellList>&\n"
                ")\n"
            )   << "Cell " << cellI << "contains face labels out of range: "
                << curCell << " Max face index = " << faces_.size()
                << abort(FatalError);
        }
    }

    // Set the primitive mesh
    initMesh(cLst);
}


void Foam::polyMesh::resetPrimitives
(
    const Xfer<pointField>& points,
    const Xfer<faceList>& faces,
    const Xfer<labelList>& owner,
    const Xfer<labelList>& neighbour,
    const labelList& patchSizes,
    const labelList& patchStarts,
    const bool validBoundary
)
{
    // Clear addressing. Keep geometric props and updateable props for mapping.
    clearAddressing(true);

    // Take over new primitive data.
    // Optimized to avoid overwriting data at all
    if (&points)
    {
        points_.transfer(points());
        bounds_ = boundBox(points_, validBoundary);
    }

    if (&faces)
    {
        faces_.transfer(faces());
    }

    if (&owner)
    {
        owner_.transfer(owner());
    }

    if (&neighbour)
    {
        neighbour_.transfer(neighbour());
    }


    // Reset patch sizes and starts
    forAll(boundary_, patchI)
    {
        boundary_[patchI] = polyPatch
        (
            boundary_[patchI],
            boundary_,
            patchI,
            patchSizes[patchI],
            patchStarts[patchI]
        );
    }


    // Flags the mesh files as being changed
    setInstance(time().timeName());

    // Check if the faces and cells are valid
    forAll(faces_, faceI)
    {
        const face& curFace = faces_[faceI];

        if (min(curFace) < 0 || max(curFace) > points_.size())
        {
            FatalErrorIn
            (
                "polyMesh::polyMesh::resetPrimitives\n"
                "(\n"
                "    const Xfer<pointField>&,\n"
                "    const Xfer<faceList>&,\n"
                "    const Xfer<labelList>& owner,\n"
                "    const Xfer<labelList>& neighbour,\n"
                "    const labelList& patchSizes,\n"
                "    const labelList& patchStarts\n"
                "    const bool validBoundary\n"
                ")\n"
            )   << "Face " << faceI << " contains vertex labels out of range: "
                << curFace << " Max point index = " << points_.size()
                << abort(FatalError);
        }
    }


    // Set the primitive mesh from the owner_, neighbour_.
    // Works out from patch end where the active faces stop.
    initMesh();


    if (validBoundary)
    {
        // Note that we assume that all the patches stay the same and are
        // correct etc. so we can already use the patches to do
        // processor-processor comms.

        // Calculate topology for the patches (processor-processor comms etc.)
        boundary_.updateMesh();

        // Calculate the geometry for the patches (transformation tensors etc.)
        boundary_.calcGeometry();

        // Warn if global empty mesh
        if
        (
            (returnReduce(nPoints(), sumOp<label>()) == 0)
         || (returnReduce(nCells(), sumOp<label>()) == 0)
        )
        {
            FatalErrorIn
            (
                "polyMesh::polyMesh::resetPrimitives\n"
                "(\n"
                "    const Xfer<pointField>&,\n"
                "    const Xfer<faceList>&,\n"
                "    const Xfer<labelList>& owner,\n"
                "    const Xfer<labelList>& neighbour,\n"
                "    const labelList& patchSizes,\n"
                "    const labelList& patchStarts\n"
                "    const bool validBoundary\n"
                ")\n"
            )   << "no points or no cells in mesh" << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polyMesh::~polyMesh()
{
    clearOut();
    resetMotion();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fileName& Foam::polyMesh::dbDir() const
{
    if (objectRegistry::dbDir() == defaultRegion)
    {
        return parent().dbDir();
    }
    else
    {
        return objectRegistry::dbDir();
    }
}


Foam::fileName Foam::polyMesh::meshDir() const
{
    return dbDir()/meshSubDir;
}


const Foam::fileName& Foam::polyMesh::pointsInstance() const
{
    return points_.instance();
}


const Foam::fileName& Foam::polyMesh::facesInstance() const
{
    return faces_.instance();
}


const Foam::Vector<Foam::label>& Foam::polyMesh::geometricD() const
{
    if (geometricD_.x() == 0)
    {
        calcDirections();
    }

    return geometricD_;
}


Foam::label Foam::polyMesh::nGeometricD() const
{
    return cmptSum(geometricD() + Vector<label>::one)/2;
}


const Foam::Vector<Foam::label>& Foam::polyMesh::solutionD() const
{
    if (solutionD_.x() == 0)
    {
        calcDirections();
    }

    return solutionD_;
}


Foam::label Foam::polyMesh::nSolutionD() const
{
    return cmptSum(solutionD() + Vector<label>::one)/2;
}


const Foam::labelList& Foam::polyMesh::tetBasePtIs() const
{
    if (tetBasePtIsPtr_.empty())
    {
        if (debug)
        {
            WarningIn("const labelList& polyMesh::tetBasePtIs() const")
                << "Tet base point indices not available.  "
                << "Forcing storage of base points."
                << endl;
        }

        tetBasePtIsPtr_.reset
        (
            new labelList
            (
                polyMeshTetDecomposition::findFaceBasePts(*this)
            )
        );
    }

    return tetBasePtIsPtr_();
}


const Foam::indexedOctree<Foam::treeDataCell>&
Foam::polyMesh::cellTree() const
{
    if (cellTreePtr_.empty())
    {
        treeBoundBox overallBb(points());

        Random rndGen(261782);

        overallBb = overallBb.extend(rndGen, 1e-4);
        overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        cellTreePtr_.reset
        (
            new indexedOctree<treeDataCell>
            (
                treeDataCell
                (
                    false,          // not cache bb
                    *this,
                    FACEDIAGTETS    // use tetDecomposition for any inside test
                ),
                overallBb,
                8,              // maxLevel
                10,             // leafsize
                5.0             // duplicity
            )
        );
    }

    return cellTreePtr_();
}


void Foam::polyMesh::addPatches
(
    const List<polyPatch*>& p,
    const bool validBoundary
)
{
    if (boundaryMesh().size())
    {
        FatalErrorIn
        (
            "void polyMesh::addPatches(const List<polyPatch*>&, const bool)"
        )   << "boundary already exists"
            << abort(FatalError);
    }

    // Reset valid directions
    geometricD_ = Vector<label>::zero;
    solutionD_ = Vector<label>::zero;

    boundary_.setSize(p.size());

    // Copy the patch pointers
    forAll(p, pI)
    {
        boundary_.set(pI, p[pI]);
    }

    // parallelData depends on the processorPatch ordering so force
    // recalculation. Problem: should really be done in removeBoundary but
    // there is some info in parallelData which might be interesting inbetween
    // removeBoundary and addPatches.
    globalMeshDataPtr_.clear();

    if (validBoundary)
    {
        // Calculate topology for the patches (processor-processor comms etc.)
        boundary_.updateMesh();

        // Calculate the geometry for the patches (transformation tensors etc.)
        boundary_.calcGeometry();

        boundary_.checkDefinition();
    }
}


void Foam::polyMesh::addZones
(
    const List<pointZone*>& pz,
    const List<faceZone*>& fz,
    const List<cellZone*>& cz
)
{
    if (pointZones().size() || faceZones().size() || cellZones().size())
    {
        FatalErrorIn
        (
            "void addZones\n"
            "(\n"
            "    const List<pointZone*>&,\n"
            "    const List<faceZone*>&,\n"
            "    const List<cellZone*>&\n"
            ")"
        )   << "point, face or cell zone already exists"
            << abort(FatalError);
    }

    // Point zones
    if (pz.size())
    {
        pointZones_.setSize(pz.size());

        // Copy the zone pointers
        forAll(pz, pI)
        {
            pointZones_.set(pI, pz[pI]);
        }

        pointZones_.writeOpt() = IOobject::AUTO_WRITE;
    }

    // Face zones
    if (fz.size())
    {
        faceZones_.setSize(fz.size());

        // Copy the zone pointers
        forAll(fz, fI)
        {
            faceZones_.set(fI, fz[fI]);
        }

        faceZones_.writeOpt() = IOobject::AUTO_WRITE;
    }

    // Cell zones
    if (cz.size())
    {
        cellZones_.setSize(cz.size());

        // Copy the zone pointers
        forAll(cz, cI)
        {
            cellZones_.set(cI, cz[cI]);
        }

        cellZones_.writeOpt() = IOobject::AUTO_WRITE;
    }
}


const Foam::pointField& Foam::polyMesh::points() const
{
    if (clearedPrimitives_)
    {
        FatalErrorIn("const pointField& polyMesh::points() const")
            << "points deallocated"
            << abort(FatalError);
    }

    return points_;
}


bool Foam::polyMesh::upToDatePoints(const regIOobject& io) const
{
    return io.upToDate(points_);
}


void Foam::polyMesh::setUpToDatePoints(regIOobject& io) const
{
    io.eventNo() = points_.eventNo();
}


const Foam::faceList& Foam::polyMesh::faces() const
{
    if (clearedPrimitives_)
    {
        FatalErrorIn("const faceList& polyMesh::faces() const")
            << "faces deallocated"
            << abort(FatalError);
    }

    return faces_;
}


const Foam::labelList& Foam::polyMesh::faceOwner() const
{
    return owner_;
}


const Foam::labelList& Foam::polyMesh::faceNeighbour() const
{
    return neighbour_;
}


const Foam::pointField& Foam::polyMesh::oldPoints() const
{
    if (oldPointsPtr_.empty())
    {
        if (debug)
        {
            WarningIn("const pointField& polyMesh::oldPoints() const")
                << "Old points not available.  Forcing storage of old points"
                << endl;
        }

        oldPointsPtr_.reset(new pointField(points_));
        curMotionTimeIndex_ = time().timeIndex();
    }

    return oldPointsPtr_();
}


Foam::tmp<Foam::scalarField> Foam::polyMesh::movePoints
(
    const pointField& newPoints
)
{
    if (debug)
    {
        Info<< "tmp<scalarField> polyMesh::movePoints(const pointField&) : "
            << " Moving points for time " << time().value()
            << " index " << time().timeIndex() << endl;
    }

    moving(true);

    // Pick up old points
    if (curMotionTimeIndex_ != time().timeIndex())
    {
        // Mesh motion in the new time step
        oldPointsPtr_.clear();
        oldPointsPtr_.reset(new pointField(points_));
        curMotionTimeIndex_ = time().timeIndex();
    }

    points_ = newPoints;

    bool moveError = false;
    if (debug)
    {
        // Check mesh motion
        if (checkMeshMotion(points_, true))
        {
            moveError = true;

            Info<< "tmp<scalarField> polyMesh::movePoints"
                << "(const pointField&) : "
                << "Moving the mesh with given points will "
                << "invalidate the mesh." << nl
                << "Mesh motion should not be executed." << endl;
        }
    }

    points_.writeOpt() = IOobject::AUTO_WRITE;
    points_.instance() = time().timeName();
    points_.eventNo() = getEvent();

    tmp<scalarField> sweptVols = primitiveMesh::movePoints
    (
        points_,
        oldPoints()
    );

    // Adjust parallel shared points
    if (globalMeshDataPtr_.valid())
    {
        globalMeshDataPtr_().movePoints(points_);
    }

    // Force recalculation of all geometric data with new points

    bounds_ = boundBox(points_);
    boundary_.movePoints(points_);

    pointZones_.movePoints(points_);
    faceZones_.movePoints(points_);
    cellZones_.movePoints(points_);

    // Reset valid directions (could change with rotation)
    geometricD_ = Vector<label>::zero;
    solutionD_ = Vector<label>::zero;

    meshObject::movePoints<polyMesh>(*this);
    meshObject::movePoints<pointMesh>(*this);

    const_cast<Time&>(time()).functionObjects().movePoints(*this);


    if (debug && moveError)
    {
        // Write mesh to ease debugging. Note we want to avoid calling
        // e.g. fvMesh::write since meshPhi not yet complete.
        polyMesh::write();
    }

    return sweptVols;
}


// Reset motion by deleting old points
void Foam::polyMesh::resetMotion() const
{
    curMotionTimeIndex_ = 0;
    oldPointsPtr_.clear();
}


// Return parallel info
const Foam::globalMeshData& Foam::polyMesh::globalData() const
{
    if (globalMeshDataPtr_.empty())
    {
        if (debug)
        {
            Pout<< "polyMesh::globalData() const : "
                << "Constructing parallelData from processor topology"
                << endl;
        }
        // Construct globalMeshData using processorPatch information only.
        globalMeshDataPtr_.reset(new globalMeshData(*this));
    }

    return globalMeshDataPtr_();
}


Foam::label Foam::polyMesh::comm() const
{
    return comm_;
}


Foam::label& Foam::polyMesh::comm()
{
    return comm_;
}


void Foam::polyMesh::removeFiles(const fileName& instanceDir) const
{
    fileName meshFilesPath = thisDb().time().path()/instanceDir/meshDir();

    rm(meshFilesPath/"points");
    rm(meshFilesPath/"faces");
    rm(meshFilesPath/"owner");
    rm(meshFilesPath/"neighbour");
    rm(meshFilesPath/"cells");
    rm(meshFilesPath/"boundary");
    rm(meshFilesPath/"pointZones");
    rm(meshFilesPath/"faceZones");
    rm(meshFilesPath/"cellZones");
    rm(meshFilesPath/"meshModifiers");
    rm(meshFilesPath/"parallelData");

    // remove subdirectories
    if (isDir(meshFilesPath/"sets"))
    {
        rmDir(meshFilesPath/"sets");
    }
}


void Foam::polyMesh::removeFiles() const
{
    removeFiles(instance());
}


void Foam::polyMesh::findCellFacePt
(
    const point& pt,
    label& cellI,
    label& tetFaceI,
    label& tetPtI
) const
{
    cellI = -1;
    tetFaceI = -1;
    tetPtI = -1;

    const indexedOctree<treeDataCell>& tree = cellTree();

    // Find nearest cell to the point

    pointIndexHit info = tree.findNearest(pt, sqr(GREAT));

    if (info.hit())
    {
        label nearestCellI = tree.shapes().cellLabels()[info.index()];

        // Check the nearest cell to see if the point is inside.
        findTetFacePt(nearestCellI, pt, tetFaceI, tetPtI);

        if (tetFaceI != -1)
        {
            // Point was in the nearest cell

            cellI = nearestCellI;

            return;
        }
        else
        {
            // Check the other possible cells that the point may be in

            labelList testCells = tree.findIndices(pt);

            forAll(testCells, pCI)
            {
                label testCellI = tree.shapes().cellLabels()[testCells[pCI]];

                if (testCellI == nearestCellI)
                {
                    // Don't retest the nearest cell

                    continue;
                }

                // Check the test cell to see if the point is inside.
                findTetFacePt(testCellI, pt, tetFaceI, tetPtI);

                if (tetFaceI != -1)
                {
                    // Point was in the test cell

                    cellI = testCellI;

                    return;
                }
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::polyMesh::findCellFacePt"
            "("
                "const point&, "
                "label&, "
                "label&, "
                "label&"
            ") const"
        )   << "Did not find nearest cell in search tree."
            << abort(FatalError);
    }
}


void Foam::polyMesh::findTetFacePt
(
    const label cellI,
    const point& pt,
    label& tetFaceI,
    label& tetPtI
) const
{
    const polyMesh& mesh = *this;

    tetIndices tet(polyMeshTetDecomposition::findTet(mesh, cellI, pt));
    tetFaceI = tet.face();
    tetPtI = tet.tetPt();
}


bool Foam::polyMesh::pointInCell
(
    const point& p,
    label cellI,
    const cellRepresentation decompMode
) const
{
    switch (decompMode)
    {
        case FACEPLANES:
        {
            return primitiveMesh::pointInCell(p, cellI);
        }
        break;

        case FACECENTRETETS:
        {
            // only test that point is on inside of plane defined by cell face
            // triangles
            const cell& cFaces = cells()[cellI];

            forAll(cFaces, cFaceI)
            {
                label faceI = cFaces[cFaceI];
                const face& f = faces_[faceI];
                const point& fc = faceCentres()[faceI];
                bool isOwn = (owner_[faceI] == cellI);

                forAll(f, fp)
                {
                    label pointI;
                    label nextPointI;

                    if (isOwn)
                    {
                        pointI = f[fp];
                        nextPointI = f.nextLabel(fp);
                    }
                    else
                    {
                        pointI = f.nextLabel(fp);
                        nextPointI = f[fp];
                    }

                    triPointRef faceTri
                    (
                        points()[pointI],
                        points()[nextPointI],
                        fc
                    );

                    vector proj = p - faceTri.centre();

                    if ((faceTri.normal() & proj) > 0)
                    {
                        return false;
                    }
                }
            }
            return true;
        }
        break;

        case FACEDIAGTETS:
        {
            // only test that point is on inside of plane defined by cell face
            // triangles
            const cell& cFaces = cells()[cellI];

            forAll(cFaces, cFaceI)
            {
                label faceI = cFaces[cFaceI];
                const face& f = faces_[faceI];

                for (label tetPtI = 1; tetPtI < f.size() - 1; tetPtI++)
                {
                    // Get tetIndices of face triangle
                    tetIndices faceTetIs
                    (
                        polyMeshTetDecomposition::triangleTetIndices
                        (
                            *this,
                            faceI,
                            cellI,
                            tetPtI
                        )
                    );

                    triPointRef faceTri = faceTetIs.faceTri(*this);

                    vector proj = p - faceTri.centre();

                    if ((faceTri.normal() & proj) > 0)
                    {
                        return false;
                    }
                }
            }

            return true;
        }
        break;
    }
    return false;
}


Foam::label Foam::polyMesh::findCell
(
    const point& location,
    const cellRepresentation decompMode
) const
{
    if (Pstream::parRun() && decompMode == FACEDIAGTETS)
    {
        // Force construction of face-diagonal decomposition before testing
        // for zero cells. If parallel running a local domain might have zero
        // cells so never construct the face-diagonal decomposition (which
        // uses parallel transfers)
        (void)tetBasePtIs();
    }

    if (nCells() == 0)
    {
        return -1;
    }

    // Find the nearest cell centre to this location
    label cellI = findNearestCell(location);

    // If point is in the nearest cell return
    if (pointInCell(location, cellI, decompMode))
    {
        return cellI;
    }
    else // point is not in the nearest cell so search all cells
    {
        for (label cellI = 0; cellI < nCells(); cellI++)
        {
            if (pointInCell(location, cellI, decompMode))
            {
                return cellI;
            }
        }
        return -1;
    }
}


// ************************************************************************* //
