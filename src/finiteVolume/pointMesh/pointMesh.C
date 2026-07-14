/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "pointMesh.H"
#include "globalMeshData.H"
#include "pointFields.H"
#include "SubField.H"
#include "facePointPatch.H"
#include "MapGeometricFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointMesh, 0);
}

const Foam::HashSet<Foam::word> Foam::pointMesh::geometryFields;

const Foam::HashSet<Foam::word> Foam::pointMesh::curGeometryFields;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointMesh::pointMesh(const polyMesh& pMesh)
:
    DemandDrivenMeshObject<polyMesh, PermanentMeshObject, pointMesh>(pMesh),
    mesh_(pMesh),
    boundary_(*this, pMesh.boundary())
{
    if (debug)
    {
        Pout<< "pointMesh::pointMesh(const polyMesh&): "
            << "Constructing from polyMesh " << pMesh.name()
            << endl;
    }

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointMesh::~pointMesh()
{
    if (debug)
    {
        Pout<< "~pointMesh::pointMesh()"
            << endl;
        error::printStack(Pout);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::DimensionedField<Foam::vector, Foam::pointMesh>&
Foam::pointMesh::C() const
{
    if (!pointsPtr_.valid())
    {
        pointsPtr_ = new SlicedDimensionedField<vector, pointMesh>
        (
            IOobject
            (
                "points",
                time().name(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                true
            ),
            *this,
            dimensions::length,
            mesh_.points()
        );
    }

    return *pointsPtr_;
}


bool Foam::pointMesh::movePoints()
{
    if (debug)
    {
        Pout<< "pointMesh::movePoints(const pointField&): "
            << "Moving points." << endl;
    }

    boundary_.movePoints(mesh_.points());

    // Move the points of all pointMesh meshObjects
    meshObjects::movePoints<pointMesh>(const_cast<polyMesh&>(mesh_));

    return true;
}


void Foam::pointMesh::topoChange(const polyTopoChangeMap& map)
{
    if (debug)
    {
        Pout<< "pointMesh::topoChange(const polyTopoChangeMap&): "
            << "Topology change." << endl;
    }

    pointsPtr_.clear();
    boundary_.topoChange();

    meshObjects::topoChange<pointMesh>(const_cast<polyMesh&>(mesh_), map);
}


void Foam::pointMesh::mapMesh(const polyMeshMap& map)
{
    if (debug)
    {
        Pout<< "pointMesh::mapMesh(const polyMeshMap&): "
            << "Mesh mapping." << endl;
    }

    pointsPtr_.clear();
    boundary_.topoChange();

    meshObjects::mapMesh<pointMesh>(const_cast<polyMesh&>(mesh_), map);
}


void Foam::pointMesh::swap(polyMesh& otherMesh)
{
    if (debug)
    {
        Pout<< "pointMesh::swap(const polyMesh&): "
            << "Mesh swapping." << endl;
    }

    boundary_.reset();

    meshObjects::clearUpto
    <
        pointMesh,
        DeletableMeshObject,
        TopoChangeableMeshObject
    >
    (
        const_cast<polyMesh&>(mesh_)
    );

    meshObjects::clearUpto
    <
        pointMesh,
        DeletableMeshObject,
        TopoChangeableMeshObject
    >
    (
        otherMesh
    );
}


void Foam::pointMesh::distribute(const polyDistributionMap& map)
{
    if (debug)
    {
        Pout<< "pointMesh::distribute(const polyDistributionMap&): "
            << "Distribute." << endl;
    }

    pointsPtr_.clear();
    boundary_.topoChange();

    meshObjects::distribute<pointMesh>(const_cast<polyMesh&>(mesh_), map);
}


void Foam::pointMesh::reorderPatches
(
    const labelUList& newToOld,
    const bool validBoundary
)
{
    if (debug)
    {
        Pout<< "pointMesh::reorderPatches( const labelUList&, const bool): "
            << "Updating for reordered patches." << endl;
    }

    boundary_.shuffle(newToOld, validBoundary);
    meshObjects::reorderPatches<pointMesh>
    (
        const_cast<polyMesh&>(mesh_),
        newToOld,
        validBoundary
    );

    #define ReorderPatchFieldsType(Type, nullArg)                              \
        ReorderPatchFields<PointField<Type>>                                   \
        (                                                                      \
            const_cast<objectRegistry&>(db()),                                 \
            newToOld                                                           \
        );
    FOR_ALL_FIELD_TYPES(ReorderPatchFieldsType);
    #undef ReorderPatchFieldsType
}


void Foam::pointMesh::addPatch(const label patchi)
{
    if (debug)
    {
        Pout<< "pointMesh::addPatch(const label): "
            << "Adding patch at " << patchi << endl;
    }

    const polyBoundaryMesh& pbm = poly().boundary();

    if (pbm.size() != boundary_.size())
    {
        FatalErrorInFunction << "Problem :"
            << " pointBoundaryMesh size :" << boundary_.size()
            << " polyBoundaryMesh size :" << pbm.size()
            << exit(FatalError);
    }

    boundary_.set(patchi, facePointPatch::New(pbm[patchi], boundary_).ptr());
    meshObjects::addPatch<pointMesh>(const_cast<polyMesh&>(mesh_), patchi);

    #define AddPatchFieldsType(Type, nullArg)                                  \
        AddPatchFields<PointField<Type>>                                       \
        (                                                                      \
            const_cast<objectRegistry&>(db()),                                 \
            patchi,                                                            \
            calculatedPointPatchField<scalar>::typeName                        \
        );
    FOR_ALL_FIELD_TYPES(AddPatchFieldsType);
    #undef ReorderPatchFieldsType
}


void Foam::pointMesh::clear()
{
    if (debug)
    {
        Pout<< "pointMesh::clear(): "
            << "Clear all but permanent registered meshObjects" << endl;
    }

    // Clear all but permanent registered meshObjects
    meshObjects::clear<pointMesh, DeletableMeshObject>
    (
        const_cast<polyMesh&>(mesh_)
    );
}


void Foam::pointMesh::reset()
{
    if (debug)
    {
        Pout<< "pointMesh::reset(): "
            << "Mesh reset." << endl;
    }

    boundary_.reset();
    meshObjects::reset<pointMesh>(const_cast<polyMesh&>(mesh_));
}


// ************************************************************************* //
