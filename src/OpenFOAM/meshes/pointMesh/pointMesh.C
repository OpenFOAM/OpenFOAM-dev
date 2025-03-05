/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    DemandDrivenMeshObject<polyMesh, RepatchableMeshObject, pointMesh>(pMesh),
    GeoMesh<polyMesh>(pMesh),
    boundary_(*this, pMesh.boundaryMesh())
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

bool Foam::pointMesh::movePoints()
{
    if (debug)
    {
        Pout<< "pointMesh::movePoints(const pointField&): "
            << "Moving points." << endl;
    }

    boundary_.movePoints(GeoMesh<polyMesh>::mesh_.points());

    return true;
}


void Foam::pointMesh::topoChange(const polyTopoChangeMap& map)
{
    if (debug)
    {
        Pout<< "pointMesh::topoChange(const polyTopoChangeMap&): "
            << "Topology change." << endl;
    }
    boundary_.topoChange();
}


void Foam::pointMesh::mapMesh(const polyMeshMap& map)
{
    if (debug)
    {
        Pout<< "pointMesh::mapMesh(const polyMeshMap&): "
            << "Mesh mapping." << endl;
    }
    boundary_.topoChange();
}


void Foam::pointMesh::distribute(const polyDistributionMap& map)
{
    if (debug)
    {
        Pout<< "pointMesh::distribute(const polyDistributionMap&): "
            << "Distribute." << endl;
    }
    boundary_.topoChange();
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

    #define ReorderPatchFieldsType(Type, nullArg)                              \
        ReorderPatchFields<PointField<Type>>                                   \
        (                                                                      \
            const_cast<objectRegistry&>(thisDb()),                             \
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

    const polyBoundaryMesh& pbm = mesh().boundaryMesh();

    if (pbm.size() != boundary_.size())
    {
        FatalErrorInFunction << "Problem :"
            << " pointBoundaryMesh size :" << boundary_.size()
            << " polyBoundaryMesh size :" << pbm.size()
            << exit(FatalError);
    }

    boundary_.set(patchi, facePointPatch::New(pbm[patchi], boundary_).ptr());

    #define AddPatchFieldsType(Type, nullArg)                                  \
        AddPatchFields<PointField<Type>>                                       \
        (                                                                      \
            const_cast<objectRegistry&>(thisDb()),                             \
            patchi,                                                            \
            calculatedPointPatchField<scalar>::typeName                        \
        );
    FOR_ALL_FIELD_TYPES(AddPatchFieldsType);
    #undef ReorderPatchFieldsType
}


void Foam::pointMesh::reset()
{
    if (debug)
    {
        Pout<< "pointMesh::reset(): "
            << "Mesh reset." << endl;
    }
    boundary_.reset();
}


// ************************************************************************* //
