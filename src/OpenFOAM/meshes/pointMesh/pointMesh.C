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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointMesh::pointMesh(const polyMesh& pMesh)
:
    MeshObject<polyMesh, Foam::PatchMeshObject, pointMesh>(pMesh),
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

void Foam::pointMesh::reset(const bool validBoundary)
{
    const polyMesh& pm = operator()();
    if (debug)
    {
        Pout<< "pointMesh::reset(const bool validBoundary): "
            << "Resetting from polyMesh " << pm.name() << endl;
    }

    boundary_.reset(pm.boundaryMesh());
    if (validBoundary)
    {
        // Calculate the geometry for the patches (transformation tensors etc.)
        boundary_.calcGeometry();
    }
}


void Foam::pointMesh::updateMesh(const mapPolyMesh& map)
{
    if (debug)
    {
        Pout<< "pointMesh::updateMesh(const mapPolyMesh&): "
            << "Updating for topology changes." << endl;
        Pout<< endl;
    }
    boundary_.updateMesh();
}


void Foam::pointMesh::distribute(const mapDistributePolyMesh& map)
{
    if (debug)
    {
        Pout<< "pointMesh::distribute(const mapDistributePolyMesh&): "
            << "Distribute." << endl;
        Pout<< endl;
    }
    boundary_.updateMesh();
}


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
        Pout<< endl;
    }

    boundary_.shuffle(newToOld, validBoundary);

    objectRegistry& db = const_cast<objectRegistry&>(thisDb());
    ReorderPatchFields<pointScalarField>(db, newToOld);
    ReorderPatchFields<pointVectorField>(db, newToOld);
    ReorderPatchFields<pointSphericalTensorField>(db, newToOld);
    ReorderPatchFields<pointSymmTensorField>(db, newToOld);
    ReorderPatchFields<pointTensorField>(db, newToOld);
}


void Foam::pointMesh::addPatch(const label patchi)
{
    if (debug)
    {
        Pout<< "pointMesh::addPatch(const label): "
            << "Adding patch at " << patchi << endl;
        Pout<< endl;
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

    objectRegistry& db = const_cast<objectRegistry&>(thisDb());
    const dictionary d;
    const word patchFieldType("calculated");

    AddPatchFields<pointScalarField>(db, patchi, d, patchFieldType, Zero);
    AddPatchFields<pointVectorField>(db, patchi, d, patchFieldType, Zero);
    AddPatchFields<pointSphericalTensorField>
    (
        db,
        patchi,
        d,
        patchFieldType,
        Zero
    );
    AddPatchFields<pointSymmTensorField>(db, patchi, d, patchFieldType, Zero);
    AddPatchFields<pointTensorField>(db, patchi, d, patchFieldType, Zero);
}


// ************************************************************************* //
