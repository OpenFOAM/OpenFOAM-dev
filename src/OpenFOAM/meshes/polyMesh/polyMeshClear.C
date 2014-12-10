/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
#include "primitiveMesh.H"
#include "globalMeshData.H"
#include "MeshObject.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "pointMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::polyMesh::removeBoundary()
{
    if (debug)
    {
        Info<< "void polyMesh::removeBoundary(): "
            << "Removing boundary patches."
            << endl;
    }

    // Remove the point zones
    boundary_.clear();
    boundary_.setSize(0);

    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyMesh::clearGeom()
{
    if (debug)
    {
        Info<< "void polyMesh::clearGeom() : "
            << "clearing geometric data"
            << endl;
    }

    // Clear all geometric mesh objects
    meshObject::clear<pointMesh, GeometricMeshObject>(*this);
    meshObject::clear<polyMesh, GeometricMeshObject>(*this);

    primitiveMesh::clearGeom();

    boundary_.clearGeom();

    // Reset valid directions (could change with rotation)
    geometricD_ = Vector<label>::zero;
    solutionD_ = Vector<label>::zero;

    // Remove the stored tet base points
    tetBasePtIsPtr_.clear();
    // Remove the cell tree
    cellTreePtr_.clear();
}


void Foam::polyMesh::clearAdditionalGeom()
{
    if (debug)
    {
        Info<< "void polyMesh::clearAdditionalGeom() : "
            << "clearing additional geometric data"
            << endl;
    }

    // Remove the stored tet base points
    tetBasePtIsPtr_.clear();
    // Remove the cell tree
    cellTreePtr_.clear();
}


void Foam::polyMesh::clearAddressing(const bool isMeshUpdate)
{
    if (debug)
    {
        Info<< "void polyMesh::clearAddressing() : "
            << "clearing topology  isMeshUpdate:" << isMeshUpdate
            << endl;
    }

    if (isMeshUpdate)
    {
        // Part of a mesh update. Keep meshObjects that have an updateMesh
        // callback
        meshObject::clearUpto
        <
            pointMesh,
            TopologicalMeshObject,
            UpdateableMeshObject
        >
        (
            *this
        );
        meshObject::clearUpto
        <
            polyMesh,
            TopologicalMeshObject,
            UpdateableMeshObject
        >
        (
            *this
        );
    }
    else
    {
        meshObject::clear<pointMesh, TopologicalMeshObject>(*this);
        meshObject::clear<polyMesh, TopologicalMeshObject>(*this);
    }

    primitiveMesh::clearAddressing();

    // parallelData depends on the processorPatch ordering so force
    // recalculation
    globalMeshDataPtr_.clear();

    // Reset valid directions
    geometricD_ = Vector<label>::zero;
    solutionD_ = Vector<label>::zero;

    // Update zones
    pointZones_.clearAddressing();
    faceZones_.clearAddressing();
    cellZones_.clearAddressing();

    // Remove the stored tet base points
    tetBasePtIsPtr_.clear();

    // Remove the cell tree
    cellTreePtr_.clear();
}


void Foam::polyMesh::clearPrimitives()
{
    resetMotion();

    points_.setSize(0);
    faces_.setSize(0);
    owner_.setSize(0);
    neighbour_.setSize(0);

    clearedPrimitives_ = true;
}


void Foam::polyMesh::clearOut()
{
    clearGeom();
    clearAddressing();
}


void Foam::polyMesh::clearCellTree()
{
    if (debug)
    {
        Info<< "void polyMesh::clearCellTree() : "
            << "clearing cell tree"
            << endl;
    }

    cellTreePtr_.clear();
}


// ************************************************************************* //
