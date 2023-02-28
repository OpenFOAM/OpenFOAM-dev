/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

Class
    Foam::nonConformalBoundary

Description
    Mesh object that stores an all boundary patch and mapping to and from it
    and the mesh and the individual patches

SourceFiles
    nonConformalBoundary.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalBoundary_H
#define nonConformalBoundary_H

#include "DemandDrivenMeshObject.H"
#include "polyMesh.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class nonConformalCoupledPolyPatch;

/*---------------------------------------------------------------------------*\
                          Class nonConformalBoundary Declaration
\*---------------------------------------------------------------------------*/

class nonConformalBoundary
:
    public DemandDrivenMeshObject
    <
        polyMesh,
        MoveableMeshObject,
        nonConformalBoundary
    >
{
    // Private Data

        //- Primitive patch of the owner-orig boundary
        indirectPrimitivePatch ownerOrigBoundary_;

        //- A map from owner-orig boundary point to mesh point
        mutable autoPtr<labelList> meshPointOwnerOrigBoundaryPointPtr_;

        //- A map from owner-orig boundary point to mesh point
        mutable autoPtr<labelList> ownerOrigBoundaryPointMeshPointPtr_;

        //- A map from owner-orig boundary edge to mesh edge
        mutable autoPtr<labelList> ownerOrigBoundaryEdgeMeshEdgePtr_;

        //- Owner-orig boundary edges, referring to boundary points
        mutable autoPtr<edgeList> ownerOrigBoundaryEdgesPtr_;

        //- Owner-orig boundary edges, referring to mesh points
        mutable autoPtr<edgeList> ownerOrigBoundaryMeshEdgesPtr_;

        //- A map from patch point to owner-orig boundary points
        mutable PtrList<labelList> patchPointOwnerOrigBoundaryPointsPtr_;

        //- A map from patch edge to owner-orig boundary edge
        mutable PtrList<labelList> patchEdgeOwnerOrigBoundaryEdgesPtr_;

        //- Point normals for the owner-orig boundary
        mutable autoPtr<vectorField> ownerOrigBoundaryPointNormalsPtr_;

        //- Old point normals for the owner-orig boundary
        mutable autoPtr<vectorField> ownerOrigBoundaryPointNormals0Ptr_;


    // Private Member Functions

        //- Construct a boundary for the given patch set
        indirectPrimitivePatch boundary(const labelList& patches) const;

        //- Get a map from mesh point to owner-orig boundary points. Non
        //  owner-orig boundary points are indicated by the value -1.
        const labelList& meshPointOwnerOrigBoundaryPoint() const;

        //- Get point normals for the owner-orig boundary
        const vectorField& ownerOrigBoundaryPointNormals() const;

        //- Get old point normals for the owner-orig boundary
        const vectorField& ownerOrigBoundaryPointNormals0() const;


protected:

    friend class DemandDrivenMeshObject
    <
        polyMesh,
        MoveableMeshObject,
        nonConformalBoundary
    >;

    // Protected Constructors

        //- Construct from mesh
        explicit nonConformalBoundary(const polyMesh& mesh);


public:

    //- Runtime type information
    TypeName("nonConformalBoundary");


    // Constructors

        //- Disallow default bitwise copy construction
        nonConformalBoundary(const nonConformalBoundary&) = delete;


    //- Destructor
    ~nonConformalBoundary();


    // Member Functions

        //- Update for mesh motion
        virtual bool movePoints();

        //- Return a list of non-conformal non-coupled patch indices
        labelList nonConformalNonCoupledPatchIDs
        (
            const label side,
            label (nonConformalCoupledPolyPatch::*method)() const
        ) const;

        //- Return a list of the orig patch indices
        labelList allOrigPatchIDs() const;

        //- Return a list of the error patch indices
        labelList allErrorPatchIDs() const;

        //- Return a list of the owner-orig patch indices
        labelList ownerOrigPatchIDs() const;

        //- Return a list of the owner error patch indices
        labelList ownerErrorPatchIDs() const;

        //- Get a map from owner-orig boundary point to mesh point
        const labelList& ownerOrigBoundaryPointMeshPoint() const;

        //- Get a map from owner-orig boundary edge to mesh edge
        const labelList& ownerOrigBoundaryEdgeMeshEdge() const;

        //- Get the owner-orig boundary edges, referring to boundary points
        const edgeList& ownerOrigBoundaryEdges() const;

        //- Get the owner-orig boundary edges, referring to mesh points
        const edgeList& ownerOrigBoundaryMeshEdges() const;

        //- Get a map from patch point to owner-orig boundary point
        const labelList& patchPointOwnerOrigBoundaryPoints
        (
            const label patchi
        ) const;

        //- Get a map from patch edge to owner-orig boundary edge
        const labelList& patchEdgeOwnerOrigBoundaryEdges
        (
            const label patchi
        ) const;

        //- Get parallel consistent point normals for the patch
        tmp<vectorField> patchPointNormals(const label patchi) const;

        //- Get parallel consistent old-time point normals for the patch
        tmp<vectorField> patchPointNormals0(const label patchi) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const nonConformalBoundary&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
