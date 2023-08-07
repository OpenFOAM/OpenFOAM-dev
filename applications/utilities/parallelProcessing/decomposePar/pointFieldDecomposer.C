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

#include "pointFieldDecomposer.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointFieldDecomposer::patchFieldDecomposer::patchFieldDecomposer
(
    const pointPatch& completeMeshPatch,
    const pointPatch& procMeshPatch,
    const labelList& directAddr
)
:
    pointPatchFieldMapperPatchRef
    (
        completeMeshPatch,
        procMeshPatch
    ),
    directAddressing_(procMeshPatch.size(), -1),
    hasUnmapped_(false)
{
    // Create the inverse-addressing of the patch point labels.
    labelList pointMap(completeMeshPatch.boundaryMesh().mesh().size(), -1);

    const labelList& completeMeshPatchPoints = completeMeshPatch.meshPoints();

    forAll(completeMeshPatchPoints, pointi)
    {
        pointMap[completeMeshPatchPoints[pointi]] = pointi;
    }

    // Use the inverse point addressing to create the addressing table for this
    // patch
    const labelList& procMeshPatchPoints = procMeshPatch.meshPoints();

    forAll(procMeshPatchPoints, pointi)
    {
        directAddressing_[pointi] =
            pointMap[directAddr[procMeshPatchPoints[pointi]]];
    }

    // Check that all the patch point addresses are set
    if (directAddressing_.size() && min(directAddressing_) < 0)
    {
        hasUnmapped_ = true;

        FatalErrorInFunction
            << "Incomplete patch point addressing"
            << abort(FatalError);
    }
}


Foam::pointFieldDecomposer::pointFieldDecomposer
(
    const pointMesh& completeMesh,
    const PtrList<fvMesh>& procMeshes,
    const labelListList& pointProcAddressing
)
:
    completeMesh_(completeMesh),
    procMeshes_(procMeshes),
    pointProcAddressing_(pointProcAddressing),
    patchFieldDecomposers_(procMeshes_.size())
{
    forAll(procMeshes_, proci)
    {
        const pointMesh& procMesh = pointMesh::New(procMeshes_[proci]);

        patchFieldDecomposers_.set
        (
            proci,
            new PtrList<patchFieldDecomposer>(procMesh.boundary().size())
        );

        forAll(procMesh.boundary(), procPatchi)
        {
            if (procPatchi < completeMesh_.boundary().size())
            {
                patchFieldDecomposers_[proci].set
                (
                    procPatchi,
                    new patchFieldDecomposer
                    (
                        completeMesh_.boundary()[procPatchi],
                        procMesh.boundary()[procPatchi],
                        pointProcAddressing_[proci]
                    )
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointFieldDecomposer::~pointFieldDecomposer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::pointFieldDecomposer::decomposes(const IOobjectList& objects)
{
    bool result = false;

    #define DO_POINT_FIELDS_TYPE(Type, nullArg)                                \
        result = result                                                        \
         || !objects.lookupClass(PointField<Type>::typeName).empty();
    FOR_ALL_FIELD_TYPES(DO_POINT_FIELDS_TYPE)
    #undef DO_POINT_FIELDS_TYPE

    return result;
}


// ************************************************************************* //
