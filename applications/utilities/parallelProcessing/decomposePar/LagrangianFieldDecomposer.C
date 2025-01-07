/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "LagrangianFieldDecomposer.H"
#include "remote.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LagrangianFieldDecomposer::LagrangianFieldDecomposer
(
    const fvMesh& completeFvMesh,
    const PtrList<fvMesh>& procFvMeshes,
    const labelListList& faceProcAddressing,
    const labelListList& cellProcAddressing,
    const word& LagrangianName
)
:
    completeMesh_(completeFvMesh, LagrangianName),
    procMeshes_(procFvMeshes.size()),
    particleProcAddressing_(procFvMeshes.size())
{
    // Construct empty processor meshes
    forAll(procMeshes_, proci)
    {
        procMeshes_.set
        (
            proci,
            new LagrangianMesh
            (
                procFvMeshes[proci],
                LagrangianName,
                IOobject::NO_READ
            )
        );
    }

    // Create reverse cell addressing
    List<remote> completeCellProcCell(completeMesh_.mesh().nCells());
    forAll(cellProcAddressing, proci)
    {
        forAll(cellProcAddressing[proci], procCelli)
        {
            completeCellProcCell[cellProcAddressing[proci][procCelli]] =
                remote(proci, procCelli);
        }
    }

    // Create reverse face addressing
    List<remote> completeFaceOwnerProcFace(completeMesh_.mesh().nFaces());
    List<remote> completeFaceNeighbourProcFace(completeMesh_.mesh().nFaces());
    forAll(faceProcAddressing, proci)
    {
        forAll(faceProcAddressing[proci], procFacei)
        {
            const bool owner = faceProcAddressing[proci][procFacei] > 0;
            const label completeFacei =
                mag(faceProcAddressing[proci][procFacei]) - 1;

            (
                owner
              ? completeFaceOwnerProcFace
              : completeFaceNeighbourProcFace
            )[completeFacei] = remote(proci, procFacei);
        }
    }

    // Count the number of particles on each processor
    labelList procMeshSizes(procMeshes_.size(), 0);
    forAll(completeMesh_, i)
    {
        const label proci =
            completeCellProcCell[completeMesh_.celli()[i]].proci;

        procMeshSizes[proci] ++;
    }

    // Resize the addressing
    forAll(procMeshes_, proci)
    {
        particleProcAddressing_[proci].resize(procMeshSizes[proci], -1);
    }

    // Allocate processor geometry and topology
    PtrList<barycentricField> procCoordinates(procMeshes_.size());
    PtrList<labelField> procCellIndices(procMeshes_.size());
    PtrList<labelField> procFaceIndices(procMeshes_.size());
    PtrList<labelField> procFaceTriIndices(procMeshes_.size());
    forAll(procMeshes_, proci)
    {
        procCoordinates.set(proci, new barycentricField(procMeshSizes[proci]));
        procCellIndices.set(proci, new labelField(procMeshSizes[proci]));
        procFaceIndices.set(proci, new labelField(procMeshSizes[proci]));
        procFaceTriIndices.set(proci, new labelField(procMeshSizes[proci]));
    }

    // Distribute the elements to the processor meshes, and simultaneously
    // build addressing for distributing associated fields
    const faceList& completeFaces = completeMesh_.mesh().faces();
    labelList procIs(procMeshes_.size(), 0);
    forAll(completeMesh_, i)
    {
        const label completeCelli = completeMesh_.celli()[i];
        const label completeFacei = completeMesh_.facei()[i];

        const label proci = completeCellProcCell[completeCelli].proci;
        const label procCelli = completeCellProcCell[completeCelli].elementi;
        const label procFacei =
            completeFaceOwnerProcFace[completeFacei].proci == proci
          ? completeFaceOwnerProcFace[completeFacei].elementi
          : completeFaceNeighbourProcFace[completeFacei].elementi;

        particleProcAddressing_[proci][procIs[proci]] = i;

        procCoordinates[proci][procIs[proci]] = completeMesh_.coordinates()[i];
        procCellIndices[proci][procIs[proci]] = procCelli;
        procFaceIndices[proci][procIs[proci]] = procFacei;
        procFaceTriIndices[proci][procIs[proci]] = completeMesh_.faceTrii()[i];

        // If the owner has changed then the face will be numbered around in
        // the opposite direction. Change the face triangle index accordingly.
        if (faceProcAddressing[proci][procFacei] < 0)
        {
            const label completeFaceSize = completeFaces[completeFacei].size();

            procFaceTriIndices[proci][procIs[proci]] =
                completeFaceSize - 1 - procFaceTriIndices[proci][procIs[proci]];
        }

        procIs[proci] ++;
    }

    // Inject into the processor meshes
    forAll(procMeshes_, proci)
    {
        procMeshes_[proci].inject
        (
            procCoordinates[proci],
            procCellIndices[proci],
            procFaceIndices[proci],
            procFaceTriIndices[proci]
        );
    }

    // Write
    forAll(procMeshes_, proci)
    {
        procMeshes_[proci].write();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LagrangianFieldDecomposer::~LagrangianFieldDecomposer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::LagrangianFieldDecomposer::decomposes
(
    const IOobjectList& objects
)
{
    bool result = false;

    #define DECOMPOSES_LAGRANGIAN_FIELDS_TYPE(Type, nullArg)                   \
        result =                                                               \
            result                                                             \
         || decomposes<LagrangianField<Type>>(objects)                         \
         || decomposes<LagrangianDynamicField<Type>>(objects)                  \
         || decomposes<LagrangianInternalField<Type>>(objects)                 \
         || decomposes<LagrangianInternalDynamicField<Type>>(objects);
    DECOMPOSES_LAGRANGIAN_FIELDS_TYPE(label, )
    FOR_ALL_FIELD_TYPES(DECOMPOSES_LAGRANGIAN_FIELDS_TYPE)
    #undef DECOMPOSES_LAGRANGIAN_FIELDS_TYPE

    return result;
}


// ************************************************************************* //
