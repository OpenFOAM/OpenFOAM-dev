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

#include "lagrangianFieldDecomposer.H"
#include "passiveParticleCloud.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lagrangianFieldDecomposer::lagrangianFieldDecomposer
(
    const fvMesh& completeMesh,
    const PtrList<fvMesh>& procMeshes,
    const labelListList& faceProcAddressing,
    const labelListList& cellProcAddressing,
    const word& cloudName
)
:
    completeMesh_(completeMesh),
    procMeshes_(procMeshes),
    particleProcAddressing_(procMeshes_.size()),
    cloudName_(cloudName)
{
    // Create reverse cell addressing
    List<remote> cellProcCell(completeMesh_.nCells());
    forAll(cellProcAddressing, proci)
    {
        forAll(cellProcAddressing[proci], procCelli)
        {
            cellProcCell[cellProcAddressing[proci][procCelli]] =
                remote(proci, procCelli);
        }
    }

    // Create reverse face addressing
    List<remote> faceOwnerProcFace(completeMesh_.nFaces());
    List<remote> faceNeighbourProcFace(completeMesh_.nFaces());
    forAll(faceProcAddressing, proci)
    {
        forAll(faceProcAddressing[proci], procFacei)
        {
            const bool owner = faceProcAddressing[proci][procFacei] > 0;
            const label facei = mag(faceProcAddressing[proci][procFacei]) - 1;

            (owner ? faceOwnerProcFace : faceNeighbourProcFace)[facei] =
                remote(proci, procFacei);
        }
    }

    // Read the complete positions
    const passiveParticleCloud completePositions
    (
        completeMesh_,
        cloudName_,
        false
    );

    // Construct empty clouds for processor positions
    PtrList<passiveParticleCloud> procPositions(procMeshes_.size());
    forAll(procMeshes_, proci)
    {
        procPositions.set
        (
            proci,
            new passiveParticleCloud
            (
                procMeshes_[proci],
                cloudName_,
                IDLList<passiveParticle>()
            )
        );
    }

    // Count the number of particles on each processor
    labelList procNParticles(procMeshes_.size(), 0);
    forAllConstIter(passiveParticleCloud, completePositions, iter)
    {
        const passiveParticle& p = iter();
        const label proci = cellProcCell[p.cell()].proci;

        procNParticles[proci] ++;
    }

    // Resize the addressing
    forAll(procMeshes_, proci)
    {
        particleProcAddressing_[proci].resize(procNParticles[proci], -1);
    }

    // Distribute positions to the processor meshes
    label completeParticlei = 0;
    labelList procParticlei(procMeshes_.size(), 0);
    forAllConstIter(passiveParticleCloud, completePositions, iter)
    {
        const passiveParticle& p = iter();
        const label proci = cellProcCell[p.cell()].proci;
        const label procCelli = cellProcCell[p.cell()].elementi;
        const label procFacei =
            faceOwnerProcFace[p.tetFace()].proci == proci
          ? faceOwnerProcFace[p.tetFace()].elementi
          : faceNeighbourProcFace[p.tetFace()].elementi;

        particleProcAddressing_[proci][procParticlei[proci]] =
            completeParticlei;

        procPositions[proci].append
        (
            new passiveParticle
            (
                procMeshes_[proci],
                p.coordinates(),
                procCelli,
                procFacei,
                p.procTetPt
                (
                    completeMesh_,
                    procMeshes_[proci],
                    procCelli,
                    procFacei
                )
            )
        );

        completeParticlei ++;
        procParticlei[proci] ++;
    }

    // Write
    forAll(procPositions, proci)
    {
        IOPosition<passiveParticleCloud>(procPositions[proci]).write();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::lagrangianFieldDecomposer::~lagrangianFieldDecomposer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::lagrangianFieldDecomposer::decomposes(const IOobjectList& objects)
{
    bool result = false;

    #define DO_LAGRANGIAN_FIELDS_TYPE(Type, nullArg)                           \
        result = result                                                        \
         || !objects.lookupClass(IOField<Type>::typeName).empty()              \
         || !objects.lookupClass(IOField<Field<Type>>::typeName).empty()       \
         || !objects.lookupClass(CompactIOField<Field<Type>>::typeName).empty();
    DO_LAGRANGIAN_FIELDS_TYPE(label, )
    FOR_ALL_FIELD_TYPES(DO_LAGRANGIAN_FIELDS_TYPE)
    #undef DO_LAGRANGIAN_FIELDS_TYPE

    return result;
}


// ************************************************************************* //
