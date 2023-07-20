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

#include "lagrangianFieldReconstructor.H"
#include "passiveParticleCloud.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lagrangianFieldReconstructor::lagrangianFieldReconstructor
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
    cloudName_(cloudName)
{
    // Construct and empty cloud for the complete positions
    passiveParticleCloud completePositions
    (
        completeMesh_,
        cloudName_,
        IDLList<passiveParticle>()
    );

    forAll(procMeshes_, proci)
    {
        // Read the processor positions
        Cloud<passiveParticle> procPositions
        (
            procMeshes_[proci],
            cloudName_,
            false
        );

        // Combine the processor's positions into the complete cloud
        forAllConstIter(Cloud<passiveParticle>, procPositions, iter)
        {
            const passiveParticle& p = iter();
            const label completeCelli = cellProcAddressing[proci][p.cell()];
            const label completeFacei =
                mag(faceProcAddressing[proci][p.tetFace()]) - 1;

            completePositions.append
            (
                new passiveParticle
                (
                    completeMesh_,
                    p.coordinates(),
                    completeCelli,
                    completeFacei,
                    p.procTetPt
                    (
                        procMeshes_[proci],
                        completeMesh_,
                        completeCelli,
                        completeFacei
                    )
                )
            );
        }
    }

    // Write
    IOPosition<Cloud<passiveParticle>>(completePositions).write();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::lagrangianFieldReconstructor::~lagrangianFieldReconstructor()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::lagrangianFieldReconstructor::reconstructs
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
)
{
    bool result = false;

    #define DO_LAGRANGIAN_FIELDS_TYPE(Type, nullArg)                           \
        result = result                                                        \
         || reconstructs<IOField<Type>>(objects, selectedFields)               \
         || reconstructs<IOField<Field<Type>>>(objects, selectedFields)        \
         || reconstructs<CompactIOField<Field<Type>>>(objects, selectedFields);
    DO_LAGRANGIAN_FIELDS_TYPE(label, )
    FOR_ALL_FIELD_TYPES(DO_LAGRANGIAN_FIELDS_TYPE)
    #undef DO_LAGRANGIAN_FIELDS_TYPE

    return result;
}


// ************************************************************************* //
