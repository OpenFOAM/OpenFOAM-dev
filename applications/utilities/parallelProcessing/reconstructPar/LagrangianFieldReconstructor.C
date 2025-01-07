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

#include "LagrangianFieldReconstructor.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LagrangianFieldReconstructor::LagrangianFieldReconstructor
(
    const fvMesh& completeFvMesh,
    const PtrList<fvMesh>& procFvMeshes,
    const labelListList& faceProcAddressing,
    const labelListList& cellProcAddressing,
    const word& LagrangianName
)
:
    completeMesh_(completeFvMesh, LagrangianName, IOobject::NO_READ),
    procMeshes_(procFvMeshes.size())
{
    // Read the processor meshes
    forAll(procMeshes_, proci)
    {
        procMeshes_.set
        (
            proci,
            new LagrangianMesh(procFvMeshes[proci], LagrangianName)
        );
    }

    // Determine the size of the complete mesh
    label completeSize = 0;
    forAll(procMeshes_, proci)
    {
        completeSize += procMeshes_[proci].size();
    }

    // Construct complete geometry and topology
    barycentricField completeCoordinates(completeSize);
    labelField completeCellIndices(completeSize, -1);
    labelField completeFaceIndices(completeSize, -1);
    labelField completeFaceTriIndices(completeSize, -1);
    label i0 = 0;
    forAll(procMeshes_, proci)
    {
        const label procSize = procMeshes_[proci].size();

        const labelField& procCelli = procMeshes_[proci].celli();
        const labelField& procFacei = procMeshes_[proci].facei();

        SubList<barycentric>(completeCoordinates, procSize, i0) =
            procMeshes_[proci].coordinates();
        SubList<label>(completeCellIndices, procSize, i0) =
            UIndirectList<label>(cellProcAddressing[proci], procCelli)();
        SubList<label>(completeFaceIndices, procSize, i0) =
            UIndirectList<label>(faceProcAddressing[proci], procFacei)();
        SubList<label>(completeFaceTriIndices, procSize, i0) =
            procMeshes_[proci].faceTrii();

        i0 += procMeshes_[proci].size();
    }

    // If the owner has changed then the face will be numbered around in
    // the opposite direction. Change the face triangle index accordingly.
    const faceList& completeFaces = completeMesh_.mesh().faces();
    forAll(completeCoordinates, i)
    {
        if (completeFaceIndices[i] < 0)
        {
            const label completeFacei = - completeFaceIndices[i] - 1;
            const label completeFaceSize = completeFaces[completeFacei].size();

            completeFaceTriIndices[i] =
                completeFaceSize - 1 - completeFaceTriIndices[i];
        }
    }

    // Remove the turning information to recover the actual face indices
    completeFaceIndices = mag(completeFaceIndices) - 1;

    // Inject into the complete mesh
    completeMesh_.inject
    (
        completeCoordinates,
        completeCellIndices,
        completeFaceIndices,
        completeFaceTriIndices
    );

    // Write
    completeMesh_.write();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LagrangianFieldReconstructor::~LagrangianFieldReconstructor()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::LagrangianFieldReconstructor::reconstructs
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
)
{
    bool result = false;

    #define RECONSTRUCTS_LAGRANGIAN_FIELDS_TYPE(Type, nullArg)                 \
        result =                                                               \
            result                                                             \
         || reconstructs<LagrangianField<Type>>                                \
            (                                                                  \
                objects,                                                       \
                selectedFields                                                 \
            )                                                                  \
         || reconstructs<LagrangianDynamicField<Type>>                         \
            (                                                                  \
                objects,                                                       \
                selectedFields                                                 \
            )                                                                  \
         || reconstructs<LagrangianInternalField<Type>>                        \
            (                                                                  \
                objects,                                                       \
                selectedFields                                                 \
            )                                                                  \
         || reconstructs<LagrangianInternalDynamicField<Type>>                 \
            (                                                                  \
                objects,                                                       \
                selectedFields                                                 \
            );
    RECONSTRUCTS_LAGRANGIAN_FIELDS_TYPE(label, )
    FOR_ALL_FIELD_TYPES(RECONSTRUCTS_LAGRANGIAN_FIELDS_TYPE)
    #undef RECONSTRUCTS_LAGRANGIAN_FIELDS_TYPE

    return result;
}


// ************************************************************************* //
