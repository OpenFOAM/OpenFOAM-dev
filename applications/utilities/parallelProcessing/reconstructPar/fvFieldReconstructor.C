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

#include "fvFieldReconstructor.H"
#include "processorCyclicFvPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::fvFieldReconstructor::completePatchID
(
    const label proci,
    const label procPatchi
) const
{
    const fvPatch& procPatch = procMeshes_[proci].boundary()[procPatchi];

    if (procPatchi < completeMesh_.boundary().size())
    {
        return procPatchi;
    }
    else if (isA<processorCyclicFvPatch>(procPatch))
    {
        return refCast<const processorCyclicFvPatch>(procPatch).referPatchID();
    }
    else
    {
        return -1;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvFieldReconstructor::fvFieldReconstructor
(
    const fvMesh& completeMesh,
    const PtrList<fvMesh>& procMeshes,
    const labelListList& faceProcAddressing,
    const labelListList& cellProcAddressing,
    const PtrList<surfaceLabelField::Boundary>& faceProcAddressingBf
)
:
    completeMesh_(completeMesh),
    procMeshes_(procMeshes),
    faceProcAddressing_(faceProcAddressing),
    cellProcAddressing_(cellProcAddressing),
    faceProcAddressingBf_(faceProcAddressingBf)
{
    forAll(procMeshes_, proci)
    {
        const fvMesh& procMesh = procMeshes_[proci];

        if
        (
            faceProcAddressing[proci].size() != procMesh.nFaces()
         || cellProcAddressing[proci].size() != procMesh.nCells()
        )
        {
            FatalErrorInFunction
                << "Size of maps does not correspond to size of mesh"
                << " for processor " << proci << endl
                << "faceProcAddressing : " << faceProcAddressing[proci].size()
                << " nFaces : " << procMesh.nFaces() << endl
                << "cellProcAddressing : " << cellProcAddressing[proci].size()
                << " nCell : " << procMesh.nCells() << endl
                << " nFaces : " << procMesh.boundary().size()
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvFieldReconstructor::reconstructs
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
)
{
    bool result = false;

    #define DO_FV_FIELDS_TYPE(Type, nullArg)                                   \
        result = result                                                        \
         || reconstructs<VolField<Type>::Internal>(objects, selectedFields)    \
         || reconstructs<VolField<Type>>(objects, selectedFields)              \
         || reconstructs<SurfaceField<Type>>(objects, selectedFields);
    FOR_ALL_FIELD_TYPES(DO_FV_FIELDS_TYPE)
    #undef DO_FV_FIELDS_TYPE

    return result;
}


// ************************************************************************* //
