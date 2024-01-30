/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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

#include "nonConformalMappedFvPatchBase.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalMappedFvPatchBase, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::nonConformalMappedFvPatchBase::fromNeighbour
(
    const Field<Type>& nbrFld
) const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    return
        map
        (
            nbrMesh().polyFacesBf()[nbrFvPatch().index()],
            nbrFld,
            mesh.polyFacesBf()[patch().index()]
        );
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::nonConformalMappedFvPatchBase::toNeighbour
(
    const Field<Type>& fld
) const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    return
        map
        (
            mesh.polyFacesBf()[patch().index()],
            fld,
            nbrMesh().polyFacesBf()[nbrFvPatch().index()]
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nonConformalMappedFvPatchBase::nonConformalMappedFvPatchBase
(
    const fvPatch& patch
)
:
    mappedFvPatchBaseBase(patch),
    mapper_(refCast<const nonConformalMappedPatchBase>(patch.patch()))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalMappedFvPatchBase::~nonConformalMappedFvPatchBase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

IMPLEMENT_MAPPED_FV_PATCH_BASE_FROM_AND_TO_NEIGHBOUR
(
    label,
    nonConformalMappedFvPatchBase
);


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_MAPPED_FV_PATCH_BASE_FROM_AND_TO_NEIGHBOUR,
    nonConformalMappedFvPatchBase
);


// ************************************************************************* //
