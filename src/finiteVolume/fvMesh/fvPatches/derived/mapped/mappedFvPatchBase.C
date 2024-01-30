/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "mappedFvPatchBase.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedFvPatchBase, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::mappedFvPatchBase::fromNeighbour
(
    const Field<Type>& nbrFld
) const
{
    return mapper_.fromNeighbour(nbrFld);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::mappedFvPatchBase::toNeighbour
(
    const Field<Type>& fld
) const
{
    return mapper_.toNeighbour(fld);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mappedFvPatchBase::mappedFvPatchBase(const fvPatch& patch)
:
    mappedFvPatchBaseBase(patch),
    mapper_(refCast<const mappedPatchBase>(patch.patch()))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedFvPatchBase::~mappedFvPatchBase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_MAPPED_FV_PATCH_BASE_FROM_AND_TO_NEIGHBOUR,
    mappedFvPatchBase
);


// ************************************************************************* //
