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

Namespace
    Foam::fvMeshStitcherTools

Description
    Collection of free functions utilised by the stitching process

SourceFiles
    fvMeshStitcherTools.C
    fvMeshStitcherToolsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshStitcherTools_H
#define fvMeshStitcherTools_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshStitcherTools
{

//- Map a field with an (optional) addressing offset
template<class Type>
tmp<Field<Type>> fieldMap
(
    const Field<Type>& f,
    const labelUList& addr,
    const label addr0 = 0
);

//- Tmp variant of above
template<class Type>
tmp<Field<Type>> fieldMap
(
    const tmp<Field<Type>>& f,
    const labelUList& addr,
    const label addr0 = 0
);

//- Reverse map a field with an (optional) addressing offset, initialising the
//  value to zero and summing repeated indices
template<class Type>
tmp<Field<Type>> fieldRMapSum
(
    const Field<Type>& f,
    const label size,
    const labelUList& addr,
    const label addr0 = 0
);

//- Tmp variant of above
template<class Type>
tmp<Field<Type>> fieldRMapSum
(
    const tmp<Field<Type>>& f,
    const label size,
    const labelUList& addr,
    const label addr0 = 0
);

//- Alias for surface boundary fields to reduce verbosity of method
//  definitions below
template<class Type>
using SurfaceFieldBoundary =
    GeometricBoundaryField<Type, fvsPatchField, surfaceMesh>;

//- Return the total non-conformal area associated with each original face
tmp<SurfaceFieldBoundary<scalar>> origNcMagSfb(const fvMesh& mesh);

//- Extract the non-conformal parts of the boundary field and store it on the
//  conformal faces
template<class Type>
tmp<SurfaceFieldBoundary<Type>> conformedNcBoundaryField
(
    const SurfaceFieldBoundary<Type>& fieldb
);

//- Extract the original parts of the boundary field and store it
//  on the conformal faces
template<class Type>
tmp<SurfaceFieldBoundary<Type>> conformedOrigBoundaryField
(
    const SurfaceFieldBoundary<Type>& fieldb
);

//- Combine non-conformal and original parts of the boundary field from the
//  conformal faces to construct the complete non-conformal boundary field
template<class Type>
tmp<SurfaceFieldBoundary<Type>> unconformedBoundaryField
(
    const SurfaceFieldBoundary<Type>& ncFieldb,
    const SurfaceFieldBoundary<Type>& origFieldb
);

//- Synchronise the boundary field by combining corresponding
//  values across couplings with the given weightings
template<class Type>
tmp<SurfaceFieldBoundary<Type>> synchronisedBoundaryField
(
    const SurfaceFieldBoundary<Type>& fieldb,
    const bool flip,
    const scalar ownerWeight,
    const scalar neighbourWeight
);

//- Synchronise the boundary field by combining corresponding
//  values across couplings with equal weightings
template<class Type>
tmp<SurfaceFieldBoundary<Type>> synchronisedBoundaryField
(
    const SurfaceFieldBoundary<Type>& fieldb
);

} // End namespace patchToPatchTools
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMeshStitcherToolsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
