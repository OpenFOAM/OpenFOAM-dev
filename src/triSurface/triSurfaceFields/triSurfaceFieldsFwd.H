/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#ifndef triSurfaceFieldsFwd_H
#define triSurfaceFieldsFwd_H

#include "DimensionedFieldFwd.H"
#include "fieldTypes.H"
#include "triSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef DimensionedField<label, triSurface, Field> triSurfaceLabelField;
typedef DimensionedField<scalar, triSurface, Field>
    triSurfaceScalarField;
typedef DimensionedField<vector, triSurface, Field>
    triSurfaceVectorField;
typedef DimensionedField<sphericalTensor, triSurface, Field>
    triSurfaceSphericalTensorField;
typedef DimensionedField<symmTensor, triSurface, Field>
    triSurfaceSymmTensorField;
typedef DimensionedField<tensor, triSurface, Field>
    triSurfaceTensorField;

typedef DimensionedField<label, triSurfacePointMesh, Field>
    triSurfacePointLabelField;
typedef DimensionedField<scalar, triSurfacePointMesh, Field>
    triSurfacePointScalarField;
typedef DimensionedField<vector, triSurfacePointMesh, Field>
    triSurfacePointVectorField;
typedef DimensionedField<sphericalTensor, triSurfacePointMesh, Field>
    triSurfacePointSphericalTensorField;
typedef DimensionedField<symmTensor, triSurfacePointMesh, Field>
    triSurfacePointSymmTensorField;
typedef DimensionedField<tensor, triSurfacePointMesh, Field>
    triSurfacePointTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
