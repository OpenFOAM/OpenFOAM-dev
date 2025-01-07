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

#ifndef LagrangianFieldsFwd_H
#define LagrangianFieldsFwd_H

#include "DynamicFieldFwd.H"
#include "DimensionedFieldFwd.H"
#include "GeometricFieldFwd.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class LagrangianMesh;

template<class Type>
class LagrangianPatchField;

template<class Type>
class LagrangianFieldSource;

template<class Type>
using LagrangianField = GeometricField<Type, LagrangianMesh>;

template<class Type>
using LagrangianInternalField = DimensionedField<Type, LagrangianMesh>;

template<class Type>
using LagrangianPrimitiveDynamicField = DynamicField<Type>;

template<class Type>
using LagrangianDynamicField =
    GeometricField<Type, LagrangianMesh, LagrangianPrimitiveDynamicField>;

template<class Type>
using LagrangianInternalDynamicField =
    DimensionedField<Type, LagrangianMesh, LagrangianPrimitiveDynamicField>;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef LagrangianInternalField<label> LagrangianLabelInternalField;
typedef LagrangianInternalField<scalar> LagrangianScalarInternalField;
typedef LagrangianInternalField<vector> LagrangianVectorInternalField;
typedef LagrangianInternalField<sphericalTensor>
    LagrangianSphericalTensorInternalField;
typedef LagrangianInternalField<symmTensor> LagrangianSymmTensorInternalField;
typedef LagrangianInternalField<tensor> LagrangianTensorInternalField;

typedef LagrangianField<label> LagrangianLabelField;
typedef LagrangianField<scalar> LagrangianScalarField;
typedef LagrangianField<vector> LagrangianVectorField;
typedef LagrangianField<sphericalTensor> LagrangianSphericalTensorField;
typedef LagrangianField<symmTensor> LagrangianSymmTensorField;
typedef LagrangianField<tensor> LagrangianTensorField;

typedef LagrangianInternalDynamicField<label>
    LagrangianLabelInternalDynamicField;
typedef LagrangianInternalDynamicField<scalar>
    LagrangianScalarInternalDynamicField;
typedef LagrangianInternalDynamicField<vector>
    LagrangianVectorInternalDynamicField;
typedef LagrangianInternalDynamicField<sphericalTensor>
    LagrangianSphericalTensorInternalDynamicField;
typedef LagrangianInternalDynamicField<symmTensor>
    LagrangianSymmTensorInternalDynamicField;
typedef LagrangianInternalDynamicField<tensor>
    LagrangianTensorInternalDynamicField;

typedef LagrangianDynamicField<label> LagrangianLabelDynamicField;
typedef LagrangianDynamicField<scalar> LagrangianScalarDynamicField;
typedef LagrangianDynamicField<vector> LagrangianVectorDynamicField;
typedef LagrangianDynamicField<sphericalTensor>
    LagrangianSphericalTensorDynamicField;
typedef LagrangianDynamicField<symmTensor> LagrangianSymmTensorDynamicField;
typedef LagrangianDynamicField<tensor> LagrangianTensorDynamicField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
