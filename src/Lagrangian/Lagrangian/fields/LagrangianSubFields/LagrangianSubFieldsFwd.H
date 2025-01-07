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

#ifndef LagrangianSubFieldsFwd_H
#define LagrangianSubFieldsFwd_H

#include "DimensionedFieldFwd.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class LagrangianSubMesh;

template<class>
class SubField;

template<class Type, template<class> class PrimitiveField = Field>
using LagrangianSubField =
    DimensionedField<Type, LagrangianSubMesh, PrimitiveField>;

template<class Type>
using LagrangianSubSubField =
    DimensionedField<Type, LagrangianSubMesh, SubField>;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef LagrangianSubField<label> LagrangianSubLabelField;
typedef LagrangianSubField<scalar> LagrangianSubScalarField;
typedef LagrangianSubField<vector> LagrangianSubVectorField;
typedef LagrangianSubField<sphericalTensor> LagrangianSubSphericalTensorField;
typedef LagrangianSubField<symmTensor> LagrangianSubSymmTensorField;
typedef LagrangianSubField<tensor> LagrangianSubTensorField;

typedef LagrangianSubSubField<label> LagrangianSubLabelSubField;
typedef LagrangianSubSubField<scalar> LagrangianSubScalarSubField;
typedef LagrangianSubSubField<vector> LagrangianSubVectorSubField;
typedef LagrangianSubSubField<sphericalTensor>
    LagrangianSubSphericalTensorSubField;
typedef LagrangianSubSubField<symmTensor> LagrangianSubSymmTensorSubField;
typedef LagrangianSubSubField<tensor> LagrangianSubTensorSubField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
