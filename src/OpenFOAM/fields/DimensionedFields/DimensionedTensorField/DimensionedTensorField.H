/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

InClass
    Foam::DimensionedTensorField

Description
    Tensor specific part of the implementation of DimensionedField.

SourceFiles
    DimensionedTensorField.C

\*---------------------------------------------------------------------------*/

#ifndef DimensionedTensorField_H
#define DimensionedTensorField_H

#include "DimensionedSphericalTensorField.H"
#include "DimensionedSymmTensorField.H"
#include "tensor.H"

#define TEMPLATE template<class GeoMesh>
#include "DimensionedFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_FUNCTION(scalar, tensor, tr, transform)
UNARY_FUNCTION(sphericalTensor, tensor, sph, transform)
UNARY_FUNCTION(symmTensor, tensor, symm, transform)
UNARY_FUNCTION(symmTensor, tensor, twoSymm, transform)
UNARY_FUNCTION(tensor, tensor, skew, transform)
UNARY_FUNCTION(tensor, tensor, dev, transform)
UNARY_FUNCTION(tensor, tensor, dev2, transform)
UNARY_FUNCTION(scalar, tensor, det, transform)
UNARY_FUNCTION(tensor, tensor, cof, cof)
UNARY_FUNCTION(tensor, tensor, inv, inv)
UNARY_FUNCTION(vector, tensor, eigenValues, transform)
UNARY_FUNCTION(tensor, tensor, eigenVectors, sign)

UNARY_FUNCTION(vector, symmTensor, eigenValues, transform)
UNARY_FUNCTION(symmTensor, symmTensor, eigenVectors, sign)


// * * * * * * * * * * * * * * * global operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(vector, tensor, *, hdual, transform)
UNARY_OPERATOR(tensor, vector, *, hdual, transform)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

#ifdef NoRepository
    #include "DimensionedTensorField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
