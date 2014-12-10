/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::sphericalTensorFieldField

Description
    Specialisation of FieldField<Field, T> for sphericalTensor.

SourceFiles
    sphericalTensorFieldField.C

\*---------------------------------------------------------------------------*/

#ifndef sphericalTensorFieldField_H
#define sphericalTensorFieldField_H

#include "FieldField.H"
#include "sphericalTensor.H"

#define TEMPLATE template<template<class> class Field>
#include "FieldFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_FUNCTION(scalar, sphericalTensor, tr)
UNARY_FUNCTION(sphericalTensor, sphericalTensor, sph)
UNARY_FUNCTION(scalar, sphericalTensor, det)
UNARY_FUNCTION(sphericalTensor, sphericalTensor, inv)


// * * * * * * * * * * * * * * * global operators  * * * * * * * * * * * * * //

BINARY_OPERATOR(sphericalTensor, scalar, sphericalTensor, /, divide)
BINARY_TYPE_OPERATOR(sphericalTensor, scalar, sphericalTensor, /, divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

#ifdef NoRepository
#   include "sphericalTensorFieldField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
