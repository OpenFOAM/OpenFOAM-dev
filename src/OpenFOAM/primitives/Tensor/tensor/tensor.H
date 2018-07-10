/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Typedef
    Foam::tensor

Description
    Tensor of scalars.

SourceFiles
    tensor.C

\*---------------------------------------------------------------------------*/

#ifndef tensor_H
#define tensor_H

#include "Tensor.H"
#include "vector.H"
#include "sphericalTensor.H"
#include "symmTensor.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef Tensor<scalar> tensor;

vector eigenValues(const tensor& T);
vector eigenVector
(
    const tensor& T,
    const scalar lambda,
    const vector& direction1,
    const vector& direction2
);
tensor eigenVectors(const tensor& T, const vector& lambdas);
tensor eigenVectors(const tensor& T);

vector eigenValues(const symmTensor& T);
vector eigenVector
(
    const symmTensor& T,
    const scalar lambda,
    const vector& direction1,
    const vector& direction2
);
tensor eigenVectors(const symmTensor& T, const vector& lambdas);
tensor eigenVectors(const symmTensor& T);

//- Data associated with tensor type are contiguous
template<>
inline bool contiguous<tensor>() {return true;}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
