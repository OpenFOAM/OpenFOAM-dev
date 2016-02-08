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

\*---------------------------------------------------------------------------*/

#include "Constant.H"
#include "PolynomialEntry.H"
#include "Sine.H"
#include "CSV.H"
#include "Table.H"
#include "TableFile.H"

#include "label.H"
#include "scalar.H"
#include "vector.H"
#include "sphericalTensor.H"
#include "symmTensor.H"
#include "tensor.H"

namespace Foam
{
    makeFunction1(label);
    makeFunction1Type(Constant, label);
    // Polynomial functions and interpolation do evaluate to label
    // Instead evaluate a scalar and convert to label as appropriate

    makeFunction1(scalar);
    makeFunction1Type(Constant, scalar);
    makeFunction1Type(Polynomial, scalar);
    makeFunction1Type(Sine, scalar);
    makeFunction1Type(CSV, scalar);
    makeFunction1Type(Table, scalar);
    makeFunction1Type(TableFile, scalar);

    makeFunction1(vector);
    makeFunction1Type(Constant, vector);
    makeFunction1Type(Polynomial, vector);
    makeFunction1Type(Sine, vector);
    makeFunction1Type(CSV, vector);
    makeFunction1Type(Table, vector);
    makeFunction1Type(TableFile, vector);

    makeFunction1(sphericalTensor);
    makeFunction1Type(Constant, sphericalTensor);
    makeFunction1Type(Polynomial, sphericalTensor);
    makeFunction1Type(Sine, sphericalTensor);
    makeFunction1Type(CSV, sphericalTensor);
    makeFunction1Type(Table, sphericalTensor);
    makeFunction1Type(TableFile, sphericalTensor);

    makeFunction1(symmTensor);
    makeFunction1Type(Constant, symmTensor);
    makeFunction1Type(Polynomial, symmTensor);
    makeFunction1Type(Sine, symmTensor);
    makeFunction1Type(CSV, symmTensor);
    makeFunction1Type(Table, symmTensor);
    makeFunction1Type(TableFile, symmTensor);

    makeFunction1(tensor);
    makeFunction1Type(Constant, tensor);
    makeFunction1Type(Polynomial, tensor);
    makeFunction1Type(Sine, tensor);
    makeFunction1Type(CSV, tensor);
    makeFunction1Type(Table, tensor);
    makeFunction1Type(TableFile, tensor);
}


// ************************************************************************* //
