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
    makeDataEntry(label);
    makeDataEntryType(Constant, label);
    // Polynomial functions and interpolation do evaluate to label
    // Instead evaluate a scalar and convert to label as appropriate

    makeDataEntry(scalar);
    makeDataEntryType(Constant, scalar);
    makeDataEntryType(Polynomial, scalar);
    makeDataEntryType(Sine, scalar);
    makeDataEntryType(CSV, scalar);
    makeDataEntryType(Table, scalar);
    makeDataEntryType(TableFile, scalar);

    makeDataEntry(vector);
    makeDataEntryType(Constant, vector);
    makeDataEntryType(Polynomial, vector);
    makeDataEntryType(Sine, vector);
    makeDataEntryType(CSV, vector);
    makeDataEntryType(Table, vector);
    makeDataEntryType(TableFile, vector);

    makeDataEntry(sphericalTensor);
    makeDataEntryType(Constant, sphericalTensor);
    makeDataEntryType(Polynomial, sphericalTensor);
    makeDataEntryType(Sine, sphericalTensor);
    makeDataEntryType(CSV, sphericalTensor);
    makeDataEntryType(Table, sphericalTensor);
    makeDataEntryType(TableFile, sphericalTensor);

    makeDataEntry(symmTensor);
    makeDataEntryType(Constant, symmTensor);
    makeDataEntryType(Polynomial, symmTensor);
    makeDataEntryType(Sine, symmTensor);
    makeDataEntryType(CSV, symmTensor);
    makeDataEntryType(Table, symmTensor);
    makeDataEntryType(TableFile, symmTensor);

    makeDataEntry(tensor);
    makeDataEntryType(Constant, tensor);
    makeDataEntryType(Polynomial, tensor);
    makeDataEntryType(Sine, tensor);
    makeDataEntryType(CSV, tensor);
    makeDataEntryType(Table, tensor);
    makeDataEntryType(TableFile, tensor);
}


// ************************************************************************* //
