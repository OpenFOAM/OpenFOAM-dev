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

#include "tensorField.H"
#include "symmTensorField.H"
#include "transformField.H"

#define TEMPLATE
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

UNARY_FUNCTION(scalar, tensor, tr)
UNARY_FUNCTION(sphericalTensor, tensor, sph)
UNARY_FUNCTION(symmTensor, tensor, symm)
UNARY_FUNCTION(symmTensor, tensor, twoSymm)
UNARY_FUNCTION(tensor, tensor, skew)
UNARY_FUNCTION(tensor, tensor, dev)
UNARY_FUNCTION(tensor, tensor, dev2)
UNARY_FUNCTION(scalar, tensor, det)
UNARY_FUNCTION(tensor, tensor, cof)
UNARY_FUNCTION(tensor, tensor, inv)

void inv
(
    Field<tensor>& tf,
    const UList<tensor>& tf1,
    const Vector<label>& solutionD
)
{
    if (tf.empty())
    {
        return;
    }

    if (solutionD.x() == -1 || solutionD.y() == -1 || solutionD.z() == -1)
    {
        tensorField tf1Plus(tf1);

        if (solutionD.x() == -1)
        {
            tf1Plus += tensor(1,0,0,0,0,0,0,0,0);
        }

        if (solutionD.y() == -1)
        {
            tf1Plus += tensor(0,0,0,0,1,0,0,0,0);
        }

        if (solutionD.z() == -1)
        {
            tf1Plus += tensor(0,0,0,0,0,0,0,0,1);
        }

        TFOR_ALL_F_OP_FUNC_F(tensor, tf, =, inv, tensor, tf1Plus)

        if (solutionD.x() == -1)
        {
            tf -= tensor(1,0,0,0,0,0,0,0,0);
        }

        if (solutionD.y() == -1)
        {
            tf -= tensor(0,0,0,0,1,0,0,0,0);
        }

        if (solutionD.z() == -1)
        {
            tf -= tensor(0,0,0,0,0,0,0,0,1);
        }
    }
    else
    {
        TFOR_ALL_F_OP_FUNC_F(tensor, tf, =, inv, tensor, tf1)
    }
}

tmp<tensorField> inv
(
    const Field<tensor>& tf,
    const Vector<label>& solutionD
)
{
    tmp<tensorField> result(new tensorField(tf.size()));
    inv(result.ref(), tf);
    return result;
}

tmp<tensorField> inv
(
    const SubField<tensor>& tf,
    const Vector<label>& solutionD
)
{
    tmp<tensorField> result(new tensorField(tf.size()));
    inv(result.ref(), tf);
    return result;
}

tmp<tensorField> inv
(
    const tmp<tensorField>& tf,
    const Vector<label>& solutionD
)
{
    tmp<tensorField> tRes = New(tf);
    inv(tRes.ref(), tf());
    tf.clear();
    return tRes;
}

UNARY_FUNCTION(vector, tensor, eigenValues)
UNARY_FUNCTION(tensor, tensor, eigenVectors)

UNARY_FUNCTION(vector, symmTensor, eigenValues)
UNARY_FUNCTION(tensor, symmTensor, eigenVectors)


template<>
tmp<Field<tensor>> transformFieldMask<tensor>
(
    const symmTensorField& stf
)
{
    tmp<tensorField> tRes(new tensorField(stf.size()));
    tensorField& res = tRes.ref();
    TFOR_ALL_F_OP_F(tensor, res, =, symmTensor, stf)
    return tRes;
}

template<>
tmp<Field<tensor>> transformFieldMask<tensor>
(
    const tmp<symmTensorField>& tstf
)
{
    tmp<Field<tensor>> ret = transformFieldMask<tensor>(tstf());
    tstf.clear();
    return ret;
}


// * * * * * * * * * * * * * * * global operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(vector, tensor, *, hdual)
UNARY_OPERATOR(tensor, vector, *, hdual)

BINARY_OPERATOR(vector, vector, tensor, /, divide)
BINARY_TYPE_OPERATOR(vector, vector, tensor, /, divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
