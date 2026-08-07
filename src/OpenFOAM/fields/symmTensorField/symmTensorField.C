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

#include "symmTensorField.H"
#include "transformField.H"

#define TEMPLATE
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

UNARY_FUNCTION(symmTensor, vector, sqr)
UNARY_FUNCTION(symmTensor, symmTensor, innerSqr)

UNARY_FUNCTION(scalar, symmTensor, tr)
UNARY_FUNCTION(sphericalTensor, symmTensor, sph)
UNARY_FUNCTION(symmTensor, symmTensor, symm)
UNARY_FUNCTION(symmTensor, symmTensor, twoSymm)
UNARY_FUNCTION(symmTensor, symmTensor, dev)
UNARY_FUNCTION(symmTensor, symmTensor, dev2)
UNARY_FUNCTION(scalar, symmTensor, det)
UNARY_FUNCTION(symmTensor, symmTensor, cof)
UNARY_FUNCTION(symmTensor, symmTensor, inv)

void inv
(
    Field<symmTensor>& tf,
    const UList<symmTensor>& tf1,
    const Vector<label>& solutionD
)
{
    if (tf.empty())
    {
        return;
    }

    if (solutionD.x() == -1 || solutionD.y() == -1 || solutionD.z() == -1)
    {
        symmTensorField tf1Plus(tf1);

        if (solutionD.x() == -1)
        {
            tf1Plus += symmTensor(1,0,0,0,0,0);
        }

        if (solutionD.y() == -1)
        {
            tf1Plus += symmTensor(0,0,0,1,0,0);
        }

        if (solutionD.z() == -1)
        {
            tf1Plus += symmTensor(0,0,0,0,0,1);
        }

        TFOR_ALL_F_OP_FUNC_F(symmTensor, tf, =, inv, symmTensor, tf1Plus)

        if (solutionD.x() == -1)
        {
            tf -= symmTensor(1,0,0,0,0,0);
        }

        if (solutionD.y() == -1)
        {
            tf -= symmTensor(0,0,0,1,0,0);
        }

        if (solutionD.z() == -1)
        {
            tf -= symmTensor(0,0,0,0,0,1);
        }
    }
    else
    {
        TFOR_ALL_F_OP_FUNC_F(symmTensor, tf, =, inv, symmTensor, tf1)
    }
}

tmp<symmTensorField> inv
(
    const Field<symmTensor>& tf,
    const Vector<label>& solutionD
)
{
    tmp<symmTensorField> result(new symmTensorField(tf.size()));
    inv(result.ref(), tf, solutionD);
    return result;
}

tmp<symmTensorField> inv
(
    const SubField<symmTensor>& tf,
    const Vector<label>& solutionD
)
{
    tmp<symmTensorField> result(new symmTensorField(tf.size()));
    inv(result.ref(), tf, solutionD);
    return result;
}

tmp<symmTensorField> inv
(
    const tmp<symmTensorField>& tf,
    const Vector<label>& solutionD
)
{
    tmp<symmTensorField> tRes = New(tf);
    inv(tRes.ref(), tf(), solutionD);
    tf.clear();
    return tRes;
}


template<>
tmp<Field<symmTensor>> transformFieldMask<symmTensor>
(
    const tensorField& tf
)
{
    return symm(tf);
}

template<>
tmp<Field<symmTensor>> transformFieldMask<symmTensor>
(
    const tmp<tensorField>& ttf
)
{
    tmp<Field<symmTensor>> ret = transformFieldMask<symmTensor>(ttf());
    ttf.clear();
    return ret;
}


template<>
tmp<Field<symmTensor>> transformFieldMask<symmTensor>
(
    const symmTensorField& stf
)
{
    return stf;
}

template<>
tmp<Field<symmTensor>> transformFieldMask<symmTensor>
(
    const tmp<symmTensorField>& tstf
)
{
    return tstf;
}


// * * * * * * * * * * * * * * * global operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(vector, symmTensor, *, hdual)

BINARY_OPERATOR(tensor, symmTensor, symmTensor, &, dot)
BINARY_TYPE_OPERATOR(tensor, symmTensor, symmTensor, &, dot)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
