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

\*---------------------------------------------------------------------------*/

#include "transformField.H"
#include "FieldM.H"
#include "diagTensor.H"

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

void Foam::transform
(
    vectorField& rtf,
    const quaternion& q,
    const vectorField& tf
)
{
    tensor t = q.R();
    TFOR_ALL_F_OP_FUNC_S_F(vector, rtf, =, transform, tensor, t, vector, tf)
}


Foam::tmp<Foam::vectorField> Foam::transform
(
    const quaternion& q,
    const vectorField& tf
)
{
    tmp<vectorField > tranf(new vectorField(tf.size()));
    transform(tranf.ref(), q, tf);
    return tranf;
}


Foam::tmp<Foam::vectorField> Foam::transform
(
    const quaternion& q,
    const tmp<vectorField>& ttf
)
{
    tmp<vectorField > tranf = New(ttf);
    transform(tranf.ref(), q, ttf());
    ttf.clear();
    return tranf;
}


void Foam::transformPoints
(
    vectorField& rtf,
    const septernion& tr,
    const vectorField& tf
)
{
    vector T = tr.t();

    // Check if any translation
    if (mag(T) > vSmall)
    {
        TFOR_ALL_F_OP_F_OP_S(vector, rtf, =, vector, tf, -, vector, T);
    }
    else
    {
        rtf = tf;
    }

    // Check if any rotation
    if (mag(tr.r().R() - I) > small)
    {
        transform(rtf, tr.r(), rtf);
    }
}


Foam::tmp<Foam::vectorField> Foam::transformPoints
(
    const septernion& tr,
    const vectorField& tf
)
{
    tmp<vectorField > tranf(new vectorField(tf.size()));
    transformPoints(tranf.ref(), tr, tf);
    return tranf;
}


Foam::tmp<Foam::vectorField> Foam::transformPoints
(
    const septernion& tr,
    const tmp<vectorField>& ttf
)
{
    tmp<vectorField > tranf = New(ttf);
    transformPoints(tranf.ref(), tr, ttf());
    ttf.clear();
    return tranf;
}


// ************************************************************************* //
