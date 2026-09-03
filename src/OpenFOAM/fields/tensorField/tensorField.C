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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

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

        forAll(tf, i)
        {
            tf[i] = inv(tf1Plus[i]);
        }

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
        forAll(tf, i)
        {
            tf[i] = inv(tf1[i]);
        }
    }
}

tmp<tensorField> inv
(
    const Field<tensor>& tf,
    const Vector<label>& solutionD
)
{
    tmp<tensorField> result(new tensorField(tf.size()));
    inv(result.ref(), tf, solutionD);
    return result;
}

tmp<tensorField> inv
(
    const SubField<tensor>& tf,
    const Vector<label>& solutionD
)
{
    tmp<tensorField> result(new tensorField(tf.size()));
    inv(result.ref(), tf, solutionD);
    return result;
}

tmp<tensorField> inv
(
    const tmp<tensorField>& tf,
    const Vector<label>& solutionD
)
{
    tmp<tensorField> tRes = tensorField::New(tf, false);
    inv(tRes.ref(), tf(), solutionD);
    tf.clear();
    return tRes;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
