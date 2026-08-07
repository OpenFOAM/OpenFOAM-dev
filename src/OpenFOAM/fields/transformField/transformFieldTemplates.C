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

#include "transformField.H"
#include "FieldM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

template<class Type>
List<Type> transform(const tensor& t, const UList<Type>& list)
{
    List<Type> newList(list.size());

    forAll(list, i)
    {
        newList[i] = transform(t, list[i]);
    }

    return newList;
}


template<class Type>
void transform
(
    Field<Type>& rtf,
    const tensorField& trf,
    const Field<Type>& tf
)
{
    TFOR_ALL_F_OP_FUNC_F_F
    (
        Type, rtf, =, transform, tensor, trf, Type, tf
    )
}


template<class Type>
void transform
(
    Field<Type>& rtf,
    const tensor& t,
    const Field<Type>& tf
)
{
    TFOR_ALL_F_OP_FUNC_S_F(Type, rtf, =, transform, tensor, t, Type, tf)
}


template<class Type1, class Type2>
tmp<Field<Type1>> transformFieldMask(const Field<Type2>& f)
{
    return f;
}

template<class Type1, class Type2>
tmp<Field<Type1>> transformFieldMask(const tmp<Field<Type2>>& tf)
{
    return tmp<Field<Type1>>(tf.ptr());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
