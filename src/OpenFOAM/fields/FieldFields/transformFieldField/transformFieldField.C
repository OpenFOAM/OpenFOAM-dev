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

Description
    Spatial transformation functions for FieldField.

\*---------------------------------------------------------------------------*/

#include "transformFieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

template<template<class> class Field, class Type>
void transform
(
    FieldField<Field, Type>& rtf,
    const FieldField<Field, tensor>& trf,
    const FieldField<Field, Type>& tf
)
{
    forAll(rtf, i)
    {
        transform(rtf[i], trf[i], tf[i]);
    }
}


template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> transform
(
    const FieldField<Field, tensor>& trf,
    const FieldField<Field, Type>& tf
)
{
    tmp<FieldField<Field, Type>> tranf
    (
        FieldField<Field, Type>::NewCalculatedType(tf)
    );
    transform(tranf(), trf, tf);
    return tranf;
}


template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> transform
(
    const FieldField<Field, tensor>& trf,
    const tmp<FieldField<Field, Type>>& ttf
)
{
    tmp<FieldField<Field, Type>> tranf(ttf.ptr());
    transform(tranf(), trf, tranf());
    return tranf;
}


template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> transform
(
    const tmp<FieldField<Field, tensor>>& ttrf,
    const FieldField<Field, Type>& tf
)
{
    tmp<FieldField<Field, Type>> tranf
    (
        FieldField<Field, Type>::NewCalculatedType(tf)
    );
    transform(tranf(), ttrf(), tf);
    ttrf.clear();
    return tranf;
}


template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> transform
(
    const tmp<FieldField<Field, tensor>>& ttrf,
    const tmp<FieldField<Field, Type>>& ttf
)
{
    tmp<FieldField<Field, Type>> tranf(ttf.ptr());
    transform(tranf(), ttrf(), tranf());
    ttrf.clear();
    return tranf;
}


template<template<class> class Field, class Type>
void transform
(
    FieldField<Field, Type>& rtf,
    const tensor& t,
    const FieldField<Field, Type>& tf
)
{
    forAll(rtf, i)
    {
        transform(rtf[i], t, tf[i]);
    }
}


template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> transform
(
    const tensor& t,
    const FieldField<Field, Type>& tf
)
{
    tmp<FieldField<Field, Type>> tranf
    (
        FieldField<Field, Type>::NewCalculatedType(tf)
    );
    transform(tranf(), t, tf);
    return tranf;
}


template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> transform
(
    const tensor& t,
    const tmp<FieldField<Field, Type>>& ttf
)
{
    tmp<FieldField<Field, Type>> tranf(ttf.ptr());
    transform(tranf(), t, tranf());
    return tranf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
