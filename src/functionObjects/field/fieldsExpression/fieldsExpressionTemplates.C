/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2026 OpenFOAM Foundation
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

#include "fieldsExpression.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<template<class> class GeoField, class Op, class A, class B, class>
bool Foam::functionObjects::fieldsExpression::opAndStore
(
    const GeoField<A>& a,
    const GeoField<B>& b
)
{
    store(resultName_, eval(Op::value(a, b)));

    return true;
}


template<template<class> class GeoField, class Op, class ... Args>
bool Foam::functionObjects::fieldsExpression::opAndStore(const Args& ...)
{
    return false;
}


template<template<class> class GeoField, class Op, class A, class B>
bool Foam::functionObjects::fieldsExpression::foldAB(const label i)
{
    if
    (
        i == 0
     && foundObject<GeoField<A>>(fieldNames_[0])
    )
    {
        clearObject(resultName_);
        store
        (
            resultName_,
            lookupObject<GeoField<A>>(fieldNames_[0]).clone()
        );

        return true;
    }

    if
    (
        i > 0
     && foundObject<GeoField<A>>(resultName_)
     && foundObject<GeoField<B>>(fieldNames_[i])
    )
    {
        tmp<GeoField<A>> a = lookupObject<GeoField<A>>(resultName_).clone();
        const GeoField<B>& b = lookupObject<GeoField<B>>(fieldNames_[i]);

        clearObject(resultName_);
        return opAndStore<GeoField, Op>(a(), b);
    }

    return false;
}


template<template<class> class GeoField, class Op, class A>
bool Foam::functionObjects::fieldsExpression::foldA(const label i)
{
    bool success = false;

    #define processType(Type, none) \
        success = success || foldAB<GeoField, Op, A, Type>(i);
    FOR_ALL_FIELD_TYPES(processType);
    #undef processType

    return success;
}


template<template<class> class GeoField, class Op>
bool Foam::functionObjects::fieldsExpression::fold(const label i)
{
    bool success = false;

    #define processType(Type, none) \
        success = success || foldA<GeoField, Op, Type>(i);
    FOR_ALL_FIELD_TYPES(processType);
    #undef processType

    return success;
}


template<template<class> class GeoField, class Op>
bool Foam::functionObjects::fieldsExpression::calcGeoFieldOp()
{
    forAll(fieldNames_, i)
    {
        if (!fold<GeoField, Op>(i))
        {
            return false;
        }
    }

    return true;
}


template<class Op>
bool Foam::functionObjects::fieldsExpression::calcOp()
{
    return
        calcGeoFieldOp<VolField, Op>()
     || calcGeoFieldOp<VolInternalField, Op>()
     || calcGeoFieldOp<SurfaceField, Op>();
}


// ************************************************************************* //
