/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "DimensionedField.H"
#include "dimensionedType.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define checkFieldAssignment(df1, df2)                                         \
                                                                               \
    if                                                                         \
    (                                                                          \
        static_cast<const regIOobject*>(&df1)                                  \
     == static_cast<const regIOobject*>(&df2)                                  \
    )                                                                          \
    {                                                                          \
        FatalErrorInFunction                                                   \
            << "attempted assignment to self for field "                       \
            << (df1).name() << abort(FatalError);                              \
    }


#define checkFieldOperation(df1, df2, op)                                      \
                                                                               \
    if (&(df1).mesh() != &(df2).mesh())                                        \
    {                                                                          \
        FatalErrorInFunction                                                   \
            << "different mesh for fields "                                    \
            << (df1).name() << " and " << (df2).name()                         \
            << " during operation " <<  op                                     \
            << abort(FatalError);                                              \
    }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class GeoMesh, template<class> class PrimitiveField>
DimensionedField<Type, GeoMesh, PrimitiveField>::DimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& dims,
    const PrimitiveField<Type>& field
)
:
    regIOobject(io),
    PrimitiveField<Type>(field),
    OldTimeField<DimensionedField>(this->time().timeIndex()),
    mesh_(mesh),
    dimensions_(dims)
{
    if (field.size() && field.size() != GeoMesh::size(mesh))
    {
        FatalErrorInFunction
            << "size of field = " << field.size()
            << " is not the same as the size of mesh = "
            << GeoMesh::size(mesh)
            << abort(FatalError);
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
DimensionedField<Type, GeoMesh, PrimitiveField>::DimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& dims,
    const bool checkIOFlags
)
:
    regIOobject(io),
    PrimitiveField<Type>(GeoMesh::size(mesh)),
    OldTimeField<DimensionedField>(this->time().timeIndex()),
    mesh_(mesh),
    dimensions_(dims)
{
    // Expand dynamic primitive fields to their full size
    PrimitiveField<Type>::setSize(GeoMesh::size(mesh));

    if (checkIOFlags)
    {
        readIfPresent();
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
DimensionedField<Type, GeoMesh, PrimitiveField>::DimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensioned<Type>& dt,
    const bool checkIOFlags
)
:
    regIOobject(io),
    PrimitiveField<Type>(GeoMesh::size(mesh), dt.value()),
    OldTimeField<DimensionedField>(this->time().timeIndex()),
    mesh_(mesh),
    dimensions_(dt.dimensions())
{
    if (checkIOFlags)
    {
        readIfPresent();
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
DimensionedField<Type, GeoMesh, PrimitiveField>::DimensionedField
(
    const DimensionedField<Type, GeoMesh, PrimitiveField>& df
)
:
    regIOobject(df),
    PrimitiveField<Type>(df),
    OldTimeField<DimensionedField>(df),
    mesh_(df.mesh_),
    dimensions_(df.dimensions_)
{}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
DimensionedField<Type, GeoMesh, PrimitiveField>::DimensionedField
(
    DimensionedField<Type, GeoMesh, PrimitiveField>&& df
)
:
    regIOobject(move(df)),
    PrimitiveField<Type>(move(df)),
    OldTimeField<DimensionedField>(move(df)),
    mesh_(df.mesh_),
    dimensions_(move(df.dimensions_))
{}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
DimensionedField<Type, GeoMesh, PrimitiveField>::DimensionedField
(
    const DimensionedField<Type, GeoMesh, PrimitiveField2>& df
)
:
    regIOobject(df),
    PrimitiveField<Type>(df),
    OldTimeField<DimensionedField>(this->time().timeIndex()),
    mesh_(df.mesh()),
    dimensions_(df.dimensions())
{}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
DimensionedField<Type, GeoMesh, PrimitiveField>::DimensionedField
(
    DimensionedField<Type, GeoMesh, PrimitiveField2>& df,
    const bool reuse
)
:
    regIOobject(df, reuse),
    PrimitiveField<Type>(df, reuse),
    OldTimeField<DimensionedField>(this->time().timeIndex()),
    mesh_(df.mesh()),
    dimensions_(df.dimensions())
{}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
DimensionedField<Type, GeoMesh, PrimitiveField>::DimensionedField
(
    const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>& tdf
)
:
    regIOobject(tdf(), tdf.isTmp()),
    PrimitiveField<Type>
    (
        const_cast<DimensionedField<Type, GeoMesh, PrimitiveField>&>(tdf()),
        tdf.isTmp()
    ),
    OldTimeField<DimensionedField>(this->time().timeIndex()),
    mesh_(tdf().mesh_),
    dimensions_(tdf().dimensions_)
{
    tdf.clear();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
DimensionedField<Type, GeoMesh, PrimitiveField>::DimensionedField
(
    const IOobject& io,
    const DimensionedField<Type, GeoMesh, PrimitiveField2>& df,
    const bool checkIOFlags
)
:
    regIOobject(io),
    PrimitiveField<Type>(df),
    OldTimeField<DimensionedField>(this->time().timeIndex()),
    mesh_(df.mesh_),
    dimensions_(df.dimensions_)
{
    if (!checkIOFlags || !readIfPresent())
    {
        copyOldTimes(io, df);
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
DimensionedField<Type, GeoMesh, PrimitiveField>::DimensionedField
(
    const IOobject& io,
    DimensionedField<Type, GeoMesh, PrimitiveField2>& df,
    const bool reuse,
    const bool checkIOFlags
)
:
    regIOobject(io, df),
    PrimitiveField<Type>(df, reuse),
    OldTimeField<DimensionedField>(this->time().timeIndex()),
    mesh_(df.mesh_),
    dimensions_(df.dimensions_)
{
    if (checkIOFlags)
    {
        readIfPresent();
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
DimensionedField<Type, GeoMesh, PrimitiveField>::DimensionedField
(
    const IOobject& io,
    const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>& tdf,
    const bool checkIOFlags
)
:
    regIOobject(io),
    PrimitiveField<Type>
    (
        const_cast<DimensionedField<Type, GeoMesh, PrimitiveField>&>(tdf()),
        tdf.isTmp()
    ),
    OldTimeField<DimensionedField>(this->time().timeIndex()),
    mesh_(tdf().mesh_),
    dimensions_(tdf().dimensions_)
{
    tdf.clear();

    if (checkIOFlags)
    {
        readIfPresent();
    }
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
DimensionedField<Type, GeoMesh, PrimitiveField>::DimensionedField
(
    const word& newName,
    const DimensionedField<Type, GeoMesh, PrimitiveField2>& df
)
:
    regIOobject(newName, df, newName != df.name()),
    PrimitiveField<Type>(df),
    OldTimeField<DimensionedField>(this->time().timeIndex()),
    mesh_(df.mesh_),
    dimensions_(df.dimensions_)
{
    copyOldTimes(newName, df);
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
DimensionedField<Type, GeoMesh, PrimitiveField>::DimensionedField
(
    const word& newName,
    DimensionedField<Type, GeoMesh, PrimitiveField2>& df,
    const bool reuse
)
:
    regIOobject(newName, df, true),
    PrimitiveField<Type>(df, reuse),
    OldTimeField<DimensionedField>(this->time().timeIndex()),
    mesh_(df.mesh_),
    dimensions_(df.dimensions_)
{}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
DimensionedField<Type, GeoMesh, PrimitiveField>::DimensionedField
(
    const word& newName,
    const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>& tdf
)
:
    regIOobject(newName, tdf(), true),
    PrimitiveField<Type>
    (
        const_cast<DimensionedField<Type, GeoMesh, PrimitiveField>&>(tdf()),
        tdf.isTmp()
    ),
    OldTimeField<DimensionedField>(this->time().timeIndex()),
    mesh_(tdf().mesh_),
    dimensions_(tdf().dimensions_)
{
    tdf.clear();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>
DimensionedField<Type, GeoMesh, PrimitiveField>::clone() const
{
    return tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>
    (
        new DimensionedField<Type, GeoMesh, PrimitiveField>(*this)
    );
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh, PrimitiveField>>
DimensionedField<Type, GeoMesh, PrimitiveField>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensionSet& ds,
    const PrimitiveField<Type>& field
)
{
    const bool cacheTmp = mesh.thisDb().cacheTemporaryObject(name);

    return tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>
    (
        new DimensionedField<Type, GeoMesh, PrimitiveField>
        (
            IOobject
            (
                name,
                mesh.thisDb().time().name(),
                mesh.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                cacheTmp
            ),
            mesh,
            ds,
            field
        ),
        cacheTmp
    );
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh, PrimitiveField>>
DimensionedField<Type, GeoMesh, PrimitiveField>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensionSet& ds
)
{
    const bool cacheTmp = mesh.thisDb().cacheTemporaryObject(name);

    return tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>
    (
        new DimensionedField<Type, GeoMesh, PrimitiveField>
        (
            IOobject
            (
                name,
                mesh.thisDb().time().name(),
                mesh.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                cacheTmp
            ),
            mesh,
            ds,
            false
        ),
        cacheTmp
    );
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh, PrimitiveField>>
DimensionedField<Type, GeoMesh, PrimitiveField>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensioned<Type>& dt
)
{
    const bool cacheTmp = mesh.thisDb().cacheTemporaryObject(name);

    return tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>
    (
        new DimensionedField<Type, GeoMesh, PrimitiveField>
        (
            IOobject
            (
                name,
                mesh.thisDb().time().name(),
                mesh.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                cacheTmp
            ),
            mesh,
            dt,
            false
        ),
        cacheTmp
    );
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh, PrimitiveField>>
DimensionedField<Type, GeoMesh, PrimitiveField>::New
(
    const word& newName,
    const DimensionedField<Type, GeoMesh, PrimitiveField2>& df
)
{
    const bool cacheTmp = df.db().cacheTemporaryObject(newName);

    return tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>
    (
        new DimensionedField<Type, GeoMesh, PrimitiveField>
        (
            IOobject
            (
                newName,
                df.instance(),
                df.local(),
                df.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                cacheTmp
            ),
            df
        ),
        cacheTmp
    );
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh, PrimitiveField>>
DimensionedField<Type, GeoMesh, PrimitiveField>::New
(
    const word& newName,
    const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>& tdf
)
{
    const bool cacheTmp = tdf().db().cacheTemporaryObject(newName);

    return tmp<DimensionedField<Type, GeoMesh, Field>>
    (
        new DimensionedField<Type, GeoMesh, Field>
        (
            IOobject
            (
                newName,
                tdf().instance(),
                tdf().local(),
                tdf().db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                cacheTmp
            ),
            tdf
        ),
        cacheTmp
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, class GeoMesh, template<class> class PrimitiveField>
DimensionedField<Type, GeoMesh, PrimitiveField>::~DimensionedField()
{
    db().cacheTemporaryObject(*this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMesh, template<class> class PrimitiveField>
PrimitiveField<Type>&
Foam::DimensionedField<Type, GeoMesh, PrimitiveField>::primitiveFieldRef()
{
    this->setUpToDate();
    storeOldTimes();
    return *this;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
tmp
<
    DimensionedField
    <
        typename DimensionedField<Type, GeoMesh, PrimitiveField>::cmptType,
        GeoMesh,
        Field
    >
>
DimensionedField<Type, GeoMesh, PrimitiveField>::component
(
    const direction d
) const
{
    tmp<DimensionedField<cmptType, GeoMesh, Field>> result
    (
        DimensionedField<cmptType, GeoMesh, Field>::New
        (
            name() + ".component(" + ::Foam::name(d) + ')',
            mesh_,
            dimensions_
        )
    );

    Foam::component(result.ref(), *this, d);

    return result;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
void DimensionedField<Type, GeoMesh, PrimitiveField>::replace
(
    const direction d,
    const DimensionedField
    <
        typename DimensionedField<Type, GeoMesh, PrimitiveField>::cmptType,
        GeoMesh,
        PrimitiveField2
    >& df
)
{
    PrimitiveField<Type>::replace(d, df);
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
void DimensionedField<Type, GeoMesh, PrimitiveField>::replace
(
    const direction d,
    const tmp
    <
        DimensionedField
        <
            typename DimensionedField<Type, GeoMesh, PrimitiveField>::cmptType,
            GeoMesh,
            PrimitiveField2
        >
    >& tdf
)
{
    replace(d, tdf());
    tdf.clear();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<Type, GeoMesh, Field>>
DimensionedField<Type, GeoMesh, PrimitiveField>::T() const
{
    tmp<DimensionedField<Type, GeoMesh, Field>> result
    (
        DimensionedField<Type, GeoMesh, Field>::New
        (
            name() + ".T()",
            mesh_,
            dimensions_
        )
    );

    Foam::T(result.ref(), *this);

    return result;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
dimensioned<Type>
DimensionedField<Type, GeoMesh, PrimitiveField>::average() const
{
    dimensioned<Type> Average
    (
        this->name() + ".average()",
        this->dimensions(),
        gAverage(primitiveField())
    );

    return Average;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
dimensioned<Type>
DimensionedField<Type, GeoMesh, PrimitiveField>::weightedAverage
(
    const DimensionedField<scalar, GeoMesh, PrimitiveField2>& weightField
) const
{
    return
    (
        dimensioned<Type>
        (
            this->name() + ".weightedAverage(weights)",
            this->dimensions(),
            gSum(weightField*primitiveField())/gSum(weightField)
        )
    );
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
dimensioned<Type>
DimensionedField<Type, GeoMesh, PrimitiveField>::weightedAverage
(
    const tmp<DimensionedField<scalar, GeoMesh, PrimitiveField2>>& tweightField
) const
{
    dimensioned<Type> wa = weightedAverage(tweightField());
    tweightField.clear();
    return wa;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
void DimensionedField<Type, GeoMesh, PrimitiveField>::reset
(
    const DimensionedField<Type, GeoMesh, PrimitiveField2>& df
)
{
    checkFieldAssignment(*this, df);

    dimensions_ = df.dimensions();
    PrimitiveField<Type>::operator=(df.primitiveField());
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void DimensionedField<Type, GeoMesh, PrimitiveField>::reset
(
    const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>& tdf
)
{
    const DimensionedField<Type, GeoMesh, PrimitiveField>& df = tdf();

    checkFieldAssignment(*this, df);

    dimensions_ = df.dimensions();

    if (tdf.isTmp())
    {
        PrimitiveField<Type>::transfer(tdf.ref());
    }
    else
    {
        PrimitiveField<Type>::operator=(df.primitiveField());
    }

    tdf.clear();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
void DimensionedField<Type, GeoMesh, PrimitiveField>::reset
(
    const tmp<DimensionedField<Type, GeoMesh, PrimitiveField2>>& tdf
)
{
    reset(tdf());

    tdf.clear();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type, class GeoMesh, template<class> class PrimitiveField>
void DimensionedField<Type, GeoMesh, PrimitiveField>::operator=
(
    const DimensionedField<Type, GeoMesh, PrimitiveField>& df
)
{
    checkFieldAssignment(*this, df);
    checkFieldOperation(*this, df, "=");

    dimensions_ = df.dimensions();
    PrimitiveField<Type>::operator=(df);
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void DimensionedField<Type, GeoMesh, PrimitiveField>::operator=
(
    DimensionedField<Type, GeoMesh, PrimitiveField>&& df
)
{
    checkFieldAssignment(*this, df);
    checkFieldOperation(*this, df, "=");

    dimensions_ = move(df.dimensions());
    PrimitiveField<Type>::operator=(move(df));
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void DimensionedField<Type, GeoMesh, PrimitiveField>::operator=
(
    const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>& tdf
)
{
    const DimensionedField<Type, GeoMesh, PrimitiveField>& df = tdf();

    checkFieldAssignment(*this, df);
    checkFieldOperation(*this, df, "=");

    dimensions_ = df.dimensions();

    if (tdf.isTmp())
    {
        primitiveFieldRef().transfer(tdf.ref());
    }
    else
    {
        primitiveFieldRef() = df.primitiveField();
    }

    tdf.clear();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
void DimensionedField<Type, GeoMesh, PrimitiveField>::operator=
(
    const DimensionedField<Type, GeoMesh, PrimitiveField2>& df
)
{
    checkFieldOperation(*this, df, "=");

    dimensions_ = df.dimensions();
    PrimitiveField<Type>::operator=(df);
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
void DimensionedField<Type, GeoMesh, PrimitiveField>::operator=
(
    const tmp<DimensionedField<Type, GeoMesh, PrimitiveField2>>& tdf
)
{
    const DimensionedField<Type, GeoMesh, PrimitiveField2>& df = tdf();

    checkFieldOperation(*this, df, "=");

    dimensions_ = df.dimensions();
    PrimitiveField<Type>::operator=(df);
    tdf.clear();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void DimensionedField<Type, GeoMesh, PrimitiveField>::operator=
(
    const dimensioned<Type>& dt
)
{
    dimensions_ = dt.dimensions();
    PrimitiveField<Type>::operator=(dt.value());
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void DimensionedField<Type, GeoMesh, PrimitiveField>::operator=(const zero&)
{
    PrimitiveField<Type>::operator=(Zero);
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
void DimensionedField<Type, GeoMesh, PrimitiveField>::operator==
(
    const DimensionedField<Type, GeoMesh, PrimitiveField2>& df
)
{
    this->operator=(df);
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
template<template<class> class PrimitiveField2>
void DimensionedField<Type, GeoMesh, PrimitiveField>::operator==
(
    const tmp<DimensionedField<Type, GeoMesh, PrimitiveField2>>& tdf
)
{
    this->operator=(tdf);
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void DimensionedField<Type, GeoMesh, PrimitiveField>::operator==
(
    const dimensioned<Type>& dt
)
{
    this->operator=(dt);
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
void DimensionedField<Type, GeoMesh, PrimitiveField>::operator==(const zero&)
{
    this->operator=(Zero);
}


#define COMPUTED_ASSIGNMENT(TYPE, op)                                          \
                                                                               \
template<class Type, class GeoMesh, template<class> class PrimitiveField>      \
template<template<class> class PrimitiveField2>                                \
void DimensionedField<Type, GeoMesh, PrimitiveField>::operator op              \
(                                                                              \
    const DimensionedField<TYPE, GeoMesh, PrimitiveField2>& df                 \
)                                                                              \
{                                                                              \
    checkFieldOperation(*this, df, #op);                                       \
                                                                               \
    dimensions_ op df.dimensions();                                            \
    PrimitiveField<Type>::operator op(df);                                     \
}                                                                              \
                                                                               \
template<class Type, class GeoMesh, template<class> class PrimitiveField>      \
template<template<class> class PrimitiveField2>                                \
void DimensionedField<Type, GeoMesh, PrimitiveField>::operator op              \
(                                                                              \
    const tmp<DimensionedField<TYPE, GeoMesh, PrimitiveField2>>& tdf           \
)                                                                              \
{                                                                              \
    operator op(tdf());                                                        \
    tdf.clear();                                                               \
}                                                                              \
                                                                               \
template<class Type, class GeoMesh, template<class> class PrimitiveField>      \
void DimensionedField<Type, GeoMesh, PrimitiveField>::operator op              \
(                                                                              \
    const dimensioned<TYPE>& dt                                                \
)                                                                              \
{                                                                              \
    dimensions_ op dt.dimensions();                                            \
    PrimitiveField<Type>::operator op(dt.value());                             \
}

COMPUTED_ASSIGNMENT(Type, +=)
COMPUTED_ASSIGNMENT(Type, -=)
COMPUTED_ASSIGNMENT(scalar, *=)
COMPUTED_ASSIGNMENT(scalar, /=)

#undef COMPUTED_ASSIGNMENT


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef checkFieldAssignment
#undef checkFieldOperation

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DimensionedFieldIO.C"
#include "DimensionedFieldFunctions.C"

// ************************************************************************* //
