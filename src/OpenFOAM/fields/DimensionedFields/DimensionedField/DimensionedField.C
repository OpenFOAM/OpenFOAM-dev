/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

// check mesh for two fields
#define checkField(df1, df2, op)                                    \
if (&(df1).mesh() != &(df2).mesh())                                 \
{                                                                   \
    FatalErrorInFunction                                            \
        << "different mesh for fields "                             \
        << (df1).name() << " and " << (df2).name()                  \
        << " during operatrion " <<  op                             \
        << abort(FatalError);                                       \
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
DimensionedField<Type, GeoMesh>::DimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& dims,
    const Field<Type>& field
)
:
    regIOobject(io),
    Field<Type>(field),
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


template<class Type, class GeoMesh>
DimensionedField<Type, GeoMesh>::DimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& dims,
    const bool checkIOFlags
)
:
    regIOobject(io),
    Field<Type>(GeoMesh::size(mesh)),
    mesh_(mesh),
    dimensions_(dims)
{
    if (checkIOFlags)
    {
        readIfPresent();
    }
}


template<class Type, class GeoMesh>
DimensionedField<Type, GeoMesh>::DimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensioned<Type>& dt,
    const bool checkIOFlags
)
:
    regIOobject(io),
    Field<Type>(GeoMesh::size(mesh), dt.value()),
    mesh_(mesh),
    dimensions_(dt.dimensions())
{
    if (checkIOFlags)
    {
        readIfPresent();
    }
}


template<class Type, class GeoMesh>
DimensionedField<Type, GeoMesh>::DimensionedField
(
    const DimensionedField<Type, GeoMesh>& df
)
:
    regIOobject(df),
    Field<Type>(df),
    mesh_(df.mesh_),
    dimensions_(df.dimensions_)
{}


template<class Type, class GeoMesh>
DimensionedField<Type, GeoMesh>::DimensionedField
(
    DimensionedField<Type, GeoMesh>& df,
    bool reuse
)
:
    regIOobject(df, reuse),
    Field<Type>(df, reuse),
    mesh_(df.mesh_),
    dimensions_(df.dimensions_)
{}


template<class Type, class GeoMesh>
DimensionedField<Type, GeoMesh>::DimensionedField
(
    DimensionedField<Type, GeoMesh>&& df
)
:
    regIOobject(move(df), true),
    Field<Type>(move(df)),
    mesh_(df.mesh_),
    dimensions_(move(df.dimensions_))
{}


template<class Type, class GeoMesh>
DimensionedField<Type, GeoMesh>::DimensionedField
(
    const tmp<DimensionedField<Type, GeoMesh>>& tdf
)
:
    regIOobject(tdf(), tdf.isTmp()),
    Field<Type>
    (
        const_cast<DimensionedField<Type, GeoMesh>&>(tdf()),
        tdf.isTmp()
    ),
    mesh_(tdf().mesh_),
    dimensions_(tdf().dimensions_)
{
    tdf.clear();
}


template<class Type, class GeoMesh>
DimensionedField<Type, GeoMesh>::DimensionedField
(
    const IOobject& io,
    const DimensionedField<Type, GeoMesh>& df
)
:
    regIOobject(io),
    Field<Type>(df),
    mesh_(df.mesh_),
    dimensions_(df.dimensions_)
{}


template<class Type, class GeoMesh>
DimensionedField<Type, GeoMesh>::DimensionedField
(
    const IOobject& io,
    DimensionedField<Type, GeoMesh>& df,
    bool reuse
)
:
    regIOobject(io, df),
    Field<Type>(df, reuse),
    mesh_(df.mesh_),
    dimensions_(df.dimensions_)
{}


template<class Type, class GeoMesh>
DimensionedField<Type, GeoMesh>::DimensionedField
(
    const word& newName,
    const DimensionedField<Type, GeoMesh>& df
)
:
    regIOobject(newName, df, newName != df.name()),
    Field<Type>(df),
    mesh_(df.mesh_),
    dimensions_(df.dimensions_)
{}


template<class Type, class GeoMesh>
DimensionedField<Type, GeoMesh>::DimensionedField
(
    const word& newName,
    DimensionedField<Type, GeoMesh>& df,
    bool reuse
)
:
    regIOobject(newName, df, true),
    Field<Type>(df, reuse),
    mesh_(df.mesh_),
    dimensions_(df.dimensions_)
{}


template<class Type, class GeoMesh>
DimensionedField<Type, GeoMesh>::DimensionedField
(
    const word& newName,
    const tmp<DimensionedField<Type, GeoMesh>>& tdf
)
:
    regIOobject(newName, tdf(), true),
    Field<Type>
    (
        const_cast<DimensionedField<Type, GeoMesh>&>(tdf()),
        tdf.isTmp()
    ),
    mesh_(tdf().mesh_),
    dimensions_(tdf().dimensions_)
{
    tdf.clear();
}


template<class Type, class GeoMesh>
tmp<DimensionedField<Type, GeoMesh>>
DimensionedField<Type, GeoMesh>::clone() const
{
    return tmp<DimensionedField<Type, GeoMesh>>
    (
        new DimensionedField<Type, GeoMesh>(*this)
    );
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
DimensionedField<Type, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensionSet& ds
)
{
    const bool cacheTmp = mesh.thisDb().cacheTemporaryObject(name);

    return tmp<DimensionedField<Type, GeoMesh>>
    (
        new DimensionedField<Type, GeoMesh>
        (
            IOobject
            (
                name,
                mesh.thisDb().time().timeName(),
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


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
DimensionedField<Type, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensioned<Type>& dt
)
{
    const bool cacheTmp = mesh.thisDb().cacheTemporaryObject(name);

    return tmp<DimensionedField<Type, GeoMesh>>
    (
        new DimensionedField<Type, GeoMesh>
        (
            IOobject
            (
                name,
                mesh.thisDb().time().timeName(),
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


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
DimensionedField<Type, GeoMesh>::New
(
    const word& newName,
    const DimensionedField<Type, GeoMesh>& df
)
{
    const bool cacheTmp = df.db().cacheTemporaryObject(newName);

    return tmp<DimensionedField<Type, GeoMesh>>
    (
        new DimensionedField<Type, GeoMesh>
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


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
DimensionedField<Type, GeoMesh>::New
(
    const word& newName,
    const tmp<DimensionedField<Type, GeoMesh>>& tdf
)
{
    const bool cacheTmp = tdf().db().cacheTemporaryObject(newName);

    return tmp<DimensionedField<Type, GeoMesh>>
    (
        new DimensionedField<Type, GeoMesh>
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

template<class Type, class GeoMesh>
DimensionedField<Type, GeoMesh>::~DimensionedField()
{
    db().cacheTemporaryObject(*this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
tmp
<
    DimensionedField
        <typename DimensionedField<Type, GeoMesh>::cmptType, GeoMesh>
>
DimensionedField<Type, GeoMesh>::component
(
    const direction d
) const
{
    tmp<DimensionedField<cmptType, GeoMesh>> result
    (
        DimensionedField<cmptType, GeoMesh>::New
        (
            name() + ".component(" + ::Foam::name(d) + ')',
            mesh_,
            dimensions_
        )
    );

    Foam::component(result.ref(), *this, d);

    return result;
}


template<class Type, class GeoMesh>
void DimensionedField<Type, GeoMesh>::replace
(
    const direction d,
    const DimensionedField
        <typename DimensionedField<Type, GeoMesh>::cmptType, GeoMesh>& df
)
{
    Field<Type>::replace(d, df);
}


template<class Type, class GeoMesh>
void DimensionedField<Type, GeoMesh>::replace
(
    const direction d,
    const tmp
    <
        DimensionedField
            <typename DimensionedField<Type, GeoMesh>::cmptType, GeoMesh>
    >& tdf
)
{
    replace(d, tdf());
    tdf.clear();
}


template<class Type, class GeoMesh>
tmp<DimensionedField<Type, GeoMesh>>
DimensionedField<Type, GeoMesh>::T() const
{
    tmp<DimensionedField<Type, GeoMesh>> result
    (
        DimensionedField<Type, GeoMesh>::New
        (
            name() + ".T()",
            mesh_,
            dimensions_
        )
    );

    Foam::T(result.ref(), *this);

    return result;
}


template<class Type, class GeoMesh>
dimensioned<Type> DimensionedField<Type, GeoMesh>::average() const
{
    dimensioned<Type> Average
    (
        this->name() + ".average()",
        this->dimensions(),
        gAverage(field())
    );

    return Average;
}


template<class Type, class GeoMesh>
dimensioned<Type> DimensionedField<Type, GeoMesh>::weightedAverage
(
    const DimensionedField<scalar, GeoMesh>& weightField
) const
{
    return
    (
        dimensioned<Type>
        (
            this->name() + ".weightedAverage(weights)",
            this->dimensions(),
            gSum(weightField*field())/gSum(weightField)
        )
    );
}


template<class Type, class GeoMesh>
dimensioned<Type> DimensionedField<Type, GeoMesh>::weightedAverage
(
    const tmp<DimensionedField<scalar, GeoMesh>>& tweightField
) const
{
    dimensioned<Type> wa = weightedAverage(tweightField());
    tweightField.clear();
    return wa;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
void DimensionedField<Type, GeoMesh>::operator=
(
    const DimensionedField<Type, GeoMesh>& df
)
{
    // Check for assignment to self
    if (this == &df)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    checkField(*this, df, "=");

    dimensions_ = df.dimensions();
    Field<Type>::operator=(df);
}


template<class Type, class GeoMesh>
void DimensionedField<Type, GeoMesh>::operator=
(
    DimensionedField<Type, GeoMesh>&& df
)
{
    // Check for assignment to self
    if (this == &df)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    checkField(*this, df, "=");

    dimensions_ = move(df.dimensions());
    Field<Type>::operator=(move(df));
}


template<class Type, class GeoMesh>
void DimensionedField<Type, GeoMesh>::operator=
(
    const tmp<DimensionedField<Type, GeoMesh>>& tdf
)
{
    const DimensionedField<Type, GeoMesh>& df = tdf();

    // Check for assignment to self
    if (this == &df)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    checkField(*this, df, "=");

    dimensions_ = df.dimensions();
    this->transfer(const_cast<DimensionedField<Type, GeoMesh>&>(df));
    tdf.clear();
}


template<class Type, class GeoMesh>
void DimensionedField<Type, GeoMesh>::operator=
(
    const dimensioned<Type>& dt
)
{
    dimensions_ = dt.dimensions();
    Field<Type>::operator=(dt.value());
}


template<class Type, class GeoMesh>
void DimensionedField<Type, GeoMesh>::operator=(const zero&)
{
    Field<Type>::operator=(Zero);
}


#define COMPUTED_ASSIGNMENT(TYPE, op)                                          \
                                                                               \
template<class Type, class GeoMesh>                                            \
void DimensionedField<Type, GeoMesh>::operator op                              \
(                                                                              \
    const DimensionedField<TYPE, GeoMesh>& df                                  \
)                                                                              \
{                                                                              \
    checkField(*this, df, #op);                                                \
                                                                               \
    dimensions_ op df.dimensions();                                            \
    Field<Type>::operator op(df);                                              \
}                                                                              \
                                                                               \
template<class Type, class GeoMesh>                                            \
void DimensionedField<Type, GeoMesh>::operator op                              \
(                                                                              \
    const tmp<DimensionedField<TYPE, GeoMesh>>& tdf                            \
)                                                                              \
{                                                                              \
    operator op(tdf());                                                        \
    tdf.clear();                                                               \
}                                                                              \
                                                                               \
template<class Type, class GeoMesh>                                            \
void DimensionedField<Type, GeoMesh>::operator op                              \
(                                                                              \
    const dimensioned<TYPE>& dt                                                \
)                                                                              \
{                                                                              \
    dimensions_ op dt.dimensions();                                            \
    Field<Type>::operator op(dt.value());                                      \
}

COMPUTED_ASSIGNMENT(Type, +=)
COMPUTED_ASSIGNMENT(Type, -=)
COMPUTED_ASSIGNMENT(scalar, *=)
COMPUTED_ASSIGNMENT(scalar, /=)

#undef COMPUTED_ASSIGNMENT


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef checkField

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DimensionedFieldIO.C"
#include "DimensionedFieldFunctions.C"

// ************************************************************************* //
