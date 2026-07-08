/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "DimensionedFunction1.H"

// * * * * * * * * * * * * * * Private Constructors  * * * * * * * * * * * * //

template<class Type>
Foam::DimensionedFunction1<Type>::DimensionedFunction1
(
    const word& name,
    const Function1s::unitSets& units,
    const dictionary& dict
)
:
    xDimensions_(units.x.dimensions()),
    valueDimensions_(units.value.dimensions()),
    function_(Function1<Type>::New(name, units, dict))
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
template
<
    class GeoMesh,
    template<class> class ResultPrimitiveField,
    template<class> class XPrimitiveField
>
void Foam::DimensionedFunction1<Type>::value
(
    DimensionedField<scalar, GeoMesh, ResultPrimitiveField>& result,
    const DimensionedField<scalar, GeoMesh, XPrimitiveField>& x
) const
{
    x.dimensions() = xDimensions_;

    result.dimensions() = valueDimensions_;

    result.primitiveFieldRef() = function_->value(x.primitiveField());
}


template<class Type>
template
<
    class GeoMesh,
    template<class> class ResultPrimitiveField,
    template<class> class XPrimitiveField
>
void Foam::DimensionedFunction1<Type>::derivative
(
    DimensionedField<scalar, GeoMesh, ResultPrimitiveField>& result,
    const DimensionedField<scalar, GeoMesh, XPrimitiveField>& x
) const
{
    x.dimensions() = xDimensions_;

    result.dimensions() = valueDimensions_/xDimensions_;

    result.primitiveFieldRef() = function_->derivative(x.primitiveField());
}


template<class Type>
template
<
    class GeoMesh,
    template<class> class ResultPrimitiveField,
    template<class> class XPrimitiveField
>
void Foam::DimensionedFunction1<Type>::value
(
    GeometricField<scalar, GeoMesh, ResultPrimitiveField>& result,
    const GeometricField<scalar, GeoMesh, XPrimitiveField>& x
) const
{
    value(result.internalFieldRef(), x.internalField());

    typename GeometricField<Type, GeoMesh, ResultPrimitiveField>::Boundary&
        bresult = result.boundaryFieldRef();

    const typename GeometricField<scalar, GeoMesh, XPrimitiveField>::Boundary&
        bx = x.boundaryField();

    forAll(bresult, patchi)
    {
        bresult[patchi] = function_->value(bx[patchi]);
    }
}


template<class Type>
template
<
    class GeoMesh,
    template<class> class ResultPrimitiveField,
    template<class> class XPrimitiveField
>
void Foam::DimensionedFunction1<Type>::derivative
(
    GeometricField<scalar, GeoMesh, ResultPrimitiveField>& result,
    const GeometricField<scalar, GeoMesh, XPrimitiveField>& x
) const
{
    derivative(result.internalFieldRef(), x.internalField());

    typename GeometricField<Type, GeoMesh, ResultPrimitiveField>::Boundary&
        bresult = result.boundaryFieldRef();

    const typename GeometricField<scalar, GeoMesh, XPrimitiveField>::Boundary&
        bx = x.boundaryField();

    forAll(bresult, patchi)
    {
        bresult[patchi] = function_->derivative(bx[patchi]);
    }
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::DimensionedFunction1<Type>>
Foam::DimensionedFunction1<Type>::New
(
    const word& name,
    const Function1s::unitSets& units,
    const dictionary& dict
)
{
    return
        autoPtr<DimensionedFunction1<Type>>
        (
            new DimensionedFunction1(name, units, dict)
        );
}


template<class Type>
Foam::autoPtr<Foam::DimensionedFunction1<Type>>
Foam::DimensionedFunction1<Type>::New
(
    const word& name,
    const unitSet& xUnits,
    const unitSet& valueUnits,
    const dictionary& dict
)
{
    return
        autoPtr<DimensionedFunction1<Type>>
        (
            new DimensionedFunction1(name, {xUnits, valueUnits}, dict)
        );
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::DimensionedFunction1<Type>::value
(
    const scalarField& x
) const
{
    return function_->value(x);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::DimensionedFunction1<Type>::derivative
(
    const scalarField& x
) const
{
    return function_->derivative(x);
}


template<class Type>
template<class GeoMesh, template<class> class XPrimitiveField>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedFunction1<Type>::value
(
    const DimensionedField<scalar, GeoMesh, XPrimitiveField>& x
) const
{
    tmp<DimensionedField<Type, GeoMesh>> tresult
    (
        DimensionedField<Type, GeoMesh>::New
        (
            function_->name() + "(" + x.name() + ')',
            x.mesh(),
            valueDimensions_
        )
    );

    value(tresult.ref(), x);

    return tresult;
}


template<class Type>
template<class GeoMesh, template<class> class XPrimitiveField>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedFunction1<Type>::derivative
(
    const DimensionedField<scalar, GeoMesh, XPrimitiveField>& x
) const
{
    tmp<DimensionedField<Type, GeoMesh>> tresult
    (
        DimensionedField<Type, GeoMesh>::New
        (
            function_->name() + "'(" + x.name() + ')',
            x.mesh(),
            valueDimensions_/xDimensions_
        )
    );

    derivative(tresult.ref(), x);

    return tresult;
}


template<class Type>
template<class GeoMesh, template<class> class XPrimitiveField>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedFunction1<Type>::value
(
    const tmp<DimensionedField<scalar, GeoMesh, XPrimitiveField>>& tx
) const
{
    tmp<DimensionedField<Type, GeoMesh>> tresult = value(tx());

    tx.clear();

    return tresult;
}


template<class Type>
template<class GeoMesh, template<class> class XPrimitiveField>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedFunction1<Type>::derivative
(
    const tmp<DimensionedField<scalar, GeoMesh, XPrimitiveField>>& tx
) const
{
    tmp<DimensionedField<Type, GeoMesh>> tresult = derivative(tx());

    tx.clear();

    return tresult;
}


template<class Type>
template<class GeoMesh, template<class> class XPrimitiveField>
Foam::tmp<Foam::GeometricField<Type, GeoMesh>>
Foam::DimensionedFunction1<Type>::value
(
    const GeometricField<scalar, GeoMesh, XPrimitiveField>& x
) const
{
    tmp<GeometricField<Type, GeoMesh>> tresult
    (
        GeometricField<Type, GeoMesh>::New
        (
            function_->name() + "(" + x.name() + ')',
            x.mesh(),
            valueDimensions_
        )
    );

    value(tresult.ref(), x);

    return tresult;
}


template<class Type>
template<class GeoMesh, template<class> class XPrimitiveField>
Foam::tmp<Foam::GeometricField<Type, GeoMesh>>
Foam::DimensionedFunction1<Type>::derivative
(
    const GeometricField<scalar, GeoMesh, XPrimitiveField>& x
) const
{
    tmp<GeometricField<Type, GeoMesh>> tresult
    (
        GeometricField<Type, GeoMesh>::New
        (
            function_->name() + "'(" + x.name() + ')',
            x.mesh(),
            valueDimensions_/xDimensions_
        )
    );

    derivative(tresult.ref(), x);

    return tresult;
}


template<class Type>
template<class GeoMesh, template<class> class XPrimitiveField>
Foam::tmp<Foam::GeometricField<Type, GeoMesh>>
Foam::DimensionedFunction1<Type>::value
(
    const tmp<GeometricField<scalar, GeoMesh, XPrimitiveField>>& tx
) const
{
    tmp<GeometricField<Type, GeoMesh>> tresult = value(tx());

    tx.clear();

    return tresult;
}


template<class Type>
template<class GeoMesh, template<class> class XPrimitiveField>
Foam::tmp<Foam::GeometricField<Type, GeoMesh>>
Foam::DimensionedFunction1<Type>::derivative
(
    const tmp<GeometricField<scalar, GeoMesh, XPrimitiveField>>& tx
) const
{
    tmp<GeometricField<Type, GeoMesh>> tresult = derivative(tx());

    tx.clear();

    return tresult;
}


// ************************************************************************* //
