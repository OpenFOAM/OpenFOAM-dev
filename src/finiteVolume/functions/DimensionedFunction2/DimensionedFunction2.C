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

#include "DimensionedFunction2.H"

// * * * * * * * * * * * * * * Private Constructors  * * * * * * * * * * * * //

template<class Type>
Foam::DimensionedFunction2<Type>::DimensionedFunction2
(
    const word& name,
    const Function2s::unitSets& units,
    const dictionary& dict
)
:
    xDimensions_(units.x.dimensions()),
    yDimensions_(units.y.dimensions()),
    valueDimensions_(units.value.dimensions()),
    function_(Function2<Type>::New(name, units, dict))
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
template
<
    class GeoMesh,
    template<class> class ResultPrimitiveField,
    template<class> class XPrimitiveField,
    template<class> class YPrimitiveField
>
void Foam::DimensionedFunction2<Type>::value
(
    DimensionedField<scalar, GeoMesh, ResultPrimitiveField>& result,
    const DimensionedField<scalar, GeoMesh, XPrimitiveField>& x,
    const DimensionedField<scalar, GeoMesh, YPrimitiveField>& y
) const
{
    x.dimensions() = xDimensions_;

    y.dimensions() = yDimensions_;

    result.primitiveFieldRef() =
        function_->value(x.primitiveField(), y.primitiveField());
}


template<class Type>
template
<
    class GeoMesh,
    template<class> class ResultPrimitiveField,
    template<class> class XPrimitiveField,
    template<class> class YPrimitiveField
>
void Foam::DimensionedFunction2<Type>::value
(
    GeometricField<scalar, GeoMesh, ResultPrimitiveField>& result,
    const GeometricField<scalar, GeoMesh, XPrimitiveField>& x,
    const GeometricField<scalar, GeoMesh, YPrimitiveField>& y
) const
{
    value(result.internalFieldRef(), x.internalField(), y.internalField());

    typename GeometricField<Type, GeoMesh>::Boundary&
        bresult = result.boundaryFieldRef();

    const typename GeometricField<scalar, GeoMesh, XPrimitiveField>::Boundary&
        bx = x.boundaryField();

    const typename GeometricField<scalar, GeoMesh, YPrimitiveField>::Boundary&
        by = y.boundaryField();

    forAll(bresult, patchi)
    {
        bresult[patchi] = function_->value(bx[patchi], by[patchi]);
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::DimensionedFunction2<Type>>
Foam::DimensionedFunction2<Type>::New
(
    const word& name,
    const Function2s::unitSets& units,
    const dictionary& dict
)
{
    return
        autoPtr<DimensionedFunction2<Type>>
        (
            new DimensionedFunction2(name, units, dict)
        );
}


template<class Type>
Foam::autoPtr<Foam::DimensionedFunction2<Type>>
Foam::DimensionedFunction2<Type>::New
(
    const word& name,
    const unitSet& xUnits,
    const unitSet& yUnits,
    const unitSet& valueUnits,
    const dictionary& dict
)
{
    return
        autoPtr<DimensionedFunction2<Type>>
        (
            new DimensionedFunction2(name, {xUnits, yUnits, valueUnits}, dict)
        );
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::DimensionedFunction2<Type>::value
(
    const scalarField& x,
    const scalarField& y
) const
{
    return function_->value(x, y);
}


template<class Type>
template
<
    class GeoMesh,
    template<class> class XPrimitiveField,
    template<class> class YPrimitiveField
>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedFunction2<Type>::value
(
    const DimensionedField<scalar, GeoMesh, XPrimitiveField>& x,
    const DimensionedField<scalar, GeoMesh, YPrimitiveField>& y
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

    value(tresult.ref(), x, y);

    return tresult;
}


template<class Type>
template
<
    class GeoMesh,
    template<class> class XPrimitiveField,
    template<class> class YPrimitiveField
>
Foam::tmp<Foam::GeometricField<Type, GeoMesh>>
Foam::DimensionedFunction2<Type>::value
(
    const GeometricField<scalar, GeoMesh, XPrimitiveField>& x,
    const GeometricField<scalar, GeoMesh, YPrimitiveField>& y
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

    value(tresult.ref(), x, y);

    return tresult;
}


// ************************************************************************* //
