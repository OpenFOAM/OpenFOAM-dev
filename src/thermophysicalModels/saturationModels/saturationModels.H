/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

Namespace
    Foam::saturationModels

Description
    Namespace containing common functionality for saturationPressureModel and
    saturationTemperatureModel

\*---------------------------------------------------------------------------*/

#ifndef saturationModels_H
#define saturationModels_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Structure to store a dimensioned coefficient of the saturation model
struct coefficient
:
    public dimensionedScalar
{
    //- Inherit Constructors
    using dimensionedScalar::dimensionedScalar;

    //- Cast to a scalar
    operator scalar() const
    {
        return dimensionedScalar::value();
    }
};

//- Construct a scalar field from a uniform value
inline tmp<scalarField> evaluate
(
    const scalarField& argument,
    const word& name,
    const dimensionedScalar& value
)
{
    return tmp<scalarField>(new scalarField(argument.size(), value.value()));
}

//- Construct a volScalarField::Internal or a volScalarField field from a
//  uniform value
template<class FieldType>
inline tmp<FieldType> evaluate
(
    const FieldType& argument,
    const word& name,
    const dimensionedScalar& value
)
{
    return FieldType::New(name, argument.mesh(), value);
}

//- Use a model's scalarField method to construct a volScalarField::Internal
template<class Model>
inline tmp<volScalarField::Internal> evaluate
(
    const volScalarField::Internal& p,
    const word& name,
    const dimensionSet& dimensions,
    const Model& model,
    tmp<scalarField> (Model::*method)(const scalarField&) const
)
{
    tmp<volScalarField::Internal> tResult
    (
        volScalarField::Internal::New(name, p.mesh(), dimensions)
    );

    volScalarField::Internal& result = tResult.ref();

    result.primitiveFieldRef() = (model.*method)(p.primitiveField());

    return tResult;
}

//- Use a model's scalarField method to construct a volScalarField
template<class Model>
inline tmp<volScalarField> evaluate
(
    const volScalarField& p,
    const word& name,
    const dimensionSet& dimensions,
    const Model& model,
    tmp<scalarField> (Model::*method)(const scalarField&) const
)
{
    tmp<volScalarField> tResult
    (
        volScalarField::New(name, p.mesh(), dimensions)
    );

    volScalarField& result = tResult.ref();

    result.primitiveFieldRef() = (model.*method)(p.primitiveField());

    volScalarField::Boundary& resultBf = result.boundaryFieldRef();

    forAll(result.boundaryField(), patchi)
    {
        resultBf[patchi] = (model.*method)(p.boundaryField()[patchi]);
    }

    return tResult;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace saturationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
