/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2021 OpenFOAM Foundation
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

#include "strainRateFunction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{
namespace generalisedNewtonianViscosityModels
{
    defineTypeNameAndDebug(strainRateFunction, 0);

    addToRunTimeSelectionTable
    (
        generalisedNewtonianViscosityModel,
        strainRateFunction,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarModels::generalisedNewtonianViscosityModels::strainRateFunction::
strainRateFunction
(
    const dictionary& viscosityProperties,
    const Foam::viscosity& viscosity,
    const volVectorField& U
)
:
    strainRateViscosityModel(viscosityProperties, viscosity, U),
    strainRateFunction_
    (
        Function1<scalar>::New
        (
            "function",
            viscosityProperties.optionalSubDict(typeName + "Coeffs")
        )
    )
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::laminarModels::generalisedNewtonianViscosityModels::
strainRateFunction::read
(
    const dictionary& viscosityProperties
)
{
    strainRateViscosityModel::read(viscosityProperties);

    strainRateFunction_.clear();
    strainRateFunction_ = Function1<scalar>::New
    (
        "function",
        viscosityProperties.optionalSubDict
        (
            typeName + "Coeffs"
        )
    );

    return true;
}


Foam::tmp<Foam::volScalarField>
Foam::laminarModels::generalisedNewtonianViscosityModels::strainRateFunction::
nu
(
    const volScalarField& nu0,
    const volScalarField& strainRate
) const
{
    tmp<volScalarField> tnu
    (
        volScalarField::New
        (
            IOobject::groupName(type() + ":nu", nu0.group()),
            nu0.mesh(),
            dimensionedScalar(dimViscosity, 0)
        )
    );

    tnu.ref().primitiveFieldRef() = strainRateFunction_->value(strainRate);

    volScalarField::Boundary& nuBf = tnu.ref().boundaryFieldRef();
    const volScalarField::Boundary& sigmaBf = strainRate.boundaryField();

    forAll(nuBf, patchi)
    {
        nuBf[patchi] = strainRateFunction_->value(sigmaBf[patchi]);
    }

    return tnu;
}


// ************************************************************************* //
