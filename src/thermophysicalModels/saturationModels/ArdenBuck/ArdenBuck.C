/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2022 OpenFOAM Foundation
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

#include "ArdenBuck.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{
    defineTypeNameAndDebug(ArdenBuck, 0);
    addToRunTimeSelectionTable(saturationPressureModel, ArdenBuck, dictionary);
}
}


static const Foam::dimensionedScalar zeroC
(
    "zeroC",
    Foam::dimTemperature,
    273.15
);


static const Foam::dimensionedScalar A("A", Foam::dimPressure, 611.21);
static const Foam::dimensionedScalar B("B", Foam::dimless, 18.678);
static const Foam::dimensionedScalar C("C", Foam::dimTemperature, 234.5);
static const Foam::dimensionedScalar D("D", Foam::dimTemperature, 257.14);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FieldType>
Foam::tmp<FieldType>
Foam::saturationModels::ArdenBuck::xByTC(const FieldType& TC) const
{
    return (B - TC/C)/(D + TC);
}


template<class FieldType>
Foam::tmp<FieldType>
Foam::saturationModels::ArdenBuck::pSat(const FieldType& T) const
{
    const FieldType TC(T - zeroC);

    return A*exp(TC*xByTC(TC));
}


template<class FieldType>
Foam::tmp<FieldType>
Foam::saturationModels::ArdenBuck::pSatPrime(const FieldType& T) const
{
    const FieldType TC(T - zeroC);
    const FieldType x(xByTC(TC));

    return A*exp(TC*x)*(D*x - TC/C)/(D + TC);
}


template<class FieldType>
Foam::tmp<FieldType>
Foam::saturationModels::ArdenBuck::lnPSat(const FieldType& T) const
{
    const FieldType TC(T - zeroC);

    return log(A.value()) + TC*xByTC(TC);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationModels::ArdenBuck::ArdenBuck(const dictionary& dict)
:
    saturationPressureModel()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationModels::ArdenBuck::~ArdenBuck()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

IMPLEMENT_PSAT(ArdenBuck, volScalarField::Internal);


IMPLEMENT_PSAT(ArdenBuck, volScalarField);


// ************************************************************************* //
