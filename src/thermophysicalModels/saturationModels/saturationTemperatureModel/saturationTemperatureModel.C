/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2025 OpenFOAM Foundation
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

#include "saturationModels.H"
#include "saturationTemperatureModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FieldType>
Foam::tmp<FieldType> Foam::saturationTemperatureModel::Tsat
(
    const FieldType& p
) const
{
    return
        saturationModels::evaluate
        (
            p,
            "Tsat",
            dimTemperature,
            *this,
            &saturationTemperatureModel::Tsat
        );
}


template<class FieldType>
Foam::tmp<FieldType> Foam::saturationTemperatureModel::TsatPrime
(
    const FieldType& p
) const
{
    return
        saturationModels::evaluate
        (
            p,
            "TsatPrime",
            dimTemperature/dimPressure,
            *this,
            &saturationTemperatureModel::TsatPrime
        );
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(saturationTemperatureModel, 0);
    defineRunTimeSelectionTable(saturationTemperatureModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationTemperatureModel::saturationTemperatureModel()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationTemperatureModel::~saturationTemperatureModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

IMPLEMENT_TSAT(saturationTemperatureModel, volScalarField::Internal);


IMPLEMENT_TSAT(saturationTemperatureModel, volScalarField);


// ************************************************************************* //
