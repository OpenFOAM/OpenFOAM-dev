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
#include "saturationPressureModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(saturationPressureModel, 0);
    defineRunTimeSelectionTable(saturationPressureModel, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FieldType>
Foam::tmp<FieldType> Foam::saturationPressureModel::pSat
(
    const FieldType& T
) const
{
    return
        saturationModels::evaluate
        (
            T,
            "pSat",
            dimPressure,
            *this,
            &saturationPressureModel::pSat
        );
}


template<class FieldType>
Foam::tmp<FieldType> Foam::saturationPressureModel::pSatPrime
(
    const FieldType& T
) const
{
    return
        saturationModels::evaluate
        (
            T,
            "pSatPrime",
            dimPressure/dimTemperature,
            *this,
            &saturationPressureModel::pSatPrime
        );
}


template<class FieldType>
Foam::tmp<FieldType> Foam::saturationPressureModel::lnPSat
(
    const FieldType& T
) const
{
    return
        saturationModels::evaluate
        (
            T,
            "lnPSat",
            dimless,
            *this,
            &saturationPressureModel::lnPSat
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationPressureModel::saturationPressureModel()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationPressureModel::~saturationPressureModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

IMPLEMENT_PSAT(saturationPressureModel, volScalarField::Internal);


IMPLEMENT_PSAT(saturationPressureModel, volScalarField);


// ************************************************************************* //
