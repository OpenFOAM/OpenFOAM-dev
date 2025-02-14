/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2025 OpenFOAM Foundation
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

#include "function1Temperature.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{
    defineTypeNameAndDebug(function1Temperature, 0);
    addToRunTimeSelectionTable
    (
        saturationTemperatureModel,
        function1Temperature,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationModels::function1Temperature::function1Temperature
(
    const dictionary& dict
)
:
    saturationTemperatureModel(),
    function_
    (
        Function1<scalar>::New("function", dimPressure, dimTemperature, dict)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationModels::function1Temperature::~function1Temperature()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::saturationModels::function1Temperature::Tsat
(
    const scalarField& p
) const
{
    return function_->value(p);
}


Foam::tmp<Foam::scalarField>
Foam::saturationModels::function1Temperature::TsatPrime
(
    const scalarField& p
) const
{
    const scalar dp(rootSmall*gAverage(p));

    return
        (
            function_->value(p + dp/2)
          - function_->value(p - dp/2)
        )/dp;
}


// ************************************************************************* //
