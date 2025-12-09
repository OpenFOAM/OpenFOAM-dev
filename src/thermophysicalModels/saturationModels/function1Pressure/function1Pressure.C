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

\*---------------------------------------------------------------------------*/

#include "function1Pressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{
    defineTypeNameAndDebug(function1Pressure, 0);
    addToRunTimeSelectionTable
    (
        saturationPressureModel,
        function1Pressure,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationModels::function1Pressure::function1Pressure
(
    const dictionary& dict
)
:
    saturationPressureModel(),
    function_
    (
        Function1<scalar>::New("function", dimTemperature, dimPressure, dict)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationModels::function1Pressure::~function1Pressure()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::saturationModels::function1Pressure::pSat
(
    const scalarField& T
) const
{
    return function_->value(T);
}


Foam::tmp<Foam::scalarField>
Foam::saturationModels::function1Pressure::pSatPrime
(
    const scalarField& T
) const
{
    const scalar dT(rootSmall);

    return
        (
            function_->value(T + dT/2)
          - function_->value(T - dT/2)
        )/dT;
}


Foam::tmp<Foam::scalarField>
Foam::saturationModels::function1Pressure::lnPSat(const scalarField& T) const
{
    return log(pSat(T));
}


// ************************************************************************* //
