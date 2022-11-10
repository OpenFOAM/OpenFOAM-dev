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

#include "constantTemperature.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{
    defineTypeNameAndDebug(constantTemperature, 0);
    addToRunTimeSelectionTable
    (
        saturationTemperatureModel,
        constantTemperature,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FieldType>
Foam::tmp<FieldType>
Foam::saturationModels::constantTemperature::Tsat(const FieldType& p) const
{
    return FieldType::New("Tsat", p.mesh(), Tsat_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationModels::constantTemperature::constantTemperature
(
    const dictionary& dict
)
:
    saturationTemperatureModel(),
    Tsat_("value", dimTemperature, dict)
{}


Foam::saturationModels::constantTemperature::constantTemperature
(
    const dimensionedScalar& Tsat
)
:
    saturationTemperatureModel(),
    Tsat_(Tsat)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationModels::constantTemperature::
~constantTemperature()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

IMPLEMENT_TSAT(constantTemperature, volScalarField::Internal);


IMPLEMENT_TSAT(constantTemperature, volScalarField);


// ************************************************************************* //
