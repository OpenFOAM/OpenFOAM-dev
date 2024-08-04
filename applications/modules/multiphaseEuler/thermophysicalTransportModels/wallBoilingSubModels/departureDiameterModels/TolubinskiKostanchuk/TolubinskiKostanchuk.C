/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2024 OpenFOAM Foundation
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

#include "TolubinskiKostanchuk.H"
#include "wallBoilingModelsCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace departureDiameterModels
{
    defineTypeNameAndDebug(TolubinskiKostanchuk, 0);
    addToRunTimeSelectionTable
    (
        departureDiameterModel,
        TolubinskiKostanchuk,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ScalarFieldType>
Foam::tmp<ScalarFieldType>
Foam::wallBoilingModels::departureDiameterModels::TolubinskiKostanchuk::
calculate
(
    const ScalarFieldType& Tl,
    const ScalarFieldType& Tsatw
) const
{
    auto dRef = coefficient<ScalarFieldType>::value(dRef_);
    auto dMax = coefficient<ScalarFieldType>::value(dMax_);
    auto dMin = coefficient<ScalarFieldType>::value(dMin_);

    const dimensionedScalar T45_(dimTemperature, 45);

    auto T45 = coefficient<ScalarFieldType>::value(T45_);

    return max(min(dRef*exp(- (Tsatw - Tl)/T45), dMax), dMin);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::departureDiameterModels::TolubinskiKostanchuk::
TolubinskiKostanchuk
(
    const dictionary& dict
)
:
    departureDiameterModel(),
    dRef_("dRef", dimLength, dict, 6e-4),
    dMax_("dMax", dimLength, dict, 0.0014),
    dMin_("dMin", dimLength, dict, 1e-6)
{}


Foam::wallBoilingModels::departureDiameterModels::TolubinskiKostanchuk::
TolubinskiKostanchuk
(
    const TolubinskiKostanchuk& model
)
:
    departureDiameterModel(),
    dRef_(model.dRef_),
    dMax_(model.dMax_),
    dMin_(model.dMin_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallBoilingModels::departureDiameterModels::TolubinskiKostanchuk::
~TolubinskiKostanchuk()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::departureDiameterModels::TolubinskiKostanchuk::
dDeparture
(
    const phaseModel& liquid,
    const phaseModel& vapour,
    const label patchi,
    const scalarField& Tl,
    const scalarField& Tsatw,
    const scalarField& L
) const
{
    return calculate(Tl, Tsatw);
}


Foam::tmp<Foam::volScalarField>
Foam::wallBoilingModels::departureDiameterModels::
TolubinskiKostanchuk::dDeparture
(
    const phaseModel& liquid,
    const phaseModel& vapour,
    const phaseModel& solid,
    const volScalarField& Tf,
    const volScalarField& Tsatw,
    const volScalarField& L
) const
{
    return calculate(liquid.thermo().T(), Tsatw);
}


void Foam::wallBoilingModels::departureDiameterModels::TolubinskiKostanchuk::
write
(
    Ostream& os
) const
{
    departureDiameterModel::write(os);
    writeEntry(os, "dRef", dRef_);
    writeEntry(os, "dMax", dMax_);
    writeEntry(os, "dMin", dMin_);
}


// ************************************************************************* //
