/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2023 OpenFOAM Foundation
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

#include "Cole.H"
#include "wallBoilingModelsCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace departureFrequencyModels
{
    defineTypeNameAndDebug(Cole, 0);
    addToRunTimeSelectionTable
    (
        departureFrequencyModel,
        Cole,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ScalarFieldType>
Foam::tmp<ScalarFieldType>
Foam::wallBoilingModels::departureFrequencyModels::Cole::calculate
(
    const fvMesh& mesh,
    const ScalarFieldType& dDep,
    const ScalarFieldType& rhoLiquid,
    const ScalarFieldType& rhoVapour
) const
{
    auto g = coefficient<ScalarFieldType>::value
    (
        mesh.lookupObject<uniformDimensionedVectorField>("g")
    );

    const dimensionedScalar dRhoMin_(dimDensity, 0.1);
    auto dRhoMin = coefficient<ScalarFieldType>::value(dRhoMin_);

    return
        sqrt
        (
            4*mag(g)*max(rhoLiquid - rhoVapour, dRhoMin)/(3*dDep*rhoLiquid)
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::departureFrequencyModels::Cole::Cole
(
    const dictionary& dict
)
:
    departureFrequencyModel()
{}


Foam::wallBoilingModels::departureFrequencyModels::Cole::Cole
(
    const Cole& model
)
:
    departureFrequencyModel(model)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallBoilingModels::departureFrequencyModels::Cole::~Cole()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::departureFrequencyModels::Cole::fDeparture
(
    const phaseModel& liquid,
    const phaseModel& vapour,
    const label patchi,
    const scalarField& Tl,
    const scalarField& Tsatw,
    const scalarField& L,
    const scalarField& dDep
) const
{
    return
        calculate
        (
            liquid.mesh(),
            dDep,
            liquid.thermo().rho(patchi)(),
            vapour.thermo().rho(patchi)()
        );
}


Foam::tmp<Foam::volScalarField>
Foam::wallBoilingModels::departureFrequencyModels::Cole::fDeparture
(
    const phaseModel& liquid,
    const phaseModel& vapour,
    const phaseModel& solid,
    const volScalarField& Tl,
    const volScalarField& Tsatw,
    const volScalarField& L,
    const volScalarField& dDep
) const
{
    return
        calculate
        (
            liquid.mesh(),
            dDep,
            liquid.rho(),
            vapour.rho()
        );
}


// ************************************************************************* //
