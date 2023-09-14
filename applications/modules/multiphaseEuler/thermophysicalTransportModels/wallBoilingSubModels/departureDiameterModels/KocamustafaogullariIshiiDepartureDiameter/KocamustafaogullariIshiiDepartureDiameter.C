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

#include "KocamustafaogullariIshiiDepartureDiameter.H"
#include "wallBoilingModelsCoefficient.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace departureDiameterModels
{
    defineTypeNameAndDebug(KocamustafaogullariIshiiDepartureDiameter, 0);
    addToRunTimeSelectionTable
    (
        departureDiameterModel,
        KocamustafaogullariIshiiDepartureDiameter,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ScalarFieldType>
Foam::tmp<ScalarFieldType>
Foam::wallBoilingModels::departureDiameterModels::
KocamustafaogullariIshiiDepartureDiameter::calculate
(
    const fvMesh& mesh,
    const ScalarFieldType& Tl,
    const ScalarFieldType& Tsatw,
    const ScalarFieldType& L,
    const ScalarFieldType& rhoLiquid,
    const ScalarFieldType& rhoVapour,
    const ScalarFieldType& sigma
) const
{
    auto g = coefficient<ScalarFieldType>::value
    (
        mesh.lookupObject<uniformDimensionedVectorField>("g")
    );

    auto phi = coefficient<ScalarFieldType>::value(phi_);

    return
        0.0012*pow((rhoLiquid - rhoVapour)/rhoVapour, 0.9)*0.0208*phi
       *sqrt(sigma/(mag(g)*(rhoLiquid - rhoVapour)));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::departureDiameterModels::
KocamustafaogullariIshiiDepartureDiameter::
KocamustafaogullariIshiiDepartureDiameter
(
    const dictionary& dict
)
:
    departureDiameterModel(),
    phi_("phi", dimless, dict)
{}


Foam::wallBoilingModels::departureDiameterModels::
KocamustafaogullariIshiiDepartureDiameter::
KocamustafaogullariIshiiDepartureDiameter
(
    const KocamustafaogullariIshiiDepartureDiameter& model
)
:
    departureDiameterModel(),
    phi_(model.phi_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallBoilingModels::departureDiameterModels::
KocamustafaogullariIshiiDepartureDiameter::
~KocamustafaogullariIshiiDepartureDiameter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::departureDiameterModels::
KocamustafaogullariIshiiDepartureDiameter::dDeparture
(
    const phaseModel& liquid,
    const phaseModel& vapour,
    const label patchi,
    const scalarField& Tl,
    const scalarField& Tsatw,
    const scalarField& L
) const
{

    return
        calculate
        (
            liquid.mesh(),
            Tl,
            Tsatw,
            L,
            static_cast<const scalarField&>
            (
                liquid.rho().boundaryField()[patchi]
            ),
            static_cast<const scalarField&>
            (
                vapour.rho().boundaryField()[patchi]
            ),
            liquid.fluid().sigma(phaseInterfaceKey(liquid, vapour), patchi)()
        );
}


Foam::tmp<Foam::volScalarField>
Foam::wallBoilingModels::departureDiameterModels::
KocamustafaogullariIshiiDepartureDiameter::dDeparture
(
    const phaseModel& liquid,
    const phaseModel& vapour,
    const phaseModel& solid,
    const volScalarField& Tf,
    const volScalarField& Tsatw,
    const volScalarField& L
) const
{
    return calculate
    (
        liquid.mesh(),
        liquid.thermo().T(),
        Tsatw,
        L,
        liquid.rho(),
        vapour.rho(),
        liquid.fluid().sigma(phaseInterfaceKey(liquid, vapour))()
    );
}


void Foam::wallBoilingModels::departureDiameterModels::
KocamustafaogullariIshiiDepartureDiameter::write(Ostream& os) const
{
    departureDiameterModel::write(os);
    writeEntry(os, "phi", phi_);
}


// ************************************************************************* //
