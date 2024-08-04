/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2024 OpenFOAM Foundation
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

#include "KocamustafaogullariIshiiDepartureFrequency.H"
#include "wallBoilingModelsCoefficient.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace departureFrequencyModels
{
    defineTypeNameAndDebug(KocamustafaogullariIshiiDepartureFrequency, 0);
    addToRunTimeSelectionTable
    (
        departureFrequencyModel,
        KocamustafaogullariIshiiDepartureFrequency,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ScalarFieldType>
Foam::tmp<ScalarFieldType>
Foam::wallBoilingModels::departureFrequencyModels::
KocamustafaogullariIshiiDepartureFrequency::calculate
(
    const fvMesh& mesh,
    const ScalarFieldType& dDep,
    const ScalarFieldType& rhoLiquid,
    const ScalarFieldType& rhoVapour,
    const ScalarFieldType& sigma
) const
{
    auto g = coefficient<ScalarFieldType>::value
    (
        mesh.lookupObject<uniformDimensionedVectorField>("g")
    );

    auto Cf = coefficient<ScalarFieldType>::value(Cf_);

    return
        (Cf/dDep)*pow025(sigma*mag(g)*(rhoLiquid - rhoVapour)/sqr(rhoLiquid));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::departureFrequencyModels::
KocamustafaogullariIshiiDepartureFrequency::
KocamustafaogullariIshiiDepartureFrequency
(
    const dictionary& dict
)
:
    departureFrequencyModel(),
    Cf_("Cf", dimless, dict, 1.18)
{}


Foam::wallBoilingModels::departureFrequencyModels::
KocamustafaogullariIshiiDepartureFrequency::
KocamustafaogullariIshiiDepartureFrequency
(
    const KocamustafaogullariIshiiDepartureFrequency& model
)
:
    departureFrequencyModel(model),
    Cf_(model.Cf_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallBoilingModels::departureFrequencyModels::
KocamustafaogullariIshiiDepartureFrequency::
~KocamustafaogullariIshiiDepartureFrequency()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::departureFrequencyModels::
KocamustafaogullariIshiiDepartureFrequency::fDeparture
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
Foam::wallBoilingModels::departureFrequencyModels::
KocamustafaogullariIshiiDepartureFrequency::fDeparture
(
    const phaseModel& liquid,
    const phaseModel& vapour,
    const phaseModel& solid,
    const volScalarField& Tf,
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
            vapour.rho(),
            liquid.fluid().sigma(phaseInterfaceKey(liquid, vapour))()
        );
}


void Foam::wallBoilingModels::departureFrequencyModels::
KocamustafaogullariIshiiDepartureFrequency::write(Ostream& os) const
{
    departureFrequencyModel::write(os);
    writeKeyword(os, "Cf") << Cf_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //
