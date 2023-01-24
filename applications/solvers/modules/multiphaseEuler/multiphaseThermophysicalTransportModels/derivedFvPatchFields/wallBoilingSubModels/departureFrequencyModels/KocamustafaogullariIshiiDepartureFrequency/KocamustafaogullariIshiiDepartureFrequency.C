/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2023 OpenFOAM Foundation
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
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "compressibleMomentumTransportModels.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "phaseSystem.H"

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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::departureFrequencyModels::
KocamustafaogullariIshiiDepartureFrequency::
KocamustafaogullariIshiiDepartureFrequency(const dictionary& dict)
:
    departureFrequencyModel(),
    Cf_(dict.lookupOrDefault<scalar>("Cf", 1.18))
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
    // Gravitational acceleration
    const uniformDimensionedVectorField& g =
        liquid.mesh().lookupObject<uniformDimensionedVectorField>("g");

    const scalarField rhoLiquid(liquid.thermo().rho(patchi));
    const scalarField rhoVapor(min(vapour.thermo().rho(patchi), rhoLiquid));

    const tmp<volScalarField> tsigma
    (
        liquid.fluid().sigma(phaseInterface(liquid, vapour))
    );

    const volScalarField& sigma = tsigma();
    const fvPatchScalarField& sigmaw = sigma.boundaryField()[patchi];

    return (Cf_/dDep)*pow025
    (
        sigmaw*mag(g.value())*(rhoLiquid - rhoVapor)/sqr(rhoLiquid)
    );
}


void Foam::wallBoilingModels::departureFrequencyModels::
KocamustafaogullariIshiiDepartureFrequency::write(Ostream& os) const
{
    departureFrequencyModel::write(os);
    writeKeyword(os, "Cf") << Cf_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //
