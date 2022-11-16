/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::departureDiameterModels::
KocamustafaogullariIshiiDepartureDiameter::
KocamustafaogullariIshiiDepartureDiameter
(
    const dictionary& dict
)
:
    departureDiameterModel(),
    phi_(dict.lookup<scalar>("phi"))
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
    const phaseModel& vapor,
    const label patchi,
    const scalarField& Tl,
    const scalarField& Tsatw,
    const scalarField& L
) const
{
    // Gravitational acceleration
    const uniformDimensionedVectorField& g =
        liquid.mesh().lookupObject<uniformDimensionedVectorField>("g");

    const scalarField rhoLiquid(liquid.thermo().rho(patchi));
    const scalarField rhoVapor(vapor.thermo().rho(patchi));

    const scalarField rhoM((rhoLiquid - rhoVapor)/rhoVapor);

    const scalarField sigmaw
    (
        liquid.fluid().sigma(phaseInterface(liquid, vapor), patchi)
    );

    return
        0.0012*pow(rhoM, 0.9)*0.0208*phi_
       *sqrt(sigmaw/(mag(g.value())*(rhoLiquid - rhoVapor)));
}


void Foam::wallBoilingModels::departureDiameterModels::
KocamustafaogullariIshiiDepartureDiameter::write(Ostream& os) const
{
    departureDiameterModel::write(os);
    writeEntry(os, "phi", phi_);
}


// ************************************************************************* //
