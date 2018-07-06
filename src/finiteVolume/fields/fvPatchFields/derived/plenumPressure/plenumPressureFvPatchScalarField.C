/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "plenumPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plenumPressureFvPatchScalarField::plenumPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    gamma_(1.4),
    R_(287.04),
    supplyMassFlowRate_(1.0),
    supplyTotalTemperature_(300.0),
    plenumVolume_(1.0),
    plenumDensity_(1.0),
    plenumDensityOld_(1.0),
    plenumTemperature_(300.0),
    plenumTemperatureOld_(300.0),
    rho_(1.0),
    hasRho_(false),
    inletAreaRatio_(1.0),
    inletDischargeCoefficient_(1.0),
    timeScale_(0.0),
    timeIndex_(-1),
    phiName_("phi"),
    UName_("U")
{}


Foam::plenumPressureFvPatchScalarField::plenumPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    gamma_(readScalar(dict.lookup("gamma"))),
    R_(readScalar(dict.lookup("R"))),
    supplyMassFlowRate_(readScalar(dict.lookup("supplyMassFlowRate"))),
    supplyTotalTemperature_
    (
        readScalar(dict.lookup("supplyTotalTemperature"))
    ),
    plenumVolume_(readScalar(dict.lookup("plenumVolume"))),
    plenumDensity_(readScalar(dict.lookup("plenumDensity"))),
    plenumTemperature_(readScalar(dict.lookup("plenumTemperature"))),
    rho_(1.0),
    hasRho_(false),
    inletAreaRatio_(readScalar(dict.lookup("inletAreaRatio"))),
    inletDischargeCoefficient_
    (
        readScalar(dict.lookup("inletDischargeCoefficient"))
    ),
    timeScale_(dict.lookupOrDefault<scalar>("timeScale", 0.0)),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    UName_(dict.lookupOrDefault<word>("U", "U"))
{
    if (dict.found("rho"))
    {
        rho_ = readScalar(dict.lookup("rho"));
        hasRho_ = true;
    }
}


Foam::plenumPressureFvPatchScalarField::plenumPressureFvPatchScalarField
(
    const plenumPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    gamma_(ptf.gamma_),
    R_(ptf.R_),
    supplyMassFlowRate_(ptf.supplyMassFlowRate_),
    supplyTotalTemperature_(ptf.supplyTotalTemperature_),
    plenumVolume_(ptf.plenumVolume_),
    plenumDensity_(ptf.plenumDensity_),
    plenumTemperature_(ptf.plenumTemperature_),
    rho_(ptf.rho_),
    hasRho_(ptf.hasRho_),
    inletAreaRatio_(ptf.inletAreaRatio_),
    inletDischargeCoefficient_(ptf.inletDischargeCoefficient_),
    timeScale_(ptf.timeScale_),
    phiName_(ptf.phiName_),
    UName_(ptf.UName_)
{}


Foam::plenumPressureFvPatchScalarField::plenumPressureFvPatchScalarField
(
    const plenumPressureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    gamma_(tppsf.gamma_),
    R_(tppsf.R_),
    supplyMassFlowRate_(tppsf.supplyMassFlowRate_),
    supplyTotalTemperature_(tppsf.supplyTotalTemperature_),
    plenumVolume_(tppsf.plenumVolume_),
    plenumDensity_(tppsf.plenumDensity_),
    plenumTemperature_(tppsf.plenumTemperature_),
    rho_(tppsf.rho_),
    hasRho_(tppsf.hasRho_),
    inletAreaRatio_(tppsf.inletAreaRatio_),
    inletDischargeCoefficient_(tppsf.inletDischargeCoefficient_),
    timeScale_(tppsf.timeScale_),
    phiName_(tppsf.phiName_),
    UName_(tppsf.UName_)
{}


Foam::plenumPressureFvPatchScalarField::plenumPressureFvPatchScalarField
(
    const plenumPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    gamma_(tppsf.gamma_),
    R_(tppsf.R_),
    supplyMassFlowRate_(tppsf.supplyMassFlowRate_),
    supplyTotalTemperature_(tppsf.supplyTotalTemperature_),
    plenumVolume_(tppsf.plenumVolume_),
    plenumDensity_(tppsf.plenumDensity_),
    plenumTemperature_(tppsf.plenumTemperature_),
    rho_(tppsf.rho_),
    hasRho_(tppsf.hasRho_),
    inletAreaRatio_(tppsf.inletAreaRatio_),
    inletDischargeCoefficient_(tppsf.inletDischargeCoefficient_),
    timeScale_(tppsf.timeScale_),
    phiName_(tppsf.phiName_),
    UName_(tppsf.UName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::plenumPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Patch properties
    const fvPatchField<scalar>& p = *this;
    const fvPatchField<scalar>& p_old =
        db().lookupObject<volScalarField>
        (
            internalField().name()
        ).oldTime().boundaryField()[patch().index()];
    const fvPatchField<vector>& U =
        patch().lookupPatchField<volVectorField, vector>(UName_);
    const fvsPatchField<scalar>& phi =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    // Get the timestep
    const scalar dt = db().time().deltaTValue();

    // Check if operating at a new time index and update the old-time properties
    // if so
    if (timeIndex_ != db().time().timeIndex())
    {
        timeIndex_ = db().time().timeIndex();
        plenumDensityOld_ = plenumDensity_;
        plenumTemperatureOld_ = plenumTemperature_;
    }

    // Calculate the current mass flow rate
    scalar massFlowRate(1.0);
    if (phi.internalField().dimensions() == dimVelocity*dimArea)
    {
        if (hasRho_)
        {
            massFlowRate = - gSum(rho_*phi);
        }
        else
        {
            FatalErrorInFunction
                << "The density must be specified when using a volumetric flux."
                << exit(FatalError);
        }
    }
    else if
    (
        phi.internalField().dimensions()
     == dimDensity*dimVelocity*dimArea
    )
    {
        if (hasRho_)
        {
            FatalErrorInFunction
                << "The density must be not specified when using a mass flux."
                << exit(FatalError);
        }
        else
        {
            massFlowRate = - gSum(phi);
        }
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of phi are not correct"
            << "\n    on patch " << patch().name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath() << nl
            << exit(FatalError);
    }

    // Calcaulate the specific heats
    const scalar cv = R_/(gamma_ - 1), cp = R_*gamma_/(gamma_ - 1);

    // Calculate the new plenum properties
    plenumDensity_ =
        plenumDensityOld_
      + (dt/plenumVolume_)*(supplyMassFlowRate_ - massFlowRate);
    plenumTemperature_ =
        plenumTemperatureOld_
      + (dt/(plenumDensity_*cv*plenumVolume_))
      * (
            supplyMassFlowRate_
           *(cp*supplyTotalTemperature_ - cv*plenumTemperature_)
          - massFlowRate*R_*plenumTemperature_
        );
    const scalar plenumPressure = plenumDensity_*R_*plenumTemperature_;

    // Squared velocity magnitude at exit of channels
    const scalarField U_e(magSqr(U/inletAreaRatio_));

    // Exit temperature to plenum temperature ratio
    const scalarField r
    (
        1.0 - (gamma_ - 1.0)*U_e/(2.0*gamma_*R_*plenumTemperature_)
    );

    // Quadratic coefficient (others not needed as b = +1.0 and c = -1.0)
    const scalarField a
    (
        (1.0 - r)/(r*r*inletDischargeCoefficient_*inletDischargeCoefficient_)
    );

    // Isentropic exit temperature to plenum temperature ratio
    const scalarField s(2.0/(1.0 + sqrt(1.0 + 4.0*a)));

    // Exit pressure to plenum pressure ratio
    const scalarField t(pow(s, gamma_/(gamma_ - 1.0)));

    // Limit to prevent outflow
    const scalarField p_new
    (
        (1.0 - pos0(phi))*t*plenumPressure + pos0(phi)*max(p, plenumPressure)
    );

    // Relaxation fraction
    const scalar oneByFraction = timeScale_/dt;
    const scalar fraction = oneByFraction < 1.0 ? 1.0 : 1.0/oneByFraction;

    // Set the new value
    operator==((1.0 - fraction)*p_old + fraction*p_new);
    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::plenumPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("gamma") << gamma_
        << token::END_STATEMENT << nl;
    os.writeKeyword("R") << R_
        << token::END_STATEMENT << nl;
    os.writeKeyword("supplyMassFlowRate") << supplyMassFlowRate_
        << token::END_STATEMENT << nl;
    os.writeKeyword("supplyTotalTemperature") << supplyTotalTemperature_
        << token::END_STATEMENT << nl;
    os.writeKeyword("plenumVolume") << plenumVolume_
        << token::END_STATEMENT << nl;
    os.writeKeyword("plenumDensity") << plenumDensity_
        << token::END_STATEMENT << nl;
    os.writeKeyword("plenumTemperature") << plenumTemperature_
        << token::END_STATEMENT << nl;
    if (hasRho_)
    {
        os.writeKeyword("rho") << rho_
            << token::END_STATEMENT << nl;
    }
    os.writeKeyword("inletAreaRatio") << inletAreaRatio_
        << token::END_STATEMENT << nl;
    os.writeKeyword("inletDischargeCoefficient") << inletDischargeCoefficient_
        << token::END_STATEMENT << nl;
    writeEntryIfDifferent<scalar>(os, "timeScale", 0.0, timeScale_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        plenumPressureFvPatchScalarField
    );
}

// ************************************************************************* //
