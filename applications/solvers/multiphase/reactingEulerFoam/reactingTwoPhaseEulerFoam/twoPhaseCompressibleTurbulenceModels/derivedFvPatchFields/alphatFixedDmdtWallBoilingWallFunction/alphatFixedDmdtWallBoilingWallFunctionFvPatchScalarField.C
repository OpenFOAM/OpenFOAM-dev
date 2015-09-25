/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

#include "alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField.H"
#include "compressibleTurbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "twoPhaseSystem.H"
#include "phaseSystem.H"
#include "ThermalPhaseChangePhaseSystem.H"
#include "MomentumTransferPhaseSystem.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::maxExp_
    = 50.0;
scalar alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::tolerance_
    = 0.01;
label alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::maxIters_
    = 10;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorIn
        (
            "alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField"
            "::checkType()"
        )
            << "Patch type for patch " << patch().name() << " must be wall\n"
            << "Current patch type is " << patch().type() << nl
            << exit(FatalError);
    }
}


tmp<scalarField>
alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::Psmooth
(
    const scalarField& Prat
) const
{
    return 9.24*(pow(Prat, 0.75) - 1.0)*(1.0 + 0.28*exp(-0.007*Prat));
}


tmp<scalarField>
alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::yPlusTherm
(
    const scalarField& P,
    const scalarField& Prat
) const
{

    tmp<scalarField> typtf(new scalarField(this->size()));
    scalarField& yptf = typtf();

    forAll(yptf, faceI)
    {
        scalar ypt = 11.0;

        for (int i=0; i<maxIters_; i++)
        {
            scalar f = ypt - (log(E_*ypt)/kappa_ + P[faceI])/Prat[faceI];
            scalar df = 1.0 - 1.0/(ypt*kappa_*Prat[faceI]);
            scalar yptNew = ypt - f/df;

            if (yptNew < VSMALL)
            {
                yptf[faceI] =  0;
            }
            else if (mag(yptNew - ypt) < tolerance_)
            {
                yptf[faceI] = yptNew;
            }
            else
            {
                ypt = yptNew;
            }
        }

        yptf[faceI] = ypt;
    }

    return typtf;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::
alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphatPhaseChangeWallFunctionFvPatchScalarField(p, iF),
    Prt_(0.85),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
    fixedDmdt_(0.0)
{
    checkType();
}


alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::
alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphatPhaseChangeWallFunctionFvPatchScalarField(p, iF, dict),
    Prt_(dict.lookupOrDefault<scalar>("Prt", 0.85)),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8)),
    fixedDmdt_(dict.lookupOrDefault<scalar>("fixedDmdt", 0.0))
{}


alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::
alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
(
    const alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphatPhaseChangeWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    Prt_(ptf.Prt_),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_),
    fixedDmdt_(ptf.fixedDmdt_)
{}


alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::
alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
(
    const alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField& awfpsf
)
:
    alphatPhaseChangeWallFunctionFvPatchScalarField(awfpsf),
    Prt_(awfpsf.Prt_),
    Cmu_(awfpsf.Cmu_),
    kappa_(awfpsf.kappa_),
    E_(awfpsf.E_),
    fixedDmdt_(awfpsf.fixedDmdt_)
{}


alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::
alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
(
    const alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField& awfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphatPhaseChangeWallFunctionFvPatchScalarField(awfpsf, iF),
    Prt_(awfpsf.Prt_),
    Cmu_(awfpsf.Cmu_),
    kappa_(awfpsf.kappa_),
    E_(awfpsf.E_),
    fixedDmdt_(awfpsf.fixedDmdt_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Lookup the fluid model
    const ThermalPhaseChangePhaseSystem
    <
        MomentumTransferPhaseSystem<twoPhaseSystem>
    >& fluid =
        refCast
        <
            const ThermalPhaseChangePhaseSystem
            <
                MomentumTransferPhaseSystem<twoPhaseSystem>
            >
        >
        (
            db().lookupObject<phaseSystem>("phaseProperties")
        );

    const phaseModel& liquid
    (
        fluid.phase1().name() == dimensionedInternalField().group()
      ? fluid.phase1()
      : fluid.phase2()
    );

    const label patchi = patch().index();

    // Retrieve turbulence properties from model
    const compressibleTurbulenceModel& turbModel =
        db().lookupObject<compressibleTurbulenceModel>
        (
            IOobject::groupName
            (
                compressibleTurbulenceModel::propertiesName,
                dimensionedInternalField().group()
            )
        );

    const scalar Cmu25 = pow025(Cmu_);

    const scalarField& y = turbModel.y()[patchi];

    const tmp<scalarField> tmuw = turbModel.mu(patchi);
    const scalarField& muw = tmuw();

    const tmp<scalarField> talphaw = liquid.thermo().alpha(patchi);
    const scalarField& alphaw = talphaw();

    scalarField& alphatw = *this;

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    const fvPatchScalarField& kw = k.boundaryField()[patchi];

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));
    const scalarField magGradUw(mag(Uw.snGrad()));

    const fvPatchScalarField& rhow = turbModel.rho().boundaryField()[patchi];
    const fvPatchScalarField& hew =
        liquid.thermo().he().boundaryField()[patchi];

    const fvPatchScalarField& Tw =
        liquid.thermo().T().boundaryField()[patchi];

    scalarField Tp(Tw.patchInternalField());

    // Heat flux [W/m2] - lagging alphatw
    const scalarField qDot
    (
        (alphatw + alphaw)*hew.snGrad()
    );

    scalarField uTau(Cmu25*sqrt(kw));

    scalarField yPlus(uTau*y/(muw/rhow));

    scalarField Pr(muw/alphaw);

    // Molecular-to-turbulent Prandtl number ratio
    scalarField Prat(Pr/Prt_);

    // Thermal sublayer thickness
    scalarField P(this->Psmooth(Prat));

    scalarField yPlusTherm(this->yPlusTherm(P, Prat));

    scalarField alphatConv(this->size(), 0.0);

    // Populate boundary values
    forAll(alphatw, faceI)
    {
        // Evaluate new effective thermal diffusivity
        scalar alphaEff = 0.0;
        if (yPlus[faceI] < yPlusTherm[faceI])
        {
            scalar A = qDot[faceI]*rhow[faceI]*uTau[faceI]*y[faceI];
            scalar B = qDot[faceI]*Pr[faceI]*yPlus[faceI];
            scalar C = Pr[faceI]*0.5*rhow[faceI]*uTau[faceI]*sqr(magUp[faceI]);
            alphaEff = A/(B + C + VSMALL);
        }
        else
        {
            scalar A = qDot[faceI]*rhow[faceI]*uTau[faceI]*y[faceI];
            scalar B =
                qDot[faceI]*Prt_*(1.0/kappa_*log(E_*yPlus[faceI]) + P[faceI]);
            scalar magUc =
                uTau[faceI]/kappa_*log(E_*yPlusTherm[faceI]) - mag(Uw[faceI]);
            scalar C =
                0.5*rhow[faceI]*uTau[faceI]
               *(Prt_*sqr(magUp[faceI]) + (Pr[faceI] - Prt_)*sqr(magUc));
            alphaEff = A/(B + C + VSMALL);
        }

        // Update convective heat transfer turbulent thermal diffusivity
        alphatConv[faceI] = max(0.0, alphaEff - alphaw[faceI]);
    }

    dmdt_ = fixedDmdt_;

    operator==(alphatConv);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("Prt") << Prt_ << token::END_STATEMENT << nl;
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
    os.writeKeyword("fixedDmdt") << fixedDmdt_ << token::END_STATEMENT << nl;
    dmdt_.writeEntry("dmdt", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
