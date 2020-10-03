/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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

#include "alphatPhaseJayatillekeWallFunctionFvPatchScalarField.H"
#include "phaseSystem.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar alphatPhaseJayatillekeWallFunctionFvPatchScalarField::maxExp_  = 50.0;
scalar alphatPhaseJayatillekeWallFunctionFvPatchScalarField::tolerance_ = 0.01;
label alphatPhaseJayatillekeWallFunctionFvPatchScalarField::maxIters_ = 10;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

tmp<scalarField> alphatPhaseJayatillekeWallFunctionFvPatchScalarField::Psmooth
(
    const scalarField& Prat
) const
{
    return 9.24*(pow(Prat, 0.75) - 1.0)*(1.0 + 0.28*exp(-0.007*Prat));
}


tmp<scalarField>
alphatPhaseJayatillekeWallFunctionFvPatchScalarField::yPlusTherm
(
    const nutWallFunctionFvPatchScalarField& nutw,
    const scalarField& P,
    const scalarField& Prat
) const
{
    tmp<scalarField> typsf(new scalarField(this->size()));
    scalarField& ypsf = typsf.ref();

    forAll(ypsf, facei)
    {
        scalar ypt = 11.0;

        for (int i=0; i<maxIters_; i++)
        {
            const scalar f =
                ypt - (log(nutw.E()*ypt)/nutw.kappa() + P[facei])/Prat[facei];
            const scalar df = 1 - 1.0/(ypt*nutw.kappa()*Prat[facei]);
            const scalar yptNew = ypt - f/df;

            if (yptNew < vSmall)
            {
                ypsf[facei] = 0;
                break;
            }
            else if (mag(yptNew - ypt) < tolerance_)
            {
                ypsf[facei] = yptNew;
            }
            else
            {
                ypt = yptNew;
            }
        }

        ypsf[facei] = ypt;
    }

    return typsf;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatPhaseJayatillekeWallFunctionFvPatchScalarField::
alphatPhaseJayatillekeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Prt_(0.85)
{}


alphatPhaseJayatillekeWallFunctionFvPatchScalarField::
alphatPhaseJayatillekeWallFunctionFvPatchScalarField
(
    const alphatPhaseJayatillekeWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Prt_(ptf.Prt_)
{}


alphatPhaseJayatillekeWallFunctionFvPatchScalarField::
alphatPhaseJayatillekeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    Prt_(dict.lookupOrDefault<scalar>("Prt", 0.85))
{}


alphatPhaseJayatillekeWallFunctionFvPatchScalarField::
alphatPhaseJayatillekeWallFunctionFvPatchScalarField
(
    const alphatPhaseJayatillekeWallFunctionFvPatchScalarField& awfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(awfpsf, iF),
    Prt_(awfpsf.Prt_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField>
alphatPhaseJayatillekeWallFunctionFvPatchScalarField::calcAlphat
(
    const scalarField& prevAlphat
) const
{
    // Lookup the fluid model
    const phaseSystem& fluid =
        db().lookupObject<phaseSystem>("phaseProperties");

    const phaseModel& phase
    (
        fluid.phases()[internalField().group()]
    );

    const label patchi = patch().index();

    // Retrieve turbulence properties from model
    const phaseCompressibleMomentumTransportModel& turbModel =
        db().lookupObject<phaseCompressibleMomentumTransportModel>
        (
            IOobject::groupName(momentumTransportModel::typeName, phase.name())
        );

    const nutWallFunctionFvPatchScalarField& nutw =
        nutWallFunctionFvPatchScalarField::nutw(turbModel, patchi);

    const scalar Cmu25 = pow025(nutw.Cmu());

    const scalarField& y = turbModel.y()[patchi];

    const tmp<scalarField> tmuw = turbModel.mu(patchi);
    const scalarField& muw = tmuw();

    const tmp<scalarField> talphaw = phase.thermo().alpha(patchi);
    const scalarField& alphaw = talphaw();

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    const fvPatchScalarField& kw = k.boundaryField()[patchi];

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));
    const scalarField magGradUw(mag(Uw.snGrad()));

    const fvPatchScalarField& rhow = turbModel.rho().boundaryField()[patchi];
    const fvPatchScalarField& hew =
        phase.thermo().he().boundaryField()[patchi];

    const fvPatchScalarField& Tw =
        phase.thermo().T().boundaryField()[patchi];

    const scalarField Tp(Tw.patchInternalField());

    // Heat flux [W/m^2] - lagging alphatw
    const scalarField qDot
    (
        (prevAlphat + alphaw)*hew.snGrad()
    );

    const scalarField uTau(Cmu25*sqrt(kw));

    const scalarField yPlus(uTau*y/(muw/rhow));

    const scalarField Pr(muw/alphaw);

    // Molecular-to-turbulent Prandtl number ratio
    const scalarField Prat(Pr/Prt_);

    // Thermal sublayer thickness
    const scalarField P(this->Psmooth(Prat));

    const scalarField yPlusTherm(this->yPlusTherm(nutw, P, Prat));

    tmp<scalarField> talphatConv(new scalarField(this->size()));
    scalarField& alphatConv = talphatConv.ref();

    // Populate boundary values
    forAll(alphatConv, facei)
    {
        // Evaluate new effective thermal diffusivity
        scalar alphaEff = 0.0;
        if (yPlus[facei] < yPlusTherm[facei])
        {
            const scalar A = qDot[facei]*rhow[facei]*uTau[facei]*y[facei];

            const scalar B = qDot[facei]*Pr[facei]*yPlus[facei];

            const scalar C =
                Pr[facei]*0.5*rhow[facei]*uTau[facei]*sqr(magUp[facei]);

            alphaEff = A/(B + C + vSmall);
        }
        else
        {
            const scalar A = qDot[facei]*rhow[facei]*uTau[facei]*y[facei];

            const scalar B =
                qDot[facei]*Prt_
               *(1.0/nutw.kappa()*log(nutw.E()*yPlus[facei]) + P[facei]);

            const scalar magUc =
                uTau[facei]/nutw.kappa()
               *log(nutw.E()*yPlusTherm[facei]) - mag(Uw[facei]);

            const scalar C =
                0.5*rhow[facei]*uTau[facei]
               *(Prt_*sqr(magUp[facei]) + (Pr[facei] - Prt_)*sqr(magUc));

            alphaEff = A/(B + C + vSmall);
        }

        // Update convective heat transfer turbulent thermal diffusivity
        alphatConv[facei] = max(0.0, alphaEff - alphaw[facei]);
    }

    return talphatConv;
}


void alphatPhaseJayatillekeWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    operator==(calcAlphat(*this));

    fixedValueFvPatchScalarField::updateCoeffs();
}


void alphatPhaseJayatillekeWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    writeEntry(os, "Prt", Prt_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphatPhaseJayatillekeWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
