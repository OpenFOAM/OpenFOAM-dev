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

#include "alphatWallBoilingWallFunctionFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"

#include "twoPhaseSystem.H"
#include "phaseSystem.H"
#include "ThermalPhaseChangePhaseSystem.H"
#include "MomentumTransferPhaseSystem.H"
#include "compressibleTurbulenceModel.H"
#include "ThermalDiffusivity.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "saturationModel.H"
#include "wallFvPatch.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatWallBoilingWallFunctionFvPatchScalarField::
alphatWallBoilingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(p, iF),
    relax_(0.5),
    AbyV_(p.size(), 0.0),
    alphatConv_(p.size(), 0.0)
{
    AbyV_ = this->patch().magSf();
    forAll(AbyV_,facei)
    {
        label faceCelli = this->patch().faceCells()[facei];
        AbyV_[facei] /= iF.mesh().V()[faceCelli];
    }
}


alphatWallBoilingWallFunctionFvPatchScalarField::
alphatWallBoilingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(p, iF, dict),
    relax_(dict.lookupOrDefault<scalar>("relax", 0.5)),
    AbyV_(p.size(), 0.0),
    alphatConv_(p.size(), 0.0)
{
    if (dict.found("alphatConv"))
    {
        alphatConv_ = scalarField("alphatConv", dict, p.size());
    }

    AbyV_ = this->patch().magSf();
    forAll(AbyV_,facei)
    {
        label faceCelli = this->patch().faceCells()[facei];
        AbyV_[facei] /= iF.mesh().V()[faceCelli];
    }
}


alphatWallBoilingWallFunctionFvPatchScalarField::
alphatWallBoilingWallFunctionFvPatchScalarField
(
    const alphatWallBoilingWallFunctionFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField
    (
        psf,
        p,
        iF,
        mapper
    ),
    relax_(psf.relax_),
    AbyV_(psf.AbyV_),
    alphatConv_(psf.alphatConv_, mapper)
{}


alphatWallBoilingWallFunctionFvPatchScalarField::
alphatWallBoilingWallFunctionFvPatchScalarField
(
    const alphatWallBoilingWallFunctionFvPatchScalarField& psf
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(psf),
    relax_(psf.relax_),
    AbyV_(psf.AbyV_),
    alphatConv_(psf.alphatConv_)
{}


alphatWallBoilingWallFunctionFvPatchScalarField::
alphatWallBoilingWallFunctionFvPatchScalarField
(
    const alphatWallBoilingWallFunctionFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(psf, iF),
    relax_(psf.relax_),
    AbyV_(psf.AbyV_),
    alphatConv_(psf.alphatConv_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void alphatWallBoilingWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Lookup the fluid model
    const ThermalPhaseChangePhaseSystem<MomentumTransferPhaseSystem
    <
        twoPhaseSystem>
    >&
        fluid
      = refCast
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

    const phaseModel& vapor(fluid.otherPhase(liquid));

    const label patchi = patch().index();

    // Retrieve turbulence properties from model
    const phaseCompressibleTurbulenceModel& turbModel = liquid.turbulence();

    const tmp<scalarField> tnutw = turbModel.nut(patchi);

    const scalar Cmu25(pow025(Cmu_));

    const scalarField& y = turbModel.y()[patchi];

    const tmp<scalarField> tmuw = turbModel.mu(patchi);
    const scalarField& muw = tmuw();

    const tmp<scalarField> talphaw = liquid.thermo().alpha(patchi);
    const scalarField& alphaw = talphaw();

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
    const scalarField Tc(Tw.patchInternalField());

    scalarField uTau(Cmu25*sqrt(kw));

    scalarField yPlus(uTau*y/(muw/rhow));

    scalarField Pr(muw/alphaw);

    // Molecular-to-turbulent Prandtl number ratio
    scalarField Prat(Pr/Prt_);

    // Thermal sublayer thickness
    scalarField P(this->Psmooth(Prat));

    scalarField yPlusTherm(this->yPlusTherm(P, Prat));

    const scalarField rhoc(rhow.patchInternalField());

    tmp<volScalarField> trhoVapor = vapor.thermo().rho();
    const volScalarField& rhoVapor = trhoVapor();
    const fvPatchScalarField& rhoVaporw =
       rhoVapor.boundaryField()[patchi];
    const scalarField rhoVaporp(rhoVaporw.patchInternalField());

    tmp<volScalarField> tCp = liquid.thermo().Cp();
    const volScalarField& Cp = tCp();
    const fvPatchScalarField& Cpw = Cp.boundaryField()[patchi];

    // Saturation temperature
    const tmp<volScalarField> tTsat =
        fluid.saturation().Tsat(liquid.thermo().p());
    const volScalarField& Tsat = tTsat();
    const fvPatchScalarField& Tsatw(Tsat.boundaryField()[patchi]);
    const scalarField Tsatc(Tsatw.patchInternalField());

    // Gravitational acceleration
    const uniformDimensionedVectorField& g =
        db().lookupObject<uniformDimensionedVectorField>("g");

    const fvPatchScalarField& pw =
        liquid.thermo().p().boundaryField()[patchi];

    const scalarField L
    (
        vapor.thermo().he(pw,Tsatc,patchi)-hew.patchInternalField()
    );

    // Liquid temperature at y+=250 is estimated from logarithmic
    // thermal wall function (Koncar, Krepper & Egorov, 2005)
    scalarField Tplus_y250(Prt_*(log(E_*250)/kappa_ + P));
    scalarField Tplus(Prt_*(log(E_*yPlus)/kappa_ + P));
    scalarField Tl(Tw - (Tplus_y250/Tplus)*(Tw - Tc));
    Tl = max(Tc - 40, min(Tc, Tl));

    // Nucleation site density:
    // Reformulation of Lemmert & Chawla (Egorov & Menter, 2004)
    const scalarField N
    (
        0.8*9.922e5*pow(max((Tw - Tsatw)/10, scalar(0)), 1.805)
    );

    // Bubble departure diameter:
    // Tolubinski and Kostanchuk (1970)
    const scalarField Tsub(max(Tsatw - Tl, scalar(0)));
    const scalarField Ddep
    (
        max(min(0.0006*exp(-Tsub/45), scalar(0.0014)), scalar(1e-6))
    );

    // Bubble departure frequency:
    // Cole (1960)
    const scalarField F
    (
        sqrt
        (
            4*mag(g).value()*(max(rhoc - rhoVaporp, scalar(0.1)))/(3*Ddep*rhow)
        )
    );

    // Area fractions:

    // Del Valle & Kenning (1985)
    const scalarField Ja(rhoc*Cpw*Tsub/(rhoVaporp*L));
    const scalarField Al(4.8*exp(-Ja/80));

    // Liquid phase fraction at the wall
    const scalarField liquidw(liquid.boundaryField()[patchi]);

    // Damp boiling at high void fractions.
    const scalarField W(min(liquidw/0.2, scalar(0.1)));

    const scalarField A2(W*min(M_PI*sqr(Ddep)*N*Al/4, scalar(1)));
    const scalarField A1(max(1 - A2, scalar(1e-4)));
    const scalarField A2E(W*min(M_PI*sqr(Ddep)*N*Al/4, scalar(5)));

    // Wall evaporation heat flux [kg/s3 = J/m2s]
    const scalarField Qe((1.0/6.0)*A2E*Ddep*rhoVaporw*F*L);

    // Volumetric mass source in the near wall cell due to the wall boiling
    dmdt_ = (1 - relax_)*dmdt_ + relax_*Qe*AbyV_/L;

    // Volumetric source in the near wall cell due to the wall boiling
    mDotL_ = dmdt_*L;

    // Quenching heat transfer coefficient
    const scalarField hQ
    (
        2*(alphaw*Cpw)*F*sqrt((0.8/F)/(M_PI*alphaw/rhow))
    );

    // Quenching heat flux
    const scalarField Qq(A2*hQ*max(Tw - Tl, scalar(0)));

    // Convective heat flux
    alphatConv_ = calcAlphat(alphatConv_);
    //const scalarField Qc(A1*(alphatConv_ + alphaw)*hew.snGrad());

    // Effective thermal diffusivity that corresponds to the calculated
    // convective, quenching and evaporative heat fluxes

    operator==
    (
        A1*alphatConv_ + (Qq + Qe)/max(liquidw*hew.snGrad(), scalar(1e-16))
    );

    if(debug)
    {
        Info<< "  L: " << gMin(L) << " - " << gMax(L) << endl;
        Info<< "  Tl: " << gMin(Tl) << " - " << gMax(Tl) << endl;
        Info<< "  N: " << gMin(N) << " - " << gMax(N) << endl;
        Info<< "  Ddep: " << gMin(Ddep) << " - " << gMax(Ddep) << endl;
        Info<< "  F: " << gMin(F) << " - " << gMax(F) << endl;
        Info<< "  Al: " << gMin(Al) << " - " << gMax(Al) << endl;
        Info<< "  A1: " << gMin(A1) << " - " << gMax(A1) << endl;
        Info<< "  A2: " << gMin(A2) << " - " << gMax(A2) << endl;
        Info<< "  A2E: " << gMin(A2E) << " - " << gMax(A2E) << endl;
        Info<< "  dmdtW: " << gMin(dmdt_) << " - " << gMax(dmdt_) << endl;
        const scalarField Qc(A1*(alphatConv_ + alphaw)*hew.snGrad());
        Info<< "  Qc: " << gMin(Qc) << " - " << gMax(Qc) << endl;
        Info<< "  Qq: " << gMin(Qq) << " - " << gMax(Qq) << endl;
        Info<< "  Qe: " << gMin(Qe) << " - " << gMax(Qe) << endl;
        const scalarField QEff(liquidw*(*this + alphaw)*hew.snGrad());
        Info<< "  QEff: " << gMin(QEff) << " - " << gMax(QEff) << endl;
        Info<< "  alphat: " << gMin(*this) << " - " << gMax(*this) << endl;
        Info<< "  alphatConv: " << gMin(alphatConv_)
            << " - " << gMax(alphatConv_) << endl;
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void alphatWallBoilingWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("relax") << relax_ << token::END_STATEMENT << nl;
    dmdt_.writeEntry("dmdt", os);
    alphatConv_.writeEntry("alphatConv", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphatWallBoilingWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
