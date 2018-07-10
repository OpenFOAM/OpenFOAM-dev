/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "phaseSystem.H"
#include "compressibleTurbulenceModel.H"
#include "ThermalDiffusivity.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "saturationModel.H"
#include "wallFvPatch.H"
#include "uniformDimensionedFields.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* Foam::NamedEnum
<
    Foam::compressible::
    alphatWallBoilingWallFunctionFvPatchScalarField::phaseType,
    2
>::names[] =
{
    "vapor",
    "liquid"
};

const Foam::NamedEnum
<
    Foam::compressible::
    alphatWallBoilingWallFunctionFvPatchScalarField::phaseType,
    2
>
Foam::compressible::
alphatWallBoilingWallFunctionFvPatchScalarField::phaseTypeNames_;


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
    otherPhaseName_("vapor"),
    phaseType_(liquidPhase),
    relax_(0.5),
    AbyV_(p.size(), 0),
    alphatConv_(p.size(), 0),
    dDep_(p.size(), 1e-5),
    qq_(p.size(), 0),
    partitioningModel_(nullptr),
    nucleationSiteModel_(nullptr),
    departureDiamModel_(nullptr),
    departureFreqModel_(nullptr)
{
    AbyV_ = this->patch().magSf();
    forAll(AbyV_, facei)
    {
        const label faceCelli = this->patch().faceCells()[facei];
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
    otherPhaseName_(dict.lookup("otherPhase")),
    phaseType_(phaseTypeNames_.read(dict.lookup("phaseType"))),
    relax_(dict.lookupOrDefault<scalar>("relax", 0.5)),
    AbyV_(p.size(), 0),
    alphatConv_(p.size(), 0),
    dDep_(p.size(), 1e-5),
    qq_(p.size(), 0),
    partitioningModel_(nullptr),
    nucleationSiteModel_(nullptr),
    departureDiamModel_(nullptr),
    departureFreqModel_(nullptr)
{

    // Check that otherPhaseName != this phase
    if (internalField().group() == otherPhaseName_)
    {
        FatalErrorInFunction
            << "otherPhase should be the name of the vapor phase that "
            << "corresponds to the liquid base of vice versa" << nl
            << "This phase: " << internalField().group() << nl
            << "otherPhase: " << otherPhaseName_
            << abort(FatalError);
    }

    switch (phaseType_)
    {
        case vaporPhase:
        {
            partitioningModel_ =
                wallBoilingModels::partitioningModel::New
                (
                    dict.subDict("partitioningModel")
                );

            dmdt_ = 0;

            break;
        }
        case liquidPhase:
        {
            partitioningModel_ =
                wallBoilingModels::partitioningModel::New
                (
                    dict.subDict("partitioningModel")
                );

            nucleationSiteModel_ =
                wallBoilingModels::nucleationSiteModel::New
                (
                    dict.subDict("nucleationSiteModel")
                );

            departureDiamModel_ =
                wallBoilingModels::departureDiameterModel::New
                (
                    dict.subDict("departureDiamModel")
                );

            departureFreqModel_ =
                wallBoilingModels::departureFrequencyModel::New
                (
                    dict.subDict("departureFreqModel")
                );

            if (dict.found("dDep"))
            {
                dDep_ = scalarField("dDep", dict, p.size());
            }

            if (dict.found("qQuenching"))
            {
                qq_ = scalarField("qQuenching", dict, p.size());
            }

            break;
        }
    }

    if (dict.found("alphatConv"))
    {
        alphatConv_ = scalarField("alphatConv", dict, p.size());
    }

    AbyV_ = this->patch().magSf();
    forAll(AbyV_, facei)
    {
        const label faceCelli = this->patch().faceCells()[facei];
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
    otherPhaseName_(psf.otherPhaseName_),
    phaseType_(psf.phaseType_),
    relax_(psf.relax_),
    AbyV_(psf.AbyV_),
    alphatConv_(psf.alphatConv_, mapper),
    dDep_(psf.dDep_, mapper),
    qq_(psf.qq_, mapper),
    partitioningModel_(psf.partitioningModel_),
    nucleationSiteModel_(psf.nucleationSiteModel_),
    departureDiamModel_(psf.departureDiamModel_),
    departureFreqModel_(psf.departureFreqModel_)
{}


alphatWallBoilingWallFunctionFvPatchScalarField::
alphatWallBoilingWallFunctionFvPatchScalarField
(
    const alphatWallBoilingWallFunctionFvPatchScalarField& psf
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(psf),
    otherPhaseName_(psf.otherPhaseName_),
    phaseType_(psf.phaseType_),
    relax_(psf.relax_),
    AbyV_(psf.AbyV_),
    alphatConv_(psf.alphatConv_),
    dDep_(psf.dDep_),
    qq_(psf.qq_),
    partitioningModel_(psf.partitioningModel_),
    nucleationSiteModel_(psf.nucleationSiteModel_),
    departureDiamModel_(psf.departureDiamModel_),
    departureFreqModel_(psf.departureFreqModel_)
{}


alphatWallBoilingWallFunctionFvPatchScalarField::
alphatWallBoilingWallFunctionFvPatchScalarField
(
    const alphatWallBoilingWallFunctionFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(psf, iF),
    otherPhaseName_(psf.otherPhaseName_),
    phaseType_(psf.phaseType_),
    relax_(psf.relax_),
    AbyV_(psf.AbyV_),
    alphatConv_(psf.alphatConv_),
    dDep_(psf.dDep_),
    qq_(psf.qq_),
    partitioningModel_(psf.partitioningModel_),
    nucleationSiteModel_(psf.nucleationSiteModel_),
    departureDiamModel_(psf.departureDiamModel_),
    departureFreqModel_(psf.departureFreqModel_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool alphatWallBoilingWallFunctionFvPatchScalarField::
activePhasePair(const phasePairKey& phasePair) const
{
    if (phasePair == phasePairKey(otherPhaseName_, internalField().group()))
    {
        return true;
    }
    else
    {
        return false;
    }
}

const scalarField& alphatWallBoilingWallFunctionFvPatchScalarField::
dmdt(const phasePairKey& phasePair) const
{
    if (activePhasePair(phasePair))
    {
        return dmdt_;
    }
    else
    {
        FatalErrorInFunction
            << " dmdt requested for invalid phasePair!"
            << abort(FatalError);

        return dmdt_;
    }
}

const scalarField& alphatWallBoilingWallFunctionFvPatchScalarField::
mDotL(const phasePairKey& phasePair) const
{
    if (activePhasePair(phasePair))
    {
        return mDotL_;
    }
    else
    {
        FatalErrorInFunction
            << " mDotL requested for invalid phasePair!"
            << abort(FatalError);

        return mDotL_;
    }
}

void alphatWallBoilingWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Check that partitioningModel has been constructed
    if (!partitioningModel_.valid())
    {
        FatalErrorInFunction
            << "partitioningModel has not been constructed!"
            << abort(FatalError);
    }

    // Lookup the fluid model
    const phaseSystem& fluid =
        refCast<const phaseSystem>
        (
            db().lookupObject<phaseSystem>("phaseProperties")
        );

    const saturationModel& satModel =
        db().lookupObject<saturationModel>("saturationModel");

    const label patchi = patch().index();

    switch (phaseType_)
    {
        case vaporPhase:
        {
            const phaseModel& vapor
            (
                fluid.phases()[internalField().group()]
            );

            // Vapor Liquid phase fraction at the wall
            const scalarField vaporw(vapor.boundaryField()[patchi]);

            // NOTE! Assumes 1-thisPhase for liquid fraction in
            // multiphase simulations
            const scalarField fLiquid
            (
                partitioningModel_->fLiquid(1-vaporw)
            );

            operator==
            (
                calcAlphat(*this)*(1 - fLiquid)/max(vaporw, scalar(1e-8))
            );
            break;
        }
        case liquidPhase:
        {
            // Check that nucleationSiteModel has been constructed
            if (!nucleationSiteModel_.valid())
            {
                FatalErrorInFunction
                    << "nucleationSiteModel has not been constructed!"
                    << abort(FatalError);
            }

            // Check that departureDiameterModel has been constructed
            if (!departureDiamModel_.valid())
            {
                FatalErrorInFunction
                    << "departureDiameterModel has not been constructed!"
                    << abort(FatalError);
            }

            // Check that nucleationSiteModel has been constructed
            if (!departureFreqModel_.valid())
            {
                FatalErrorInFunction
                    << "departureFrequencyModel has not been constructed!"
                    << abort(FatalError);
            }

            const phaseModel& liquid
            (
                fluid.phases()[internalField().group()]
            );

            const phaseModel& vapor(fluid.phases()[otherPhaseName_]);

            // Retrieve turbulence properties from models
            const phaseCompressibleTurbulenceModel& turbModel =
                db().lookupObject<phaseCompressibleTurbulenceModel>
                (
                    IOobject::groupName
                    (
                        turbulenceModel::propertiesName,
                        liquid.name()
                    )
                );
            const phaseCompressibleTurbulenceModel& vaporTurbModel =
                db().lookupObject<phaseCompressibleTurbulenceModel>
                (
                    IOobject::groupName
                    (
                        turbulenceModel::propertiesName,
                        vapor.name()
                    )
                );

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

            const fvPatchVectorField& Uw =
                turbModel.U().boundaryField()[patchi];
            const scalarField magUp(mag(Uw.patchInternalField() - Uw));
            const scalarField magGradUw(mag(Uw.snGrad()));

            const fvPatchScalarField& rhow =
                turbModel.rho().boundaryField()[patchi];
            const fvPatchScalarField& hew =
            liquid.thermo().he().boundaryField()[patchi];

            const fvPatchScalarField& Tw =
                liquid.thermo().T().boundaryField()[patchi];
            const scalarField Tc(Tw.patchInternalField());

            const scalarField uTau(Cmu25*sqrt(kw));

            const scalarField yPlus(uTau*y/(muw/rhow));

            const scalarField Pr(muw/alphaw);

            // Molecular-to-turbulent Prandtl number ratio
            const scalarField Prat(Pr/Prt_);

            // Thermal sublayer thickness
            const scalarField P(this->Psmooth(Prat));

            const scalarField yPlusTherm(this->yPlusTherm(P, Prat));

            const fvPatchScalarField& rhoLiquidw =
                turbModel.rho().boundaryField()[patchi];

            const fvPatchScalarField& rhoVaporw =
                vaporTurbModel.rho().boundaryField()[patchi];

            tmp<volScalarField> tCp = liquid.thermo().Cp();
            const volScalarField& Cp = tCp();
            const fvPatchScalarField& Cpw = Cp.boundaryField()[patchi];

            // Saturation temperature
            const tmp<volScalarField> tTsat =
                satModel.Tsat(liquid.thermo().p());

            const volScalarField& Tsat = tTsat();
            const fvPatchScalarField& Tsatw(Tsat.boundaryField()[patchi]);
            const scalarField Tsatc(Tsatw.patchInternalField());

            const fvPatchScalarField& pw =
                liquid.thermo().p().boundaryField()[patchi];

            const scalarField L
            (
                vapor.thermo().he(pw, Tsatc, patchi) - hew.patchInternalField()
            );

            // Liquid phase fraction at the wall
            const scalarField liquidw(liquid.boundaryField()[patchi]);

            const scalarField fLiquid(partitioningModel_->fLiquid(liquidw));

            // Convective thermal diffusivity
            alphatConv_ = calcAlphat(alphatConv_);

            for (label i=0; i<10; i++)
            {
                // Liquid temperature at y+=250 is estimated from logarithmic
                // thermal wall function (Koncar, Krepper & Egorov, 2005)
                const scalarField Tplus_y250(Prt_*(log(E_*250)/kappa_ + P));

                const scalarField Tplus(Prt_*(log(E_*yPlus)/kappa_ + P));
                scalarField Tl(Tw - (Tplus_y250/Tplus)*(Tw - Tc));
                Tl = max(Tc - 40, Tl);

                // Nucleation site density:
                const scalarField N
                (
                    nucleationSiteModel_->N
                    (
                        liquid,
                        vapor,
                        patchi,
                        Tl,
                        Tsatw,
                        L
                    )
                );

                // Bubble departure diameter:
                dDep_ = departureDiamModel_->dDeparture
                (
                    liquid,
                    vapor,
                    patchi,
                    Tl,
                    Tsatw,
                    L
                );

                // Bubble departure frequency:
                const scalarField fDep
                (
                    departureFreqModel_->fDeparture
                    (
                        liquid,
                        vapor,
                        patchi,
                        dDep_
                    )
                );

                // Area fractions:

                // Del Valle & Kenning (1985)
                const scalarField Ja
                (
                    rhoLiquidw*Cpw*(Tsatw - Tl)/(rhoVaporw*L)
                );

                const scalarField Al
                (
                    fLiquid*4.8*exp( min(-Ja/80, log(vGreat)))
                );

                const scalarField A2(min(pi*sqr(dDep_)*N*Al/4, scalar(1)));
                const scalarField A1(max(1 - A2, scalar(1e-4)));
                const scalarField A2E(min(pi*sqr(dDep_)*N*Al/4, scalar(5)));

                // Volumetric mass source in the near wall cell due to the
                // wall boiling
                dmdt_ =
                    (1 - relax_)*dmdt_
                  + relax_*(1.0/6.0)*A2E*dDep_*rhoVaporw*fDep*AbyV_;

                // Volumetric source in the near wall cell due to the wall
                // boiling
                mDotL_ = dmdt_*L;

                // Quenching heat transfer coefficient
                const scalarField hQ
                (
                    2*(alphaw*Cpw)*fDep*sqrt((0.8/fDep)/(pi*alphaw/rhow))
                );

                // Quenching heat flux
                qq_ = (A2*hQ*max(Tw - Tl, scalar(0)));

                // Effective thermal diffusivity that corresponds to the
                // calculated convective, quenching and evaporative heat fluxes

                operator==
                (
                    (
                        A1*alphatConv_
                      + (qq_ + qe())/max(hew.snGrad(), scalar(1e-16))
                    )
                   /max(liquidw, scalar(1e-8))
                );

                scalarField TsupPrev(max((Tw - Tsatw), scalar(0)));
                const_cast<fvPatchScalarField&>(Tw).evaluate();
                scalarField TsupNew(max((Tw - Tsatw), scalar(0)));

                scalar maxErr(max(mag(TsupPrev - TsupNew)));

                if (debug)
                {
                    const scalarField qc
                    (
                        fLiquid*A1*(alphatConv_ + alphaw)*hew.snGrad()
                    );

                    const scalarField qEff
                    (
                        liquidw*(*this + alphaw)*hew.snGrad()
                    );

                    Info<< "  L: " << gMin(L) << " - " << gMax(L) << endl;
                    Info<< "  Tl: " << gMin(Tl) << " - " << gMax(Tl) << endl;
                    Info<< "  N: " << gMin(N) << " - " << gMax(N) << endl;
                    Info<< "  dDep_: " << gMin(dDep_) << " - "
                        << gMax(dDep_) << endl;
                    Info<< "  fDep: " << gMin(fDep) << " - "
                        << gMax(fDep) << endl;
                    Info<< "  Al: " << gMin(Al) << " - " << gMax(Al) << endl;
                    Info<< "  A1: " << gMin(A1) << " - " << gMax(A1) << endl;
                    Info<< "  A2: " << gMin(A2) << " - " << gMax(A2) << endl;
                    Info<< "  A2E: " << gMin(A2E) << " - "
                        << gMax(A2E) << endl;
                    Info<< "  dmdtW: " << gMin(dmdt_) << " - "
                        << gMax(dmdt_) << endl;
                    Info<< "  qc: " << gMin(qc) << " - " << gMax(qc) << endl;
                    Info<< "  qq: " << gMin(fLiquid*qq()) << " - "
                        << gMax(fLiquid*qq()) << endl;
                    Info<< "  qe: " << gMin(fLiquid*qe()) << " - "
                        << gMax(fLiquid*qe()) << endl;
                    Info<< "  qEff: " << gMin(qEff) << " - "
                        << gMax(qEff) << endl;
                    Info<< "  alphat: " << gMin(*this) << " - "
                        << gMax(*this) << endl;
                    Info<< "  alphatConv: " << gMin(alphatConv_)
                        << " - " << gMax(alphatConv_) << endl;
                }

                if (maxErr < 1e-1)
                {
                    if (i > 0)
                    {
                        Info<< "Wall boiling wall function iterations: "
                            << i + 1 << endl;
                    }
                    break;
                }

            }
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown phase type. Valid types are: "
                << phaseTypeNames_ << nl << exit(FatalError);
        }
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void alphatWallBoilingWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);

    os.writeKeyword("phaseType") << phaseTypeNames_[phaseType_]
        << token::END_STATEMENT << nl;

    os.writeKeyword("relax") << relax_ << token::END_STATEMENT << nl;

    switch (phaseType_)
    {
        case vaporPhase:
        {
            os.writeKeyword("partitioningModel") << nl;
            os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
            partitioningModel_->write(os);
            os << decrIndent << indent << token::END_BLOCK << nl;
            break;
        }
        case liquidPhase:
        {
            os.writeKeyword("partitioningModel") << nl;
            os << indent << token::BEGIN_BLOCK << incrIndent << nl;
            partitioningModel_->write(os);
            os << decrIndent << indent << token::END_BLOCK << nl;

            os.writeKeyword("nucleationSiteModel") << nl;
            os << indent << token::BEGIN_BLOCK << incrIndent << nl;
            nucleationSiteModel_->write(os);
            os << decrIndent << indent << token::END_BLOCK << nl;

            os.writeKeyword("departureDiamModel") << nl;
            os << indent << token::BEGIN_BLOCK << incrIndent << nl;
            departureDiamModel_->write(os);
            os << decrIndent << indent << token::END_BLOCK << nl;

            os.writeKeyword("departureFreqModel") << nl;
            os << indent << token::BEGIN_BLOCK << incrIndent << nl;
            departureFreqModel_->write(os);
            os << decrIndent << indent << token::END_BLOCK << nl;

            break;
        }
    }

    os.writeKeyword("otherPhase") << otherPhaseName_ << token::END_STATEMENT
    << nl;
    dmdt_.writeEntry("dmdt", os);
    dDep_.writeEntry("dDep", os);
    qq_.writeEntry("qQuenching", os);
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
