/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2022 OpenFOAM Foundation
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
#include "phaseSystem.H"
#include "heatTransferPhaseSystem.H"
#include "compressibleMomentumTransportModels.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "interfaceSaturationTemperatureModel.H"
#include "rhoMulticomponentThermo.H"
#include "addToRunTimeSelectionTable.H"

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
    alphatPhaseChangeWallFunctionFvPatchScalarField(p, iF),
    phaseType_(liquidPhase),
    AbyV_(p.size(), 0),
    alphatConv_(p.size(), 0),
    dDep_(p.size(), 1e-5),
    fDep_(p.size(), 0),
    N_(p.size(), 0),
    fLiquid_(p.size(), 0),
    qq_(p.size(), 0),
    qe_(p.size(), 0),
    tau_(0.8),
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
    alphatPhaseChangeWallFunctionFvPatchScalarField(p, iF, dict),
    phaseType_(phaseTypeNames_.read(dict.lookup("phaseType"))),
    AbyV_(p.size(), 0),
    alphatConv_(p.size(), 0),
    dDep_(p.size(), 1e-5),
    fDep_(p.size(), 0),
    N_(p.size(), 0),
    fLiquid_(p.size(), 0),
    qq_(p.size(), 0),
    qe_(p.size(), 0),
    tau_(dict.lookupOrDefault<scalar>("bubbleWaitingTimeRatio", 0.8)),
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
            << "corresponds to the liquid base or vice versa" << nl
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

            dmdtf_ = 0;

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

            if (dict.found("dDeparture"))
            {
                dDep_ = scalarField("dDeparture", dict, p.size());
            }

            if (dict.found("depFrequency"))
            {
                fDep_ = scalarField("depFrequency", dict, p.size());
            }

            if (dict.found("nucSiteDensity"))
            {
                N_ = scalarField("nucSiteDensity", dict, p.size());
            }

            if (dict.found("wallLiquidFraction"))
            {
                fLiquid_ = scalarField("wallLiquidFraction", dict, p.size());
            }

            if (dict.found("qQuenching"))
            {
                qq_ = scalarField("qQuenching", dict, p.size());
            }

            if (dict.found("qEvaporation"))
            {
                qq_ = scalarField("qEvaporation", dict, p.size());
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
    alphatPhaseChangeWallFunctionFvPatchScalarField
    (
        psf,
        p,
        iF,
        mapper
    ),
    phaseType_(psf.phaseType_),
    AbyV_(mapper(psf.AbyV_)),
    alphatConv_(mapper(psf.alphatConv_)),
    dDep_(mapper(psf.dDep_)),
    fDep_(mapper(psf.fDep_)),
    N_(mapper(psf.N_)),
    fLiquid_(mapper(psf.fLiquid_)),
    qq_(mapper(psf.qq_)),
    qe_(mapper(psf.qe_)),
    tau_(psf.tau_),
    partitioningModel_(psf.partitioningModel_, false),
    nucleationSiteModel_(psf.nucleationSiteModel_, false),
    departureDiamModel_(psf.departureDiamModel_, false),
    departureFreqModel_(psf.departureFreqModel_, false)
{}


alphatWallBoilingWallFunctionFvPatchScalarField::
alphatWallBoilingWallFunctionFvPatchScalarField
(
    const alphatWallBoilingWallFunctionFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphatPhaseChangeWallFunctionFvPatchScalarField(psf, iF),
    phaseType_(psf.phaseType_),
    AbyV_(psf.AbyV_),
    alphatConv_(psf.alphatConv_),
    dDep_(psf.dDep_),
    fDep_(psf.fDep_),
    N_(psf.N_),
    fLiquid_(psf.fLiquid_),
    qq_(psf.qq_),
    qe_(psf.qe_),
    tau_(psf.tau_),
    partitioningModel_(psf.partitioningModel_, false),
    nucleationSiteModel_(psf.nucleationSiteModel_, false),
    departureDiamModel_(psf.departureDiamModel_, false),
    departureFreqModel_(psf.departureFreqModel_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void alphatWallBoilingWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    alphatPhaseChangeWallFunctionFvPatchScalarField::autoMap(m);

    m(AbyV_, AbyV_);
    m(alphatConv_, alphatConv_);
    m(dDep_, dDep_);
    m(fDep_, fDep_);
    m(N_, N_);
    m(fLiquid_, fLiquid_);
    m(qq_, qq_);
    m(qe_, qe_);
}


void alphatWallBoilingWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    alphatPhaseChangeWallFunctionFvPatchScalarField::rmap(ptf, addr);

    const alphatWallBoilingWallFunctionFvPatchScalarField& tiptf =
        refCast<const alphatWallBoilingWallFunctionFvPatchScalarField>(ptf);

    AbyV_.rmap(tiptf.AbyV_, addr);
    alphatConv_.rmap(tiptf.alphatConv_, addr);
    dDep_.rmap(tiptf.dDep_, addr);
    fDep_.rmap(tiptf.fDep_, addr);
    N_.rmap(tiptf.N_, addr);
    fLiquid_.rmap(tiptf.fLiquid_, addr);
    qq_.rmap(tiptf.qq_, addr);
    qe_.rmap(tiptf.qe_, addr);
}


void alphatWallBoilingWallFunctionFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    alphatPhaseChangeWallFunctionFvPatchScalarField::reset(ptf);

    const alphatWallBoilingWallFunctionFvPatchScalarField& tiptf =
        refCast<const alphatWallBoilingWallFunctionFvPatchScalarField>(ptf);

    AbyV_.reset(tiptf.AbyV_);
    alphatConv_.reset(tiptf.alphatConv_);
    dDep_.reset(tiptf.dDep_);
    fDep_.reset(tiptf.fDep_);
    N_.reset(tiptf.N_);
    fLiquid_.reset(tiptf.fLiquid_);
    qq_.reset(tiptf.qq_);
    qe_.reset(tiptf.qe_);
}


void alphatWallBoilingWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Lookup the fluid model
    const phaseSystem& fluid =
        db().lookupObject<phaseSystem>(phaseSystem::propertiesName);

    const word volatileSpecie(fluid.lookupOrDefault<word>("volatile", "none"));

    const label patchi = patch().index();

    switch (phaseType_)
    {
        case vaporPhase:
        {
            const phaseModel& vapor = fluid.phases()[internalField().group()];

            // Vapor phase fraction at the wall
            const scalarField& vaporw = vapor.boundaryField()[patchi];

            // Partitioning
            // NOTE! Assumes that there is only only one liquid phase and all
            // other phases are vapor

            const phaseModel& liquid = fluid.phases()[otherPhaseName_];
            const scalarField& liquidw = liquid.boundaryField()[patchi];
            fLiquid_ = partitioningModel_->fLiquid(liquidw);

            operator==
            (
                calcAlphat(*this)*(vaporw/(1 - liquidw + small) )
               *(1 - fLiquid_)/max(vaporw, scalar(1e-8))
            );
            break;
        }
        case liquidPhase:
        {
            const phaseModel& liquid = fluid.phases()[internalField().group()];
            const phaseModel& vapor = fluid.phases()[otherPhaseName_];

            const phaseInterface interface(vapor, liquid);

            if
            (
                fluid.foundInterfacialModel
                <
                    interfaceSaturationTemperatureModel
                >(interface)
            )
            {
                // Retrieve turbulence properties from models
                const phaseCompressible::momentumTransportModel& turbModel
                  = db().lookupType<phaseCompressible::momentumTransportModel>
                    (
                        liquid.name()
                    );
                const phaseCompressible::momentumTransportModel& vaporTurbModel
                  = db().lookupType<phaseCompressible::momentumTransportModel>
                    (
                        vapor.name()
                    );

                const nutWallFunctionFvPatchScalarField& nutw =
                    nutWallFunctionFvPatchScalarField::nutw(turbModel, patchi);

                const scalar Cmu25(pow025(nutw.Cmu()));

                const scalarField& y = turbModel.y()[patchi];

                const tmp<scalarField> tnuw = turbModel.nu(patchi);
                const scalarField& nuw = tnuw();

                const rhoThermo& lThermo = liquid.thermo();

                const tmp<scalarField> talphaw
                (
                    lThermo.kappa().boundaryField()[patchi]
                   /lThermo.Cp().boundaryField()[patchi]
                );
                const scalarField& alphaw = talphaw();

                const tmp<volScalarField> tk = turbModel.k();
                const volScalarField& k = tk();
                const fvPatchScalarField& kw = k.boundaryField()[patchi];

                const fvPatchVectorField& Uw =
                    turbModel.U().boundaryField()[patchi];
                const scalarField magUp(mag(Uw.patchInternalField() - Uw));
                const scalarField magGradUw(mag(Uw.snGrad()));

                const fvPatchScalarField& rhoLiquidw =
                    turbModel.rho().boundaryField()[patchi];

                const fvPatchScalarField& rhoVaporw =
                    vaporTurbModel.rho().boundaryField()[patchi];

                const fvPatchScalarField& hew =
                    lThermo.he().boundaryField()[patchi];

                const fvPatchScalarField& Tw =
                    lThermo.T().boundaryField()[patchi];

                const scalarField Tc(Tw.patchInternalField());

                const scalarField uTau(Cmu25*sqrt(kw));

                const scalarField yPlus(uTau*y/nuw);

                const scalarField Pr(rhoLiquidw*nuw/alphaw);

                // Molecular-to-turbulent Prandtl number ratio
                const scalarField Prat(Pr/Prt_);

                // Thermal sublayer thickness
                const scalarField P(this->Psmooth(Prat));

                const scalarField yPlusTherm(this->yPlusTherm(nutw, P, Prat));

                const scalarField Cpw(lThermo.Cp(Tw, patchi));

                // Saturation temperature
                const interfaceSaturationTemperatureModel& satModel =
                    fluid.lookupInterfacialModel
                    <
                        interfaceSaturationTemperatureModel
                    >(interface);
                const tmp<volScalarField> tTsat = satModel.Tsat(lThermo.p());
                const volScalarField& Tsat = tTsat();
                const fvPatchScalarField& Tsatw(Tsat.boundaryField()[patchi]);

                // Latent heat
                const scalarField L
                (
                    volatileSpecie != "none"
                  ? -refCast<const heatTransferPhaseSystem>(fluid)
                    .Li
                     (
                         interface,
                         volatileSpecie,
                         dmdtf_,
                         Tsat,
                         patch().faceCells(),
                         heatTransferPhaseSystem::latentHeatScheme::upwind
                     )
                  : -refCast<const heatTransferPhaseSystem>(fluid)
                    .L
                     (
                         interface,
                         dmdtf_,
                         Tsat,
                         patch().faceCells(),
                         heatTransferPhaseSystem::latentHeatScheme::upwind
                     )
                );

                // Liquid phase fraction at the wall
                const scalarField liquidw(liquid.boundaryField()[patchi]);

                // Partitioning
                fLiquid_ = partitioningModel_->fLiquid(liquidw);

                // Convective thermal diffusivity
                alphatConv_ = calcAlphat(alphatConv_);

                label maxIter(10);
                for (label i=0; i<maxIter; i++)
                {
                    // Liquid temperature at y+=250 is estimated from
                    // logarithmic thermal wall function (Koncar, Krepper &
                    // Egorov, 2005)
                    const scalarField Tplus_y250
                    (
                        Prt_*(log(nutw.E()*250)/nutw.kappa() + P)
                    );

                    const scalarField Tplus
                    (
                        Prt_*(log(nutw.E()*yPlus)/nutw.kappa() + P)
                    );

                    const scalarField Tl
                    (
                        max
                        (
                            Tc - 40,
                            Tw - (Tplus_y250/Tplus)*(Tw - Tc)
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
                    fDep_ =
                        departureFreqModel_->fDeparture
                        (
                            liquid,
                            vapor,
                            patchi,
                            Tl,
                            Tsatw,
                            L,
                            dDep_
                        );

                    // Nucleation site density:
                    N_ =
                        nucleationSiteModel_->N
                        (
                            liquid,
                            vapor,
                            patchi,
                            Tl,
                            Tsatw,
                            L,
                            dDep_,
                            fDep_
                        );

                    // Area fractions:

                    // Del Valle & Kenning (1985)
                    const scalarField Ja
                    (
                        rhoLiquidw*Cpw*(Tsatw - Tl)/(rhoVaporw*L)
                    );

                    const scalarField Al
                    (
                        fLiquid_*4.8*exp(min(-Ja/80, log(vGreat)))
                    );

                    scalarField A2(min(pi*sqr(dDep_)*N_*Al/4, scalar(1)));
                    const scalarField A1(max(1 - A2, scalar(1e-4)));
                    scalarField A2E(min(pi*sqr(dDep_)*N_*Al/4, scalar(5)));

                    if (volatileSpecie != "none" && !liquid.pure())
                    {
                        const volScalarField& Yvolatile =
                            liquid.Y(volatileSpecie);
                        A2E *= Yvolatile.boundaryField()[patchi];
                        A2 *= Yvolatile.boundaryField()[patchi];
                    }

                    // Volumetric mass source in the near wall cell due to the
                    // wall boiling
                    dmdtf_ =
                        (1 - relax_)*dmdtf_
                      + relax_*(1.0/6.0)*A2E*dDep_*rhoVaporw*fDep_*AbyV_;

                    // Quenching heat transfer coefficient
                    const scalarField hQ
                    (
                        2*(alphaw*Cpw)*fDep_
                       *sqrt
                        (
                            (tau_/max(fDep_, small))/(pi*alphaw/rhoLiquidw)
                        )
                    );

                    // Quenching heat flux
                    qq_ =
                        (1 - relax_)*qq_
                      + relax_*(A2*hQ*max(Tw - Tl, scalar(0)));

                    // Evaporation heat flux
                    qe_ = dmdtf_*L/AbyV_;

                    // Effective thermal diffusivity that corresponds to the
                    // calculated convective, quenching and evaporative heat
                    // fluxes

                    operator==
                    (
                        (
                            A1*alphatConv_
                          + (qq_ + qe_)/max(hew.snGrad(), scalar(1e-16))
                        )
                       /max(liquidw, scalar(1e-8))
                    );

                    const scalarField TsupPrev(max((Tw - Tsatw), scalar(0)));
                    const_cast<fvPatchScalarField&>(Tw).evaluate();
                    const scalarField TsupNew(max((Tw - Tsatw), scalar(0)));

                    const scalar maxErr(gMax(mag(TsupPrev - TsupNew)));

                    if (debug)
                    {
                        const scalarField qc
                        (
                            fLiquid_*A1*(alphatConv_ + alphaw)*hew.snGrad()
                        );

                        const scalarField qEff
                        (
                            liquidw*(*this + alphaw)*hew.snGrad()
                        );

                        Info<< "  L: " << gMin(L) << " - " << gMax(L) << endl;
                        Info<< "  Tl: " << gMin(Tl) << " - " << gMax(Tl)
                            << endl;
                        Info<< "  N: " << gMin(N_) << " - " << gMax(N_) << endl;
                        Info<< "  dDep_: " << gMin(dDep_) << " - "
                            << gMax(dDep_) << endl;
                        Info<< "  fDep: " << gMin(fDep_) << " - "
                            << gMax(fDep_) << endl;
                        Info<< "  Al: " << gMin(Al) << " - " << gMax(Al)
                            << endl;
                        Info<< "  A1: " << gMin(A1) << " - " << gMax(A1)
                            << endl;
                        Info<< "  A2: " << gMin(A2) << " - " << gMax(A2)
                            << endl;
                        Info<< "  A2E: " << gMin(A2E) << " - "
                            << gMax(A2E) << endl;
                        Info<< "  dmdtW: " << gMin(dmdtf_) << " - "
                            << gMax(dmdtf_) << endl;
                        Info<< "  qc: " << gMin(qc) << " - " << gMax(qc)
                            << endl;
                        Info<< "  qq: " << gMin(fLiquid_*qq_) << " - "
                            << gMax(fLiquid_*qq_) << endl;
                        Info<< "  qe: " << gMin(fLiquid_*qe_) << " - "
                            << gMax(fLiquid_*qe_) << endl;
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
                    else if (i == (maxIter - 1))
                    {
                        Info<< "Maximum number of wall boiling wall function "
                            << "iterations (" << maxIter << ") reached." << endl
                            << "Maximum change in wall temperature on last "
                            << "iteration: " << maxErr << endl;
                    }

                }
                break;
            }
            else
            {
                Info<< "Saturation model for interface " << interface.name()
                    << " not found. Wall boiling disabled." << endl;

                operator== (alphatConv_);
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
    alphatPhaseChangeWallFunctionFvPatchScalarField::write(os);

    writeEntry(os, "phaseType", phaseTypeNames_[phaseType_]);
    writeEntry(os, "bubbleWaitingTimeRatio", tau_);
    writeEntry(os, "alphatConv", alphatConv_);
    writeEntry(os, "dDeparture", dDep_);
    writeEntry(os, "depFrequency", fDep_);
    writeEntry(os, "nucSiteDensity", N_);
    writeEntry(os, "wallLiquidFraction", fLiquid_);
    writeEntry(os, "qQuenching", qq_);
    writeEntry(os, "qEvaporative", qe_);

    switch (phaseType_)
    {
        case vaporPhase:
        {
            writeKeyword(os, "partitioningModel") << nl;
            os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
            partitioningModel_->write(os);
            os << decrIndent << indent << token::END_BLOCK << nl;
            break;
        }
        case liquidPhase:
        {
            writeKeyword(os, "partitioningModel") << nl;
            os << indent << token::BEGIN_BLOCK << incrIndent << nl;
            partitioningModel_->write(os);
            os << decrIndent << indent << token::END_BLOCK << nl;

            writeKeyword(os, "nucleationSiteModel") << nl;
            os << indent << token::BEGIN_BLOCK << incrIndent << nl;
            nucleationSiteModel_->write(os);
            os << decrIndent << indent << token::END_BLOCK << nl;

            writeKeyword(os, "departureDiamModel") << nl;
            os << indent << token::BEGIN_BLOCK << incrIndent << nl;
            departureDiamModel_->write(os);
            os << decrIndent << indent << token::END_BLOCK << nl;

            writeKeyword(os, "departureFreqModel") << nl;
            os << indent << token::BEGIN_BLOCK << incrIndent << nl;
            departureFreqModel_->write(os);
            os << decrIndent << indent << token::END_BLOCK << nl;

            break;
        }
    }
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
