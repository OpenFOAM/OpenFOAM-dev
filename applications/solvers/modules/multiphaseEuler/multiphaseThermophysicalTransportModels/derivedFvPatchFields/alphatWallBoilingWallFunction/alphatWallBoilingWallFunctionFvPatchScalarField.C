/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
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
#include "alphatJayatillekeWallFunctionFvPatchScalarField.H"
#include "fluidThermophysicalTransportModel.H"
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
    Foam::compressible::alphatWallBoilingWallFunctionFvPatchScalarField::
        phaseType,
    2
>::names[] = {"vapor", "liquid"};

const Foam::NamedEnum
<
    Foam::compressible::alphatWallBoilingWallFunctionFvPatchScalarField::
        phaseType,
    2
>
Foam::compressible::alphatWallBoilingWallFunctionFvPatchScalarField::
    phaseTypeNames_;


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
    fixedValueFvPatchScalarField(p, iF),
    alphatPhaseChangeWallFunctionBase(),

    phaseType_(liquidPhase),
    useLiquidTemperatureWallFunction_(true),
    relax_(1),
    Prt_(0.85),
    tau_(0.8),

    partitioningModel_(nullptr),
    nucleationSiteModel_(nullptr),
    departureDiameterModel_(nullptr),
    departureFrequencyModel_(nullptr),

    fLiquid_(p.size(), 0),
    dDep_(p.size(), 1e-5),
    fDep_(p.size(), 0),
    N_(p.size(), 0),
    qq_(p.size(), 0),
    qe_(p.size(), 0),
    dmdtf_(p.size(), 0)
{}


alphatWallBoilingWallFunctionFvPatchScalarField::
alphatWallBoilingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    alphatPhaseChangeWallFunctionBase(p, iF, dict),

    phaseType_(phaseTypeNames_.read(dict.lookup("phaseType"))),
    useLiquidTemperatureWallFunction_
    (
        dict.lookupOrDefault<Switch>("useLiquidTemperatureWallFunction", true)
    ),
    relax_(dict.lookupOrDefault<scalar>("relax", 1)),
    Prt_(dict.lookupOrDefault<scalar>("Prt", 0.85)),
    tau_(dict.lookupOrDefault<scalar>("bubbleWaitingTimeRatio", 0.8)),

    partitioningModel_(nullptr),
    nucleationSiteModel_(nullptr),
    departureDiameterModel_(nullptr),
    departureFrequencyModel_(nullptr),

    fLiquid_(p.size(), 0),
    dDep_(p.size(), 1e-5),
    fDep_(p.size(), 0),
    N_(p.size(), 0),
    qq_(p.size(), 0),
    qe_(p.size(), 0),
    dmdtf_(p.size(), 0)
{
    // Sub-Models
    partitioningModel_ =
        wallBoilingModels::partitioningModel::New
        (
            dict.subDict("partitioningModel")
        );

    if (phaseType_ == liquidPhase)
    {
        nucleationSiteModel_ =
            wallBoilingModels::nucleationSiteModel::New
            (
                dict.subDict("nucleationSiteModel")
            );
        departureDiameterModel_ =
            wallBoilingModels::departureDiameterModel::New
            (
                dict.subDict("departureDiamModel")
            );
        departureFrequencyModel_ =
            wallBoilingModels::departureFrequencyModel::New
            (
                dict.subDict("departureFreqModel")
            );
    }

    // State
    if (dict.found("wallLiquidFraction"))
    {
        fLiquid_ = scalarField("wallLiquidFraction", dict, p.size());
    }

    if (phaseType_ == liquidPhase)
    {
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
        if (dict.found("qQuenching"))
        {
            qq_ = scalarField("qQuenching", dict, p.size());
        }
        if (dict.found("qEvaporation"))
        {
            qe_ = scalarField("qEvaporation", dict, p.size());
        }
        if (dict.found("dmdtf"))
        {
            dmdtf_ = scalarField("dmdtf", dict, p.size());
        }
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
    fixedValueFvPatchScalarField(psf, p, iF, mapper),
    alphatPhaseChangeWallFunctionBase(psf),

    phaseType_(psf.phaseType_),
    useLiquidTemperatureWallFunction_(psf.useLiquidTemperatureWallFunction_),
    relax_(psf.relax_),
    Prt_(psf.Prt_),
    tau_(psf.tau_),

    partitioningModel_(psf.partitioningModel_, false),
    nucleationSiteModel_(psf.nucleationSiteModel_, false),
    departureDiameterModel_(psf.departureDiameterModel_, false),
    departureFrequencyModel_(psf.departureFrequencyModel_, false),

    fLiquid_(mapper(psf.fLiquid_)),
    dDep_(mapper(psf.dDep_)),
    fDep_(mapper(psf.fDep_)),
    N_(mapper(psf.N_)),
    qq_(mapper(psf.qq_)),
    qe_(mapper(psf.qe_)),
    dmdtf_(mapper(psf.dmdtf_))
{}


alphatWallBoilingWallFunctionFvPatchScalarField::
alphatWallBoilingWallFunctionFvPatchScalarField
(
    const alphatWallBoilingWallFunctionFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(psf, iF),
    alphatPhaseChangeWallFunctionBase(psf),

    phaseType_(psf.phaseType_),
    useLiquidTemperatureWallFunction_(psf.useLiquidTemperatureWallFunction_),
    relax_(psf.relax_),
    Prt_(psf.Prt_),
    tau_(psf.tau_),

    partitioningModel_(psf.partitioningModel_, false),
    nucleationSiteModel_(psf.nucleationSiteModel_, false),
    departureDiameterModel_(psf.departureDiameterModel_, false),
    departureFrequencyModel_(psf.departureFrequencyModel_, false),

    fLiquid_(psf.fLiquid_),
    dDep_(psf.dDep_),
    fDep_(psf.fDep_),
    N_(psf.N_),
    qq_(psf.qq_),
    qe_(psf.qe_),
    dmdtf_(psf.dmdtf_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void alphatWallBoilingWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);

    m(fLiquid_, fLiquid_);
    m(dDep_, dDep_);
    m(fDep_, fDep_);
    m(N_, N_);
    m(qq_, qq_);
    m(qe_, qe_);
    m(dmdtf_, dmdtf_);
}


void alphatWallBoilingWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const alphatWallBoilingWallFunctionFvPatchScalarField& tiptf =
        refCast<const alphatWallBoilingWallFunctionFvPatchScalarField>(ptf);

    fLiquid_.rmap(tiptf.fLiquid_, addr);
    dDep_.rmap(tiptf.dDep_, addr);
    fDep_.rmap(tiptf.fDep_, addr);
    N_.rmap(tiptf.N_, addr);
    qq_.rmap(tiptf.qq_, addr);
    qe_.rmap(tiptf.qe_, addr);
    dmdtf_.rmap(tiptf.dmdtf_, addr);
}


void alphatWallBoilingWallFunctionFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    fixedValueFvPatchScalarField::reset(ptf);

    const alphatWallBoilingWallFunctionFvPatchScalarField& tiptf =
        refCast<const alphatWallBoilingWallFunctionFvPatchScalarField>(ptf);

    fLiquid_.reset(tiptf.fLiquid_);
    dDep_.reset(tiptf.dDep_);
    fDep_.reset(tiptf.fDep_);
    N_.reset(tiptf.N_);
    qq_.reset(tiptf.qq_);
    qe_.reset(tiptf.qe_);
    dmdtf_.reset(tiptf.dmdtf_);
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

            // Vapor thermophysical transport model
            const fluidThermophysicalTransportModel& vaporTtm =
                db().lookupType<fluidThermophysicalTransportModel>
                (
                    vapor.name()
                );

            // Vapor phase fraction at the wall
            const scalarField& vaporw = vapor.boundaryField()[patchi];

            // Partitioning. Note: Assumes that there is only only one liquid
            // phase and all other phases are vapor.
            const phaseModel& liquid = fluid.phases()[otherPhaseName_];
            const scalarField& liquidw = liquid.boundaryField()[patchi];
            fLiquid_ = partitioningModel_->fLiquid(liquidw);

            // Vapour thermal diffusivity
            const scalarField alphatConv
            (
                alphatJayatillekeWallFunctionFvPatchScalarField::alphat
                (
                    vaporTtm,
                    Prt_,
                    patch().index()
                )
            );

            operator==
            (
                alphatConv*(vaporw/(1 - liquidw + small) )
               *(1 - fLiquid_)/max(vaporw, scalar(1e-8))
            );

            break;
        }

        case liquidPhase:
        {
            const phaseModel& liquid = fluid.phases()[internalField().group()];
            const phaseModel& vapor = fluid.phases()[otherPhaseName_];

            const phaseInterface interface(vapor, liquid);

            // Liquid thermophysical and momentum transport models
            const fluidThermophysicalTransportModel& liquidTtm =
                db().lookupType<fluidThermophysicalTransportModel>
                (
                    liquid.name()
                );
            const compressibleMomentumTransportModel& liquidMtm =
                liquidTtm.momentumTransport();

            // Convective thermal diffusivity
            const scalarField alphatConv
            (
                alphatJayatillekeWallFunctionFvPatchScalarField::alphat
                (
                    liquidTtm,
                    Prt_,
                    patch().index()
                )
            );

            // Quit if no saturation temperature model exists for this
            // interface. Saturation modelling is considered to be the
            // fundamental enabling model for thermal mass transfers.
            if
            (
               !fluid.foundInterfacialModel
                <
                    interfaceSaturationTemperatureModel
                >
                (interface)
            )
            {
                Info<< "Saturation model for interface " << interface.name()
                    << " not found. Wall boiling disabled." << endl;

                operator==(alphatConv);

                break;
            }

            // Lookup and calculate turbulent and thermal properties
            const nutWallFunctionFvPatchScalarField& nutw =
                nutWallFunctionFvPatchScalarField::nutw(liquidMtm, patchi);

            const scalar Cmu25(pow025(nutw.Cmu()));

            const scalarField& y = liquidMtm.y()[patchi];

            const tmp<scalarField> tnuw = liquidMtm.nu(patchi);
            const scalarField& nuw = tnuw();

            const tmp<scalarField> talphaw
            (
                liquid.thermo().kappa().boundaryField()[patchi]
               /liquid.thermo().Cp().boundaryField()[patchi]
            );
            const scalarField& alphaw = talphaw();

            const tmp<volScalarField> tk = liquidMtm.k();
            const volScalarField& k = tk();
            const fvPatchScalarField& kw = k.boundaryField()[patchi];

            const fvPatchVectorField& Uw =
                liquidMtm.U().boundaryField()[patchi];
            const scalarField magUp(mag(Uw.patchInternalField() - Uw));
            const scalarField magGradUw(mag(Uw.snGrad()));

            const tmp<scalarField> trhoLiquidw = liquid.thermo().rho(patchi);
            const scalarField rhoLiquidw = trhoLiquidw();

            const tmp<scalarField> trhoVaporw = vapor.thermo().rho(patchi);
            const scalarField rhoVaporw = trhoVaporw();

            const fvPatchScalarField& hew =
                liquid.thermo().he().boundaryField()[patchi];

            const fvPatchScalarField& Tw =
                liquid.thermo().T().boundaryField()[patchi];

            const scalarField Tc(Tw.patchInternalField());

            const scalarField uTau(Cmu25*sqrt(kw));

            const scalarField yPlus(uTau*y/nuw);

            const scalarField Pr(rhoLiquidw*nuw/alphaw);

            // Molecular-to-turbulent Prandtl number ratio
            const scalarField Prat(Pr/Prt_);

            // Thermal sublayer thickness
            const scalarField P
            (
                alphatJayatillekeWallFunctionFvPatchScalarField::P(Prat)
            );
            const scalarField yPlusTherm
            (
                alphatJayatillekeWallFunctionFvPatchScalarField::yPlusTherm
                (
                    nutw,
                    P,
                    Prat
                )
            );

            const scalarField Cpw(liquid.thermo().Cp(Tw, patchi));

            // Saturation temperature
            const interfaceSaturationTemperatureModel& satModel =
                fluid.lookupInterfacialModel
                <
                    interfaceSaturationTemperatureModel
                >(interface);
            const tmp<volScalarField> tTsat =
                satModel.Tsat(liquid.thermo().p());
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

            // Iterative solution for the wall temperature
            label maxIter(10);
            for (label i=0; i<maxIter; i++)
            {
                scalarField Tl(Tc);

                if (useLiquidTemperatureWallFunction_)
                {
                    // Liquid temperature at y+=250 is estimated from the
                    // logarithmic thermal wall function of Koncar, Krepper
                    // & Egorov (2005)
                    const scalarField TyPlus250
                    (
                        Prt_*(log(nutw.E()*250)/nutw.kappa() + P)
                    );

                    const scalarField TyPlus
                    (
                        Prt_
                       *(
                            log(nutw.E()*max(yPlus, scalar(11)))
                           /nutw.kappa()
                          + P
                        )
                    );

                    Tl = Tw - (TyPlus250/TyPlus)*(Tw - Tc);
                }

                // Bubble departure diameter
                dDep_ =
                    departureDiameterModel_->dDeparture
                    (
                        liquid,
                        vapor,
                        patchi,
                        Tl,
                        Tsatw,
                        L
                    );

                // Bubble departure frequency
                fDep_ =
                    departureFrequencyModel_->fDeparture
                    (
                        liquid,
                        vapor,
                        patchi,
                        Tl,
                        Tsatw,
                        L,
                        dDep_
                    );

                // Nucleation site density
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

                // Patch area by neighbouring cell volume ratio
                const scalarField AbyV
                (
                    patch().magSf()
                   /scalarField
                    (
                        patch().boundaryMesh().mesh().V(),
                        patch().faceCells()
                    )
                );

                // Volumetric mass source in the near wall cell due to the
                // wall boiling
                dmdtf_ =
                    (1 - relax_)*dmdtf_
                  + relax_*(1.0/6.0)*A2E*dDep_*rhoVaporw*fDep_*AbyV;

                // Quenching heat transfer coefficient
                const scalarField hQ
                (
                    2*alphaw*Cpw*fDep_
                   *sqrt((tau_/max(fDep_, small))/(pi*alphaw/rhoLiquidw))
                );

                // Quenching heat flux
                qq_ =
                    (1 - relax_)*qq_
                  + relax_*(A2*hQ*max(Tw - Tl, scalar(0)));

                // Evaporation heat flux
                qe_ = dmdtf_*L/AbyV;

                // Set an effective thermal diffusivity that corresponds to the
                // calculated convective, quenching and evaporative heat fluxes
                operator==
                (
                    (
                        A1*alphatConv
                      + (qq_ + qe_)/max(hew.snGrad(), scalar(1e-16))
                    )
                   /max(liquidw, scalar(1e-8))
                );

                // Evaluate the temperature condition and estimate the
                // remaining residual error
                const scalarField TsupPrev(max((Tw - Tsatw), scalar(0)));
                const_cast<fvPatchScalarField&>(Tw).evaluate();
                const scalarField TsupNew(max((Tw - Tsatw), scalar(0)));
                const scalar maxErr(gMax(mag(TsupPrev - TsupNew)));

                if (maxErr < 1e-1)
                {
                    if (i > 0)
                    {
                        Info<< "Wall boiling wall function iterations: "
                            << i + 1 << endl;
                    }

                    break;
                }

                if (i == maxIter - 1)
                {
                    Info<< "Maximum number of wall boiling wall function "
                        << "iterations (" << maxIter << ") reached." << endl
                        << "Maximum change in wall temperature on last "
                        << "iteration: " << maxErr << endl;
                }
            }

            break;
        }
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void alphatWallBoilingWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchScalarField::write(os);
    alphatPhaseChangeWallFunctionBase::write(os);

    // Controls
    writeEntry(os, "phaseType", phaseTypeNames_[phaseType_]);
    writeEntry
    (
        os,
        "useLiquidTemperatureWallFunction",
        useLiquidTemperatureWallFunction_
    );
    writeEntry(os, "relax", relax_);
    writeEntry(os, "Prt", Prt_);
    writeEntry(os, "bubbleWaitingTimeRatio", tau_);

    // Sub-models
    writeKeyword(os, "partitioningModel") << nl;
    os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
    partitioningModel_->write(os);
    os  << decrIndent << indent << token::END_BLOCK << nl;

    if (phaseType_ == liquidPhase)
    {
        writeKeyword(os, "partitioningModel") << nl;
        os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
        partitioningModel_->write(os);
        os  << decrIndent << indent << token::END_BLOCK << nl;

        writeKeyword(os, "nucleationSiteModel") << nl;
        os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
        nucleationSiteModel_->write(os);
        os  << decrIndent << indent << token::END_BLOCK << nl;

        writeKeyword(os, "departureDiamModel") << nl;
        os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
        departureDiameterModel_->write(os);
        os  << decrIndent << indent << token::END_BLOCK << nl;

        writeKeyword(os, "departureFreqModel") << nl;
        os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
        departureFrequencyModel_->write(os);
        os  << decrIndent << indent << token::END_BLOCK << nl;
    }

    // State
    writeEntry(os, "wallLiquidFraction", fLiquid_);

    if (phaseType_ == liquidPhase)
    {
        writeEntry(os, "dDeparture", dDep_);
        writeEntry(os, "depFrequency", fDep_);
        writeEntry(os, "nucSiteDensity", N_);
        writeEntry(os, "qQuenching", qq_);
        writeEntry(os, "qEvaporative", qe_);
        writeEntry(os, "dmdtf", dmdtf_);
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
