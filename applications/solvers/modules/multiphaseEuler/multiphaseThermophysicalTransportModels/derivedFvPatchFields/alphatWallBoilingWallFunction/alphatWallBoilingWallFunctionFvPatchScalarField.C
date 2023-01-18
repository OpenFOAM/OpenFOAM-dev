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

#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "mixedFvPatchFields.H"

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

struct alphatWallBoilingWallFunctionFvPatchScalarField::properties
{
    // Data

        //- Wall function field
        const alphatWallBoilingWallFunctionFvPatchScalarField& field;

        //- Phase
        const phaseModel& phase;

        //- Other phase
        const phaseModel& otherPhase;

        //- Volume fraction
        const scalarField& alphaw;

        //- Other volume fraction
        const scalarField& otherAlphaw;

        //- Interface
        const phaseInterface interface;

        //- Phase thermophysical transport model
        const fluidThermophysicalTransportModel& ttm;

        //- Phase convective turbulent thermal diffusivity
        const scalarField alphatConv;


    //- Constructor
    properties
    (
        const alphatWallBoilingWallFunctionFvPatchScalarField& field,
        const phaseModel& phase,
        const phaseModel& otherPhase
    )
    :
        field(field),
        phase(phase),
        otherPhase(otherPhase),
        alphaw(phase.boundaryField()[patchi()]),
        otherAlphaw(otherPhase.boundaryField()[patchi()]),
        interface(phase, otherPhase),
        ttm
        (
            field.db().lookupType<fluidThermophysicalTransportModel>
            (
                phase.name()
            )
        ),
        alphatConv
        (
            alphatJayatillekeWallFunctionFvPatchScalarField::alphat
            (
                ttm,
                field.Prt_,
                patchi()
            )
        )
    {}


    // Member Fuctions

        //- Patch
        inline const fvPatch& patch() const
        {
            return field.patch();
        }

        //- Patch index
        inline label patchi() const
        {
            return patch().index();
        }
};


struct alphatWallBoilingWallFunctionFvPatchScalarField::boilingLiquidProperties
:
    public alphatWallBoilingWallFunctionFvPatchScalarField::properties
{
    // Data

        //- Name of the volatile specie
        const word volatileSpecie;

        //- Patch area by neighbouring cell volume ratio
        const scalarField AbyV;

        //- Liquid density
        const tmp<scalarField> trhoLiquidw;
        const scalarField& rhoLiquidw;

        //- Vapour density
        const tmp<scalarField> trhoVapourw;
        const scalarField& rhoVapourw;

        //- Liquid heat capacity
        const scalarField& Cpw;

        //- Liquid laminar kinematic viscosity
        const tmp<scalarField> tnuw;
        const scalarField& nuw;

        //- Liquid laminar thermal diffusivity
        const scalarField kappaByCp;

        //- Liquid viscosity wall function
        const nutWallFunctionFvPatchScalarField& nutw;

        //- Dimensionless wall distance
        const scalarField yPlus;

        //- Smoothing function
        const scalarField P;

        //- Cell temperature
        const scalarField Tc;

        //- Saturation temperature
        const scalarField Tsat;

        //- Latent heat
        const scalarField L;


    //- Constructor
    boilingLiquidProperties
    (
        const alphatWallBoilingWallFunctionFvPatchScalarField& field,
        const phaseModel& liquid,
        const phaseModel& vapour
    )
    :
        properties(field, liquid, vapour),
        volatileSpecie
        (
            liquid.fluid().lookupOrDefault<word>("volatile", "none")
        ),
        AbyV
        (
            patch().magSf()
           /scalarField
            (
                patch().boundaryMesh().mesh().V(),
                patch().faceCells()
            )
        ),
        trhoLiquidw(liquid.thermo().rho(patchi())),
        rhoLiquidw(trhoLiquidw()),
        trhoVapourw(vapour.thermo().rho(patchi())),
        rhoVapourw(trhoVapourw()),
        Cpw(liquid.thermo().Cp().boundaryField()[patchi()]),
        tnuw(liquid.thermo().nu(patchi())),
        nuw(tnuw()),
        kappaByCp
        (
            liquid.thermo().kappa().boundaryField()[patchi()]
           /liquid.thermo().Cp().boundaryField()[patchi()]
        ),
        nutw
        (
            nutWallFunctionFvPatchScalarField::nutw
            (
                ttm.momentumTransport(),
                patchi()
            )
        ),
        yPlus
        (
            pow025(nutw.Cmu())
           *sqrt(ttm.momentumTransport().k()().boundaryField()[patchi()])
           *ttm.momentumTransport().y()[patchi()]
           /nuw
        ),
        P
        (
            alphatJayatillekeWallFunctionFvPatchScalarField::P
            (
                rhoLiquidw*nuw/kappaByCp/field.Prt_
            )
        ),
        Tc
        (
            liquid.thermo().T().boundaryField()[patchi()].patchInternalField()
        ),
        Tsat
        (
            liquid.fluid().lookupInterfacialModel
            <
                interfaceSaturationTemperatureModel
            >
            (interface)
           .Tsat(liquid.thermo().p())()
           .boundaryField()[patchi()]
        ),
        L
        (
            volatileSpecie != "none"
         ? -refCast<const heatTransferPhaseSystem>(liquid.fluid())
           .Li
            (
                interface,
                volatileSpecie,
                scalarField(patch().size(), +1),
                Tsat,
                patch().faceCells(),
                heatTransferPhaseSystem::latentHeatScheme::upwind
            )
         : -refCast<const heatTransferPhaseSystem>(liquid.fluid())
           .L
            (
                interface,
                scalarField(patch().size(), +1),
                Tsat,
                patch().faceCells(),
                heatTransferPhaseSystem::latentHeatScheme::upwind
            )
        )
    {}
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

tmp<scalarField>
alphatWallBoilingWallFunctionFvPatchScalarField::calcBoiling
(
    const boilingLiquidProperties& props,
    const scalarField& Tw,
    scalarField& dDep,
    scalarField& fDep,
    scalarField& N,
    scalarField& qq,
    scalarField& qe,
    scalarField& dmdtf
) const
{
    scalarField Tl;
    if (!useLiquidTemperatureWallFunction_)
    {
        Tl = props.Tc;
    }
    else
    {
        // Liquid temperature at y+=250 is estimated from the
        // logarithmic thermal wall function of Koncar, Krepper
        // & Egorov (2005)
        const scalarField TyPlus250
        (
            Prt_
           *(
                log(props.nutw.E()*250)
               /props.nutw.kappa()
              + props.P
           )
        );

        const scalarField TyPlus
        (
            Prt_
           *(
                log(props.nutw.E()*max(props.yPlus, scalar(11)))
               /props.nutw.kappa()
              + props.P
            )
        );

        Tl = Tw - (TyPlus250/TyPlus)*(Tw - props.Tc);
    }

    // Bubble departure diameter
    dDep =
        departureDiameterModel_->dDeparture
        (
            props.phase,
            props.otherPhase,
            patch().index(),
            Tl,
            props.Tsat,
            props.L
        );

    // Bubble departure frequency
    fDep =
        departureFrequencyModel_->fDeparture
        (
            props.phase,
            props.otherPhase,
            patch().index(),
            Tl,
            props.Tsat,
            props.L,
            dDep
        );

    // Nucleation site density
    N =
        nucleationSiteModel_->N
        (
            props.phase,
            props.otherPhase,
            patch().index(),
            Tl,
            props.Tsat,
            props.L,
            dDep,
            fDep
        );

    // Del Valle & Kenning (1985)
    const scalarField Ja
    (
        props.rhoLiquidw
       *props.Cpw
       *(props.Tsat - Tl)
       /(props.rhoVapourw*props.L)
    );

    const scalarField Al
    (
        fLiquid_*4.8*exp(min(-Ja/80, log(vGreat)))
    );

    scalarField A2(min(pi*sqr(dDep)*N*Al/4, scalar(1)));
    const scalarField A1(max(1 - A2, scalar(1e-4)));
    scalarField A2E(min(pi*sqr(dDep)*N*Al/4, scalar(5)));

    if (props.volatileSpecie != "none" && !props.phase.pure())
    {
        const scalarField& Yvolatile =
            props.phase
           .Y(props.volatileSpecie)
           .boundaryField()[patch().index()];
        A2E *= Yvolatile;
        A2 *= Yvolatile;
    }

    // Volumetric mass source in the near wall cell due to the
    // wall boiling
    dmdtf = (1.0/6.0)*A2E*dDep*props.rhoVapourw*fDep*props.AbyV;

    // Quenching heat transfer coefficient
    const scalarField hQ
    (
        2*props.kappaByCp*props.Cpw*fDep
       *sqrt((tau_/max(fDep, small))/(pi*props.kappaByCp/props.rhoLiquidw))
    );

    // Quenching heat flux
    qq = A2*hQ*max(Tw - Tl, scalar(0));

    // Evaporation heat flux
    qe = dmdtf*props.L/props.AbyV;

    // Return total sum of convective, quenching and evaporative heat fluxes
    const scalarField gradTw
    (
        patch().deltaCoeffs()*max(Tw - props.Tc, small*props.Tc)
    );
    return A1*props.alphatConv*props.Cpw*gradTw + qq_ + qe_;
}


tmp<scalarField>
alphatWallBoilingWallFunctionFvPatchScalarField::calcBoiling
(
    const boilingLiquidProperties& props,
    const scalarField& Tw
) const
{
    scalarField dDep(dDep_);
    scalarField fDep(fDep_);
    scalarField N(N_);
    scalarField qq(qq_);
    scalarField qe(qe_);
    scalarField dmdtf(dmdtf_);

    return calcBoiling(props, Tw, dDep, fDep, N, qq, qe, dmdtf);
}


tmp<scalarField>
alphatWallBoilingWallFunctionFvPatchScalarField::evaluateBoiling
(
    const boilingLiquidProperties& props,
    const scalarField& Tw
)
{
    return calcBoiling(props, Tw, dDep_, fDep_, N_, qq_, qe_, dmdtf_);
}


const fvPatchScalarField&
alphatWallBoilingWallFunctionFvPatchScalarField::getTemperaturePatchField
(
    const boilingLiquidProperties& props,
    scalarField& isFixed,
    scalarField& h,
    scalarField& hTaPlusQa
) const
{
    isFixed.setSize(patch().size());
    h.setSize(patch().size());
    hTaPlusQa.setSize(patch().size());

    const fvPatchScalarField& Tw =
        props.phase.thermo().T().boundaryField()[patch().index()];

    if (isA<fixedValueFvPatchScalarField>(Tw))
    {
        isFixed = 1;
        h = rootVGreat;
        hTaPlusQa = rootVGreat*Tw;
    }
    else if (isA<zeroGradientFvPatchScalarField>(Tw))
    {
        isFixed = 0;
        h = 0;
        hTaPlusQa = 0;
    }
    else if (isA<fixedGradientFvPatchScalarField>(Tw))
    {
        const fixedGradientFvPatchScalarField& Twm =
            refCast<const fixedGradientFvPatchScalarField>(Tw);

        isFixed = 0;
        h = 0;
        hTaPlusQa = (*this)*props.Cpw*Twm.gradient();
    }
    else if (isA<mixedFvPatchScalarField>(Tw))
    {
        const mixedFvPatchScalarField& Twm =
            refCast<const mixedFvPatchScalarField>(Tw);

        isFixed = pos(Twm.valueFraction() - 1 + rootSmall);
        h =
            Twm.valueFraction()
           /max(1 - Twm.valueFraction(), rootVSmall)
           *(*this)*props.Cpw*patch().deltaCoeffs();
        hTaPlusQa =
            h*Twm.refValue()
          + (*this)*props.Cpw*Twm.refGrad();
    }
    else
    {
        FatalErrorInFunction
            << "Temperature boundary condition type not recognised"
            << exit(FatalError);
    }

    return Tw;
}


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
    tolerance_(rootSmall),

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
    tolerance_(dict.lookupOrDefault<scalar>("tolerance", rootSmall)),

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
        if (dict.found("qEvaporative"))
        {
            qe_ = scalarField("qEvaporative", dict, p.size());
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
    tolerance_(psf.tolerance_),

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
    tolerance_(psf.tolerance_),

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

    // Lookup the fluid model and the phases
    const phaseSystem& fluid =
        db().lookupObject<phaseSystem>(phaseSystem::propertiesName);

    switch (phaseType_)
    {
        case vaporPhase:
        {
            const phaseModel& vapor = fluid.phases()[internalField().group()];
            const phaseModel& liquid = fluid.phases()[otherPhaseName_];

            // Construct boiling properties
            const properties props(*this, vapor, liquid);

            // Partitioning. Note: Assumes that there is only only one liquid
            // phase and all other phases are vapor.
            fLiquid_ = partitioningModel_->fLiquid(props.otherAlphaw);

            operator==
            (
                (1 - fLiquid_)
               /max(1 - props.otherAlphaw, rootSmall)
               *props.alphatConv
            );

            break;
        }

        case liquidPhase:
        {
            const phaseModel& liquid = fluid.phases()[internalField().group()];
            const phaseModel& vapor = fluid.phases()[otherPhaseName_];

            // Boiling is enabled by the presence of saturation temperature
            // modelling. This is consistent with interfacial thermal phase
            // changes.
            if
            (
               !fluid.foundInterfacialModel
                <
                    interfaceSaturationTemperatureModel
                >(phaseInterface(liquid, vapor))
            )
            {
                // Construct non-boiling properties
                const properties props(*this, liquid, vapor);

                Info<< "Saturation model for interface "
                    << props.interface.name()
                    << " not found. Wall boiling disabled." << endl;

                operator==(props.alphatConv);
            }
            else
            {
                // Construct boiling properties
                const boilingLiquidProperties props(*this, liquid, vapor);

                // Partitioning. Note: Assumes that there is only only one
                // liquid phase and all other phases are vapor.
                fLiquid_ = partitioningModel_->fLiquid(props.alphaw);

                // Get the temperature boundary condition and extract its
                // physical parameters
                scalarField TwpfIsFixed, TwpfH, TwpfHTaPlusQa;
                const fvPatchScalarField& Twpf =
                    getTemperaturePatchField
                    (
                        props,
                        TwpfIsFixed,
                        TwpfH,
                        TwpfHTaPlusQa
                    );

                // Define the residual. This should be monotonic in Tw.
                auto R = [&](const scalarField& Tw)
                {
                    return calcBoiling(props, Tw) - TwpfHTaPlusQa + TwpfH*Tw;
                };

                // Solve using interval bisection. Boiling cannot occur below
                // the saturation temperature, so that is taken to be the lower
                // bound. The upper bound is harder to define. For now, take
                // twice the current superheat as the maximum. The solution is
                // likely to be below this value. If it is not, then the
                // iteration will converge to this upper limit, creating a new
                // superheat of twice the current value. Subsequent time-steps
                // will then double this superheat again and again until it
                // does eventually bound the solution.
                const scalarField isBoiling(neg(R(props.Tsat)));
                scalarField Tw0(props.Tsat);
                scalarField Tw1
                (
                    max
                    (
                        Twpf + (Twpf - props.Tsat),
                        props.Tsat*(1 + sqrt(tolerance_))
                    )
                );
                scalar e =
                    gMax((1 - TwpfIsFixed)*isBoiling*(Tw1 - Tw0)/(Tw0 + Tw1));
                for (; e > tolerance_; e /= 2)
                {
                    const scalarField TwM((Tw0 + Tw1)/2);
                    const scalarField rM(R(TwM));
                    Tw0 = pos(rM)*Tw0 + neg0(rM)*TwM;
                    Tw1 = pos(rM)*TwM + neg0(rM)*Tw1;
                }
                const scalarField Tw
                (
                    TwpfIsFixed*Twpf + (1 - TwpfIsFixed)*(Tw0 + Tw1)/2
                );

                // Use solution to re-evaluate the boiling and set the thermal
                // diffusivity to recover the calculated heat flux
                const scalarField gradTw
                (
                    patch().deltaCoeffs()*max(Tw - props.Tc, small*props.Tc)
                );

                const scalarField q(evaluateBoiling(props, Tw));

                operator==
                (
                    isBoiling*q/props.Cpw/gradTw/max(props.alphaw, rootSmall)
                  + (1 - isBoiling)*props.alphatConv
                );
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
    writeEntry(os, "tolerance", tolerance_);

    // Parameters
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
