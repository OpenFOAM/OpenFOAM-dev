/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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

#include "wallBoiling.H"

#include "multiphaseEuler.H"

#include "fluidThermophysicalTransportModel.H"

#include "saturationTemperatureModel.H"
#include "partitioningModel.H"
#include "nucleationSiteModel.H"
#include "departureDiameterModel.H"
#include "departureFrequencyModel.H"

#include "alphatBoilingWallFunctionFvPatchScalarField.H"
#include "alphatJayatillekeWallFunctionFvPatchScalarField.H"
#include "wallBoilingPhaseChangeRateFvPatchScalarField.H"
#include "zeroFixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

#include "addToRunTimeSelectionTable.H"

/*---------------------------------------------------------------------------*\
                Class wallBoiling::laggedProperties Declaration
\*---------------------------------------------------------------------------*/

struct Foam::fv::wallBoiling::laggedProperties
{
    //- Typedef to shorten the name of the Jayatilleke wall function
    typedef
        compressible::alphatJayatillekeWallFunctionFvPatchScalarField
        alphatJayatillekeWallFunction;

    //- Wall function field
    const wallBoiling& model;

    //- Patch index
    const label patchi;

    //- Liquid volume fraction
    const scalarField& alphaLiquid;

    //- Vapour volume fraction
    const scalarField& alphaVapour;

    //- Liquid thermophysical transport model
    const fluidThermophysicalTransportModel& ttmLiquid;

    //- Vapour thermophysical transport model
    const fluidThermophysicalTransportModel& ttmVapour;

    //- Liquid convective turbulent thermal diffusivity
    const scalarField alphatConvLiquid;

    //- Phase convective turbulent thermal diffusivity
    const scalarField alphatConvVapour;

    //- Patch area by neighbouring cell volume ratio
    const scalarField AbyV;

    //- Liquid density
    private: const tmp<scalarField> trhoLiquid;
    public: const scalarField& rhoLiquid;

    //- Vapour density
    private: const tmp<scalarField> trhoVapour;
    public: const scalarField& rhoVapour;

    //- Liquid heat capacity
    const scalarField& CpLiquid;

    //- Liquid laminar kinematic viscosity
    private: const tmp<scalarField> tnuLiquid;
    public: const scalarField& nuLiquid;

    //- Liquid laminar thermal diffusivity
    const scalarField kappaByCpLiquid;

    //- Liquid viscosity wall function
    const nutWallFunctionFvPatchScalarField& nutLiquid;

    //- Dimensionless wall distance
    const scalarField yPlusLiquid;

    //- Smoothing function
    const scalarField P;

    //- Cell temperature
    const scalarField TcLiquid;

    //- Saturation temperature
    const scalarField Tsat;

    //- Latent heat
    const scalarField L;

    //- Patch
    inline const fvPatch& patch() const
    {
        return model.mesh().boundary()[patchi];
    }

    //- Constructor
    laggedProperties
    (
        const wallBoiling& model,
        const label patchi
    )
    :
        model(model),
        patchi(patchi),
        alphaLiquid(model.liquid_.boundaryField()[patchi]),
        alphaVapour(model.vapour_.boundaryField()[patchi]),
        ttmLiquid
        (
            model.mesh().lookupType<fluidThermophysicalTransportModel>
            (
                model.liquid_.name()
            )
        ),
        ttmVapour
        (
            model.mesh().lookupType<fluidThermophysicalTransportModel>
            (
                model.vapour_.name()
            )
        ),
        alphatConvLiquid
        (
            alphatJayatillekeWallFunction::alphat(ttmLiquid, model.Prt_, patchi)
        ),
        alphatConvVapour
        (
            alphatJayatillekeWallFunction::alphat(ttmVapour, model.Prt_, patchi)
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
        trhoLiquid(model.liquid_.thermo().rho(patchi)),
        rhoLiquid(trhoLiquid()),
        trhoVapour(model.vapour_.thermo().rho(patchi)),
        rhoVapour(trhoVapour()),
        CpLiquid(model.liquid_.thermo().Cp().boundaryField()[patchi]),
        tnuLiquid(model.liquid_.fluidThermo().nu(patchi)),
        nuLiquid(tnuLiquid()),
        kappaByCpLiquid
        (
            model.liquid_.thermo().kappa().boundaryField()[patchi]/CpLiquid
        ),
        nutLiquid
        (
            nutWallFunctionFvPatchScalarField::nutw
            (
                ttmLiquid.momentumTransport(),
                patchi
            )
        ),
        yPlusLiquid
        (
            pow025(nutLiquid.Cmu())
           *sqrt(ttmLiquid.momentumTransport().k()().boundaryField()[patchi])
           *ttmLiquid.momentumTransport().yb()[patchi]
           /nuLiquid
        ),
        P
        (
            alphatJayatillekeWallFunction::P
            (
                rhoLiquid*nuLiquid/kappaByCpLiquid/model.Prt_
            )
        ),
        TcLiquid
        (
            model.liquid_.thermo().T().boundaryField()[patchi]
           .patchInternalField()
        ),
        Tsat
        (
            model.saturationModelPtr_->Tsat
            (
                static_cast<const scalarField&>
                (
                    model.liquid_.fluidThermo().p().boundaryField()[patchi]
                )
            )
        ),
        L
        (
            model.L
            (
                patchi,
                model.liquid_.thermo().T().boundaryField()[patchi]
            )
        )
    {}
};


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(wallBoiling, 0);
    addToRunTimeSelectionTable(fvModel, wallBoiling, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::wallBoiling::readCoeffs(const dictionary& dict)
{
    saturationModelPtr_.reset
    (
        saturationTemperatureModel::New
        (
            "saturationTemperature",
            dict
        ).ptr()
    );

    partitioningModel_.reset
    (
        wallBoilingModels::partitioningModel::New
        (
            dict.subDict("partitioningModel")
        ).ptr()
    );

    nucleationSiteModel_.reset
    (
        wallBoilingModels::nucleationSiteModel::New
        (
            dict.subDict("nucleationSiteModel")
        ).ptr()
    );

    departureDiameterModel_.reset
    (
        wallBoilingModels::departureDiameterModel::New
        (
            dict.subDict("departureDiameterModel")
        ).ptr()
    );

    departureFrequencyModel_.reset
    (
        wallBoilingModels::departureFrequencyModel::New
        (
            dict.subDict("departureFrequencyModel")
        ).ptr()
    );

    tolerance_ = dict.lookupOrDefault<scalar>("tolerance", unitless, rootSmall);

    liquidTemperatureWallFunction_ =
        dict.lookupOrDefault<Switch>("useLiquidTemperatureWallFunction", true);

    Prt_ = dict.lookupOrDefault<scalar>("Prt", dimless, 0.85);

    bubbleWaitingTimeRatio_ =
        dict.lookupOrDefault<scalar>("bubbleWaitingTimeRatio", dimless, 0.8);
}


Foam::wordList Foam::fv::wallBoiling::mDotBoundaryTypes() const
{
    wordList boundaryTypes
    (
        mesh().boundary().size(),
        zeroFixedValueFvPatchScalarField::typeName
    );

    forAll(alphatLiquid_.boundaryField(), patchi)
    {
        const bool liquidIsBoiling =
            isA<alphatBoilingWallFunctionFvPatchScalarField>
            (
                alphatLiquid_.boundaryField()[patchi]
            );
        const bool vapourIsBoiling =
            isA<alphatBoilingWallFunctionFvPatchScalarField>
            (
                alphatVapour_.boundaryField()[patchi]
            );

        if (liquidIsBoiling != vapourIsBoiling)
        {
            FatalErrorInFunction
                << "The field "
                << (liquidIsBoiling ? alphatLiquid_ : alphatVapour_).name()
                << " has a boiling wall function on patch "
                << mesh().boundary()[patchi].name() << " but "
                << (vapourIsBoiling ? alphatLiquid_ : alphatVapour_).name()
                << " does not" << exit(FatalError);
        }

        if (liquidIsBoiling)
        {
            boundaryTypes[patchi] =
                wallBoilingPhaseChangeRateFvPatchScalarField::typeName;
        }
    }

    return boundaryTypes;
}


const Foam::fvPatchScalarField&
Foam::fv::wallBoiling::getLiquidTemperaturePatchField
(
    const laggedProperties& lagProps,
    scalarField& isFixed,
    scalarField& h,
    scalarField& hTaPlusQa
) const
{
    isFixed.setSize(lagProps.patch().size());
    h.setSize(lagProps.patch().size());
    hTaPlusQa.setSize(lagProps.patch().size());

    const fvPatchScalarField& T =
        liquid_.thermo().T().boundaryField()[lagProps.patchi];

    const fvPatchScalarField& alphat =
        alphatLiquid_.boundaryField()[lagProps.patchi];

    if (isA<fixedValueFvPatchScalarField>(T))
    {
        isFixed = 1;
        h = rootVGreat;
        hTaPlusQa = rootVGreat*T;
    }
    else if (isA<zeroGradientFvPatchScalarField>(T))
    {
        isFixed = 0;
        h = 0;
        hTaPlusQa = 0;
    }
    else if (isA<fixedGradientFvPatchScalarField>(T))
    {
        const fixedGradientFvPatchScalarField& Tfg =
            refCast<const fixedGradientFvPatchScalarField>(T);

        isFixed = 0;
        h = 0;
        hTaPlusQa = alphat*lagProps.CpLiquid*Tfg.gradient();
    }
    else if (isA<mixedFvPatchScalarField>(T))
    {
        const mixedFvPatchScalarField& Tm =
            refCast<const mixedFvPatchScalarField>(T);

        isFixed = pos(Tm.valueFraction() - 1 + rootSmall);
        h =
            Tm.valueFraction()
           /max(1 - Tm.valueFraction(), rootVSmall)
           *alphat
           *lagProps.CpLiquid
           *lagProps.patch().deltaCoeffs();
        hTaPlusQa =
            h*Tm.refValue()
          + alphat*lagProps.CpLiquid*Tm.refGrad();
    }
    else
    {
        FatalErrorInFunction
            << "Temperature boundary condition type not recognised"
            << exit(FatalError);
    }

    return T;
}


Foam::tmp<Foam::scalarField> Foam::fv::wallBoiling::calcBoiling
(
    const laggedProperties& lagProps,
    const scalarField& TLiquid,
    const scalarField& wetFraction,
    scalarField& dDeparture,
    scalarField& fDeparture,
    scalarField& nucleationSiteDensity,
    scalarField& qQuenching,
    scalarField& qEvaporative,
    scalarField& mDot
) const
{
    using constant::mathematical::pi;

    scalarField Tl;
    if (!liquidTemperatureWallFunction_)
    {
        Tl = lagProps.TcLiquid;
    }
    else
    {
        // Liquid temperature at y+=250 is estimated from the logarithmic
        // thermal wall function of Koncar, Krepper & Egorov (2005)
        const scalarField TyPlus250
        (
            Prt_
           *(
                log(lagProps.nutLiquid.E()*250)
               /lagProps.nutLiquid.kappa()
              + lagProps.P
           )
        );

        const scalarField TyPlus
        (
            Prt_
           *(
                log
                (
                    lagProps.nutLiquid.E()
                   *max(lagProps.yPlusLiquid, scalar(11))
                )
               /lagProps.nutLiquid.kappa()
              + lagProps.P
            )
        );

        Tl = TLiquid - (TyPlus250/TyPlus)*(TLiquid - lagProps.TcLiquid);
    }

    // Bubble departure diameter
    dDeparture =
        departureDiameterModel_->dDeparture
        (
            liquid_,
            vapour_,
            lagProps.patchi,
            Tl,
            lagProps.Tsat,
            lagProps.L
        );

    // Bubble departure frequency
    fDeparture =
        departureFrequencyModel_->fDeparture
        (
            liquid_,
            vapour_,
            lagProps.patchi,
            Tl,
            lagProps.Tsat,
            lagProps.L,
            dDeparture
        );

    // Nucleation site density
    nucleationSiteDensity =
        nucleationSiteModel_->nucleationSiteDensity
        (
            liquid_,
            vapour_,
            lagProps.patchi,
            Tl,
            lagProps.Tsat,
            lagProps.L,
            dDeparture,
            fDeparture
        );

    // Del Valle & Kenning (1985)
    const scalarField Ja
    (
        lagProps.rhoLiquid
       *lagProps.CpLiquid
       *max(lagProps.Tsat - Tl, scalar(0))
       /(lagProps.rhoVapour*lagProps.L)
    );

    const scalarField Al
    (
        wetFraction*4.8*exp(min(-Ja/80, log(vGreat)))
    );

    scalarField A2
    (
        min(pi*sqr(dDeparture)*nucleationSiteDensity*Al/4, scalar(1))
    );
    const scalarField A1(max(1 - A2, scalar(1e-4)));
    scalarField A2E
    (
        min(pi*sqr(dDeparture)*nucleationSiteDensity*Al/4, scalar(5))
    );

    // Volumetric mass source in the near wall cell due to the
    // wall boiling
    mDot = (1.0/6.0)*A2E*dDeparture*lagProps.rhoVapour*fDeparture*lagProps.AbyV;

    // Quenching heat transfer coefficient
    const scalarField hQ
    (
        2*lagProps.kappaByCpLiquid*lagProps.CpLiquid*fDeparture
       *sqrt
        (
            bubbleWaitingTimeRatio_
           /max(fDeparture, small)
           /(pi*lagProps.kappaByCpLiquid/lagProps.rhoLiquid)
        )
    );

    // Quenching heat flux
    qQuenching = A2*hQ*max(TLiquid - Tl, scalar(0));

    // Evaporation heat flux
    qEvaporative = mDot*lagProps.L/lagProps.AbyV;

    // Return total sum of convective, quenching and evaporative heat fluxes
    const scalarField gradT
    (
        lagProps.patch().deltaCoeffs()
       *max(TLiquid - lagProps.TcLiquid, small*lagProps.TcLiquid)
    );
    return
        A1*lagProps.alphatConvLiquid*lagProps.CpLiquid*gradT
      + qQuenching
      + qEvaporative;
}


Foam::tmp<Foam::scalarField> Foam::fv::wallBoiling::calcBoiling
(
    const wallBoilingPhaseChangeRateFvPatchScalarField& mDot,
    const laggedProperties& lagProps,
    const scalarField& TLiquid
) const
{
    scalarField dDeparture(mDot.dDeparture_);
    scalarField fDeparture(mDot.fDeparture_);
    scalarField nucleationSiteDensity(mDot.nucleationSiteDensity_);
    scalarField qQuenching(mDot.qQuenching_);
    scalarField qEvaporative(mDot.qEvaporative_);
    scalarField mDotCopy(mDot);

    return
        calcBoiling
        (
            lagProps,
            TLiquid,
            mDot.wetFraction_,
            dDeparture,
            fDeparture,
            nucleationSiteDensity,
            qQuenching,
            qEvaporative,
            mDotCopy
        );
}


Foam::tmp<Foam::scalarField> Foam::fv::wallBoiling::evaluateBoiling
(
    wallBoilingPhaseChangeRateFvPatchScalarField& mDot,
    const laggedProperties& lagProps,
    const scalarField& TLiquid
) const
{
    return
        calcBoiling
        (
            lagProps,
            TLiquid,
            mDot.wetFraction_,
            mDot.dDeparture_,
            mDot.fDeparture_,
            mDot.nucleationSiteDensity_,
            mDot.qQuenching_,
            mDot.qEvaporative_,
            mDot
        );
}


void Foam::fv::wallBoiling::correctMDot() const
{
    Info<< type() << ": " << name() << endl << incrIndent;

    //- Reset the phase-change rates in all the near-wall cells
    forAll(mDot_.boundaryField(), patchi)
    {
        if (!isBoiling(patchi)) continue;

        const labelUList& faceCells = mesh().boundary()[patchi].faceCells();
        forAll(faceCells, i)
        {
            mDot_[faceCells[i]] = scalar(0);
        }
    }

    // Loop the boiling patches, evaluate the model, and sum the phase change
    // rates into the adjacent cells
    forAll(mDot_.boundaryField(), patchi)
    {
        if (!isBoiling(patchi)) continue;

        // Access the wall-boiling phase-change patch field for this patch
        wallBoilingPhaseChangeRateFvPatchScalarField& mDot = mDotPfRef(patchi);

        // Construct lagged properties
        const laggedProperties lagProps(*this, patchi);

        // Evaluate the wetted fraction
        mDot.wetFraction_ =
            partitioningModel_->wetFraction(lagProps.alphaLiquid);

        // Set the vapour turbulent thermal diffusivity
        mDot.alphatVapour_ =
            (1 - mDot.wetFraction_)
           /max(1 - lagProps.alphaLiquid, rootSmall)
           *lagProps.alphatConvVapour;

        // Get the temperature boundary condition and extract its
        // physical parameters
        scalarField TLiquidIsFixed, TLiquidH, TLiquidHTaPlusQa;
        const fvPatchScalarField& TLiquid =
            getLiquidTemperaturePatchField
            (
                lagProps,
                TLiquidIsFixed,
                TLiquidH,
                TLiquidHTaPlusQa
            );

        // Define the residual. This should be monotonic in T.
        auto R = [&](const scalarField& T)
        {
            return
                calcBoiling(mDot, lagProps, T)
              - TLiquidHTaPlusQa
              + TLiquidH*T;
        };

        // Solve using interval bisection. Boiling cannot occur below the
        // saturation temperature, so that is taken to be the lower bound. The
        // upper bound is harder to define. For now, take twice the current
        // superheat as the maximum. The solution is likely to be below this
        // value. If it is not, then the iteration will converge to this upper
        // limit, creating a new superheat of twice the current value.
        // Subsequent time-steps will then double this superheat again and
        // again until it does eventually bound the solution.
        const scalarField isBoiling(neg(R(lagProps.Tsat)));
        scalarField T0(lagProps.Tsat);
        scalarField T1
        (
            max
            (
                TLiquid + (TLiquid - lagProps.Tsat),
                lagProps.Tsat*(1 + sqrt(tolerance_))
            )
        );
        scalar e =
            gMax((1 - TLiquidIsFixed)*isBoiling*(T1 - T0)/(T0 + T1));
        for (; e > tolerance_; e /= 2)
        {
            const scalarField Tm((T0 + T1)/2);
            const scalarField rM(R(Tm));
            T0 = pos(rM)*T0 + neg0(rM)*Tm;
            T1 = pos(rM)*Tm + neg0(rM)*T1;
        }
        const scalarField T
        (
            TLiquidIsFixed*TLiquid + (1 - TLiquidIsFixed)*(T0 + T1)/2
        );

        // Use solution to re-evaluate the boiling and set the thermal
        // diffusivity to recover the calculated heat flux
        const scalarField gradT
        (
            lagProps.patch().deltaCoeffs()
           *max(T - lagProps.TcLiquid, small*lagProps.TcLiquid)
        );

        const scalarField q(evaluateBoiling(mDot, lagProps, T));

        const scalarField alphatBoilingLiquid
        (
            q/lagProps.CpLiquid/gradT/max(lagProps.alphaLiquid, rootSmall)
        );

        mDot.alphatLiquid_ =
            isBoiling*alphatBoilingLiquid
          + (1 - isBoiling)*lagProps.alphatConvLiquid;

        Info<< indent << mesh().boundary()[patchi].name()
            << ": min/mean/max mDot = " << gMin(mDot) << '/'
            << gAverage(mDot) << '/' << gMax(mDot) << endl;

        // Sum the phase change rate into the internal field
        const labelUList& faceCells = mesh().boundary()[patchi].faceCells();
        forAll(faceCells, i)
        {
            mDot_[faceCells[i]] += mDot[i];
        }
    }

    Info<< decrIndent;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::wallBoiling::wallBoiling
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    phaseChange(name, modelType, mesh, dict, wordList()),
    nucleation(),
    fluid_(mesh().lookupObject<phaseSystem>(phaseSystem::propertiesName)),
    liquid_(fluid_.phases()[phaseNames().first()]),
    vapour_(fluid_.phases()[phaseNames().second()]),
    alphatLiquid_
    (
        mesh().lookupObject<volScalarField>
        (
            IOobject::groupName("alphat", liquid_.name())
        )
    ),
    alphatVapour_
    (
        mesh().lookupObject<volScalarField>
        (
            IOobject::groupName("alphat", vapour_.name())
        )
    ),
    p_rgh_
    (
        mesh().lookupObject<solvers::multiphaseEuler>(solver::typeName).p_rgh
    ),
    tolerance_(NaN),
    liquidTemperatureWallFunction_(true),
    Prt_(NaN),
    bubbleWaitingTimeRatio_(NaN),
    saturationModelPtr_(nullptr),
    partitioningModel_(nullptr),
    nucleationSiteModel_(nullptr),
    departureDiameterModel_(nullptr),
    departureFrequencyModel_(nullptr),
    pressureEquationIndex_(-1),
    mDot_
    (
        IOobject
        (
            name + ":mDot",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimDensity/dimTime, scalar(0)),
        mDotBoundaryTypes()
    )
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::wallBoiling::isBoiling(const label patchi) const
{
    return
        isA<wallBoilingPhaseChangeRateFvPatchScalarField>
        (
            mDot_.boundaryFieldRef()[patchi]
        );
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::wallBoiling::Lfraction() const
{
    // Put all the latent heat into the liquid
    return
        volScalarField::Internal::New
        (
            name() + ":Lfraction",
            mesh(),
            dimensionedScalar(dimless, scalar(0))
        );
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::wallBoiling::d() const
{
    tmp<volScalarField> tdVapour = vapour_.d();
    const volScalarField::Internal& dVapour = tdVapour();

    const volScalarField::Internal mask
    (
        neg(mag(mDot_) - dimensionedScalar(dimDensity/dimTime, rootVSmall))
    );

    tmp<volScalarField::Internal> td =
        volScalarField::Internal::New(name() + ":d", mask*dVapour);
    volScalarField::Internal& d = td.ref();

    forAll(mDot_.boundaryField(), patchi)
    {
        if (!isBoiling(patchi)) continue;

        // Access the wall-boiling phase-change patch field for this patch
        const wallBoilingPhaseChangeRateFvPatchScalarField& mDot =
            mDotPf(patchi);

        // Average the diameter into the internal field
        const labelUList& faceCells = mesh().boundary()[patchi].faceCells();
        forAll(faceCells, i)
        {
            if (mask[faceCells[i]]) continue;

            d[faceCells[i]] += mDot.dDeparture_[i]*mDot[i]/mDot_[faceCells[i]];
        }
    }

    return td;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::wallBoiling::nDot() const
{
    tmp<volScalarField::Internal> tnDot =
        volScalarField::Internal::New
        (
            name() + ":nDot",
            mesh(),
            dimensionedScalar(inv(dimTime), scalar(0))
        );
    volScalarField::Internal& nDot = tnDot.ref();

    forAll(mDot_.boundaryField(), patchi)
    {
        if (!isBoiling(patchi)) continue;

        // Access the wall-boiling phase-change patch field for this patch
        const wallBoilingPhaseChangeRateFvPatchScalarField& mDot =
            mDotPf(patchi);

        // Sum the frequency into the internal field
        const labelUList& faceCells = mesh().boundary()[patchi].faceCells();
        forAll(faceCells, i)
        {
            nDot[faceCells[i]] += mDot.fDeparture_[i];
        }
    }

    return tnDot;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::wallBoiling::tau() const
{
    NotImplemented;
    return tmp<volScalarField::Internal>(nullptr);
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::wallBoiling::mDot() const
{
    return mDot_.internalField();
}


const Foam::wallBoilingPhaseChangeRateFvPatchScalarField&
Foam::fv::wallBoiling::mDotPf(const label patchi) const
{
    if (!isBoiling(patchi))
    {
        FatalErrorInFunction
            << "Patch " << mesh().boundary()[patchi].name()
            << " is not boiling" << exit(FatalError);
    }

    return
        refCast<const wallBoilingPhaseChangeRateFvPatchScalarField>
        (
            mDot_.boundaryField()[patchi]
        );
}


Foam::wallBoilingPhaseChangeRateFvPatchScalarField&
Foam::fv::wallBoiling::mDotPfRef(const label patchi) const
{
    if (!isBoiling(patchi))
    {
        FatalErrorInFunction
            << "Patch " << mesh().boundary()[patchi].name()
            << " is not boiling" << exit(FatalError);
    }

    return
        refCast<wallBoilingPhaseChangeRateFvPatchScalarField>
        (
            mDot_.boundaryFieldRef()[patchi]
        );
}


void Foam::fv::wallBoiling::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn
) const
{
    // Pressure equation (i.e., continuity, linearised in the pressure)
    if
    (
        (&alpha == &liquid_ || &alpha == &vapour_)
     && (&rho == &liquid_.rho() || &rho == &vapour_.rho())
     && &eqn.psi() == &p_rgh_
    )
    {
        // Ensure that the source is up-to date if this is the first call in
        // the current phase loop
        if (pressureEquationIndex_ % 2 == 0) correctMDot();
        pressureEquationIndex_ ++;
    }

    // Let the base class add the actual source
    massTransfer::addSup(alpha, rho, eqn);
}


void Foam::fv::wallBoiling::correct()
{
    // Reset the p_rgh equation solution counter
    pressureEquationIndex_ = 0;

    // Correct the total phase change rate
    correctMDot();
}


bool Foam::fv::wallBoiling::read(const dictionary& dict)
{
    if (phaseChange::read(dict))
    {
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
