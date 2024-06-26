/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2024 OpenFOAM Foundation
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

#include "wallBoilingHeatTransfer.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"
#include "interfaceSaturationTemperatureModel.H"
#include "diameterModel.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(wallBoilingHeatTransfer, 0);
    addToRunTimeSelectionTable
    (
        heatTransferModel,
        wallBoilingHeatTransfer,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::wallBoilingHeatTransfer::wallBoilingHeatTransfer
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    heatTransferModel(dict, interface, registerObject),
    interface_
    (
        interface.modelCast<wallBoilingHeatTransfer, dispersedPhaseInterface>()
    ),
    otherInterface_
    (
        interface
       .modelCast<wallBoilingHeatTransfer, sidedPhaseInterface>()
       .otherInterface()
    ),
    vapourPhaseName_(dict.lookup("vapourPhase")),
    heatTransferModel_
    (
        heatTransferModel::New
        (
            dict.subDict("heatTransferModel"),
            interface,
            false,
            false
        )
    ),
    relax_(dict.lookupOrDefault<scalar>("relax", scalar(1))),
    partitioningModel_(nullptr),
    nucleationSiteModel_(nullptr),
    departureDiamModel_(nullptr),
    departureFreqModel_(nullptr),
    wetFraction_
    (
        IOobject
        (
            IOobject::groupName("fWallBoiling", interface_.name()),
            interface_.mesh().time().name(),
            interface_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        interface_.mesh(),
        dimensionedScalar(dimless,1)
    ),
    dDep_
    (
        IOobject
        (
            IOobject::groupName("departureDiameter", interface_.name()),
            interface_.mesh().time().name(),
            interface_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        interface_.mesh(),
        dimensionedScalar(dimLength, 1e-4)
    ),
    fDep_
    (
        IOobject
        (
            IOobject::groupName("departureFrequency", interface_.name()),
            interface_.mesh().time().name(),
            interface_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        interface_.mesh(),
        dimensionedScalar(inv(dimTime), 0)
    ),
    nucleationSiteDensity_
    (
        IOobject
        (
            IOobject::groupName("nucleationSites", interface_.name()),
            interface_.mesh().time().name(),
            interface_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        interface_.mesh(),
        dimensionedScalar(dimless/dimArea, 0)
    ),
    dmdtf_
    (
        IOobject
        (
            IOobject::groupName(typedName("dmdtf"), interface_.name()),
            interface_.mesh().time().name(),
            interface_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        interface_.mesh(),
        dimensionedScalar(dimDensity/dimTime, 0)
    ),
    qq_
    (
        IOobject
        (
            IOobject::groupName("qq", interface_.name()),
            interface_.mesh().time().name(),
            interface_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        interface_.mesh(),
        dimensionedScalar(dimEnergy/dimTime/dimVolume, 0)
    ),
    Tsurface_
    (
        IOobject
        (
            IOobject::groupName("Tsurface", interface_.name()),
            interface_.mesh().time().name(),
            interface_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        interface_.mesh(),
        dimensionedScalar(dimTemperature, 0)
    ),
    K_
    (
        IOobject
        (
            IOobject::groupName(typedName("K"), interface_.name()),
            interface_.mesh().time().name(),
            interface_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        heatTransferModel_->K()
    )
{
    partitioningModel_ =
        Foam::wallBoilingModels::partitioningModel::New
        (
            dict.subDict("partitioningModel")
        );

    nucleationSiteModel_ =
        Foam::wallBoilingModels::nucleationSiteModel::New
        (
            dict.subDict("nucleationSiteModel")
        );

    departureDiamModel_ =
        Foam::wallBoilingModels::departureDiameterModel::New
        (
            dict.subDict("departureDiameterModel")
        );

    departureFreqModel_ =
        Foam::wallBoilingModels::departureFrequencyModel::New
        (
            dict.subDict("departureFrequencyModel")
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatTransferModels::wallBoilingHeatTransfer::
~wallBoilingHeatTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::heatTransferModels::wallBoilingHeatTransfer::K
(
    const scalar residualAlpha
) const
{
    const phaseSystem& fluid = interface_.fluid();

    const phaseModel& liquid = interface_.continuous();
    const phaseModel& vapour = fluid.phases()[vapourPhaseName_];
    const phaseModel& solid = interface_.dispersed();

    const rhoThermo& lThermo = liquid.thermo();
    const rhoThermo& vThermo = vapour.thermo();
    const rhoThermo& sThermo = solid.thermo();

    // Estimate the surface temperature from the surrounding temperature and
    // heat transfer coefficients. Note that a lagged value of K is used in
    // this calculation.
    {
        static const dimensionedScalar smallK(dimK, rootVSmall);

        const heatTransferModel& otherModel =
            fluid.lookupInterfacialModel<heatTransferModel>(otherInterface_());
        const volScalarField otherK(otherModel.K());

        Tsurface_ =
            (K_*lThermo.T() + otherK*sThermo.T())
           /max(K_ + otherK, smallK);
    }

    const volScalarField Tsat
    (
        fluid.lookupInterfacialModel
        <
            interfaceSaturationTemperatureModel
        >
        (phaseInterface(liquid, vapour))
       .Tsat(liquid.fluidThermo().p())()
    );

    const volScalarField L
    (
        vThermo.ha(liquid.fluidThermo().p(), Tsat) - lThermo.ha()
    );

    // Wetted fraction
    wetFraction_ =
        partitioningModel_->wetFraction(liquid/max(1 - solid, small));

    // Bubble departure diameter
    dDep_ =
        departureDiamModel_->dDeparture
        (
            liquid,
            vapour,
            solid,
            Tsurface_,
            Tsat,
            L
        );

    // Bubble departure frequency
    fDep_ =
        departureFreqModel_->fDeparture
        (
            liquid,
            vapour,
            solid,
            Tsurface_,
            Tsat,
            L,
            dDep_
        );

    // Nucleation site density
    nucleationSiteDensity_ =
        nucleationSiteModel_->nucleationSiteDensity
        (
            liquid,
            vapour,
            solid,
            Tsurface_,
            Tsat,
            L,
            dDep_,
            fDep_
        );

    const tmp<volScalarField> tlrho(lThermo.rho());
    const volScalarField& lrho = tlrho();
    const tmp<volScalarField> tvrho(vThermo.rho());
    const volScalarField& vrho = tvrho();
    const volScalarField& lCp = lThermo.Cp();
    const volScalarField& lkappa = lThermo.kappa();
    //const volScalarField lalpha(lkappa/lCp);

    // Area fractions

    // Del Valle & Kenning (1985)
    const volScalarField Ja
    (
        lrho*lCp*(Tsat - min(Tsurface_, Tsat))/(lrho*L)
    );

    const volScalarField Al
    (
        wetFraction_*4.8*exp(min(-Ja/80, log(vGreat)))
    );

    volScalarField A2
    (
        min(pi*sqr(dDep_)*nucleationSiteDensity_*Al/4, scalar(1))
    );
    const scalarField A1(max(1 - A2, scalar(1e-4)));
    volScalarField A2E
    (
        min(pi*sqr(dDep_)*nucleationSiteDensity_*Al/4, scalar(5))
    );

    const volScalarField Av(solid.diameter().Av());

    // Volumetric mass source in due to the wall boiling and bulk nucleation
    dmdtf_ =
        relax_*(1.0/6.0)*A2E*dDep_*vrho*fDep_*Av
      + (1 - relax_)*dmdtf_;

    const dimensionedScalar zeroT(dimTemperature, small);
    const dimensionedScalar smallT(dimTemperature, small);
    const dimensionedScalar smallF(dimless/dimTime, small);

    // Quenching heat transfer coefficient
    const volScalarField hQ
    (
        2*lkappa*fDep_
       *sqrt((0.8/max(fDep_, smallF))/(pi*(lkappa/lCp)/lrho))
    );

    // Volumetric quenching heat flux
    qq_ =
        (1 - relax_)*qq_
      + relax_*(A2*hQ*max(Tsurface_ - lThermo.T(), zeroT))*Av;

    // Heat transfer coefficient
    K_ =
        heatTransferModel_->K(residualAlpha)
     + (dmdtf_*L + qq_)/max(Tsurface_ - liquid.thermo().T(), smallT);

    return K_;
}


bool Foam::heatTransferModels::wallBoilingHeatTransfer::
activePhaseInterface(const phaseInterfaceKey& phaseInterface) const
{
    const phaseModel& liquid = interface_.continuous();
    const phaseSystem& fluid = liquid.fluid();
    const phaseModel& vapour(fluid.phases()[vapourPhaseName_]);

    if (phaseInterface == phaseInterfaceKey(vapour, liquid))
    {
        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::heatTransferModels::wallBoilingHeatTransfer::
flipSign() const
{
    const phaseModel& liquid = interface_.continuous();
    const phaseSystem& fluid = liquid.fluid();
    const phaseModel& vapour(fluid.phases()[vapourPhaseName_]);

    return vapour.name() != phaseInterfaceKey(vapour, liquid).first();
}


const Foam::volScalarField&
Foam::heatTransferModels::wallBoilingHeatTransfer::dmdtf() const
{
    return dmdtf_;
}


bool Foam::heatTransferModels::wallBoilingHeatTransfer::
writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //
