/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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

#include "COxidationHurtMitchell.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::COxidationHurtMitchell<CloudType>::COxidationHurtMitchell
(
    const dictionary& dict,
    CloudType& owner
)
:
    SurfaceReactionModel<CloudType>(dict, owner, typeName),
    Sb_(this->coeffDict().template lookup<scalar>("Sb")),
    CsLocalId_(-1),
    ashLocalId_(-1),
    O2GlobalId_(owner.composition().carrierId("O2")),
    CO2GlobalId_(owner.composition().carrierId("CO2")),
    WC_(0.0),
    WO2_(0.0),
    HcCO2_(0.0),
    heatOfReaction_(-1.0)
{
    // Determine Cs and ash ids
    label idSolid = owner.composition().idSolid();
    CsLocalId_ = owner.composition().localId(idSolid, "C");
    ashLocalId_ = owner.composition().localId(idSolid, "ash", true);

    // Set local copies of thermo properties
    WO2_ = owner.composition().carrier().Wi(O2GlobalId_);
    const scalar WCO2 = owner.composition().carrier().Wi(CO2GlobalId_);
    WC_ = WCO2 - WO2_;

    HcCO2_ = owner.composition().carrier().Hf(CO2GlobalId_);

    const scalar YCloc = owner.composition().Y0(idSolid)[CsLocalId_];
    const scalar YSolidTot = owner.composition().YMixture0()[idSolid];
    Info<< "    C(s): particle mass fraction = " << YCloc*YSolidTot << endl;

    if (this->coeffDict().readIfPresent("heatOfReaction", heatOfReaction_))
    {
        Info<< "    Using user specified heat of reaction: "
            << heatOfReaction_ << " [J/kg]" << endl;
    }
}


template<class CloudType>
Foam::COxidationHurtMitchell<CloudType>::COxidationHurtMitchell
(
    const COxidationHurtMitchell<CloudType>& srm
)
:
    SurfaceReactionModel<CloudType>(srm),
    Sb_(srm.Sb_),
    CsLocalId_(srm.CsLocalId_),
    ashLocalId_(srm.ashLocalId_),
    O2GlobalId_(srm.O2GlobalId_),
    CO2GlobalId_(srm.CO2GlobalId_),
    WC_(srm.WC_),
    WO2_(srm.WO2_),
    HcCO2_(srm.HcCO2_),
    heatOfReaction_(srm.heatOfReaction_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::COxidationHurtMitchell<CloudType>::~COxidationHurtMitchell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::COxidationHurtMitchell<CloudType>::calculate
(
    const scalar dt,
    const label celli,
    const scalar d,
    const scalar T,
    const scalar Tc,
    const scalar pc,
    const scalar rhoc,
    const scalar mass,
    const scalarField& YGas,
    const scalarField& YLiquid,
    const scalarField& YSolid,
    const scalarField& YMixture,
    const scalar N,
    scalarField& dMassGas,
    scalarField& dMassLiquid,
    scalarField& dMassSolid,
    scalarField& dMassSRCarrier
) const
{
    const label idGas = CloudType::parcelType::GAS;
    const label idSolid = CloudType::parcelType::SLD;
    const scalar Ychar = YMixture[idSolid]*YSolid[CsLocalId_];

    // Surface combustion until combustible fraction is consumed
    if (Ychar < small)
    {
        return 0.0;
    }

    const parcelThermo& thermo = this->owner().thermo();
    const basicSpecieMixture& carrier = this->owner().composition().carrier();

    // Local mass fraction of O2 in the carrier phase
    const scalar YO2 = carrier.Y(O2GlobalId_)[celli];

    // No combustion if no oxygen present
    if (YO2 < small)
    {
        return 0.0;
    }

    // Conversion from [g/cm^2) to [kg/m^2]
    const scalar convSI = 1000.0/10000.0;

    // Universal gas constant in [kcal/mol/K]
    const scalar RRcal = 1985.877534;

    // Dry mass fraction
    scalar Ydaf = YMixture[idGas] + YMixture[idSolid];
    if (ashLocalId_ != -1)
    {
        Ydaf -= YMixture[idSolid]*YSolid[ashLocalId_];
    }

    // Char percentage
    const scalar charPrc = max(0, min(Ychar/(Ydaf + rootVSmall)*100.0, 100));

    // Particle surface area
    const scalar Ap = constant::mathematical::pi*sqr(d);

    // Far field partial pressure O2 [Pa]
    // Note: Should really use the surface partial pressure
    const scalar ppO2 = max(0.0, rhoc*YO2/WO2_*RR*Tc);

    // Activation energy [kcal/mol]
    const scalar E = -5.94 + 0.355*charPrc;

    // Pre-exponential factor [g/(cm^2.s.atm^0.5)]
    const scalar lnK1750 = 2.8 - 0.0758*charPrc;
    const scalar A = exp(lnK1750 + E/RRcal/1750.0);

    // Kinetic rate of char oxidation [g/(cm^2.s.atm^0.5)]
    const scalar Rk = A*exp(-E/(RRcal*T));

    // Molar reaction rate per unit surface area [kmol/m^2/s]
    const scalar qCsLim = mass*Ychar/(WC_*Ap*dt);
    const scalar qCs = min(convSI*Rk*Foam::sqrt(ppO2/101325.0), qCsLim);

    // Calculate the number of molar units reacted [kmol]
    const scalar dOmega = qCs*Ap*dt;

    // Add to carrier phase mass transfer
    dMassSRCarrier[O2GlobalId_] += -dOmega*Sb_*WO2_;
    dMassSRCarrier[CO2GlobalId_] += dOmega*(WC_ + Sb_*WO2_);

    // Add to particle mass transfer
    dMassSolid[CsLocalId_] += dOmega*WC_;


    // Return the heat of reaction [J]
    // note: carrier sensible enthalpy exchange handled via change in mass
    if (heatOfReaction_ < 0)
    {
        const scalar HsC = thermo.solids().properties()[CsLocalId_].Hs(T);
        return dOmega*(WC_*HsC - (WC_ + Sb_*WO2_)*HcCO2_);
    }
    else
    {
        return dOmega*WC_*heatOfReaction_;
    }
}


// ************************************************************************* //
