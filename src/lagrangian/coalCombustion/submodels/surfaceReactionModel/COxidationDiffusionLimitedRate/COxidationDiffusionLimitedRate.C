/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "COxidationDiffusionLimitedRate.H"
#include "mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::COxidationDiffusionLimitedRate<CloudType>::COxidationDiffusionLimitedRate
(
    const dictionary& dict,
    CloudType& owner
)
:
    SurfaceReactionModel<CloudType>(dict, owner, typeName),
    Sb_(readScalar(this->coeffDict().lookup("Sb"))),
    D_(readScalar(this->coeffDict().lookup("D"))),
    CsLocalId_(-1),
    O2GlobalId_(owner.composition().globalCarrierId("O2")),
    CO2GlobalId_(owner.composition().globalCarrierId("CO2")),
    WC_(0.0),
    WO2_(0.0),
    HcCO2_(0.0)
{
    // Determine Cs ids
    label idSolid = owner.composition().idSolid();
    CsLocalId_ = owner.composition().localId(idSolid, "C");

    // Set local copies of thermo properties
    WO2_ = owner.thermo().carrier().W(O2GlobalId_);
    const scalar WCO2 = owner.thermo().carrier().W(CO2GlobalId_);
    WC_ = WCO2 - WO2_;

    HcCO2_ = owner.thermo().carrier().Hc(CO2GlobalId_);

    if (Sb_ < 0)
    {
        FatalErrorIn
        (
            "COxidationDiffusionLimitedRate<CloudType>"
            "("
                "const dictionary&, "
                "CloudType&"
            ")"
        )   << "Stoichiometry of reaction, Sb, must be greater than zero" << nl
            << exit(FatalError);
    }

    const scalar YCloc = owner.composition().Y0(idSolid)[CsLocalId_];
    const scalar YSolidTot = owner.composition().YMixture0()[idSolid];
    Info<< "    C(s): particle mass fraction = " << YCloc*YSolidTot << endl;
}


template<class CloudType>
Foam::COxidationDiffusionLimitedRate<CloudType>::COxidationDiffusionLimitedRate
(
    const COxidationDiffusionLimitedRate<CloudType>& srm
)
:
    SurfaceReactionModel<CloudType>(srm),
    Sb_(srm.Sb_),
    D_(srm.D_),
    CsLocalId_(srm.CsLocalId_),
    O2GlobalId_(srm.O2GlobalId_),
    CO2GlobalId_(srm.CO2GlobalId_),
    WC_(srm.WC_),
    WO2_(srm.WO2_),
    HcCO2_(srm.HcCO2_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::COxidationDiffusionLimitedRate<CloudType>::
~COxidationDiffusionLimitedRate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::COxidationDiffusionLimitedRate<CloudType>::calculate
(
    const scalar dt,
    const label cellI,
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
    // Fraction of remaining combustible material
    const label idSolid = CloudType::parcelType::SLD;
    const scalar fComb = YMixture[idSolid]*YSolid[CsLocalId_];

    // Surface combustion active combustible fraction is consumed
    if (fComb < SMALL)
    {
        return 0.0;
    }

    const SLGThermo& thermo = this->owner().thermo();

    // Local mass fraction of O2 in the carrier phase
    const scalar YO2 = thermo.carrier().Y(O2GlobalId_)[cellI];

    // Change in C mass [kg]
    scalar dmC = 4.0*mathematical::pi*d*D_*YO2*Tc*rhoc/(Sb_*(T + Tc))*dt;

    // Limit mass transfer by availability of C
    dmC = min(mass*fComb, dmC);

    // Change in O2 mass [kg]
    const scalar dmO2 = dmC/WC_*Sb_*WO2_;

    // Mass of newly created CO2 [kg]
    const scalar dmCO2 = dmC + dmO2;

    // Update local particle C mass
    dMassSolid[CsLocalId_] += dmC;

    // Update carrier O2 and CO2 mass
    dMassSRCarrier[O2GlobalId_] -= dmO2;
    dMassSRCarrier[CO2GlobalId_] += dmCO2;

    const scalar HsC = thermo.solids().properties()[CsLocalId_].Hs(T);

    // carrier sensible enthalpy exchange handled via change in mass

    // Heat of reaction [J]
    return dmC*HsC - dmCO2*HcCO2_;
}


// ************************************************************************* //
