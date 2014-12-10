/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "LiquidEvaporation.H"
#include "specie.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::tmp<Foam::scalarField> Foam::LiquidEvaporation<CloudType>::calcXc
(
    const label cellI
) const
{
    scalarField Xc(this->owner().thermo().carrier().Y().size());

    forAll(Xc, i)
    {
        Xc[i] =
            this->owner().thermo().carrier().Y()[i][cellI]
           /this->owner().thermo().carrier().W(i);
    }

    return Xc/sum(Xc);
}


template<class CloudType>
Foam::scalar Foam::LiquidEvaporation<CloudType>::Sh
(
    const scalar Re,
    const scalar Sc
) const
{
    return 2.0 + 0.6*Foam::sqrt(Re)*cbrt(Sc);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::LiquidEvaporation<CloudType>::LiquidEvaporation
(
    const dictionary& dict,
    CloudType& owner
)
:
    PhaseChangeModel<CloudType>(dict, owner, typeName),
    liquids_(owner.thermo().liquids()),
    activeLiquids_(this->coeffDict().lookup("activeLiquids")),
    liqToCarrierMap_(activeLiquids_.size(), -1),
    liqToLiqMap_(activeLiquids_.size(), -1)
{
    if (activeLiquids_.size() == 0)
    {
        WarningIn
        (
            "Foam::LiquidEvaporation<CloudType>::LiquidEvaporation"
            "("
                "const dictionary& dict, "
                "CloudType& owner"
            ")"
        )   << "Evaporation model selected, but no active liquids defined"
            << nl << endl;
    }
    else
    {
        Info<< "Participating liquid species:" << endl;

        // Determine mapping between liquid and carrier phase species
        forAll(activeLiquids_, i)
        {
            Info<< "    " << activeLiquids_[i] << endl;
            liqToCarrierMap_[i] =
                owner.composition().globalCarrierId(activeLiquids_[i]);
        }

        // Determine mapping between model active liquids and global liquids
        const label idLiquid = owner.composition().idLiquid();
        forAll(activeLiquids_, i)
        {
            liqToLiqMap_[i] =
                owner.composition().localId(idLiquid, activeLiquids_[i]);
        }
    }
}


template<class CloudType>
Foam::LiquidEvaporation<CloudType>::LiquidEvaporation
(
    const LiquidEvaporation<CloudType>& pcm
)
:
    PhaseChangeModel<CloudType>(pcm),
    liquids_(pcm.owner().thermo().liquids()),
    activeLiquids_(pcm.activeLiquids_),
    liqToCarrierMap_(pcm.liqToCarrierMap_),
    liqToLiqMap_(pcm.liqToLiqMap_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::LiquidEvaporation<CloudType>::~LiquidEvaporation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::LiquidEvaporation<CloudType>::calculate
(
    const scalar dt,
    const label cellI,
    const scalar Re,
    const scalar Pr,
    const scalar d,
    const scalar nu,
    const scalar T,
    const scalar Ts,
    const scalar pc,
    const scalar Tc,
    const scalarField& Yl,
    scalarField& dMassPC
) const
{
    // liquid volume fraction
    const scalarField X(liquids_.X(Yl));

    // immediately evaporate mass that has reached critical condition
    if ((liquids_.Tc(X) - T) < SMALL)
    {
        if (debug)
        {
            WarningIn
            (
                "void Foam::LiquidEvaporation<CloudType>::calculate"
                "("
                    "const scalar, "
                    "const label, "
                    "const scalar, "
                    "const scalar, "
                    "const scalar, "
                    "const scalar, "
                    "const scalar, "
                    "const scalar, "
                    "const scalar, "
                    "const scalar, "
                    "const scalarField&, "
                    "scalarField&"
                ") const"
            )   << "Parcel reached critical conditions: "
                << "evaporating all avaliable mass" << endl;
        }

        forAll(activeLiquids_, i)
        {
            const label lid = liqToLiqMap_[i];
            dMassPC[lid] = GREAT;
        }

        return;
    }

    // construct carrier phase species volume fractions for cell, cellI
    const scalarField Xc(calcXc(cellI));

    // calculate mass transfer of each specie in liquid
    forAll(activeLiquids_, i)
    {
        const label gid = liqToCarrierMap_[i];
        const label lid = liqToLiqMap_[i];

        // vapour diffusivity [m2/s]
        const scalar Dab = liquids_.properties()[lid].D(pc, Ts);

        // saturation pressure for species i [pa]
        // - carrier phase pressure assumed equal to the liquid vapour pressure
        //   close to the surface
        // NOTE: if pSat > pc then particle is superheated
        // calculated evaporation rate will be greater than that of a particle
        // at boiling point, but this is not a boiling model
        const scalar pSat = liquids_.properties()[lid].pv(pc, T);

        // Schmidt number
        const scalar Sc = nu/(Dab + ROOTVSMALL);

        // Sherwood number
        const scalar Sh = this->Sh(Re, Sc);

        // mass transfer coefficient [m/s]
        const scalar kc = Sh*Dab/(d + ROOTVSMALL);

        // vapour concentration at surface [kmol/m3] at film temperature
        const scalar Cs = pSat/(specie::RR*Ts);

        // vapour concentration in bulk gas [kmol/m3] at film temperature
        const scalar Cinf = Xc[gid]*pc/(specie::RR*Ts);

        // molar flux of vapour [kmol/m2/s]
        const scalar Ni = max(kc*(Cs - Cinf), 0.0);

        // mass transfer [kg]
        dMassPC[lid] += Ni*pi*sqr(d)*liquids_.properties()[lid].W()*dt;
    }
}


template<class CloudType>
Foam::scalar Foam::LiquidEvaporation<CloudType>::dh
(
    const label idc,
    const label idl,
    const scalar p,
    const scalar T
) const
{
    scalar dh = 0;

    typedef PhaseChangeModel<CloudType> parent;
    switch (parent::enthalpyTransfer_)
    {
        case (parent::etLatentHeat):
        {
            dh = liquids_.properties()[idl].hl(p, T);
            break;
        }
        case (parent::etEnthalpyDifference):
        {
            scalar hc = this->owner().composition().carrier().Ha(idc, p, T);
            scalar hp = liquids_.properties()[idl].h(p, T);

            dh = hc - hp;
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "Foam::scalar Foam::LiquidEvaporation<CloudType>::dh"
                "("
                    "const label, "
                    "const label, "
                    "const scalar, "
                    "const scalar"
                ") const"
            )   << "Unknown enthalpyTransfer type" << abort(FatalError);
        }
    }

    return dh;
}


template<class CloudType>
Foam::scalar Foam::LiquidEvaporation<CloudType>::Tvap
(
    const scalarField& Y
) const
{
    const scalarField X(liquids_.X(Y));

    return liquids_.Tpt(X);
}


template<class CloudType>
Foam::scalar Foam::LiquidEvaporation<CloudType>::TMax
(
    const scalar p,
    const scalarField& Y
) const
{
    const scalarField X(liquids_.X(Y));

    return liquids_.pvInvert(p, X);
}


// ************************************************************************* //
