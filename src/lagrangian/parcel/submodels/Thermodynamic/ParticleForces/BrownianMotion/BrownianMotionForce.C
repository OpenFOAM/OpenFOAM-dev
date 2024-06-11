/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "BrownianMotionForce.H"
#include "mathematicalConstants.H"
#include "fundamentalConstants.H"
#include "demandDrivenData.H"
#include "momentumTransportModel.H"
#include "standardNormal.H"

using namespace Foam::constant;

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class CloudType>
Foam::tmp<Foam::volScalarField>
Foam::BrownianMotionForce<CloudType>::kModel() const
{
    const objectRegistry& obr = this->owner().mesh();

    if (obr.foundType<momentumTransportModel>(this->owner().U().group()))
    {
        const momentumTransportModel& model =
            obr.lookupType<momentumTransportModel>(this->owner().U().group());

        return model.k();
    }
    else
    {
        FatalErrorInFunction
            << "Turbulence model not found in mesh database" << nl
            << "Database objects include: " << obr.sortedToc()
            << abort(FatalError);

        return tmp<volScalarField>(nullptr);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BrownianMotionForce<CloudType>::BrownianMotionForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    lambda_(this->coeffs().template lookup<scalar>("lambda")),
    turbulence_(readBool(this->coeffs().lookup("turbulence"))),
    kPtr_(nullptr),
    ownK_(false)
{}


template<class CloudType>
Foam::BrownianMotionForce<CloudType>::BrownianMotionForce
(
    const BrownianMotionForce& bmf
)
:
    ParticleForce<CloudType>(bmf),
    lambda_(bmf.lambda_),
    turbulence_(bmf.turbulence_),
    kPtr_(nullptr),
    ownK_(false)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BrownianMotionForce<CloudType>::~BrownianMotionForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::BrownianMotionForce<CloudType>::cacheFields(const bool store)
{
    if (turbulence_)
    {
        if (store)
        {
            tmp<volScalarField> tk = kModel();
            if (tk.isTmp())
            {
                kPtr_ = tk.ptr();
                ownK_ = true;
            }
            else
            {
                kPtr_ = &tk();
                ownK_ = false;
            }
        }
        else
        {
            if (ownK_ && kPtr_)
            {
                deleteDemandDrivenData(kPtr_);
                ownK_ = false;
            }
        }
    }
}


template<class CloudType>
Foam::forceSuSp Foam::BrownianMotionForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(Zero, 0.0);

    const scalar dp = p.d();
    const scalar Tc = td.Tc();

    const scalar alpha = 2.0*lambda_/dp;
    const scalar cc = 1.0 + alpha*(1.257 + 0.4*exp(-1.1/alpha));

    // Boltzmann constant
    const scalar kb = physicoChemical::k.value();

    scalar f = 0;
    if (turbulence_)
    {
        const label celli = p.cell();
        const volScalarField& k = *kPtr_;
        const scalar kc = k[celli];
        const scalar Dp = kb*Tc*cc/(3*mathematical::pi*muc*dp);
        f = sqrt(2.0*sqr(kc)*sqr(Tc)/(Dp*dt));
    }
    else
    {
        const scalar s0 =
            216*muc*kb*Tc/(sqr(mathematical::pi)*pow5(dp)*sqr(p.rho())*cc);
        f = mass*sqrt(mathematical::pi*s0/dt);
    }

    randomGenerator& rndGen = this->owner().rndGen();
    distributions::standardNormal& stdNormal = this->owner().stdNormal();

    // To generate a cubic distribution (i.e., 3 independent directions):
    // value.Su() = f*stdNormal.sample<vector>();

    // To generate a spherical distribution:
    const scalar theta = rndGen.scalar01()*twoPi;
    const scalar u = 2*rndGen.scalar01() - 1;
    const scalar a = sqrt(1 - sqr(u));
    const vector dir(a*cos(theta), a*sin(theta), u);
    value.Su() = f*mag(stdNormal.sample())*dir;

    return value;
}


// ************************************************************************* //
