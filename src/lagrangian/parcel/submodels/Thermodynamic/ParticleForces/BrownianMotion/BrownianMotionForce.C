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
    lambda_(this->coeffs().template lookup<scalar>("lambda"))
{}


template<class CloudType>
Foam::BrownianMotionForce<CloudType>::BrownianMotionForce
(
    const BrownianMotionForce& bmf
)
:
    ParticleForce<CloudType>(bmf),
    lambda_(bmf.lambda_)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BrownianMotionForce<CloudType>::~BrownianMotionForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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

    // Equation 25
    const scalar s0 =
        216*muc*kb*Tc/(sqr(mathematical::pi)*pow5(dp)*sqr(p.rho())*cc);

    // Equation 26
    const scalar f = mass*sqrt(mathematical::pi*s0/dt);

    randomGenerator& rndGen = this->owner().rndGen();
    distributions::standardNormal& stdNormal = this->owner().stdNormal();

    // To generate a cubic distribution (i.e., 3 independent directions):
    // value.Su() = f*stdNormal.sample<vector>();

    // To generate a spherical distribution:
    const scalar theta = rndGen.scalar01()*mathematical::twoPi;
    const scalar u = 2*rndGen.scalar01() - 1;
    const scalar a = sqrt(1 - sqr(u));
    const vector dir(a*cos(theta), a*sin(theta), u);
    value.Su() = f*mag(stdNormal.sample())*dir;

    return value;
}


// ************************************************************************* //
