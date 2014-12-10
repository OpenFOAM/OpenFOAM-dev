/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
#include "demandDrivenData.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "compressible/turbulenceModel/turbulenceModel.H"

using namespace Foam::constant;

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::BrownianMotionForce<CloudType>::erfInv(const scalar y) const
{
    const scalar a = 0.147;
    scalar k = 2.0/(mathematical::pi*a) +  0.5*log(1.0 - y*y);
    scalar h = log(1.0 - y*y)/a;
    scalar x = sqrt(-k + sqrt(k*k - h));

    if (y < 0.0)
    {
        return -x;
    }
    else
    {
        return x;
    }
}


template<class CloudType>
Foam::tmp<Foam::volScalarField>
Foam::BrownianMotionForce<CloudType>::kModel() const
{
    const objectRegistry& obr = this->owner().mesh();
    const word turbName = "turbulenceModel";

    if (obr.foundObject<compressible::turbulenceModel>(turbName))
    {
        const compressible::turbulenceModel& model =
            obr.lookupObject<compressible::turbulenceModel>(turbName);
        return model.k();
    }
    else if (obr.foundObject<incompressible::turbulenceModel>(turbName))
    {
        const incompressible::turbulenceModel& model =
            obr.lookupObject<incompressible::turbulenceModel>(turbName);
        return model.k();
    }
    else
    {
        FatalErrorIn
        (
            "Foam::tmp<Foam::volScalarField>"
            "Foam::BrownianMotionForce<CloudType>::kModel() const"
        )
            << "Turbulence model not found in mesh database" << nl
            << "Database objects include: " << obr.sortedToc()
            << abort(FatalError);

        return tmp<volScalarField>(NULL);
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
    rndGen_(owner.rndGen()),
    lambda_(readScalar(this->coeffs().lookup("lambda"))),
    turbulence_(readBool(this->coeffs().lookup("turbulence"))),
    kPtr_(NULL),
    ownK_(false)
{}


template<class CloudType>
Foam::BrownianMotionForce<CloudType>::BrownianMotionForce
(
    const BrownianMotionForce& bmf
)
:
    ParticleForce<CloudType>(bmf),
    rndGen_(bmf.rndGen_),
    lambda_(bmf.lambda_),
    turbulence_(bmf.turbulence_),
    kPtr_(NULL),
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
                kPtr_ = tk.operator->();
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
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(vector::zero, 0.0);

    const scalar dp = p.d();
    const scalar Tc = p.Tc();

    const scalar eta = rndGen_.sample01<scalar>();
    const scalar alpha = 2.0*lambda_/dp;
    const scalar cc = 1.0 + alpha*(1.257 + 0.4*exp(-1.1/alpha));

    const scalar sigma = physicoChemical::sigma.value();

    scalar f = 0.0;
    if (turbulence_)
    {
        const label cellI = p.cell();
        const volScalarField& k = *kPtr_;
        const scalar kc = k[cellI];
        const scalar Dp = sigma*Tc*cc/(3*mathematical::pi*muc*dp);
        f = eta/mass*sqrt(2.0*sqr(kc)*sqr(Tc)/(Dp*dt));
    }
    else
    {
        const scalar rhoRatio = p.rho()/p.rhoc();
        const scalar s0 =
            216*muc*sigma*Tc/(sqr(mathematical::pi)*pow5(dp)*(rhoRatio)*cc);
        f = eta*sqrt(mathematical::pi*s0/dt);
    }

    const scalar sqrt2 = sqrt(2.0);
    for (label i = 0; i < 3; i++)
    {
        const scalar x = rndGen_.sample01<scalar>();
        const scalar eta = sqrt2*erfInv(2*x - 1.0);
        value.Su()[i] = mass*f*eta;
    }

    return value;
}


// ************************************************************************* //
