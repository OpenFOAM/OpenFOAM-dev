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

#include "ORourkeCollision.H"
#include "mathematicalConstants.H"
#include "SLGThermo.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::ORourkeCollision<CloudType>::collide(const scalar dt)
{
    label i = 0;
    forAllIter(typename CloudType, this->owner(), iter1)
    {
        label j = 0;
        forAllIter(typename CloudType, this->owner(), iter2)
        {
            if (j > i)
            {
                parcelType& p1 = iter1();
                parcelType& p2 = iter2();

                scalar m1 = p1.nParticle()*p1.mass();
                scalar m2 = p2.nParticle()*p2.mass();

                bool massChanged = collideParcels(dt, p1, p2, m1, m2);

                if (massChanged)
                {
                    if (m1 > ROOTVSMALL)
                    {
                        const scalarField X(liquids_.X(p1.Y()));
                        p1.rho() = liquids_.rho(p1.pc(), p1.T(), X);
                        p1.Cp() = liquids_.Cp(p1.pc(), p1.T(), X);
                        p1.sigma() = liquids_.sigma(p1.pc(), p1.T(), X);
                        p1.mu() = liquids_.mu(p1.pc(), p1.T(), X);
                        p1.d() = cbrt(6.0*m1/(p1.nParticle()*p1.rho()*pi));
                    }

                    if (m2 > ROOTVSMALL)
                    {
                        const scalarField X(liquids_.X(p2.Y()));
                        p2.rho() = liquids_.rho(p2.pc(), p2.T(), X);
                        p2.Cp() = liquids_.Cp(p2.pc(), p2.T(), X);
                        p2.sigma() = liquids_.sigma(p2.pc(), p2.T(), X);
                        p2.mu() = liquids_.mu(p2.pc(), p2.T(), X);
                        p2.d() = cbrt(6.0*m2/(p2.nParticle()*p2.rho()*pi));
                    }
                }
            }
            j++;
        }
        i++;
    }

    // remove coalesced parcels that fall below minimum mass threshold
    forAllIter(typename CloudType, this->owner(), iter)
    {
        parcelType& p = iter();
        scalar mass = p.nParticle()*p.mass();

        if (mass < this->owner().constProps().minParticleMass())
        {
            this->owner().deleteParticle(p);
        }
    }
}


template<class CloudType>
bool Foam::ORourkeCollision<CloudType>::collideParcels
(
    const scalar dt,
    parcelType& p1,
    parcelType& p2,
    scalar& m1,
    scalar& m2
)
{
    const label cell1 = p1.cell();
    const label cell2 = p2.cell();

    // check if parcels belong to same cell
    if ((cell1 != cell2) || (m1 < ROOTVSMALL) || (m2 < ROOTVSMALL))
    {
        return false;
    }

    bool coalescence = false;

    const scalar Vc = this->owner().mesh().V()[cell1];
    const scalar d1 = p1.d();
    const scalar d2 = p2.d();

    scalar magUrel = mag(p1.U() - p2.U());
    scalar sumD = d1 + d2;
    scalar nu0 = 0.25*constant::mathematical::pi*sqr(sumD)*magUrel*dt/Vc;
    scalar nMin = min(p1.nParticle(), p2.nParticle());
    scalar nu = nMin*nu0;
    scalar collProb = exp(-nu);
    scalar xx = this->owner().rndGen().template sample01<scalar>();

    // collision occurs
    if (xx > collProb)
    {
        if (d1 > d2)
        {
            coalescence = collideSorted(dt, p1, p2, m1, m2);
        }
        else
        {
            coalescence = collideSorted(dt, p2, p1, m2, m1);
        }
    }

    return coalescence;
}


template<class CloudType>
bool Foam::ORourkeCollision<CloudType>::collideSorted
(
    const scalar dt,
    parcelType& p1,
    parcelType& p2,
    scalar& m1,
    scalar& m2
)
{
    bool coalescence = false;

    const scalar sigma1 = p1.sigma();
    const scalar sigma2 = p2.sigma();

    const scalar d1 = p1.d();
    const scalar d2 = p2.d();

    const scalar T1 = p1.T();
    const scalar T2 = p2.T();

    const scalar rho1 = p1.rho();
    const scalar rho2 = p2.rho();

    const vector& U1 = p1.U();
    const vector& U2 = p2.U();

    const label& nP1 = p1.nParticle();
    const label& nP2 = p2.nParticle();


    vector URel = U1 - U2;
    scalar magURel = mag(URel);

    scalar mTot = m1 + m2;

    scalar gamma = d1/max(ROOTVSMALL, d2);
    scalar f = pow3(gamma) + 2.7*gamma - 2.4*sqr(gamma);

    // mass-averaged temperature
    scalar Tave = (T1*m1 + T2*m2)/mTot;

    // interpolate to find average surface tension
    scalar sigmaAve = sigma1 + (sigma2 - sigma1)*(Tave - T1)/(T2 - T1);

    scalar Vtot = m1/rho1 + m2/rho2;
    scalar rhoAve = mTot/Vtot;

    scalar dAve = sqrt(d1*d2);
    scalar WeColl = 0.5*rhoAve*sqr(magURel)*dAve/max(ROOTVSMALL, sigmaAve);

    scalar coalesceProb = min(1.0, 2.4*f/max(ROOTVSMALL, WeColl));

    scalar prob = this->owner().rndGen().template sample01<scalar>();

    // Coalescence
    if (prob < coalesceProb && coalescence_)
    {
        coalescence = true;

        // number of the droplets that coalesce
        scalar nProb = prob*nP2/nP1;

        // conservation of mass, momentum and energy
        scalar m1Org = m1;
        scalar m2Org = m2;

        scalar dm = nP1*nProb*m2/scalar(nP2);

        m1 += dm;
        m2 -= dm;

        p1.T() = (Tave*mTot - m2*T2)/m1;

        p1.U() = (m1*U1 + (1.0 - m2/m2Org)*m2*U2)/m1;

        p1.Y() = (m1Org*p1.Y() + dm*p2.Y())/m1;

        p2.nParticle() = m2/(rho2*p2.volume());
    }
    // Grazing collision (no coalescence)
    else
    {
        scalar gf = sqrt(prob) - sqrt(coalesceProb);
        scalar denom = 1.0 - sqrt(coalesceProb);
        if (denom < 1.0e-5)
        {
            denom = 1.0;
        }
        gf /= denom;

        // if gf negative, this means that coalescence is turned off
        // and these parcels should have coalesced
        gf = max(0.0, gf);

        // gf -> 1 => v1p -> U1 ...
        // gf -> 0 => v1p -> momentum/mTot

        vector mr = m1*U1 + m2*U2;
        vector v1p = (mr + m2*gf*URel)/mTot;
        vector v2p = (mr - m1*gf*URel)/mTot;

        if (nP1 < nP2)
        {
            p1.U() = v1p;
            p2.U() = (nP1*v2p + (nP2 - nP1)*U2)/nP2;
        }
        else
        {
            p1.U() = (nP2*v1p + (nP1 - nP2)*U1)/nP1;
            p2.U() = v2p;
        }
    }

    return coalescence;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ORourkeCollision<CloudType>::ORourkeCollision
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    StochasticCollisionModel<CloudType>(dict, owner, modelName),
    liquids_
    (
        owner.db().template lookupObject<SLGThermo>("SLGThermo").liquids()
    ),
    coalescence_(this->coeffDict().lookup("coalescence"))
{}


template<class CloudType>
Foam::ORourkeCollision<CloudType>::ORourkeCollision
(
    const ORourkeCollision<CloudType>& cm
)
:
    StochasticCollisionModel<CloudType>(cm),
    liquids_(cm.liquids_),
    coalescence_(cm.coalescence_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ORourkeCollision<CloudType>::~ORourkeCollision()
{}


// ************************************************************************* //
