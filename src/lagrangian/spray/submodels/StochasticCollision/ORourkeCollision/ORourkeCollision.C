/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
#include "SLGThermo.H"
#include "CompactListList.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::ORourkeCollision<CloudType>::collide
(
    typename CloudType::parcelType::trackingData& td,
    const scalar dt
)
{
    // Create the occupancy list for the cells
    labelList occupancy(this->owner().mesh().nCells(), 0);
    forAllIter(typename CloudType, this->owner(), iter)
    {
        occupancy[iter().cell()]++;
    }

    // Initialize the sizes of the lists of parcels in each cell
    CompactListList<parcelType*> pInCell(occupancy);

    // Reset the occupancy to use as a counter
    occupancy = 0;

    // Set the parcel pointer lists for each cell
    forAllIter(typename CloudType, this->owner(), iter)
    {
        pInCell(iter().cell(), occupancy[iter().cell()]++) = &iter();
    }

    for (label celli=0; celli<this->owner().mesh().nCells(); celli++)
    {
        UList<parcelType*> pInCelli(pInCell[celli]);

        if (pInCelli.size() >= 2)
        {
            forAll(pInCelli, i)
            {
                for (label j=i+1; j<pInCelli.size(); j++)
                {
                    parcelType& p1 = *pInCelli[i];
                    parcelType& p2 = *pInCelli[j];

                    scalar m1 = p1.nParticle()*p1.mass();
                    scalar m2 = p2.nParticle()*p2.mass();

                    bool massChanged = collideParcels(dt, p1, p2, m1, m2);

                    if (massChanged)
                    {
                        if (m1 > rootVSmall)
                        {
                            const scalarField X(liquids_.X(p1.Y()));
                            p1.setCellValues(this->owner(), td);
                            p1.rho() = liquids_.rho(td.pc(), p1.T(), X);
                            p1.Cp() = liquids_.Cp(td.pc(), p1.T(), X);
                            p1.sigma() = liquids_.sigma(td.pc(), p1.T(), X);
                            p1.mu() = liquids_.mu(td.pc(), p1.T(), X);
                            p1.d() = cbrt(6.0*m1/(p1.nParticle()*p1.rho()*pi));
                        }

                        if (m2 > rootVSmall)
                        {
                            const scalarField X(liquids_.X(p2.Y()));
                            p2.setCellValues(this->owner(), td);
                            p2.rho() = liquids_.rho(td.pc(), p2.T(), X);
                            p2.Cp() = liquids_.Cp(td.pc(), p2.T(), X);
                            p2.sigma() = liquids_.sigma(td.pc(), p2.T(), X);
                            p2.mu() = liquids_.mu(td.pc(), p2.T(), X);
                            p2.d() = cbrt(6.0*m2/(p2.nParticle()*p2.rho()*pi));
                        }
                    }
                }
            }
        }
    }

    // Remove coalesced parcels that fall below minimum mass threshold
    forAllIter(typename CloudType, this->owner(), iter)
    {
        parcelType& p = iter();
        scalar mass = p.nParticle()*p.mass();

        if (mass < this->owner().constProps().minParcelMass())
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
    // Return if parcel masses are ~0
    if ((m1 < rootVSmall) || (m2 < rootVSmall))
    {
        return false;
    }

    const scalar Vc = this->owner().mesh().V()[p1.cell()];
    const scalar d1 = p1.d();
    const scalar d2 = p2.d();

    scalar magUrel = mag(p1.U() - p2.U());
    scalar sumD = d1 + d2;
    scalar nu0 = 0.25*constant::mathematical::pi*sqr(sumD)*magUrel*dt/Vc;
    scalar nMin = min(p1.nParticle(), p2.nParticle());
    scalar nu = nMin*nu0;
    scalar collProb = exp(-nu);
    scalar xx = this->owner().rndGen().template sample01<scalar>();

    // Collision occurs
    if (xx > collProb)
    {
        if (d1 > d2)
        {
            return collideSorted(dt, p1, p2, m1, m2);
        }
        else
        {
            return collideSorted(dt, p2, p1, m2, m1);
        }
    }
    else
    {
        return false;
    }
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
    const scalar nP1 = p1.nParticle();
    const scalar nP2 = p2.nParticle();

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


    vector URel = U1 - U2;
    scalar magURel = mag(URel);

    scalar mTot = m1 + m2;

    scalar gamma = d1/max(rootVSmall, d2);
    scalar f = pow3(gamma) + 2.7*gamma - 2.4*sqr(gamma);

    // Mass-averaged temperature
    scalar Tave = (T1*m1 + T2*m2)/mTot;

    // Interpolate to find average surface tension
    scalar sigmaAve = sigma1;
    if (mag(T2 - T1) > small)
    {
        sigmaAve += (sigma2 - sigma1)*(Tave - T1)/(T2 - T1);
    }

    scalar Vtot = m1/rho1 + m2/rho2;
    scalar rhoAve = mTot/Vtot;

    scalar dAve = sqrt(d1*d2);
    scalar WeColl = 0.5*rhoAve*sqr(magURel)*dAve/max(rootVSmall, sigmaAve);

    scalar coalesceProb = min(1.0, 2.4*f/max(rootVSmall, WeColl));

    scalar prob = this->owner().rndGen().template sample01<scalar>();

    // Coalescence
    if (coalescence_ && prob < coalesceProb)
    {
        // Number of the droplets that coalesce
        scalar nProb = prob*nP2/nP1;

        // Conservation of mass, momentum and energy
        scalar m1Org = m1;
        scalar m2Org = m2;

        scalar dm = nP1*nProb*m2/scalar(nP2);

        m1 += dm;
        m2 -= dm;

        p1.T() = (Tave*mTot - m2*T2)/m1;

        p1.U() = (m1*U1 + (1.0 - m2/m2Org)*m2*U2)/m1;

        p1.Y() = (m1Org*p1.Y() + dm*p2.Y())/m1;

        p2.nParticle() = m2/(rho2*p2.volume());

        return true;
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

        // If gf negative, this means that coalescence is turned off
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

        return false;
    }
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
