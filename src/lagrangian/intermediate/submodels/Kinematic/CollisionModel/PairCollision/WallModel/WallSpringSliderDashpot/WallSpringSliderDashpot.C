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

#include "WallSpringSliderDashpot.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::WallSpringSliderDashpot<CloudType>::findMinMaxProperties
(
    scalar& rMin,
    scalar& rhoMax,
    scalar& UMagMax
) const
{
    rMin = vGreat;
    rhoMax = -vGreat;
    UMagMax = -vGreat;

    forAllConstIter(typename CloudType, this->owner(), iter)
    {
        const typename CloudType::parcelType& p = iter();

        // Finding minimum diameter to avoid excessive arithmetic

        scalar dEff = p.d();

        if (useEquivalentSize_)
        {
            dEff *= cbrt(p.nParticle()*volumeFactor_);
        }

        rMin = min(dEff, rMin);

        rhoMax = max(p.rho(), rhoMax);

        UMagMax = max
        (
            mag(p.U()) + mag(p.omega())*dEff/2,
            UMagMax
        );
    }

    // Transform the minimum diameter into minimum radius
    //     rMin = dMin/2

    rMin /= 2.0;
}


template<class CloudType>
void Foam::WallSpringSliderDashpot<CloudType>::evaluateWall
(
    typename CloudType::parcelType& p,
    const point& site,
    const WallSiteData<vector>& data,
    scalar pREff,
    scalar kN,
    bool cohesion
) const
{
    vector r_PW = p.position() - site;

    vector U_PW = p.U() - data.wallData();

    scalar r_PW_mag = mag(r_PW);

    scalar normalOverlapMag = max(pREff - r_PW_mag, 0.0);

    vector rHat_PW = r_PW/(r_PW_mag + vSmall);

    scalar etaN = alpha_*sqrt(p.mass()*kN)*pow025(normalOverlapMag);

    vector fN_PW =
        rHat_PW
       *(kN*pow(normalOverlapMag, b_) - etaN*(U_PW & rHat_PW));

    // Cohesion force, energy density multiplied by the area of wall/particle
    // overlap
    if (cohesion)
    {
        fN_PW +=
           -cohesionEnergyDensity_
           *mathematical::pi*(sqr(pREff) - sqr(r_PW_mag))
           *rHat_PW;
    }

    p.f() += fN_PW;

    vector USlip_PW =
        U_PW - (U_PW & rHat_PW)*rHat_PW
      + (p.omega() ^ (pREff*-rHat_PW));

    scalar deltaT = this->owner().mesh().time().deltaTValue();

    vector& tangentialOverlap_PW =
        p.collisionRecords().matchWallRecord(-r_PW, pREff).collisionData();

    tangentialOverlap_PW += USlip_PW*deltaT;

    scalar tangentialOverlapMag = mag(tangentialOverlap_PW);

    if (tangentialOverlapMag > vSmall)
    {
        scalar kT = 8.0*sqrt(pREff*normalOverlapMag)*Gstar_;

        scalar etaT = etaN;

        // Tangential force
        vector fT_PW;

        if (kT*tangentialOverlapMag > mu_*mag(fN_PW))
        {
            // Tangential force greater than sliding friction,
            // particle slips

            fT_PW = -mu_*mag(fN_PW)*USlip_PW/mag(USlip_PW);

            tangentialOverlap_PW = Zero;
        }
        else
        {
            fT_PW = - kT*tangentialOverlap_PW - etaT*USlip_PW;
        }

        p.f() += fT_PW;

        p.torque() += (pREff*-rHat_PW) ^ fT_PW;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::WallSpringSliderDashpot<CloudType>::WallSpringSliderDashpot
(
    const dictionary& dict,
    CloudType& cloud
)
:
    WallModel<CloudType>(dict, cloud, typeName),
    Estar_(),
    Gstar_(),
    alpha_(readScalar(this->coeffDict().lookup("alpha"))),
    b_(readScalar(this->coeffDict().lookup("b"))),
    mu_(readScalar(this->coeffDict().lookup("mu"))),
    cohesionEnergyDensity_
    (
        readScalar(this->coeffDict().lookup("cohesionEnergyDensity"))
    ),
    cohesion_(false),
    collisionResolutionSteps_
    (
        readScalar
        (
            this->coeffDict().lookup("collisionResolutionSteps")
        )
    ),
    volumeFactor_(1.0),
    useEquivalentSize_(Switch(this->coeffDict().lookup("useEquivalentSize")))
{
    if (useEquivalentSize_)
    {
        volumeFactor_ = readScalar(this->coeffDict().lookup("volumeFactor"));
    }

    scalar nu = readScalar(this->coeffDict().lookup("poissonsRatio"));

    scalar E = readScalar(this->coeffDict().lookup("youngsModulus"));

    scalar pNu = this->owner().constProps().poissonsRatio();

    scalar pE = this->owner().constProps().youngsModulus();

    Estar_ = 1/((1 - sqr(pNu))/pE + (1 - sqr(nu))/E);

    Gstar_ = 1/(2*((2 + pNu - sqr(pNu))/pE + (2 + nu - sqr(nu))/E));

    cohesion_ = (mag(cohesionEnergyDensity_) > vSmall);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::WallSpringSliderDashpot<CloudType>::~WallSpringSliderDashpot()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::WallSpringSliderDashpot<CloudType>::pREff
(
    const typename CloudType::parcelType& p
) const
{
    if (useEquivalentSize_)
    {
        return p.d()/2*cbrt(p.nParticle()*volumeFactor_);
    }
    else
    {
        return p.d()/2;
    }
}


template<class CloudType>
bool Foam::WallSpringSliderDashpot<CloudType>::controlsTimestep() const
{
    return true;
}


template<class CloudType>
Foam::label Foam::WallSpringSliderDashpot<CloudType>::nSubCycles() const
{
    if (!(this->owner().size()))
    {
        return 1;
    }

    scalar rMin;
    scalar rhoMax;
    scalar UMagMax;

    findMinMaxProperties(rMin, rhoMax, UMagMax);

    // Note:  pi^(7/5)*(5/4)^(2/5) = 5.429675
    scalar minCollisionDeltaT =
        5.429675
       *rMin
       *pow(rhoMax/(Estar_*sqrt(UMagMax) + vSmall), 0.4)
       /collisionResolutionSteps_;

    return ceil(this->owner().time().deltaTValue()/minCollisionDeltaT);
}


template<class CloudType>
void Foam::WallSpringSliderDashpot<CloudType>::evaluateWall
(
    typename CloudType::parcelType& p,
    const List<point>& flatSitePoints,
    const List<WallSiteData<vector>>& flatSiteData,
    const List<point>& sharpSitePoints,
    const List<WallSiteData<vector>>& sharpSiteData
) const
{
    scalar pREff = this->pREff(p);

    scalar kN = (4.0/3.0)*sqrt(pREff)*Estar_;

    forAll(flatSitePoints, siteI)
    {
        evaluateWall
        (
            p,
            flatSitePoints[siteI],
            flatSiteData[siteI],
            pREff,
            kN,
            cohesion_
        );
    }

    forAll(sharpSitePoints, siteI)
    {
        // Treating sharp sites like flat sites, except suppress cohesion

        evaluateWall
        (
            p,
            sharpSitePoints[siteI],
            sharpSiteData[siteI],
            pREff,
            kN,
            false
        );
    }
}


// ************************************************************************* //
