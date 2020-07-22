/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "SHF.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SHF<CloudType>::SHF
(
    const dictionary& dict,
    CloudType& owner
)
:
    BreakupModel<CloudType>(dict, owner, typeName),
    weCorrCoeff_(this->coeffDict().template lookup<scalar>("weCorrCoeff")),
    weBuCrit_(this->coeffDict().template lookup<scalar>("weBuCrit")),
    weBuBag_(this->coeffDict().template lookup<scalar>("weBuBag")),
    weBuMM_(this->coeffDict().template lookup<scalar>("weBuMM")),
    ohnCoeffCrit_(this->coeffDict().template lookup<scalar>("ohnCoeffCrit")),
    ohnCoeffBag_(this->coeffDict().template lookup<scalar>("ohnCoeffBag")),
    ohnCoeffMM_(this->coeffDict().template lookup<scalar>("ohnCoeffMM")),
    ohnExpCrit_(this->coeffDict().template lookup<scalar>("ohnExpCrit")),
    ohnExpBag_(this->coeffDict().template lookup<scalar>("ohnExpBag")),
    ohnExpMM_(this->coeffDict().template lookup<scalar>("ohnExpMM")),
    cInit_(this->coeffDict().template lookup<scalar>("Cinit")),
    c1_(this->coeffDict().template lookup<scalar>("C1")),
    c2_(this->coeffDict().template lookup<scalar>("C2")),
    c3_(this->coeffDict().template lookup<scalar>("C3")),
    cExp1_(this->coeffDict().template lookup<scalar>("Cexp1")),
    cExp2_(this->coeffDict().template lookup<scalar>("Cexp2")),
    cExp3_(this->coeffDict().template lookup<scalar>("Cexp3")),
    weConst_(this->coeffDict().template lookup<scalar>("Weconst")),
    weCrit1_(this->coeffDict().template lookup<scalar>("Wecrit1")),
    weCrit2_(this->coeffDict().template lookup<scalar>("Wecrit2")),
    coeffD_(this->coeffDict().template lookup<scalar>("CoeffD")),
    onExpD_(this->coeffDict().template lookup<scalar>("OnExpD")),
    weExpD_(this->coeffDict().template lookup<scalar>("WeExpD")),
    mu_(this->coeffDict().template lookup<scalar>("mu")),
    sigma_(this->coeffDict().template lookup<scalar>("sigma")),
    d32Coeff_(this->coeffDict().template lookup<scalar>("d32Coeff")),
    cDmaxBM_(this->coeffDict().template lookup<scalar>("cDmaxBM")),
    cDmaxS_(this->coeffDict().template lookup<scalar>("cDmaxS")),
    corePerc_(this->coeffDict().template lookup<scalar>("corePerc"))
{}


template<class CloudType>
Foam::SHF<CloudType>::SHF(const SHF<CloudType>& bum)
:
    BreakupModel<CloudType>(bum),
    weCorrCoeff_(bum.weCorrCoeff_),
    weBuCrit_(bum.weBuCrit_),
    weBuBag_(bum.weBuBag_),
    weBuMM_(bum.weBuMM_),
    ohnCoeffCrit_(bum.ohnCoeffCrit_),
    ohnCoeffBag_(bum.ohnCoeffBag_),
    ohnCoeffMM_(bum.ohnCoeffMM_),
    ohnExpCrit_(bum.ohnExpCrit_),
    ohnExpBag_(bum.ohnExpBag_),
    ohnExpMM_(bum.ohnExpMM_),
    cInit_(bum.cInit_),
    c1_(bum.c1_),
    c2_(bum.c2_),
    c3_(bum.c3_),
    cExp1_(bum.cExp1_),
    cExp2_(bum.cExp2_),
    cExp3_(bum.cExp3_),
    weConst_(bum.weConst_),
    weCrit1_(bum.weCrit1_),
    weCrit2_(bum.weCrit2_),
    coeffD_(bum.coeffD_),
    onExpD_(bum.onExpD_),
    weExpD_(bum.weExpD_),
    mu_(bum.mu_),
    sigma_(bum.sigma_),
    d32Coeff_(bum.d32Coeff_),
    cDmaxBM_(bum.cDmaxBM_),
    cDmaxS_(bum.cDmaxS_),
    corePerc_(bum.corePerc_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SHF<CloudType>::~SHF()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::SHF<CloudType>::update
(
    const scalar dt,
    const vector& g,
    scalar& d,
    scalar& tc,
    scalar& ms,
    scalar& nParticle,
    scalar& KHindex,
    scalar& y,
    scalar& yDot,
    const scalar d0,
    const scalar rho,
    const scalar mu,
    const scalar sigma,
    const vector& U,
    const scalar rhoc,
    const scalar muc,
    const vector& Urel,
    const scalar Urmag,
    const scalar tMom,
    scalar& dChild,
    scalar& massChild
)
{
    Random& rndGen = this->owner().rndGen();

    bool addChild = false;

    scalar d03 = pow3(d);
    scalar rhopi6 = rho*constant::mathematical::pi/6.0;
    scalar mass0 = nParticle*rhopi6*d03;
    scalar mass = mass0;

    scalar weGas = 0.5*rhoc*sqr(Urmag)*d/sigma;
    scalar weLiquid = 0.5*rho*sqr(Urmag)*d/sigma;

    // correct the Reynolds number. Reitz is using radius instead of diameter
    scalar reLiquid   = 0.5*Urmag*d/mu;
    scalar ohnesorge  = sqrt(weLiquid)/(reLiquid + vSmall);

    scalar weGasCorr = weGas/(1.0 + weCorrCoeff_*ohnesorge);

    // update the droplet characteristic time
    tc += dt;

    // droplet deformation characteristic rate
    scalar rChar = Urmag/d*sqrt(rhoc/rho);

    // return if the characteristic deformation rate is too low for the
    // following modelling to be calculable
    if (tc*rChar < small)
    {
        return false;
    }

    // droplet deformation characteristic time
    scalar tChar = 1/rChar;

    scalar tFirst = cInit_*tChar;

    scalar tSecond = 0;
    scalar tCharSecond = 0;

    bool bag = false;
    bool multimode = false;
    bool shear = false;
    bool success = false;


    if (weGas > weConst_)
    {
        if (weGas < weCrit1_)
        {
            tCharSecond = c1_*pow((weGas - weConst_), cExp1_);
        }
        else if (weGas >= weCrit1_ && weGas <= weCrit2_)
        {
            tCharSecond = c2_*pow((weGas - weConst_), cExp2_);
        }
        else
        {
            tCharSecond = c3_*pow((weGas - weConst_), cExp3_);
        }
    }

    scalar weC = weBuCrit_*(1.0 + ohnCoeffCrit_*pow(ohnesorge, ohnExpCrit_));
    scalar weB = weBuBag_*(1.0 + ohnCoeffBag_*pow(ohnesorge, ohnExpBag_));
    scalar weMM = weBuMM_*(1.0 + ohnCoeffMM_*pow(ohnesorge, ohnExpMM_));

    if (weGas > weC && weGas < weB)
    {
        bag = true;
    }

    if (weGas >= weB && weGas <= weMM)
    {
        multimode = true;
    }

    if (weGas > weMM)
    {
        shear = true;
    }

    tSecond = tCharSecond*tChar;

    scalar tBreakUP = tFirst + tSecond;
    if (tc > tBreakUP)
    {
        scalar d32 = coeffD_*d*pow(ohnesorge, onExpD_)*pow(weGasCorr, weExpD_);

        if (bag || multimode)
        {
            scalar d05 = d32Coeff_ * d32;

            scalar x = 0.0;
            scalar yGuess = 0.0;
            scalar dGuess = 0.0;

            while(!success)
            {
                x = cDmaxBM_*rndGen.sample01<scalar>();
                dGuess = sqr(x)*d05;
                yGuess = rndGen.sample01<scalar>();

                scalar p =
                    x
                   /(2.0*sqrt(constant::mathematical::twoPi)*sigma_)
                   *exp(-0.5*sqr((x - mu_)/sigma_));

                if (yGuess < p)
                {
                    success = true;
                }
            }

            d = dGuess;
            tc = 0.0;
        }

        if (shear)
        {
            scalar dC = weConst_*sigma/(rhoc*sqr(Urmag));
            scalar d32Red = 4.0*(d32*dC)/(5.0*dC - d32);

            scalar d05 = d32Coeff_ * d32Red;

            scalar x = 0.0;
            scalar yGuess = 0.0;
            scalar dGuess = 0.0;

            while(!success)
            {
                x = cDmaxS_*rndGen.sample01<scalar>();
                dGuess = sqr(x)*d05;
                yGuess = rndGen.sample01<scalar>();

                scalar p =
                    x
                   /(2.0*sqrt(constant::mathematical::twoPi)*sigma_)
                   *exp(-0.5*sqr((x - mu_)/sigma_));

                if (yGuess<p)
                {
                    success = true;
                }
            }

            d = dC;
            dChild = dGuess;
            massChild = corePerc_*mass0;
            mass -= massChild;

            addChild = true;
            // reset timer
            tc = 0.0;
        }

        // correct nParticle to conserve mass
        nParticle = mass/(rhopi6*pow3(d));
    }

    return addChild;
}


// ************************************************************************* //
