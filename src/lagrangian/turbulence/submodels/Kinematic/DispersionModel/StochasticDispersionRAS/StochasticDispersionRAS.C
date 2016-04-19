/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "StochasticDispersionRAS.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::StochasticDispersionRAS<CloudType>::StochasticDispersionRAS
(
    const dictionary& dict,
    CloudType& owner
)
:
    DispersionRASModel<CloudType>(dict, owner)
{}


template<class CloudType>
Foam::StochasticDispersionRAS<CloudType>::StochasticDispersionRAS
(
    const StochasticDispersionRAS<CloudType>& dm
)
:
    DispersionRASModel<CloudType>(dm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::StochasticDispersionRAS<CloudType>::~StochasticDispersionRAS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::vector Foam::StochasticDispersionRAS<CloudType>::update
(
    const scalar dt,
    const label cellI,
    const vector& U,
    const vector& Uc,
    vector& UTurb,
    scalar& tTurb
)
{
    cachedRandom& rnd = this->owner().rndGen();

    const scalar cps = 0.16432;

    const scalar k = this->kPtr_->internalField()[cellI];
    const scalar epsilon =
        this->epsilonPtr_->internalField()[cellI] + ROOTVSMALL;

    const scalar UrelMag = mag(U - Uc - UTurb);

    const scalar tTurbLoc =
        min(k/epsilon, cps*pow(k, 1.5)/epsilon/(UrelMag + SMALL));


    // Parcel is perturbed by the turbulence
    if (dt < tTurbLoc)
    {
        tTurb += dt;

        if (tTurb > tTurbLoc)
        {
            tTurb = 0;

            const scalar sigma = sqrt(2*k/3.0);

            // Calculate a random direction dir distributed uniformly
            // in spherical coordinates

            const scalar theta = rnd.sample01<scalar>()*twoPi;
            const scalar u = 2*rnd.sample01<scalar>() - 1;

            const scalar a = sqrt(1 - sqr(u));
            const vector dir(a*cos(theta), a*sin(theta), u);

            UTurb = sigma*mag(rnd.GaussNormal<scalar>())*dir;
        }
    }
    else
    {
        tTurb = GREAT;
        UTurb = Zero;
    }

    return Uc + UTurb;
}


// ************************************************************************* //
