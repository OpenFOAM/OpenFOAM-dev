/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "GradientDispersionRAS.H"
#include "demandDrivenData.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::GradientDispersionRAS<CloudType>::GradientDispersionRAS
(
    const dictionary& dict,
    CloudType& owner
)
:
    DispersionRASModel<CloudType>(dict, owner),
    gradkPtr_(nullptr),
    ownGradK_(false)
{}


template<class CloudType>
Foam::GradientDispersionRAS<CloudType>::GradientDispersionRAS
(
    const GradientDispersionRAS<CloudType>& dm
)
:
    DispersionRASModel<CloudType>(dm),
    gradkPtr_(dm.gradkPtr_),
    ownGradK_(dm.ownGradK_)
{
    dm.ownGradK_ = false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::GradientDispersionRAS<CloudType>::~GradientDispersionRAS()
{
    cacheFields(false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::GradientDispersionRAS<CloudType>::cacheFields(const bool store)
{
    DispersionRASModel<CloudType>::cacheFields(store);

    if (store)
    {
        gradkPtr_ = fvc::grad(*this->kPtr_).ptr();
        ownGradK_ = true;
    }
    else
    {
        if (ownGradK_)
        {
            deleteDemandDrivenData(gradkPtr_);
            gradkPtr_ = nullptr;
            ownGradK_ = false;
        }
    }
}


template<class CloudType>
Foam::vector Foam::GradientDispersionRAS<CloudType>::update
(
    const scalar dt,
    const label celli,
    const vector& U,
    const vector& Uc,
    vector& UTurb,
    scalar& tTurb
)
{
    Random& rnd = this->owner().rndGen();

    const scalar cps = 0.16432;

    const scalar k = this->kPtr_->primitiveField()[celli];
    const scalar epsilon =
        this->epsilonPtr_->primitiveField()[celli] + rootVSmall;
    const vector& gradk = this->gradkPtr_->primitiveField()[celli];

    const scalar UrelMag = mag(U - Uc - UTurb);

    const scalar tTurbLoc =
        min(k/epsilon, cps*pow(k, 1.5)/epsilon/(UrelMag + small));


    // Parcel is perturbed by the turbulence
    if (dt < tTurbLoc)
    {
        tTurb += dt;

        if (tTurb > tTurbLoc)
        {
            tTurb = 0.0;

            const scalar sigma = sqrt(2.0*k/3.0);
            const vector dir = -gradk/(mag(gradk) + small);

            scalar fac = 0.0;

            // In 2D calculations the -grad(k) is always
            // away from the axis of symmetry
            // This creates a 'hole' in the spray and to
            // prevent this we let fac be both negative/positive
            if (this->owner().mesh().nSolutionD() == 2)
            {
                fac = rnd.scalarNormal();
            }
            else
            {
                fac = mag(rnd.scalarNormal());
            }

            UTurb = sigma*fac*dir;
        }
    }
    else
    {
        tTurb = great;
        UTurb = Zero;
    }

    return Uc + UTurb;
}


// ************************************************************************* //
