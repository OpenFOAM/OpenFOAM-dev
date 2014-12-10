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

#include "SRFForce.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SRFForce<CloudType>::SRFForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, false),
    srfPtr_(NULL)
{}


template<class CloudType>
Foam::SRFForce<CloudType>::SRFForce
(
    const SRFForce& srff
)
:
    ParticleForce<CloudType>(srff),
    srfPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SRFForce<CloudType>::~SRFForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::SRFForce<CloudType>::cacheFields(const bool store)
{
    if (store)
    {
        const typename SRF::SRFModel& model = this->mesh().template
            lookupObject<SRF::SRFModel>("SRFProperties");
        srfPtr_ = &model;
    }
    else
    {
        srfPtr_ = NULL;
    }
}


template<class CloudType>
Foam::forceSuSp Foam::SRFForce<CloudType>::calcNonCoupled
(
    const typename CloudType::parcelType& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(vector::zero, 0.0);

    const typename SRF::SRFModel& srf = *srfPtr_;

    const vector& omega = srf.omega().value();

    const vector& r = p.position();

    // Coriolis and centrifugal acceleration terms
    value.Su() =
        mass*(1.0 - p.rhoc()/p.rho())
       *(2.0*(p.U() ^ omega) + (omega ^ (r ^ omega)));

    return value;
}


// ************************************************************************* //
