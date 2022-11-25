/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "CollidingParcel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::CollidingParcel<ParcelType>::CollidingParcel
(
    const CollidingParcel<ParcelType>& p
)
:
    ParcelType(p),
    f_(p.f_),
    angularMomentum_(p.angularMomentum_),
    torque_(p.torque_),
    collisionRecords_(p.collisionRecords_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
bool Foam::CollidingParcel<ParcelType>::move
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    typename TrackCloudType::parcelType& p =
        static_cast<typename TrackCloudType::parcelType&>(*this);

    switch (td.part())
    {
        case trackingData::tpVelocityHalfStep:
        {
            // First and last leapfrog velocity adjust part, required
            // before and after tracking and force calculation

            const Pair<scalar>& sfr = td.stepFractionRange();
            const scalar dt = (sfr.second() - sfr.first())*td.trackTime()/2;

            p.U() += dt*p.f()/p.mass();
            p.angularMomentum() += dt*p.torque();

            td.keepParticle = true;
            td.sendToProc = -1;

            break;
        }

        case trackingData::tpLinearTrack:
        {
            ParcelType::move(cloud, td);

            break;
        }

        case trackingData::tpRotationalTrack:
        {
            NotImplemented;

            break;
        }
    }

    return td.keepParticle;
}


template<class ParcelType>
void Foam::CollidingParcel<ParcelType>::transformProperties
(
    const transformer& transform
)
{
    ParcelType::transformProperties(transform);
    f_ = transform.transform(f_);
    angularMomentum_ = transform.transform(angularMomentum_);
    torque_ = transform.transform(torque_);
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "CollidingParcelIO.C"

// ************************************************************************* //
