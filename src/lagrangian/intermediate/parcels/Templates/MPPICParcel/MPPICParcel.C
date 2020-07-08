/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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

#include "MPPICParcel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::MPPICParcel<ParcelType>::MPPICParcel
(
    const MPPICParcel<ParcelType>& p
)
:
    ParcelType(p),
    UCorrect_(p.UCorrect_)
{}


template<class ParcelType>
Foam::MPPICParcel<ParcelType>::MPPICParcel
(
    const MPPICParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    UCorrect_(p.UCorrect_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
bool Foam::MPPICParcel<ParcelType>::move
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar trackTime
)
{
    typename TrackCloudType::parcelType& p =
        static_cast<typename TrackCloudType::parcelType&>(*this);

    switch (td.part())
    {
        case trackingData::tpPredictTrack:
        {
            ParcelType::move(cloud, td, trackTime);

            break;
        }
        case trackingData::tpDampingNoTrack:
        {
            p.UCorrect() =
                cloud.dampingModel().velocityCorrection(p, trackTime);

            td.keepParticle = true;
            td.switchProcessor = false;

            break;
        }
        case trackingData::tpPackingNoTrack:
        {
            p.UCorrect() =
                cloud.packingModel().velocityCorrection(p, trackTime);

            td.keepParticle = true;
            td.switchProcessor = false;

            break;
        }
        case trackingData::tpCorrectTrack:
        {
            const scalar f = p.stepFraction();
            const scalar a = p.age();

            Swap(p.U(), p.UCorrect());

            ParcelType::move(cloud, td, trackTime);

            Swap(p.U(), p.UCorrect());

            p.U() += (p.stepFraction() - f)*p.UCorrect();

            p.age() = a;

            break;
        }
    }

    return td.keepParticle;
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "MPPICParcelIO.C"

// ************************************************************************* //
