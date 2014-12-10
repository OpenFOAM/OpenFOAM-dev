/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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
template<class TrackData>
bool Foam::MPPICParcel<ParcelType>::move
(
    TrackData& td,
    const scalar trackTime
)
{
    typename TrackData::cloudType::parcelType& p =
        static_cast<typename TrackData::cloudType::parcelType&>(*this);

    switch (td.part())
    {
        case TrackData::tpLinearTrack:
        {
            ParcelType::move(td, trackTime);

            break;
        }
        case TrackData::tpDampingNoTrack:
        {
            p.UCorrect() =
                td.cloud().dampingModel().velocityCorrection(p, trackTime);

            td.keepParticle = true;

            break;
        }
        case TrackData::tpPackingNoTrack:
        {
            p.UCorrect() =
                td.cloud().packingModel().velocityCorrection(p, trackTime);

            td.keepParticle = true;

            break;
        }
        case TrackData::tpCorrectTrack:
        {
            vector U = p.U();

            scalar f = p.stepFraction();

            scalar a = p.age();

            p.U() = (1.0 - f)*p.UCorrect();

            ParcelType::move(td, trackTime);

            p.U() = U + (p.stepFraction() - f)*p.UCorrect();

            p.age() = a;

            break;
        }
    }

    return td.keepParticle;
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "MPPICParcelIO.C"

// ************************************************************************* //
