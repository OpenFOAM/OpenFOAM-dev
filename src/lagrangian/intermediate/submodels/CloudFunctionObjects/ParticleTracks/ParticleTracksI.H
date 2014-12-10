/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "ParticleTracks.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CloudType>
inline Foam::label Foam::ParticleTracks<CloudType>::trackInterval() const
{
    return trackInterval_;
}


template<class CloudType>
inline Foam::label Foam::ParticleTracks<CloudType>::maxSamples() const
{
    return maxSamples_;
}


template<class CloudType>
inline const Foam::Switch& Foam::ParticleTracks<CloudType>::resetOnWrite() const
{
    return resetOnWrite_;
}


template<class CloudType>
inline const typename Foam::ParticleTracks<CloudType>::hitTableType&
Foam::ParticleTracks<CloudType>::faceHitCounter() const
{
    return faceHitCounter_;
}


template<class CloudType>
inline const Foam::Cloud<typename CloudType::parcelType>&
Foam::ParticleTracks<CloudType>::cloud() const
{
    return cloudPtr_();
}


// ************************************************************************* //
