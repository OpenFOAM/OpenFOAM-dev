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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::ReactingMultiphaseCloud<CloudType>&
Foam::ReactingMultiphaseCloud<CloudType>::cloudCopy() const
{
    return cloudCopyPtr_();
}


template<class CloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::ReactingMultiphaseCloud<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>
inline typename CloudType::particleType::constantProperties&
Foam::ReactingMultiphaseCloud<CloudType>::constProps()
{
    return constProps_;
}


template<class CloudType>
inline const Foam::DevolatilisationModel
<
    Foam::ReactingMultiphaseCloud<CloudType>
>&
Foam::ReactingMultiphaseCloud<CloudType>::devolatilisation() const
{
    return devolatilisationModel_;
}


template<class CloudType>
inline Foam::DevolatilisationModel
<
    Foam::ReactingMultiphaseCloud<CloudType>
>&
Foam::ReactingMultiphaseCloud<CloudType>::devolatilisation()
{
    return devolatilisationModel_();
}


template<class CloudType>
inline const Foam::SurfaceReactionModel
<
    Foam::ReactingMultiphaseCloud<CloudType>
>&
Foam::ReactingMultiphaseCloud<CloudType>::surfaceReaction() const
{
    return surfaceReactionModel_;
}


template<class CloudType>
inline Foam::SurfaceReactionModel
<
    Foam::ReactingMultiphaseCloud<CloudType>
>&
Foam::ReactingMultiphaseCloud<CloudType>::surfaceReaction()
{
    return surfaceReactionModel_();
}


// ************************************************************************* //
