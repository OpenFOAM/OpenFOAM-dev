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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::MPPICParcel<ParcelType>::MPPICParcel
(
    const polyMesh& owner,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    ParcelType(owner, coordinates, celli, tetFacei, tetPti),
    UCorrect_(Zero)
{}


template<class ParcelType>
inline Foam::MPPICParcel<ParcelType>::MPPICParcel
(
    const polyMesh& owner,
    const vector& position,
    const label celli
)
:
    ParcelType(owner, position, celli),
    UCorrect_(Zero)
{}


template<class ParcelType>
inline Foam::MPPICParcel<ParcelType>::MPPICParcel
(
    const polyMesh& owner,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const scalar dTarget0,
    const vector& U0,
    const vector& UCorrect0,
    const typename ParcelType::constantProperties& constProps
)
:
    ParcelType
    (
        owner,
        coordinates,
        celli,
        tetFacei,
        tetPti,
        typeId,
        nParticle0,
        d0,
        dTarget0,
        U0,
        constProps
    ),
    UCorrect_(UCorrect0)
{}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

template<class ParcelType>
inline const Foam::vector& Foam::MPPICParcel<ParcelType>::UCorrect() const
{
    return UCorrect_;
}


template<class ParcelType>
inline Foam::vector& Foam::MPPICParcel<ParcelType>::UCorrect()
{
    return UCorrect_;
}


// ************************************************************************* //
