/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "NoInjection.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoInjection<CloudType>::NoInjection
(
    const dictionary&,
    CloudType& owner,
    const word&
)
:
    InjectionModel<CloudType>(owner)
{}


template<class CloudType>
Foam::NoInjection<CloudType>::NoInjection(const NoInjection<CloudType>& im)
:
    InjectionModel<CloudType>(im.owner_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoInjection<CloudType>::~NoInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::NoInjection<CloudType>::active() const
{
    return false;
}


template<class CloudType>
Foam::scalar Foam::NoInjection<CloudType>::timeEnd() const
{
    return 0.0;
}


template<class CloudType>
Foam::label Foam::NoInjection<CloudType>::parcelsToInject
(
    const scalar,
    const scalar
)
{
    return 0;
}


template<class CloudType>
Foam::scalar Foam::NoInjection<CloudType>::volumeToInject
(
    const scalar,
    const scalar
)
{
    return 0.0;
}


template<class CloudType>
void Foam::NoInjection<CloudType>::setPositionAndCell
(
    const label,
    const label,
    const scalar,
    vector&,
    label&,
    label&,
    label&
)
{}


template<class CloudType>
void Foam::NoInjection<CloudType>::setProperties
(
    const label,
    const label,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    // set particle velocity
    parcel.U() = vector::zero;

    // set particle diameter
    parcel.d() = 0.0;
}


template<class CloudType>
bool Foam::NoInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::NoInjection<CloudType>::validInjection(const label)
{
    return false;
}


// ************************************************************************* //
