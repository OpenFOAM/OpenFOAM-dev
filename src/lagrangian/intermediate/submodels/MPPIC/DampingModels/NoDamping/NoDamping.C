/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "NoDamping.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DampingModels::NoDamping<CloudType>::NoDamping
(
    const dictionary& dict,
    CloudType& owner
)
:
    DampingModel<CloudType>(owner)
{}


template<class CloudType>
Foam::DampingModels::NoDamping<CloudType>::NoDamping
(
    const NoDamping<CloudType>& cm
)
:
    DampingModel<CloudType>(cm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DampingModels::NoDamping<CloudType>::~NoDamping()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::vector Foam::DampingModels::NoDamping<CloudType>::velocityCorrection
(
    typename CloudType::parcelType& p,
    const scalar deltaT
) const
{
    return vector::zero;
}


template<class CloudType>
bool Foam::DampingModels::NoDamping<CloudType>::active() const
{
    return false;
}


// ************************************************************************* //
