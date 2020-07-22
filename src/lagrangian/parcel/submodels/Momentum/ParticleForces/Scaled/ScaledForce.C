/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2020 OpenFOAM Foundation
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

#include "ScaledForce.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::dictionary Foam::ScaledForce<CloudType>::modelDict
(
    const dictionary& dict
) const
{
    dictionary modelDict(dict);
    modelDict.add<word>("type", dict.lookup<word>("forceType"), true);
    return modelDict;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ScaledForce<CloudType>::ScaledForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    model_
    (
        ParticleForce<CloudType>::New
        (
            owner,
            mesh,
            modelDict(dict),
            dict.lookup<word>("forceType")
        )
    ),
    factor_(this->coeffs().template lookup<scalar>("factor"))
{}


template<class CloudType>
Foam::ScaledForce<CloudType>::ScaledForce
(
    const ScaledForce<CloudType>& df
)
:
    ParticleForce<CloudType>(df),
    model_(nullptr),
    factor_(1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ScaledForce<CloudType>::~ScaledForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::ScaledForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    return factor_*model_->calcCoupled(p, td, dt, mass, Re, muc);
}


template<class CloudType>
Foam::forceSuSp Foam::ScaledForce<CloudType>::calcNonCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    return factor_*model_->calcCoupled(p, td, dt, mass, Re, muc);
}


template<class CloudType>
Foam::scalar Foam::ScaledForce<CloudType>::massAdd
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar mass
) const
{
    return factor_*model_->massAdd(p, td, mass);
}


// ************************************************************************* //
