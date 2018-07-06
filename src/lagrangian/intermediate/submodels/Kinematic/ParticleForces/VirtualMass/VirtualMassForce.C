/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "VirtualMassForce.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VirtualMassForce<CloudType>::VirtualMassForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& forceType
)
:
    PressureGradientForce<CloudType>(owner, mesh, dict, forceType),
    Cvm_(readScalar(this->coeffs().lookup("Cvm")))
{}


template<class CloudType>
Foam::VirtualMassForce<CloudType>::VirtualMassForce
(
    const VirtualMassForce& vmf
)
:
    PressureGradientForce<CloudType>(vmf),
    Cvm_(vmf.Cvm_)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VirtualMassForce<CloudType>::~VirtualMassForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::VirtualMassForce<CloudType>::cacheFields(const bool store)
{
    PressureGradientForce<CloudType>::cacheFields(store);
}


template<class CloudType>
Foam::forceSuSp Foam::VirtualMassForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value =
        PressureGradientForce<CloudType>::calcCoupled(p, td, dt, mass, Re, muc);

    value.Su() *= Cvm_;

    return value;
}


template<class CloudType>
Foam::scalar Foam::VirtualMassForce<CloudType>::massAdd
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar mass
) const
{
    return mass*td.rhoc()/p.rho()*Cvm_;
}


// ************************************************************************* //
