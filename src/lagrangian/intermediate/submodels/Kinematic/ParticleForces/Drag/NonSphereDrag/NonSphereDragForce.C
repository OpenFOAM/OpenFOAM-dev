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

#include "NonSphereDragForce.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::NonSphereDragForce<CloudType>::CdRe(const scalar Re) const
{
    return 24.0*(1.0 + a_*pow(Re, b_)) + Re*c_/(1 + d_/(Re + ROOTVSMALL));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NonSphereDragForce<CloudType>::NonSphereDragForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    phi_(readScalar(this->coeffs().lookup("phi"))),
    a_(exp(2.3288 - 6.4581*phi_ + 2.4486*sqr(phi_))),
    b_(0.0964 + 0.5565*phi_),
    c_(exp(4.9050 - 13.8944*phi_ + 18.4222*sqr(phi_) - 10.2599*pow3(phi_))),
    d_(exp(1.4681 + 12.2584*phi_ - 20.7322*sqr(phi_) + 15.8855*pow3(phi_)))
{
    if (phi_ <= 0 || phi_ > 1)
    {
        FatalErrorIn
        (
            "NonSphereDrag<CloudType>::NonSphereDrag"
            "("
                "const dictionary&, "
                "CloudType&"
            ")"
        )   << "Ratio of surface of sphere having same volume as particle to "
            << "actual surface area of particle (phi) must be greater than 0 "
            << "and less than or equal to 1" << exit(FatalError);
    }
}


template<class CloudType>
Foam::NonSphereDragForce<CloudType>::NonSphereDragForce
(
    const NonSphereDragForce<CloudType>& df
)
:
    ParticleForce<CloudType>(df),
    phi_(df.phi_),
    a_(df.a_),
    b_(df.b_),
    c_(df.c_),
    d_(df.d_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NonSphereDragForce<CloudType>::~NonSphereDragForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::NonSphereDragForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(vector::zero, 0.0);

    value.Sp() = mass*0.75*muc*CdRe(Re)/(p.rho()*sqr(p.d()));

    return value;
}


// ************************************************************************* //
