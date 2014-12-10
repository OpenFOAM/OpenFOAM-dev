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

#include "ParticleForce.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleForce<CloudType>::ParticleForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& forceType,
    const bool readCoeffs
)
:
    owner_(owner),
    mesh_(mesh),
    coeffs_(readCoeffs ? dict : dictionary::null)
{
    if (readCoeffs && (coeffs_.dictName() != forceType))
    {
        FatalIOErrorIn
        (
            "Foam::ParticleForce<CloudType>::ParticleForce"
            "("
                "CloudType&, "
                "const fvMesh&, "
                "const dictionary&, "
                "const word&, "
                "const bool"
            ")",
            dict
        )   << "Force " << forceType << " must be specified as a dictionary"
            << exit(FatalIOError);
    }
}


template<class CloudType>
Foam::ParticleForce<CloudType>::ParticleForce(const ParticleForce& pf)
:
    owner_(pf.owner_),
    mesh_(pf.mesh_),
    coeffs_(pf.coeffs_)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleForce<CloudType>::~ParticleForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleForce<CloudType>::cacheFields(const bool store)
{}


template<class CloudType>
Foam::forceSuSp Foam::ParticleForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType&,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value;
    value.Su() = vector::zero;
    value.Sp() = 0.0;

    return value;
}


template<class CloudType>
Foam::forceSuSp Foam::ParticleForce<CloudType>::calcNonCoupled
(
    const typename CloudType::parcelType&,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value;
    value.Su() = vector::zero;
    value.Sp() = 0.0;

    return value;
}


template<class CloudType>
Foam::scalar Foam::ParticleForce<CloudType>::massAdd
(
    const typename CloudType::parcelType& p,
    const scalar mass
) const
{
    return 0.0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ParticleForceNew.C"

// ************************************************************************* //
