/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "ParamagneticForce.H"
#include "demandDrivenData.H"
#include "electromagneticConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParamagneticForce<CloudType>::ParamagneticForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    HdotGradHName_
    (
        this->coeffs().template lookupOrDefault<word>("HdotGradH", "HdotGradH")
    ),
    HdotGradHInterpPtr_(nullptr),
    magneticSusceptibility_
    (
        readScalar(this->coeffs().lookup("magneticSusceptibility"))
    )
{}


template<class CloudType>
Foam::ParamagneticForce<CloudType>::ParamagneticForce
(
    const ParamagneticForce& pf
)
:
    ParticleForce<CloudType>(pf),
    HdotGradHName_(pf.HdotGradHName_),
    HdotGradHInterpPtr_(pf.HdotGradHInterpPtr_),
    magneticSusceptibility_(pf.magneticSusceptibility_)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParamagneticForce<CloudType>::~ParamagneticForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParamagneticForce<CloudType>::cacheFields(const bool store)
{
    if (store)
    {
        const volVectorField& HdotGradH =
            this->mesh().template lookupObject<volVectorField>(HdotGradHName_);

        HdotGradHInterpPtr_ = interpolation<vector>::New
        (
            this->owner().solution().interpolationSchemes(),
            HdotGradH
        ).ptr();
    }
    else
    {
        deleteDemandDrivenData(HdotGradHInterpPtr_);
    }
}


template<class CloudType>
Foam::forceSuSp Foam::ParamagneticForce<CloudType>::calcNonCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(Zero, 0.0);

    const interpolation<vector>& HdotGradHInterp = *HdotGradHInterpPtr_;

    value.Su()=
        mass*3.0*constant::electromagnetic::mu0.value()/p.rho()
       *magneticSusceptibility_/(magneticSusceptibility_ + 3)
       *HdotGradHInterp.interpolate(p.coordinates(), p.currentTetIndices());

    return value;
}


// ************************************************************************* //
