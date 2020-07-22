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

#include "PressureGradientForce.H"
#include "fvcDdt.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PressureGradientForce<CloudType>::PressureGradientForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& forceType
)
:
    ParticleForce<CloudType>(owner, mesh, dict, forceType, true),
    UName_(this->coeffs().template lookupOrDefault<word>("U", "U")),
    DUcDtInterpPtr_(nullptr)
{}


template<class CloudType>
Foam::PressureGradientForce<CloudType>::PressureGradientForce
(
    const PressureGradientForce& pgf
)
:
    ParticleForce<CloudType>(pgf),
    UName_(pgf.UName_),
    DUcDtInterpPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PressureGradientForce<CloudType>::~PressureGradientForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::PressureGradientForce<CloudType>::cacheFields(const bool store)
{
    static word fName("DUcDt");

    bool fieldExists = this->mesh().template foundObject<volVectorField>(fName);

    if (store)
    {
        if (!fieldExists)
        {
            const volVectorField& Uc = this->mesh().template
                lookupObject<volVectorField>(UName_);

            volVectorField* DUcDtPtr = new volVectorField
            (
                fName,
                fvc::ddt(Uc) + (Uc & fvc::grad(Uc))
            );

            DUcDtPtr->store();
        }

        const volVectorField& DUcDt = this->mesh().template
            lookupObject<volVectorField>(fName);

        DUcDtInterpPtr_.reset
        (
            interpolation<vector>::New
            (
                this->owner().solution().interpolationSchemes(),
                DUcDt
            ).ptr()
        );
    }
    else
    {
        DUcDtInterpPtr_.clear();

        if (fieldExists)
        {
            const volVectorField& DUcDt = this->mesh().template
                lookupObject<volVectorField>(fName);

            const_cast<volVectorField&>(DUcDt).checkOut();
        }
    }
}


template<class CloudType>
Foam::forceSuSp Foam::PressureGradientForce<CloudType>::calcCoupled
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

    vector DUcDt =
        DUcDtInterp().interpolate(p.coordinates(), p.currentTetIndices());

    value.Su() = mass*td.rhoc()/p.rho()*DUcDt;

    return value;
}


template<class CloudType>
Foam::scalar Foam::PressureGradientForce<CloudType>::massAdd
(
    const typename CloudType::parcelType&,
    const typename CloudType::parcelType::trackingData& td,
    const scalar
) const
{
    return 0.0;
}


// ************************************************************************* //
