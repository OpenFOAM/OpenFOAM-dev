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

#include "LiftForce.H"
#include "fvcCurl.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::LiftForce<CloudType>::LiftForce::Cl
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const vector& curlUc,
    const scalar Re,
    const scalar muc
) const
{
    // dummy
    return 0.0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::LiftForce<CloudType>::LiftForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& forceType
)
:
    ParticleForce<CloudType>(owner, mesh, dict, forceType, true),
    UName_(this->coeffs().template lookupOrDefault<word>("U", "U")),
    curlUcInterpPtr_(nullptr)
{}


template<class CloudType>
Foam::LiftForce<CloudType>::LiftForce(const LiftForce& lf)
:
    ParticleForce<CloudType>(lf),
    UName_(lf.UName_),
    curlUcInterpPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::LiftForce<CloudType>::~LiftForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::LiftForce<CloudType>::cacheFields(const bool store)
{
    static word fName("curlUcDt");

    bool fieldExists = this->mesh().template foundObject<volVectorField>(fName);

    if (store)
    {
        if (!fieldExists)
        {
            const volVectorField& Uc = this->mesh().template
                lookupObject<volVectorField>(UName_);

            volVectorField* curlUcPtr =
                new volVectorField(fName, fvc::curl(Uc));

            curlUcPtr->store();
        }

        const volVectorField& curlUc = this->mesh().template
            lookupObject<volVectorField>(fName);

        curlUcInterpPtr_.reset
        (
            interpolation<vector>::New
            (
                this->owner().solution().interpolationSchemes(),
                curlUc
            ).ptr()
        );
    }
    else
    {
        curlUcInterpPtr_.clear();

        if (fieldExists)
        {
            const volVectorField& curlUc = this->mesh().template
                lookupObject<volVectorField>(fName);

            const_cast<volVectorField&>(curlUc).checkOut();
        }
    }
}


template<class CloudType>
Foam::forceSuSp Foam::LiftForce<CloudType>::calcCoupled
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

    vector curlUc =
        curlUcInterp().interpolate(p.coordinates(), p.currentTetIndices());

    scalar Cl = this->Cl(p, td, curlUc, Re, muc);

    value.Su() = mass/p.rho()*td.rhoc()*Cl*((td.Uc() - p.U())^curlUc);

    return value;
}


// ************************************************************************* //
