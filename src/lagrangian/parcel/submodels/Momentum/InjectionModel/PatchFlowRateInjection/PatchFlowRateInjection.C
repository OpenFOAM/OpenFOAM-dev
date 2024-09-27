/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "PatchFlowRateInjection.H"
#include "distribution.H"
#include "mathematicalConstants.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchFlowRateInjection<CloudType>::PatchFlowRateInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName,typeName),
    patchInjectionBase(owner.mesh(), this->coeffDict().lookup("patchName")),
    phiName_(this->coeffDict().template lookupOrDefault<word>("phi", "phi")),
    rhoName_(this->coeffDict().template lookupOrDefault<word>("rho", "rho")),
    duration_(this->readDuration(dict, owner)),
    concentration_
    (
        Function1<scalar>::New
        (
            "concentration",
            this->owner().db().time().userUnits(),
            dimless,
            this->coeffDict()
        )
    ),
    parcelConcentration_
    (
        this->coeffDict().template lookup<scalar>("parcelConcentration")
    ),
    sizeDistribution_
    (
        distribution::New
        (
            dimLength,
            this->coeffDict().subDict("sizeDistribution"),
            this->sizeSampleQ(),
            owner.rndGen().generator()
        )
    )
{}


template<class CloudType>
Foam::PatchFlowRateInjection<CloudType>::PatchFlowRateInjection
(
    const PatchFlowRateInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    patchInjectionBase(im),
    phiName_(im.phiName_),
    rhoName_(im.rhoName_),
    duration_(im.duration_),
    concentration_(im.concentration_, false),
    parcelConcentration_(im.parcelConcentration_),
    sizeDistribution_(im.sizeDistribution_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchFlowRateInjection<CloudType>::~PatchFlowRateInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::PatchFlowRateInjection<CloudType>::topoChange()
{
    patchInjectionBase::topoChange(this->owner().mesh());
}


template<class CloudType>
Foam::scalar Foam::PatchFlowRateInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
Foam::scalar Foam::PatchFlowRateInjection<CloudType>::flowRate() const
{
   const polyMesh& mesh = this->owner().mesh();

    const surfaceScalarField& phi =
        mesh.lookupObject<surfaceScalarField>(phiName_);

    const scalarField& phip = phi.boundaryField()[patchId_];

    scalar flowRateIn = 0.0;
    if (phi.dimensions() == dimVolumetricFlux)
    {
        flowRateIn = max(0.0, -sum(phip));
    }
    else
    {
        const volScalarField& rho =
            mesh.lookupObject<volScalarField>(rhoName_);
        const scalarField& rhop = rho.boundaryField()[patchId_];

        flowRateIn = max(0.0, -sum(phip/rhop));
    }

    reduce(flowRateIn, sumOp<scalar>());

    return flowRateIn;
}


template<class CloudType>
Foam::scalar Foam::PatchFlowRateInjection<CloudType>::nParcelsToInject
(
    const scalar t0,
    const scalar t1
)
{
    if (t1 >= 0 && t0 < duration_)
    {
        return
            parcelConcentration_
           *concentration_->value(0.5*(t0 + t1))
           *flowRate()
           *(min(t1, duration_) - max(t0, 0));
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::PatchFlowRateInjection<CloudType>::massToInject
(
    const scalar t0,
    const scalar t1
)
{
    if (t1 >= 0 && t0 < duration_)
    {
        return
            this->owner().constProps().rho0()
           *concentration_->value(0.5*(t0 + t1))
           *flowRate()
           *(min(t1, duration_) - max(t0, 0));
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
void Foam::PatchFlowRateInjection<CloudType>::setPositionAndCell
(
    const label,
    const label,
    const scalar,
    barycentric& coordinates,
    label& celli,
    label& tetFacei,
    label& tetPti,
    label& facei
)
{
    patchInjectionBase::setPositionAndCell
    (
        this->owner().mesh(),
        this->owner().rndGen(),
        coordinates,
        celli,
        tetFacei,
        tetPti,
        facei
    );
}


template<class CloudType>
void Foam::PatchFlowRateInjection<CloudType>::setProperties
(
    const label,
    const label,
    const scalar,
    typename CloudType::parcelType::trackingData& td,
    typename CloudType::parcelType& parcel
)
{
    // Set particle velocity to carrier velocity
    parcel.U() = this->owner().U()[parcel.cell()];

    // Set particle diameter
    parcel.d() = sizeDistribution_->sample();
}


template<class CloudType>
bool Foam::PatchFlowRateInjection<CloudType>::fullyDescribed() const
{
    return false;
}


// ************************************************************************* //
