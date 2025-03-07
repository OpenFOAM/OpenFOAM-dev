/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    volumeRatio_
    (
        this->coeffDict().found("concentration")
     || this->coeffDict().found("volumeRatio")
      ? Function1<scalar>::New
        (
            this->coeffDict().found("volumeRatio")
          ? "volumeRatio"
          : "concentration",
            this->owner().db().time().userUnits(),
            dimless,
            this->coeffDict()
        )
      : autoPtr<Function1<scalar>>()
    ),
    massRatio_
    (
        this->coeffDict().found("massRatio")
      ? Function1<scalar>::New
        (
            "massRatio",
            this->owner().db().time().userUnits(),
            dimless,
            this->coeffDict()
        )
      : autoPtr<Function1<scalar>>()
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
{
    if (volumeRatio_.valid() && massRatio_.valid())
    {
        FatalIOErrorInFunction(this->coeffDict())
            << "keywords volumeRatio (or concentration) and "
            << "massRatio both defined in dictionary "
            << this->coeffDict().name() << exit(FatalIOError);
    }

    if (!volumeRatio_.valid() && !massRatio_.valid())
    {
        FatalIOErrorInFunction(this->coeffDict())
            << "keyword volumeRatio or massRatio is "
            << "undefined in dictionary "
            << this->coeffDict().name() << exit(FatalIOError);
    }
}


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
    volumeRatio_(im.volumeRatio_, false),
    massRatio_(im.massRatio_, false),
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
Foam::scalar Foam::PatchFlowRateInjection<CloudType>::volumetricFlowRate() const
{
   const polyMesh& mesh = this->owner().mesh();

    const surfaceScalarField& phi =
        mesh.lookupObject<surfaceScalarField>(phiName_);
    const scalarField& phip = phi.boundaryField()[patchId_];

    scalar flowRateIn;
    if (phi.dimensions() == dimVolumetricFlux)
    {
        flowRateIn = max(scalar(0), -sum(phip));
    }
    else
    {
        const volScalarField& rho =
            mesh.lookupObject<volScalarField>(rhoName_);
        const scalarField& rhop = rho.boundaryField()[patchId_];

        flowRateIn = max(scalar(0), -sum(phip/rhop));
    }

    reduce(flowRateIn, sumOp<scalar>());

    return flowRateIn;
}


template<class CloudType>
Foam::scalar Foam::PatchFlowRateInjection<CloudType>::massFlowRate() const
{
   const polyMesh& mesh = this->owner().mesh();

    const surfaceScalarField& phi =
        mesh.lookupObject<surfaceScalarField>(phiName_);
    const scalarField& phip = phi.boundaryField()[patchId_];

    scalar flowRateIn;
    if (phi.dimensions() == dimVolumetricFlux)
    {
        const volScalarField& rho =
            mesh.lookupObject<volScalarField>(rhoName_);
        const scalarField& rhop = rho.boundaryField()[patchId_];

        flowRateIn = max(scalar(0), -sum(phip*rhop));
    }
    else
    {
        flowRateIn = max(scalar(0), -sum(phip));
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
        const scalar tMid = (t0 + t1)/2;

        const scalar volumeFlowRateToInject =
            volumeRatio_.valid()
          ? volumeRatio_->value(tMid)
           *volumetricFlowRate()
          : massRatio_->value(tMid)
           *massFlowRate()
           /this->owner().constProps().rho0();

        return
            parcelConcentration_
           *volumeFlowRateToInject
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
        const scalar tMid = (t0 + t1)/2;

        const scalar massFlowRateToInject =
            volumeRatio_.valid()
          ? volumeRatio_->value(tMid)
           *volumetricFlowRate()
           *this->owner().constProps().rho0()
          : massRatio_->value(tMid)
           *massFlowRate();

        return
            massFlowRateToInject
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
