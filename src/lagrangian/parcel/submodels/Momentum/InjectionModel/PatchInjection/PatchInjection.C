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

#include "PatchInjection.H"
#include "distribution.H"

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class CloudType>
void Foam::PatchInjection<CloudType>::preInject
(
    typename CloudType::parcelType::trackingData& td
)
{
    InjectionModel<CloudType>::preInject(td);

    if (U0Name_ == word::null)
    {
        U0InterpPtr_.clear();
    }
    else if (U0Name_ == this->owner().U().name())
    {
        U0InterpPtr_ = tmpNrc<interpolation<vector>>(td.UInterp());
    }
    else
    {
        U0InterpPtr_ =
            tmpNrc<interpolation<vector>>
            (
                interpolation<vector>::New
                (
                    this->owner().solution().interpolationSchemes(),
                    this->owner().mesh().template lookupObject<volVectorField>
                    (
                        U0Name_
                    )
                ).ptr()
            );
    }
}


template<class CloudType>
void Foam::PatchInjection<CloudType>::postInject
(
    const label parcelsAdded,
    const scalar massAdded,
    typename CloudType::parcelType::trackingData& td
)
{
    InjectionModel<CloudType>::postInject(parcelsAdded, massAdded, td);

    U0InterpPtr_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchInjection<CloudType>::PatchInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    patchInjectionBase(owner.mesh(), this->coeffDict().lookup("patchName")),
    duration_(this->readDuration(dict, owner)),
    massFlowRate_(this->readMassFlowRate(dict, owner, duration_)),
    parcelsPerSecond_(this->readParcelsPerSecond(dict, owner)),
    U0_(vector::uniform(NaN)),
    U0Name_(word::null),
    U0InterpPtr_(nullptr),
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
    ITstream& is = this->coeffDict().lookup("U0");

    token t(is);
    is.putBack(t);

    if (t.isWord())
    {
        U0Name_ = word(is);
    }
    else
    {
        U0_ = vector(is);
    }
}


template<class CloudType>
Foam::PatchInjection<CloudType>::PatchInjection
(
    const PatchInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    patchInjectionBase(im),
    duration_(im.duration_),
    massFlowRate_(im.massFlowRate_, false),
    parcelsPerSecond_(im.parcelsPerSecond_, false),
    U0_(im.U0_),
    U0Name_(word::null),
    U0InterpPtr_(nullptr),
    sizeDistribution_(im.sizeDistribution_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchInjection<CloudType>::~PatchInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::PatchInjection<CloudType>::topoChange()
{
    patchInjectionBase::topoChange(this->owner().mesh());
}


template<class CloudType>
Foam::scalar Foam::PatchInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
Foam::scalar Foam::PatchInjection<CloudType>::nParcelsToInject
(
    const scalar t0,
    const scalar t1
)
{
    if (t1 >= 0 && t0 < duration_)
    {
        return parcelsPerSecond_->integral(max(t0, 0), min(t1, duration_));
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::PatchInjection<CloudType>::massToInject
(
    const scalar t0,
    const scalar t1
)
{
    if (t1 >= 0 && t0 < duration_)
    {
        return massFlowRate_->integral(max(t0, 0), min(t1, duration_));
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
void Foam::PatchInjection<CloudType>::setPositionAndCell
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
void Foam::PatchInjection<CloudType>::setProperties
(
    const label,
    const label,
    const scalar,
    typename CloudType::parcelType::trackingData& td,
    typename CloudType::parcelType& parcel
)
{
    // set particle velocity
    parcel.U() =
        U0InterpPtr_.valid()
      ? U0InterpPtr_().interpolate
        (
            parcel.coordinates(),
            parcel.currentTetIndices(td.mesh)
        )
      : U0_;

    // set particle diameter
    parcel.d() = sizeDistribution_->sample();
}


template<class CloudType>
bool Foam::PatchInjection<CloudType>::fullyDescribed() const
{
    return false;
}


// ************************************************************************* //
