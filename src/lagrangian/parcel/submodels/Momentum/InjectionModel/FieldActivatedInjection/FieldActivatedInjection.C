/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "FieldActivatedInjection.H"
#include "volFields.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FieldActivatedInjection<CloudType>::FieldActivatedInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    factor_(this->coeffDict().template lookup<scalar>("factor")),
    referenceField_
    (
        owner.db().objectRegistry::template lookupObject<volScalarField>
        (
            this->coeffDict().lookup("referenceField")
        )
    ),
    thresholdField_
    (
        owner.db().objectRegistry::template lookupObject<volScalarField>
        (
            this->coeffDict().lookup("thresholdField")
        )
    ),
    positionsFile_(this->coeffDict().lookup("positionsFile")),
    positions_
    (
        IOobject
        (
            positionsFile_,
            owner.db().time().constant(),
            owner.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    injectorCoordinates_(positions_.size()),
    injectorCells_(positions_.size()),
    injectorTetFaces_(positions_.size()),
    injectorTetPts_(positions_.size()),
    massTotal_(this->readMassTotal(dict, owner)),
    nParcelsPerInjector_
    (
        this->coeffDict().template lookup<label>("parcelsPerInjector")
    ),
    nParcelsInjected_(positions_.size(), 0),
    U0_(this->coeffDict().lookup("U0")),
    diameters_(positions_.size()),
    sizeDistribution_
    (
        distribution::New
        (
            this->coeffDict().subDict("sizeDistribution"),
            owner.rndGen(),
            this->sizeSampleQ()
        )
    )
{
    // Construct parcel diameters - one per injector cell
    forAll(diameters_, i)
    {
        diameters_[i] = sizeDistribution_->sample();
    }

    topoChange();
}


template<class CloudType>
Foam::FieldActivatedInjection<CloudType>::FieldActivatedInjection
(
    const FieldActivatedInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    factor_(im.factor_),
    referenceField_(im.referenceField_),
    thresholdField_(im.thresholdField_),
    positionsFile_(im.positionsFile_),
    positions_(im.positions_),
    injectorCoordinates_(im.injectorCoordinates_),
    injectorCells_(im.injectorCells_),
    injectorTetFaces_(im.injectorTetFaces_),
    injectorTetPts_(im.injectorTetPts_),
    massTotal_(im.massTotal_),
    nParcelsPerInjector_(im.nParcelsPerInjector_),
    nParcelsInjected_(im.nParcelsInjected_),
    U0_(im.U0_),
    diameters_(im.diameters_),
    sizeDistribution_(im.sizeDistribution_().clone().ptr())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FieldActivatedInjection<CloudType>::~FieldActivatedInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::FieldActivatedInjection<CloudType>::topoChange()
{
    // Set/cache the injector cells
    forAll(positions_, i)
    {
        this->findCellAtPosition
        (
            positions_[i],
            injectorCoordinates_[i],
            injectorCells_[i],
            injectorTetFaces_[i],
            injectorTetPts_[i]
        );
    }
}


template<class CloudType>
Foam::scalar Foam::FieldActivatedInjection<CloudType>::timeEnd() const
{
    return great;
}


template<class CloudType>
Foam::label Foam::FieldActivatedInjection<CloudType>::nParcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    if (sum(nParcelsInjected_) < nParcelsPerInjector_*positions_.size())
    {
        return positions_.size();
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::FieldActivatedInjection<CloudType>::massToInject
(
    const scalar time0,
    const scalar time1
)
{
    if (sum(nParcelsInjected_) < nParcelsPerInjector_*positions_.size())
    {
        return massTotal_/nParcelsPerInjector_;
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
void Foam::FieldActivatedInjection<CloudType>::setPositionAndCell
(
    const label parceli,
    const label,
    const scalar,
    barycentric& coordinates,
    label& celli,
    label& tetFacei,
    label& tetPti,
    label& facei
)
{
    const label injectorCelli = injectorCells_[parceli];

    if
    (
        nParcelsInjected_[parceli] < nParcelsPerInjector_
     && factor_*referenceField_[injectorCelli] > thresholdField_[injectorCelli]
    )
    {
        coordinates = injectorCoordinates_[parceli];
        celli = injectorCells_[parceli];
        tetFacei = injectorTetFaces_[parceli];
        tetPti = injectorTetPts_[parceli];

        nParcelsInjected_[parceli]++;
    }
}


template<class CloudType>
void Foam::FieldActivatedInjection<CloudType>::setProperties
(
    const label parceli,
    const label,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    // set particle velocity
    parcel.U() = U0_;

    // set particle diameter
    parcel.d() = diameters_[parceli];
}


template<class CloudType>
bool Foam::FieldActivatedInjection<CloudType>::fullyDescribed() const
{
    return false;
}


// ************************************************************************* //
