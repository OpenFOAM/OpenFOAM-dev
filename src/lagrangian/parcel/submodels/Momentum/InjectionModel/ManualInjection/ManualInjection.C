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

#include "ManualInjection.H"
#include "mathematicalConstants.H"
#include "PackedBoolList.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ManualInjection<CloudType>::ManualInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
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
    diameters_(positions_.size()),
    injectorCells_(positions_.size(), -1),
    injectorTetFaces_(positions_.size(), -1),
    injectorTetPts_(positions_.size(), -1),
    U0_(this->coeffDict().lookup("U0")),
    sizeDistribution_
    (
        distributionModel::New
        (
            this->coeffDict().subDict("sizeDistribution"),
            owner.rndGen()
        )
    ),
    ignoreOutOfBounds_
    (
        this->coeffDict().lookupOrDefault("ignoreOutOfBounds", false)
    )
{
    updateMesh();

    // Construct parcel diameters
    forAll(diameters_, i)
    {
        diameters_[i] = sizeDistribution_->sample();
    }

    // Determine volume of particles to inject
    this->volumeTotal_ = sum(pow3(diameters_))*pi/6.0;
}


template<class CloudType>
Foam::ManualInjection<CloudType>::ManualInjection
(
    const ManualInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    positionsFile_(im.positionsFile_),
    positions_(im.positions_),
    diameters_(im.diameters_),
    injectorCells_(im.injectorCells_),
    injectorTetFaces_(im.injectorTetFaces_),
    injectorTetPts_(im.injectorTetPts_),
    U0_(im.U0_),
    sizeDistribution_(im.sizeDistribution_().clone().ptr()),
    ignoreOutOfBounds_(im.ignoreOutOfBounds_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ManualInjection<CloudType>::~ManualInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ManualInjection<CloudType>::updateMesh()
{
    label nRejected = 0;

    PackedBoolList keep(positions_.size(), true);

    forAll(positions_, pI)
    {
        if
        (
            !this->findCellAtPosition
            (
                injectorCells_[pI],
                injectorTetFaces_[pI],
                injectorTetPts_[pI],
                positions_[pI],
                !ignoreOutOfBounds_
            )
        )
        {
            keep[pI] = false;
            nRejected++;
        }
    }

    if (nRejected > 0)
    {
        inplaceSubset(keep, positions_);
        inplaceSubset(keep, diameters_);
        inplaceSubset(keep, injectorCells_);
        inplaceSubset(keep, injectorTetFaces_);
        inplaceSubset(keep, injectorTetPts_);

        Info<< "    " << nRejected
            << " particles ignored, out of bounds" << endl;
    }
}


template<class CloudType>
Foam::scalar Foam::ManualInjection<CloudType>::timeEnd() const
{
    // Injection is instantaneous - but allow for a finite interval to
    // avoid numerical issues when interval is zero
    return rootVSmall;
}


template<class CloudType>
Foam::label Foam::ManualInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ((0.0 >= time0) && (0.0 < time1))
    {
        return positions_.size();
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::ManualInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    // All parcels introduced at SOI
    if ((0.0 >= time0) && (0.0 < time1))
    {
        return this->volumeTotal_;
    }
    else
    {
        return 0.0;
    }
}


template<class CloudType>
void Foam::ManualInjection<CloudType>::setPositionAndCell
(
    const label parcelI,
    const label,
    const scalar,
    vector& position,
    label& cellOwner,
    label& tetFacei,
    label& tetPti
)
{
    position = positions_[parcelI];
    cellOwner = injectorCells_[parcelI];
    tetFacei = injectorTetFaces_[parcelI];
    tetPti = injectorTetPts_[parcelI];
}


template<class CloudType>
void Foam::ManualInjection<CloudType>::setProperties
(
    const label parcelI,
    const label,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    // set particle velocity
    parcel.U() = U0_;

    // set particle diameter
    parcel.d() = diameters_[parcelI];
}


template<class CloudType>
bool Foam::ManualInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::ManualInjection<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
