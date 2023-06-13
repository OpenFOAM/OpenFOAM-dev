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

#include "MomentumLookupTableInjection.H"
#include "scalarIOList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MomentumLookupTableInjection<CloudType>::MomentumLookupTableInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    inputFileName_(this->coeffDict().lookup("inputFile")),
    duration_(this->readDuration(dict, owner)),
    parcelsPerSecond_(this->readParcelsPerSecond(dict, owner)),
    randomise_(readBool(this->coeffDict().lookup("randomise"))),
    injectors_
    (
        IOobject
        (
            inputFileName_,
            owner.db().time().constant(),
            owner.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    injectorCoordinates_(0),
    injectorCells_(0),
    injectorTetFaces_(0),
    injectorTetPts_(0)
{
    // Set/cache the injector cells
    injectorCells_.setSize(injectors_.size());
    injectorTetFaces_.setSize(injectors_.size());
    injectorTetPts_.setSize(injectors_.size());

    topoChange();
}


template<class CloudType>
Foam::MomentumLookupTableInjection<CloudType>::MomentumLookupTableInjection
(
    const MomentumLookupTableInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    inputFileName_(im.inputFileName_),
    duration_(im.duration_),
    parcelsPerSecond_(im.parcelsPerSecond_),
    randomise_(im.randomise_),
    injectors_(im.injectors_),
    injectorCoordinates_(im.injectorCoordinates_),
    injectorCells_(im.injectorCells_),
    injectorTetFaces_(im.injectorTetFaces_),
    injectorTetPts_(im.injectorTetPts_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MomentumLookupTableInjection<CloudType>::~MomentumLookupTableInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::MomentumLookupTableInjection<CloudType>::topoChange()
{
    // Set/cache the injector cells
    forAll(injectors_, i)
    {
        this->findCellAtPosition
        (
            injectors_[i].x(),
            injectorCoordinates_[i],
            injectorCells_[i],
            injectorTetFaces_[i],
            injectorTetPts_[i]
        );
    }
}


template<class CloudType>
Foam::scalar Foam::MomentumLookupTableInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
Foam::label Foam::MomentumLookupTableInjection<CloudType>::nParcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    if (time0 >= 0 && time0 < duration_)
    {
        return
            floor
            (
                injectorCells_.size()
               *parcelsPerSecond_.integral(time0, time1)
            );
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::MomentumLookupTableInjection<CloudType>::massToInject
(
    const scalar time0,
    const scalar time1
)
{
    scalar mass = 0;

    if (time0 >= 0 && time0 < duration_)
    {
        forAll(injectors_, i)
        {
            mass += injectors_[i].mDot()*(time1 - time0);
        }
    }

    return mass;
}


template<class CloudType>
void Foam::MomentumLookupTableInjection<CloudType>::setPositionAndCell
(
    const label parcelI,
    const label nParcels,
    const scalar time,
    barycentric& coordinates,
    label& celli,
    label& tetFacei,
    label& tetPti,
    label& facei
)
{
    label injectorI = 0;
    if (randomise_)
    {
        Random& rnd = this->owner().rndGen();
        injectorI = rnd.sampleAB<label>(0, injectorCells_.size());
    }
    else
    {
        injectorI = int64_t(parcelI)*int64_t(injectors_.size())/nParcels;
    }

    coordinates = injectorCoordinates_[injectorI];
    celli = injectorCells_[injectorI];
    tetFacei = injectorTetFaces_[injectorI];
    tetPti = injectorTetPts_[injectorI];
}


template<class CloudType>
void Foam::MomentumLookupTableInjection<CloudType>::setProperties
(
    const label parcelI,
    const label nParcels,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    label injectorI = int64_t(parcelI)*int64_t(injectors_.size())/nParcels;

    // set particle velocity
    parcel.U() = injectors_[injectorI].U();

    // set particle diameter
    parcel.d() = injectors_[injectorI].d();

    // set particle density
    parcel.rho() = injectors_[injectorI].rho();
}


template<class CloudType>
bool Foam::MomentumLookupTableInjection<CloudType>::fullyDescribed() const
{
    return true;
}


// ************************************************************************* //
