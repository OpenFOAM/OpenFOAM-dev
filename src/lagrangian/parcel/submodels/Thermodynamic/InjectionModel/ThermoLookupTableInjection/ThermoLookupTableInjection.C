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

#include "ThermoLookupTableInjection.H"
#include "scalarIOList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermoLookupTableInjection<CloudType>::ThermoLookupTableInjection
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
    injectorCoordinates_.setSize(injectors_.size());
    injectorCells_.setSize(injectors_.size());
    injectorTetFaces_.setSize(injectors_.size());
    injectorTetPts_.setSize(injectors_.size());

    topoChange();
}


template<class CloudType>
Foam::ThermoLookupTableInjection<CloudType>::ThermoLookupTableInjection
(
    const ThermoLookupTableInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    inputFileName_(im.inputFileName_),
    duration_(im.duration_),
    parcelsPerSecond_(im.parcelsPerSecond_, false),
    randomise_(im.randomise_),
    injectors_(im.injectors_),
    injectorCoordinates_(im.injectorCoordinates_),
    injectorCells_(im.injectorCells_),
    injectorTetFaces_(im.injectorTetFaces_),
    injectorTetPts_(im.injectorTetPts_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermoLookupTableInjection<CloudType>::~ThermoLookupTableInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ThermoLookupTableInjection<CloudType>::topoChange()
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
Foam::scalar Foam::ThermoLookupTableInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
Foam::scalar Foam::ThermoLookupTableInjection<CloudType>::nParcelsToInject
(
    const scalar t0,
    const scalar t1
)
{
    if (t1 >= 0 && t0 < duration_)
    {
        return
            injectorCells_.size()
           *parcelsPerSecond_->integral(max(t0, 0), min(t1, duration_));
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::ThermoLookupTableInjection<CloudType>::massToInject
(
    const scalar t0,
    const scalar t1
)
{
    scalar mass = 0;

    if (t1 >= 0 && t0 < duration_)
    {
        forAll(injectors_, i)
        {
            mass += injectors_[i].mDot()*(min(t1, duration_) - max(t0, 0));
        }
    }

    return mass;
}


template<class CloudType>
void Foam::ThermoLookupTableInjection<CloudType>::setPositionAndCell
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
        randomGenerator& rndGen = this->owner().rndGen();
        injectorI = rndGen.sampleAB<label>(0, injectorCells_.size());
    }
    else
    {
        injectorI = parcelI*injectorCells_.size()/nParcels;
    }

    coordinates = injectorCoordinates_[injectorI];
    celli = injectorCells_[injectorI];
    tetFacei = injectorTetFaces_[injectorI];
    tetPti = injectorTetPts_[injectorI];
}


template<class CloudType>
void Foam::ThermoLookupTableInjection<CloudType>::setProperties
(
    const label parcelI,
    const label nParcels,
    const scalar,
    typename CloudType::parcelType::trackingData& td,
    typename CloudType::parcelType& parcel
)
{
    label injectorI = parcelI*injectorCells_.size()/nParcels;

    // set particle velocity
    parcel.U() = injectors_[injectorI].U();

    // set particle diameter
    parcel.d() = injectors_[injectorI].d();

    // set particle density
    parcel.rho() = injectors_[injectorI].rho();

    // set particle temperature
    parcel.T() = injectors_[injectorI].T();

    // set particle specific heat capacity
    parcel.Cp() = injectors_[injectorI].Cp();
}


template<class CloudType>
bool Foam::ThermoLookupTableInjection<CloudType>::fullyDescribed() const
{
    return true;
}


// ************************************************************************* //
