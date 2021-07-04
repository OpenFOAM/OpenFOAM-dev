/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "ReactingMultiphaseCloud.H"

#include "DevolatilisationModel.H"
#include "SurfaceReactionModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::ReactingMultiphaseCloud<CloudType>::setModels()
{
    devolatilisationModel_.reset
    (
        DevolatilisationModel<ReactingMultiphaseCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    surfaceReactionModel_.reset
    (
        SurfaceReactionModel<ReactingMultiphaseCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
}


template<class CloudType>
void Foam::ReactingMultiphaseCloud<CloudType>::cloudReset
(
    ReactingMultiphaseCloud<CloudType>& c
)
{
    CloudType::cloudReset(c);

    devolatilisationModel_.reset(c.devolatilisationModel_.ptr());
    surfaceReactionModel_.reset(c.surfaceReactionModel_.ptr());

    dMassDevolatilisation_ = c.dMassDevolatilisation_;
    dMassSurfaceReaction_ = c.dMassSurfaceReaction_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ReactingMultiphaseCloud<CloudType>::ReactingMultiphaseCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const fluidThermo& carrierThermo,
    const bool readFields
)
:
    CloudType(cloudName, rho, U, g, carrierThermo, false),
    cloudCopyPtr_(nullptr),
    constProps_(this->particleProperties()),
    devolatilisationModel_(nullptr),
    surfaceReactionModel_(nullptr),
    dMassDevolatilisation_(0.0),
    dMassSurfaceReaction_(0.0)
{
    setModels();

    if (readFields)
    {
        parcelType::readFields(*this, this->composition());
        this->deleteLostParticles();
    }

    if (this->solution().resetSourcesOnStartup())
    {
        resetSourceTerms();
    }
}


template<class CloudType>
Foam::ReactingMultiphaseCloud<CloudType>::ReactingMultiphaseCloud
(
    ReactingMultiphaseCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    cloudCopyPtr_(nullptr),
    constProps_(c.constProps_),
    devolatilisationModel_(c.devolatilisationModel_->clone()),
    surfaceReactionModel_(c.surfaceReactionModel_->clone()),
    dMassDevolatilisation_(c.dMassDevolatilisation_),
    dMassSurfaceReaction_(c.dMassSurfaceReaction_)
{}


template<class CloudType>
Foam::ReactingMultiphaseCloud<CloudType>::ReactingMultiphaseCloud
(
    const fvMesh& mesh,
    const word& name,
    const ReactingMultiphaseCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    cloudCopyPtr_(nullptr),
    constProps_(),
    devolatilisationModel_(nullptr),
    surfaceReactionModel_(nullptr),
    dMassDevolatilisation_(0.0),
    dMassSurfaceReaction_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ReactingMultiphaseCloud<CloudType>::~ReactingMultiphaseCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ReactingMultiphaseCloud<CloudType>::setParcelThermoProperties
(
    parcelType& parcel,
    const scalar lagrangianDt
)
{
    CloudType::setParcelThermoProperties(parcel, lagrangianDt);

    label idGas = this->composition().idGas();
    label idLiquid = this->composition().idLiquid();
    label idSolid = this->composition().idSolid();

    parcel.YGas() = this->composition().Y0(idGas);
    parcel.YLiquid() = this->composition().Y0(idLiquid);
    parcel.YSolid() = this->composition().Y0(idSolid);
}


template<class CloudType>
void Foam::ReactingMultiphaseCloud<CloudType>::checkParcelProperties
(
    parcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    CloudType::checkParcelProperties(parcel, lagrangianDt, fullyDescribed);

    if (fullyDescribed)
    {
        label idGas = this->composition().idGas();
        label idLiquid = this->composition().idLiquid();
        label idSolid = this->composition().idSolid();

        this->checkSuppliedComposition
        (
            parcel.YGas(),
            this->composition().Y0(idGas),
            "YGas"
        );
        this->checkSuppliedComposition
        (
            parcel.YLiquid(),
            this->composition().Y0(idLiquid),
            "YLiquid"
        );
        this->checkSuppliedComposition
        (
            parcel.YSolid(),
            this->composition().Y0(idSolid),
            "YSolid"
        );
    }
}


template<class CloudType>
void Foam::ReactingMultiphaseCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<ReactingMultiphaseCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::ReactingMultiphaseCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::ReactingMultiphaseCloud<CloudType>::resetSourceTerms()
{
    CloudType::resetSourceTerms();
}


template<class CloudType>
void Foam::ReactingMultiphaseCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::trackingData td(*this);

        this->solve(*this, td);
    }
}


template<class CloudType>
void Foam::ReactingMultiphaseCloud<CloudType>::autoMap
(
    const mapPolyMesh& mapper
)
{
    Cloud<parcelType>::autoMap(mapper);

    this->updateMesh();
}


template<class CloudType>
void Foam::ReactingMultiphaseCloud<CloudType>::info()
{
    CloudType::info();

    this->devolatilisation().info(Info);
    this->surfaceReaction().info(Info);
}


template<class CloudType>
void Foam::ReactingMultiphaseCloud<CloudType>::writeFields() const
{
    if (this->compositionModel_.valid())
    {
        CloudType::particleType::writeFields(*this, this->composition());
    }
    else
    {
        CloudType::particleType::writeFields(*this);
    }
}


// ************************************************************************* //
